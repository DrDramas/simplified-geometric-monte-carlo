#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF2.h"
#include "TF3.h"
#include "TFile.h"
#include "TMath.h"

#include "TRandom2.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TStyle.h"

#include "TApplication.h"
#include "TCanvas.h"

#include "Math/ProbFunc.h"
#include "Math/DistFunc.h"

#include <atomic>
#include <iostream>
#include <fstream>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include <map>
#include <cassert>
#include <cstdlib>
#include <ctime>

using namespace std;

const double pi = TMath::Pi();

//============================================================================
// Configuration
//============================================================================
//
// All tunable parameters are read from a plain-text config file consisting
// of "key = value" lines. Lines starting with '#' and blank lines are
// ignored; inline '#' comments after a value are also stripped.
//
// Required keys:
//   nElectrons, nProtons, halfLife, aperatureRadius,
//   TKelvin, PTorr, fillGasAtomicRadius, beamIonAtomicRadius,
//   M_fillGas, M_beamIon
//
// Required list-valued keys (comma-separated; whitespace ignored):
//   numHalfLives, x_0, y_0, R  -- all four must have equal length N;
//                                  produces N beam-spot histograms.
//   De                         -- length M; produces M*9 electron-cloud
//                                  histograms (drift times 0..8 us, with
//                                  j=0 mapped to 0.1 us to avoid sigma=0).
//
// Optional keys (defaults preserve original behavior):
//   inputFile       (default: "tripleCoincidences.root")
//   inputHistName   (default: "roughBeamDistribution")
//   outputFile      (default: "AllHistograms.root")
//   nThreads        (default: 0 -- auto-detect hardware concurrency)
//
struct Config {
  // Run-size parameters
  long   nElectrons;           // number of electrons to sample per cloud
  long   nProtons;             // number of protons to sample per beam spot
  double halfLife;             // isotope half-life [s]
  double aperatureRadius;      // entrance aperture radius [mm]

  // Systematic-uncertainty parameters (now lists for sweeping).
  // The four beam-spot lists must have the same length N; index k of each
  // defines one beam-spot preset.
  vector<int>    numHalfLives;  // per-preset number of half-lives
  vector<double> x_0;           // per-preset x-origin of beam spot [mm]
  vector<double> y_0;           // per-preset y-origin of beam spot [mm]
  vector<double> R;             // per-preset radial width           [mm]
  // De is independent: one electron-cloud family is produced per entry.
  vector<double> De;            // electron diffusion coefficients

  // Chemistry parameters (used by CalculateDiffusivity)
  double TKelvin;              // temperature [K]
  double PTorr;                // gas pressure [Torr]
  double fillGasAtomicRadius;  // fill-gas atomic radius [pm]
  double beamIonAtomicRadius;  // beam-ion atomic radius [pm]
  double M_fillGas;            // molar mass of fill gas  [g/mol]
  double M_beamIon;            // molar mass of beam ion  [g/mol]

  // File/histogram names (optional in the config; see defaults in LoadConfig)
  string inputFile;            // input ROOT file path
  string inputHistName;        // name of the input TH1 to read
  string outputFile;           // output ROOT file path

  // Threading.  0 means "auto" (std::thread::hardware_concurrency()).
  int    nThreads = 0;
};

// Trim leading/trailing whitespace in place.
static void trim(string &s) {
  const char *ws = " \t\r\n";
  size_t a = s.find_first_not_of(ws);
  if (a == string::npos) { s.clear(); return; }
  size_t b = s.find_last_not_of(ws);
  s = s.substr(a, b - a + 1);
}

// Read "key = value" pairs from `path` into a map. Aborts on failure.
static map<string,string> ReadKeyValueFile(const string &path) {
  map<string,string> kv;
  ifstream in(path.c_str());
  if (!in.is_open()) {
    cerr << "ERROR: could not open config file '" << path << "'" << endl;
    exit(1);
  }
  string line;
  int lineNo = 0;
  while (getline(in, line)) {
    ++lineNo;
    // Strip comments starting with '#'
    size_t hash = line.find('#');
    if (hash != string::npos) line.erase(hash);
    trim(line);
    if (line.empty()) continue;

    size_t eq = line.find('=');
    if (eq == string::npos) {
      cerr << "WARNING: config line " << lineNo
           << " has no '=' and is ignored: " << line << endl;
      continue;
    }
    string key = line.substr(0, eq);
    string val = line.substr(eq + 1);
    trim(key); trim(val);
    if (key.empty()) continue;
    kv[key] = val;
  }
  return kv;
}

static string RequireStr(const map<string,string> &kv, const string &key) {
  map<string,string>::const_iterator it = kv.find(key);
  if (it == kv.end()) {
    cerr << "ERROR: required config key missing: " << key << endl;
    exit(1);
  }
  return it->second;
}
static string OptionalStr(const map<string,string> &kv,
                          const string &key,
                          const string &defaultVal) {
  map<string,string>::const_iterator it = kv.find(key);
  if (it == kv.end()) return defaultVal;
  return it->second;
}
static double RequireDouble(const map<string,string> &kv, const string &key) {
  string s = RequireStr(kv, key);
  istringstream iss(s);
  double v; iss >> v;
  if (iss.fail()) {
    cerr << "ERROR: config key '" << key << "' is not a number: " << s << endl;
    exit(1);
  }
  return v;
}
static long RequireLong(const map<string,string> &kv, const string &key) {
  return static_cast<long>(RequireDouble(kv, key));
}
static int RequireInt(const map<string,string> &kv, const string &key) {
  return static_cast<int>(RequireDouble(kv, key));
}
static int OptionalInt(const map<string,string> &kv,
                       const string &key, int defaultVal) {
  map<string,string>::const_iterator it = kv.find(key);
  if (it == kv.end()) return defaultVal;
  istringstream iss(it->second);
  double v; iss >> v;
  if (iss.fail()) {
    cerr << "ERROR: config key '" << key << "' is not a number: "
         << it->second << endl;
    exit(1);
  }
  return static_cast<int>(v);
}

static vector<string> SplitCSV(const string &s) {
  vector<string> out;
  string tok;
  istringstream iss(s);
  while (getline(iss, tok, ',')) {
    trim(tok);
    if (!tok.empty()) out.push_back(tok);
  }
  return out;
}

static vector<double> RequireDoubleList(const map<string,string> &kv,
                                        const string &key) {
  string raw = RequireStr(kv, key);
  vector<string> toks = SplitCSV(raw);
  if (toks.empty()) {
    cerr << "ERROR: config key '" << key << "' has no values" << endl;
    exit(1);
  }
  vector<double> out;
  out.reserve(toks.size());
  for (size_t i = 0; i < toks.size(); ++i) {
    istringstream iss(toks[i]);
    double v; iss >> v;
    if (iss.fail()) {
      cerr << "ERROR: config key '" << key << "' has non-numeric value: "
           << toks[i] << endl;
      exit(1);
    }
    out.push_back(v);
  }
  return out;
}

static vector<int> RequireIntList(const map<string,string> &kv,
                                  const string &key) {
  vector<double> d = RequireDoubleList(kv, key);
  vector<int> out; out.reserve(d.size());
  for (size_t i = 0; i < d.size(); ++i) out.push_back(static_cast<int>(d[i]));
  return out;
}

Config LoadConfig(const string &path) {
  map<string,string> kv = ReadKeyValueFile(path);
  Config c;
  c.nElectrons          = RequireLong  (kv, "nElectrons");
  c.nProtons            = RequireLong  (kv, "nProtons");
  c.halfLife            = RequireDouble(kv, "halfLife");
  c.aperatureRadius     = RequireDouble(kv, "aperatureRadius");
  c.numHalfLives        = RequireIntList   (kv, "numHalfLives");
  c.x_0                 = RequireDoubleList(kv, "x_0");
  c.y_0                 = RequireDoubleList(kv, "y_0");
  c.R                   = RequireDoubleList(kv, "R");
  c.De                  = RequireDoubleList(kv, "De");
  const size_t N = c.numHalfLives.size();
  if (c.x_0.size() != N || c.y_0.size() != N || c.R.size() != N) {
    cerr << "ERROR: numHalfLives, x_0, y_0, R must have the same length "
         << "(got " << c.numHalfLives.size() << ", " << c.x_0.size()
         << ", " << c.y_0.size() << ", " << c.R.size() << ")" << endl;
    exit(1);
  }
  c.TKelvin             = RequireDouble(kv, "TKelvin");
  c.PTorr               = RequireDouble(kv, "PTorr");
  c.fillGasAtomicRadius = RequireDouble(kv, "fillGasAtomicRadius");
  c.beamIonAtomicRadius = RequireDouble(kv, "beamIonAtomicRadius");
  c.M_fillGas           = RequireDouble(kv, "M_fillGas");
  c.M_beamIon           = RequireDouble(kv, "M_beamIon");
  c.inputFile     = OptionalStr(kv, "inputFile",     "tripleCoincidences.root");
  c.inputHistName = OptionalStr(kv, "inputHistName", "roughBeamDistribution");
  c.outputFile    = OptionalStr(kv, "outputFile",    "AllHistograms.root");
  c.nThreads      = OptionalInt(kv, "nThreads",      0);
  return c;
}

//============================================================================
// Physics
//============================================================================

double CalculateDiffusivity(const Config &c) {
  const double A    = 1.859e-3;
  const double Patm = c.PTorr / 760.;
  const double collisionDiameter =
      0.5 * (2 * c.fillGasAtomicRadius + 2 * c.beamIonAtomicRadius) * 0.01;

  const double Da = A * pow(c.TKelvin, 1.5)
                  * TMath::Sqrt(1. / c.M_fillGas + 1. / c.M_beamIon)
                  / Patm
                  / collisionDiameter
                  / collisionDiameter;
  return Da;
}

//============================================================================
// Jobs and sampling primitives
//============================================================================
//
// The unit of parallel work is "fill one histogram".  Each job is a small
// POD that tells a worker exactly which histogram to build.  Workers own
// their own TF2 and TRandom3 so nothing touches gRandom or the formula
// cache of a TF2 that would otherwise be shared between threads.
//
static std::mutex g_ioMutex;  // serializes progress prints

struct BeamSpotJob {
  string name;
  string title;
  int    nhl;
  double x0;
  double y0;
  double Rval;
  double rMax;       // aperture + diffusion boundary, precomputed
};

struct ElectronCloudJob {
  string name;
  string title;
  double sigma_e;    // Gaussian width [mm]
  double xMin, xMax, yMin, yMax;
};

// Build a fresh TF2 Gaussian sampler.  Each thread needs its own instance
// -- TF2 caches a 2-D inverse-CDF on first GetRandom2 and mutates it on
// subsequent calls, so sharing between threads corrupts the cache.
//
// IMPORTANT: Constructing a TF2 from a string expression causes ROOT to
// invoke the Cling JIT to compile the formula at runtime.  Cling has
// process-global state and is NOT thread-safe, so having 12 workers each
// parse the formula concurrently produces garbage like "redefinition of
// ROOT_TFormula_triggerAutoParse" and eventually a segfault.
//
// The fix is to compile the formula exactly once on the main thread into
// `g_f2DTemplate` below, then have each worker copy-construct from that
// template.  Copy-construction doesn't touch Cling.
static TF2 *g_f2DTemplate = nullptr;   // owned by main thread; read-only in workers

static void BuildGaussianTF2Template() {
  if (g_f2DTemplate) return;
  g_f2DTemplate = new TF2("2Dgaus_tmpl",
                          "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[2])",
                          -50, 50, -50, 50);
  g_f2DTemplate->SetNpx(1000);
  g_f2DTemplate->SetNpy(1000);
  g_f2DTemplate->SetParameter(0, 1);
  // Force the Cling JIT to happen now, on the main thread, by evaluating
  // the formula once.  After this, TF2 copies don't trigger Cling.
  g_f2DTemplate->Eval(0.0, 0.0);
}

static TF2 *MakeGaussianTF2() {
  // Copy-construct from the pre-compiled template.  This gives each thread
  // its own inverse-CDF cache and its own parameter slots, without going
  // anywhere near the Cling interpreter.
  return new TF2(*g_f2DTemplate);
}

// Fill a beam-spot histogram.  The histogram is created detached from any
// TFile (SetDirectory(nullptr)) and ownership is returned to the caller so
// the main thread can write it out after all workers have joined.
static TH2D *FillBeamSpot(const BeamSpotJob &job,
                          long nProtons,
                          TF2 *f2D,
                          TRandom3 &gen) {
  const double sigma = job.Rval / 2.;
  f2D->SetParameter(1, job.x0);
  f2D->SetParameter(2, sigma);
  f2D->SetParameter(3, job.y0);

  TH2D *h = new TH2D(job.name.c_str(), job.title.c_str(),
                     1000, -50, 50, 1000, -50, 50);
  h->SetDirectory(nullptr);
  h->GetXaxis()->SetTitle("x-position [mm]");
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->SetTitle("y-position [mm]");
  h->GetYaxis()->CenterTitle();

  double x, y, r;
  for (long i = 0; i < nProtons; ++i) {
    f2D->GetRandom2(x, y, &gen);             // thread-local RNG
    r = TMath::Sqrt(x * x + y * y);
    if (r < job.rMax) h->Fill(x, y);
  }
  return h;
}

static TH2D *FillElectronCloud(const ElectronCloudJob &job,
                               long nElectrons,
                               TF2 *f2D,
                               TRandom3 &gen) {
  f2D->SetParameter(1, 0);
  f2D->SetParameter(2, job.sigma_e);
  f2D->SetParameter(3, 0);

  TH2D *h = new TH2D(job.name.c_str(), job.title.c_str(),
                     1000, job.xMin, job.xMax, 1000, job.yMin, job.yMax);
  h->SetDirectory(nullptr);

  double x, y;
  for (long i = 0; i < nElectrons; ++i) {
    f2D->GetRandom2(x, y, &gen);
    h->Fill(x, y);
  }
  return h;
}

// Build job lists for the whole config.
static vector<BeamSpotJob> MakeBeamSpotJobs(const Config &c, double Da) {
  vector<BeamSpotJob> jobs;
  jobs.reserve(c.numHalfLives.size());
  for (size_t k = 0; k < c.numHalfLives.size(); ++k) {
    BeamSpotJob j;
    j.nhl  = c.numHalfLives[k];
    j.x0   = c.x_0[k];
    j.y0   = c.y_0[k];
    j.Rval = c.R[k];
    const double rbd = TMath::Sqrt(4 * Da * j.nhl * c.halfLife) * 10.;
    j.rMax = c.aperatureRadius + rbd;
    const double t_s = j.nhl * c.halfLife;

    char name[64], title[128];
    snprintf(name, sizeof(name), "BeamSpot_%zu", k);
    snprintf(title, sizeof(title),
             "t=%.3f s, x_{0}=%.2f, y_{0}=%.2f, R=%.2f",
             t_s, j.x0, j.y0, j.Rval);
    j.name  = name;
    j.title = title;
    jobs.push_back(j);
  }
  return jobs;
}

static vector<ElectronCloudJob> MakeElectronCloudJobs(const Config &c) {
  vector<ElectronCloudJob> jobs;
  jobs.reserve(c.De.size() * 9);
  const double xMin = -20, xMax = 20;
  const double yMin = -20, yMax = 20;

  for (size_t k = 0; k < c.De.size(); ++k) {
    for (int driftUs = 0; driftUs <= 8; ++driftUs) {
      const double t_us = (driftUs == 0) ? 0.1 : driftUs;
      const double t    = t_us * 1e-6;
      const double sigma_e = TMath::Sqrt(4 * c.De[k] * t) * 10.;  // cm->mm

      ElectronCloudJob j;
      j.sigma_e = sigma_e;
      j.xMin = xMin; j.xMax = xMax;
      j.yMin = yMin; j.yMax = yMax;

      char name[64], title[128];
      snprintf(name, sizeof(name), "EC_De%zu_t%dus", k, driftUs);
      if (driftUs == 0) {
        snprintf(title, sizeof(title),
                 "D_{e}=%.0f, drift time = 0.1 #mus", c.De[k]);
      } else {
        snprintf(title, sizeof(title),
                 "D_{e}=%.0f, drift time = %d #mus", c.De[k], driftUs);
      }
      j.name  = name;
      j.title = title;
      jobs.push_back(j);
    }
  }
  return jobs;
}

//============================================================================
// Parallel dispatch
//============================================================================
//
// All jobs (beam-spot and electron-cloud) are merged into a single indexed
// list so one thread pool handles both.  Each worker atomically pulls the
// next index from `nextJob`, fills its histogram, and stores the pointer in
// a pre-sized results[] slot -- no ordering contention.
//
enum class JobKind { BeamSpot, ElectronCloud };

struct JobHandle {
  JobKind kind;
  size_t  localIndex;   // index into the matching (bsJobs|ecJobs) vector
};

static void WorkerLoop(
    unsigned tid,
    unsigned long seedBase,
    long nProtons,
    long nElectrons,
    const vector<BeamSpotJob>      *bsJobs,
    const vector<ElectronCloudJob> *ecJobs,
    const vector<JobHandle>        *handles,
    vector<std::unique_ptr<TH2D>>  *results,
    std::atomic<size_t>            *nextJob) {

  std::unique_ptr<TF2>      f2D(MakeGaussianTF2());
  std::unique_ptr<TRandom3> gen(new TRandom3(seedBase + 7919ul * (tid + 1)));

  while (true) {
    const size_t idx = nextJob->fetch_add(1);
    if (idx >= handles->size()) break;
    const JobHandle &hd = (*handles)[idx];

    // Deterministic per-job seed so output is reproducible for a fixed
    // seedBase regardless of which thread picked up which job.
    gen->SetSeed(seedBase + 104729ul * (idx + 1));

    TH2D *h = nullptr;
    if (hd.kind == JobKind::BeamSpot) {
      const BeamSpotJob &job = (*bsJobs)[hd.localIndex];
      {
        std::lock_guard<std::mutex> lk(g_ioMutex);
        cout << "[t" << tid << "] Beam spot " << hd.localIndex
             << ": nhl=" << job.nhl
             << " x0=" << job.x0 << " y0=" << job.y0
             << " R=" << job.Rval << endl;
      }
      h = FillBeamSpot(job, nProtons, f2D.get(), *gen);
    } else {
      const ElectronCloudJob &job = (*ecJobs)[hd.localIndex];
      {
        std::lock_guard<std::mutex> lk(g_ioMutex);
        cout << "[t" << tid << "] Electron cloud " << hd.localIndex
             << ": " << job.name << "  sigma=" << job.sigma_e << " mm"
             << endl;
      }
      h = FillElectronCloud(job, nElectrons, f2D.get(), *gen);
    }
    (*results)[idx].reset(h);
  }
}

//============================================================================
// Entry point
//============================================================================
void GetParticleDistributions(const char *configPath = "GetParticleDistributions.cfg") {
  const Config c = LoadConfig(configPath);

  // Read the input drift-time reference (serial; cheap).
  TFile *InFile = TFile::Open(c.inputFile.c_str(), "READ");
  if (!InFile || InFile->IsZombie()) {
    cerr << "ERROR: could not open " << c.inputFile << endl;
    exit(1);
  }
  TH1D *h = (TH1D*)InFile->Get(c.inputHistName.c_str());
  if (!h) {
    cerr << "ERROR: histogram '" << c.inputHistName
         << "' not found in " << c.inputFile << endl;
    exit(1);
  }
  const int    nBins    = h->GetNbinsX();
  const double binWidth = h->GetBinWidth(1);

  vector<double> cpb(nBins);
  for (int i = 1; i <= nBins; ++i) cpb[i-1] = h->GetBinContent(i);

  // Open output BEFORE creating histograms destined for it.
  TFile *OutFile = new TFile(c.outputFile.c_str(), "RECREATE");

  // Drift-times reconstruction -- serial, O(total counts), negligible.
  TH1D *DriftTimes = new TH1D("DriftTimes",
                              "Longitudinal beam distribution",
                              80, 0, 8);
  DriftTimes->GetXaxis()->SetTitle("Drift times [#mus]");
  DriftTimes->GetXaxis()->CenterTitle();
  DriftTimes->GetYaxis()->SetTitle("Counts / 0.1 #mus");
  DriftTimes->GetYaxis()->CenterTitle();

  for (int i = 0; i < nBins; ++i) {
    const double t = (i - 1) * binWidth;
    const long   n = static_cast<long>(cpb[i]);
    for (long j = 0; j < n; ++j) DriftTimes->Fill(t);
  }
  DriftTimes->Write();

  // ---- Build the parallel workload -----------------------------------------
  const double Da = CalculateDiffusivity(c);
  const vector<BeamSpotJob>      bsJobs = MakeBeamSpotJobs(c, Da);
  const vector<ElectronCloudJob> ecJobs = MakeElectronCloudJobs(c);

  vector<JobHandle> handles;
  handles.reserve(bsJobs.size() + ecJobs.size());
  for (size_t i = 0; i < bsJobs.size(); ++i)
    handles.push_back({ JobKind::BeamSpot, i });
  for (size_t i = 0; i < ecJobs.size(); ++i)
    handles.push_back({ JobKind::ElectronCloud, i });

  vector<std::unique_ptr<TH2D>> results(handles.size());
  std::atomic<size_t> nextJob{0};

  unsigned nThreads = (c.nThreads > 0)
                    ? static_cast<unsigned>(c.nThreads)
                    : std::thread::hardware_concurrency();
  if (nThreads == 0) nThreads = 1;
  if (nThreads > handles.size())
    nThreads = static_cast<unsigned>(handles.size());

  const unsigned long seedBase =
      static_cast<unsigned long>(std::time(nullptr));

  cout << "Dispatching " << handles.size()
       << " histogram jobs (" << bsJobs.size() << " beam spots + "
       << ecJobs.size() << " electron clouds) across "
       << nThreads << " thread(s)." << endl;

  // Compile the TF2 formula ONCE on the main thread.  Workers copy-construct
  // from this template, which avoids concurrent calls into the Cling JIT.
  BuildGaussianTF2Template();

  if (nThreads <= 1) {
    WorkerLoop(0, seedBase, c.nProtons, c.nElectrons,
               &bsJobs, &ecJobs, &handles, &results, &nextJob);
  } else {
    vector<std::thread> pool;
    pool.reserve(nThreads);
    for (unsigned t = 0; t < nThreads; ++t) {
      pool.emplace_back(WorkerLoop,
                        t, seedBase, c.nProtons, c.nElectrons,
                        &bsJobs, &ecJobs, &handles, &results, &nextJob);
    }
    for (auto &th : pool) th.join();
  }

  // ---- Write results serially (TFile is not thread-safe) -------------------
  OutFile->cd();
  for (auto &uh : results) {
    if (uh) {
      uh->SetDirectory(OutFile);
      uh->Write();
    }
  }

  OutFile->Close();
  InFile->Close();
  delete OutFile;
  delete InFile;
}

// Allow standalone compilation: `g++ GetParticleDistributions.cpp ... -o GetParticleDistributions`
// When loaded as a ROOT macro (via CLING), __CLING__ is defined and main() is
// skipped, so `root -l 'GetParticleDistributions.cpp("cfg.cfg")'` still works as before.
#ifndef __CLING__
int main(int argc, char** argv) {
  const char* configFile = (argc > 1) ? argv[1] : "GetParticleDistributions.cfg";
  GetParticleDistributions(configFile);
  return 0;
}
#endif
