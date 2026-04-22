// SGMC.cpp -------------------------------------------------------------
//
// Monte-Carlo simulation of the proton-detector efficiency.  Reads its
// runtime configuration (input file names, histogram names, list of proton
// lab-frame energies, etc.) from a plain-text config file so that no
// recompilation is needed to change inputs.
//
// Run as a ROOT macro:
//     root -l -b -q 'SGMC.cpp("SGMC.config")'
// or compile:
//     g++ -O2 `root-config --cflags --libs` SGMC.cpp -o SGMC
//     ./SGMC SGMC.config
// -----------------------------------------------------------------------------

#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TTree.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::vector;

// -----------------------------------------------------------------------------
//  Compile-time physics / geometry constants
// -----------------------------------------------------------------------------
namespace constants {
    constexpr double pi            = 3.14159265358979323846;
    constexpr double outerRadius   = 40.0;    // active / veto boundary [mm]
    constexpr double kaptonRadius  = 25.4;    // beam-entrance aperture [mm]
    constexpr double W             = 2.6e-5;  // energy per e-ion pair in P10 [MeV]
    constexpr double zMax          = 400.0;   // active-volume length [mm]
    constexpr double tMax          = 8e-6;    // max drift time from cathode [s]

    // Bragg-curve integration grid.
    constexpr double xMin          = 0.0;     // min track length [mm]
    constexpr double xMax          = 200.0;   // max track length [mm]
    constexpr double yMin          = 0.0;     // min dE/dx [MeV/mm]
    constexpr double yMax          = 0.1;     // max dE/dx [MeV/mm]
    constexpr double eps           = 0.01;    // Bragg-curve step [mm]
}

// -----------------------------------------------------------------------------
//  Runtime configuration, loaded from the config file
// -----------------------------------------------------------------------------
//
// The upstream GetAllHistograms job now writes one electron-cloud family per
// De value: histograms are named "{electronCloudPrefix}{k}_t{j}us" for
//   k in [0, nDeValues)
//   j in [0, nDriftTimes)             (j=0 corresponds to 0.1 us upstream)
//
// The scan over De is handled by looping over k and, for each k, running the
// full simulation using that De's nDriftTimes clouds.  The Gaussian "smear"
// that used to apply an extra per-step diffusion on top of the EC-sampled
// position has been removed -- diffusion is now baked into the cloud choice.
//
struct Config {
    string         histFile;
    vector<string> stoppingFiles        = {};       // one or more [MeV/mm] tables
    string         driftTimesHist       = "DriftTimes";
    vector<string> beamSpotHists        = { "rBeamSpot" };

    // Electron-cloud lookup: name = "{prefix}{k}_t{j}us".
    string         electronCloudPrefix  = "EC_De";
    int            nDeValues            = 1;   // number of De-indexed families
    int            nDriftTimes          = 9;   // per-family clouds (j = 0..N-1)

    int            nProtons             = 10;
    vector<int>    thresholds           = { 385 };
    double         stepSize             = 0.01;
    vector<double> labEnergies;

    // Threading.  0 means "auto" (use std::thread::hardware_concurrency()).
    // Energies are distributed across threads; each thread clones its own
    // ROOT histograms and owns its own TRandom3 so no global state is shared.
    int            nThreads             = 0;

    string         summaryFile          = "efficiencySummary.txt";
};

// -----------------------------------------------------------------------------
//  Small string helpers for config parsing
// -----------------------------------------------------------------------------
static string trim(const string& s) {
    const char* ws = " \t\r\n";
    const auto a = s.find_first_not_of(ws);
    if (a == string::npos) return "";
    const auto b = s.find_last_not_of(ws);
    return s.substr(a, b - a + 1);
}

// Parse one "key = value [value ...]" line.  Returns false on comment/blank.
static bool parseConfigLine(const string& rawLine, string& key, string& value) {
    // Strip comments.
    string line = rawLine;
    const auto hash = line.find('#');
    if (hash != string::npos) line.erase(hash);
    line = trim(line);
    if (line.empty()) return false;

    const auto eq = line.find('=');
    if (eq == string::npos) {
        throw std::runtime_error("Malformed config line (no '='): " + rawLine);
    }
    key   = trim(line.substr(0, eq));
    value = trim(line.substr(eq + 1));
    if (key.empty()) {
        throw std::runtime_error("Malformed config line (empty key): " + rawLine);
    }
    return true;
}

static vector<double> parseDoubleList(const string& value, const string& key) {
    vector<double> out;
    std::istringstream iss(value);
    double x;
    while (iss >> x) out.push_back(x);
    if (out.empty()) {
        throw std::runtime_error("Expected one or more numbers for key '" + key + "'");
    }
    return out;
}

static vector<int> parseIntList(const string& value, const string& key) {
    vector<int> out;
    std::istringstream iss(value);
    int x;
    while (iss >> x) out.push_back(x);
    if (out.empty()) {
        throw std::runtime_error("Expected one or more integers for key '" + key + "'");
    }
    return out;
}

static vector<string> parseStringList(const string& value, const string& key) {
    vector<string> out;
    std::istringstream iss(value);
    string x;
    while (iss >> x) out.push_back(x);
    if (out.empty()) {
        throw std::runtime_error("Expected one or more tokens for key '" + key + "'");
    }
    return out;
}

static Config loadConfig(const string& path) {
    ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Could not open config file: " + path);
    }

    Config cfg;
    bool haveHistFile     = false;
    bool haveStoppingFile = false;
    bool haveEnergies     = false;

    string rawLine;
    while (std::getline(in, rawLine)) {
        string key, value;
        if (!parseConfigLine(rawLine, key, value)) continue;

        if      (key == "histFile")            { cfg.histFile = value;            haveHistFile = true; }
        else if (key == "stoppingFile" ||
                 key == "stoppingFiles")       { cfg.stoppingFiles = parseStringList(value, key); haveStoppingFile = true; }
        else if (key == "driftTimesHist")      { cfg.driftTimesHist = value; }
        else if (key == "beamSpotHist" ||
                 key == "beamSpotHists")       { cfg.beamSpotHists = parseStringList(value, key); }
        else if (key == "electronCloudPrefix") { cfg.electronCloudPrefix = value; }
        else if (key == "nDeValues")           { cfg.nDeValues = std::stoi(value); }
        else if (key == "nDriftTimes")         { cfg.nDriftTimes = std::stoi(value); }
        else if (key == "nProtons")            { cfg.nProtons = std::stoi(value); }
        else if (key == "threshold" ||
                 key == "thresholds")          { cfg.thresholds = parseIntList(value, key); }
        else if (key == "stepSize")            { cfg.stepSize = std::stod(value); }
        else if (key == "labEnergies")         { cfg.labEnergies = parseDoubleList(value, key); haveEnergies = true; }
        else if (key == "nThreads")            { cfg.nThreads = std::stoi(value); }
        else if (key == "summaryFile")         { cfg.summaryFile = value; }
        else {
            cerr << "Warning: unknown config key '" << key << "' -- ignoring." << endl;
        }
    }

    if (!haveHistFile)     throw std::runtime_error("Config missing required key: histFile");
    if (!haveStoppingFile || cfg.stoppingFiles.empty())
        throw std::runtime_error("Config missing required key: stoppingFile(s)");
    if (!haveEnergies)     throw std::runtime_error("Config missing required key: labEnergies");
    if (cfg.nDeValues <= 0) {
        throw std::runtime_error("nDeValues must be positive");
    }
    if (cfg.nDriftTimes <= 0) {
        throw std::runtime_error("nDriftTimes must be positive");
    }
    return cfg;
}

static void printConfig(const Config& cfg) {
    cout << "---- SGMC configuration ----"                                     << endl
         << "  histFile            = " << cfg.histFile                            << endl
         << "  stoppingFiles       =";
    for (const auto& s : cfg.stoppingFiles) cout << " " << s;
    cout << endl
         << "  driftTimesHist      = " << cfg.driftTimesHist                      << endl
         << "  beamSpotHists       =";
    for (const auto& h : cfg.beamSpotHists) cout << " " << h;
    cout << endl
         << "  electronCloudPrefix = " << cfg.electronCloudPrefix                 << endl
         << "  nDeValues           = " << cfg.nDeValues                           << endl
         << "  nDriftTimes         = " << cfg.nDriftTimes                         << endl
         << "  nProtons            = " << cfg.nProtons                            << endl
         << "  thresholds          =";
    for (int t : cfg.thresholds) cout << " " << t;
    cout << endl
         << "  stepSize            = " << cfg.stepSize                            << endl
         << "  nThreads            = " << cfg.nThreads
         << (cfg.nThreads == 0 ? "  (auto)" : "")                                 << endl
         << "  summaryFile         = " << cfg.summaryFile                         << endl
         << "  labEnergies [MeV]   =";
    for (double e : cfg.labEnergies) cout << " " << e;
    cout << endl << "-------------------------------" << endl;
}

// -----------------------------------------------------------------------------
//  Bragg-curve tools
// -----------------------------------------------------------------------------
static TGraph* buildDepositGraph(const string& stoppingFile, double Ep) {
    TGraph* stoppingGraph = new TGraph(stoppingFile.c_str());
    if (stoppingGraph->GetN() == 0) {
        delete stoppingGraph;
        throw std::runtime_error("Stopping-power file is empty or unreadable: " + stoppingFile);
    }
    TGraph* depositGraph = new TGraph();

    double z = 0.0;
    double E = Ep;
    int    n = 0;
    while (z < constants::xMax) {
        if (E < 0.0) {
            depositGraph->SetPoint(n++, z, 0.0);
            z += constants::eps;
            continue;
        }
        const double dEdz = stoppingGraph->Eval(E);
        depositGraph->SetPoint(n++, z, dEdz);
        E -= dEdz * constants::eps;
        z += constants::eps;
    }

    depositGraph->GetXaxis()->SetRangeUser(constants::xMin, constants::xMax);
    depositGraph->GetYaxis()->SetRangeUser(constants::yMin, constants::yMax);

    delete stoppingGraph;
    return depositGraph;
}

// Trapezoidal integral of a TGraph between xmin and xmax.
static double integrate(TGraph* g, double xmin, double xmax) {
    const int N = g->GetN();
    double x1, x2, y1, y2;
    double sum = 0.0;
    for (int n = 1; n < N; ++n) {
        g->GetPoint(n - 1, x1, y1);
        g->GetPoint(n,     x2, y2);
        if (x1 < xmin) continue;
        if (x2 > xmax) break;
        sum += 0.5 * (y1 + y2) * (x2 - x1);
    }
    return sum;
}

// -----------------------------------------------------------------------------
//  Geometry / physics helpers
// -----------------------------------------------------------------------------
static bool isInActiveRegion(double x, double y) {
    return TMath::Sqrt(x * x + y * y) < constants::outerRadius;
}

// Count how many of nElectrons sampled from the electron-cloud histogram
// land outside of the active region when offset by (x, y).  Diffusion is now
// entirely encoded in the per-De, per-drift-time choice of `electronCloud`,
// so no additional Gaussian smear is applied here.
//
// The thread's TRandom3 is forwarded to TH2::GetRandom2 so that none of the
// sampling touches the global gRandom -- essential for thread safety.
static int electronDiffusionEffect(int nElectrons, TH2D* electronCloud,
                                   double x, double y, TRandom3& gen) {
    double dx, dy;
    int    veto = 0;
    for (int i = 0; i < nElectrons; ++i) {
        electronCloud->GetRandom2(dx, dy, &gen);
        if (!isInActiveRegion(x + dx, y + dy)) ++veto;
    }
    return veto;
}

// -----------------------------------------------------------------------------
//  Output
// -----------------------------------------------------------------------------

// One variant's result for a single energy.  A variant is a full tuple of
// scan parameters: beam-spot histogram, veto threshold, stopping-power file,
// and electron-cloud De-family index.
struct VariantResult {
    string beamSpotHist;
    int    threshold;
    string stoppingFile;
    int    deIndex;               // which De-indexed electron-cloud family
    int    h6;
    double efficiency;
    double efficiencyError;       // 1-sigma statistical, from sqrt(h6)
};

// All results for a single proton lab energy.
struct EnergyResult {
    double                 protonEnergy;      // [MeV]
    double                 sumEnergy;         // summed deposited energy
    double                 totalElectrons;    // expected total, from E/W
    long long              sumElectrons;      // summed ionization electrons
    vector<VariantResult>  variants;          // every variant in the scan
    VariantResult          central;           // median by efficiency
    VariantResult          lower;             // min efficiency
    VariantResult          upper;             // max efficiency
};

// Efficiency and its 1-sigma statistical uncertainty from a Poisson-distributed
// count h6 out of nProtons trials.
static void efficiencyFromCounts(int h6, int nProtons,
                                 double& eff, double& dEff) {
    eff  = static_cast<double>(h6) / nProtons;
    dEff = (h6 > 0) ? eff / TMath::Sqrt(h6) : 0.0;
}

// Pick min/max/median-by-efficiency variants out of a non-empty list.
static void summarizeVariants(vector<VariantResult> variants,
                              VariantResult& lower,
                              VariantResult& upper,
                              VariantResult& central) {
    std::sort(variants.begin(), variants.end(),
              [](const VariantResult& a, const VariantResult& b) {
                  return a.efficiency < b.efficiency;
              });
    lower   = variants.front();
    upper   = variants.back();
    central = variants[variants.size() / 2];   // upper median for even counts
}

// Short, table-friendly description of a variant.
static string variantTag(const VariantResult& v) {
    char buf[256];
    std::snprintf(buf, sizeof(buf),
                  "beamSpot=%s, threshold=%d, stoppingFile=%s, DeIdx=%d",
                  v.beamSpotHist.c_str(), v.threshold,
                  v.stoppingFile.c_str(), v.deIndex);
    return string(buf);
}

static void writeSummary(const string& fileName,
                         int nProtons,
                         const vector<EnergyResult>& results) {
    FILE* out = std::fopen(fileName.c_str(), "w");
    if (!out) {
        cerr << "Error: could not open summary file " << fileName << endl;
        return;
    }

    std::fprintf(out, "# SGMC efficiency summary\n");
    std::fprintf(out, "# nProtons per variant: %d\n", nProtons);
    std::fprintf(out, "# Scan dimensions: beamSpotHist x threshold x stoppingFile x DeIndex\n");
    std::fprintf(out, "# central value is the median across variants\n");
    std::fprintf(out, "# lower / upper are the min / max efficiency across variants\n");
    std::fprintf(out, "# errors are 1-sigma statistical, = efficiency / sqrt(h6)\n\n");

    for (const auto& r : results) {
        const int keV = static_cast<int>(std::round(r.protonEnergy * 1000.0));
        std::fprintf(out, "============================================================\n");
        std::fprintf(out, "Proton energy: %.4f MeV (%d keV)\n", r.protonEnergy, keV);
        std::fprintf(out, "  Total energy deposited : %.4f MeV\n", r.sumEnergy);
        std::fprintf(out, "  Expected electrons     : %.4f  (from E/W)\n", r.totalElectrons);
        std::fprintf(out, "  Summed electrons       : %lld\n",  r.sumElectrons);
        std::fprintf(out, "  Variants evaluated     : %zu\n",   r.variants.size());

        std::fprintf(out,
            "  Central (median) : eff = %.4f +/- %.4f   [%s, h6=%d]\n",
            r.central.efficiency, r.central.efficiencyError,
            variantTag(r.central).c_str(), r.central.h6);
        std::fprintf(out,
            "  Lower  (min)     : eff = %.4f +/- %.4f   [%s, h6=%d]\n",
            r.lower.efficiency, r.lower.efficiencyError,
            variantTag(r.lower).c_str(), r.lower.h6);
        std::fprintf(out,
            "  Upper  (max)     : eff = %.4f +/- %.4f   [%s, h6=%d]\n",
            r.upper.efficiency, r.upper.efficiencyError,
            variantTag(r.upper).c_str(), r.upper.h6);

        std::fprintf(out, "\n  All variants:\n");
        std::fprintf(out, "    %-24s %10s %-40s %8s %8s %10s %10s\n",
                     "beamSpotHist", "threshold", "stoppingFile", "DeIdx",
                     "h6", "eff", "dEff");
        for (const auto& v : r.variants) {
            std::fprintf(out, "    %-24s %10d %-40s %8d %8d %10.4f %10.4f\n",
                         v.beamSpotHist.c_str(), v.threshold,
                         v.stoppingFile.c_str(), v.deIndex,
                         v.h6, v.efficiency, v.efficiencyError);
        }
        std::fprintf(out, "\n");
    }

    std::fclose(out);
}

// -----------------------------------------------------------------------------
//  Histogram loading
// -----------------------------------------------------------------------------
template <typename T>
static T* getHist(TFile* f, const string& name) {
    T* h = dynamic_cast<T*>(f->Get(name.c_str()));
    if (!h) {
        throw std::runtime_error("Histogram '" + name + "' not found in " +
                                 f->GetName());
    }
    return h;
}

// Load every electron-cloud family.  Returns a 2-D vector indexed as
//   clouds[deIndex][driftIndex]
// where cloud names follow "{prefix}{deIndex}_t{driftIndex}us".
static vector<vector<TH2D*>> loadElectronClouds(TFile* f,
                                                const string& prefix,
                                                int nDe,
                                                int nT) {
    vector<vector<TH2D*>> clouds(nDe);
    for (int k = 0; k < nDe; ++k) {
        clouds[k].reserve(nT);
        for (int j = 0; j < nT; ++j) {
            const string name = prefix + std::to_string(k)
                              + "_t" + std::to_string(j) + "us";
            clouds[k].push_back(getHist<TH2D>(f, name));
        }
    }
    return clouds;
}

// -----------------------------------------------------------------------------
//  Per-proton simulation record
// -----------------------------------------------------------------------------
struct ProtonResult {
    int    eVeto;       // electrons vetoed over this track
    bool   hitCathode;  // true if the proton reached z = zMax
};

// Simulate nProtons protons for a given (beamSpot, energy, stoppingFile, De
// family) and return each proton's eVeto count plus cathode-hit flag, along
// with diagnostic totals.  The result is independent of the veto threshold
// (thresholds are a pure post-processing cut on eVeto).
static vector<ProtonResult> simulateProtons(
        int nProtons, double stepSize,
        TH2D* beamSpot, TH1D* driftTimes,
        const vector<TH2D*>& electronClouds,
        TGraph* depositGraph,
        TRandom3& gen,
        double& sumEnergy, long long& sumElectrons) {

    vector<ProtonResult> protons;
    protons.reserve(nProtons);
    sumEnergy    = 0.0;
    sumElectrons = 0;
    const int nClouds = static_cast<int>(electronClouds.size());

    for (int np = 0; np < nProtons; ++np) {
        double x, y;
        beamSpot->GetRandom2(x, y, &gen);
        double t = driftTimes->GetRandom(&gen);      // [us]
        int    tint = static_cast<int>(std::round(t));
        t *= 1e-6;                                   // [s]

        if (tint < 0 || tint >= nClouds) {
            cerr << "Warning: drift-time index " << tint
                 << " out of range [0, " << nClouds
                 << ") -- skipping proton." << endl;
            continue;
        }

        const double theta = TMath::ACos(1.0 - 2.0 * gen.Uniform());
        const double phi   = 2.0 * constants::pi * gen.Uniform();
        const double dx = stepSize * TMath::Sin(theta) * TMath::Cos(constants::pi + phi);
        const double dy = stepSize * TMath::Sin(theta) * TMath::Sin(constants::pi + phi);
        const double dz = stepSize * TMath::Cos(theta);

        double z = t * (constants::zMax / constants::tMax);

        int    eVeto     = 0;
        double xp        = 0.0;
        double yp        = 0.0;
        double dr        = 0.0;
        double energyTrk = 0.0;
        int    nElec     = 1;
        bool   hitCathode = false;

        while (nElec > 0) {
            z += dz;
            if (z > constants::zMax) { hitCathode = true; break; }

            xp += dx;
            yp += dy;
            t  += dz * (constants::tMax / constants::zMax);

            const double dE = integrate(depositGraph, dr, dr + stepSize);
            energyTrk += dE;
            nElec = static_cast<int>(std::round(dE / constants::W));

            eVeto += electronDiffusionEffect(nElec, electronClouds[tint],
                                             x + xp, y + yp, gen);

            dr           += stepSize;
            sumElectrons += nElec;
        }

        sumEnergy += energyTrk;
        protons.push_back({ eVeto, hitCathode });
    }

    return protons;
}

// -----------------------------------------------------------------------------
//  Per-thread workspace
// -----------------------------------------------------------------------------
//
// ROOT histograms and graphs are NOT thread-safe: they cache state
// (cumulative distributions in TH1/TH2, last interpolation index in TGraph)
// that is mutated by Eval/GetRandom*.  A worker thread must therefore own
// private clones of every object it samples from, plus its own TRandom3 so
// it never touches the global gRandom.
//
// A ThreadWorkspace is built once per thread from the shared "template"
// histograms loaded in main().  Clones are owned via unique_ptr so the
// workspace frees them when it goes out of scope.
struct ThreadWorkspace {
    std::unique_ptr<TRandom3>                        gen;
    std::unique_ptr<TH1D>                            driftTimes;
    vector<std::unique_ptr<TH2D>>                    beamSpots;     // [ib]
    vector<vector<std::unique_ptr<TH2D>>>            electronClouds;// [kDe][jT]

    // Raw back-pointers for convenience (the sim code takes raw TH*/TGraph*).
    vector<TH2D*>                                    beamSpotPtrs;
    vector<vector<TH2D*>>                            electronCloudPtrs;
};

// Clone all shared histograms into a new workspace and seed its RNG.
static ThreadWorkspace makeWorkspace(TH1D* driftTimesTmpl,
                                     const vector<TH2D*>& beamSpotsTmpl,
                                     const vector<vector<TH2D*>>& ecTmpl,
                                     unsigned long seed) {
    ThreadWorkspace ws;
    ws.gen.reset(new TRandom3(seed));

    ws.driftTimes.reset(static_cast<TH1D*>(driftTimesTmpl->Clone()));
    ws.driftTimes->SetDirectory(nullptr);   // detach from any TFile

    ws.beamSpots.reserve(beamSpotsTmpl.size());
    ws.beamSpotPtrs.reserve(beamSpotsTmpl.size());
    for (TH2D* h : beamSpotsTmpl) {
        TH2D* clone = static_cast<TH2D*>(h->Clone());
        clone->SetDirectory(nullptr);
        ws.beamSpotPtrs.push_back(clone);
        ws.beamSpots.emplace_back(clone);
    }

    ws.electronClouds.resize(ecTmpl.size());
    ws.electronCloudPtrs.resize(ecTmpl.size());
    for (size_t k = 0; k < ecTmpl.size(); ++k) {
        ws.electronClouds[k].reserve(ecTmpl[k].size());
        ws.electronCloudPtrs[k].reserve(ecTmpl[k].size());
        for (TH2D* h : ecTmpl[k]) {
            TH2D* clone = static_cast<TH2D*>(h->Clone());
            clone->SetDirectory(nullptr);
            ws.electronCloudPtrs[k].push_back(clone);
            ws.electronClouds[k].emplace_back(clone);
        }
    }
    return ws;
}

// -----------------------------------------------------------------------------
//  Per-energy job (runs the full scan grid for one proton energy)
// -----------------------------------------------------------------------------
static EnergyResult runEnergyJob(double protonEnergy,
                                 const Config& cfg,
                                 ThreadWorkspace& ws,
                                 std::mutex& ioMutex) {
    EnergyResult er;
    er.protonEnergy   = protonEnergy;
    er.totalElectrons = protonEnergy / constants::W;
    er.sumEnergy      = 0.0;
    er.sumElectrons   = 0;
    int nSimRuns = 0;

    for (const auto& stoppingFile : cfg.stoppingFiles) {
        // depositGraph depends on (stoppingFile, energy).  It's built here
        // on the worker thread so TGraph::Eval isn't shared between threads.
        TGraph* depositGraph = buildDepositGraph(stoppingFile, protonEnergy);

        for (size_t ib = 0; ib < ws.beamSpotPtrs.size(); ++ib) {
            for (int kDe = 0; kDe < cfg.nDeValues; ++kDe) {
                double    simEnergy    = 0.0;
                long long simElectrons = 0;
                const auto protons = simulateProtons(
                    cfg.nProtons, cfg.stepSize,
                    ws.beamSpotPtrs[ib], ws.driftTimes.get(),
                    ws.electronCloudPtrs[kDe],
                    depositGraph, *ws.gen,
                    simEnergy, simElectrons);

                er.sumEnergy    += simEnergy;
                er.sumElectrons += simElectrons;
                ++nSimRuns;

                for (int threshold : cfg.thresholds) {
                    int h6 = 0;
                    for (const auto& p : protons) {
                        if (!p.hitCathode && p.eVeto < threshold) ++h6;
                    }
                    VariantResult v;
                    v.beamSpotHist = cfg.beamSpotHists[ib];
                    v.threshold    = threshold;
                    v.stoppingFile = stoppingFile;
                    v.deIndex      = kDe;
                    v.h6           = h6;
                    efficiencyFromCounts(h6, cfg.nProtons,
                                         v.efficiency, v.efficiencyError);
                    er.variants.push_back(v);

                    // Serialize prints so lines from different threads don't
                    // interleave.
                    std::lock_guard<std::mutex> lk(ioMutex);
                    cout << "  E = "          << protonEnergy << " MeV"
                         << "  beamSpot = "   << v.beamSpotHist
                         << "  threshold = "  << threshold
                         << "  stopping = "   << stoppingFile
                         << "  DeIdx = "      << kDe
                         << "  h6 = "         << h6 << "/" << cfg.nProtons
                         << "  eff = "        << v.efficiency
                         << endl;
                }
            }
        }

        delete depositGraph;
    }

    if (nSimRuns > 0) {
        er.sumEnergy    /= nSimRuns;
        er.sumElectrons /= nSimRuns;
    }
    summarizeVariants(er.variants, er.lower, er.upper, er.central);
    return er;
}

// -----------------------------------------------------------------------------
//  Main simulation
// -----------------------------------------------------------------------------
int SGMC(const char* configFile = "SGMC.config") {
    std::srand(static_cast<unsigned>(std::time(nullptr)));

    Config cfg;
    try {
        cfg = loadConfig(configFile);
    } catch (const std::exception& ex) {
        cerr << "Failed to load config: " << ex.what() << endl;
        return 1;
    }

    printConfig(cfg);

    TFile* rootFile = TFile::Open(cfg.histFile.c_str(), "READ");
    if (!rootFile || rootFile->IsZombie()) {
        cerr << "Could not open ROOT file: " << cfg.histFile << endl;
        return 1;
    }

    TH1D*                  driftTimes;
    vector<TH2D*>          beamSpots;
    vector<vector<TH2D*>>  electronClouds;   // [deIndex][driftIndex]
    try {
        driftTimes = getHist<TH1D>(rootFile, cfg.driftTimesHist);
        beamSpots.reserve(cfg.beamSpotHists.size());
        for (const auto& name : cfg.beamSpotHists) {
            beamSpots.push_back(getHist<TH2D>(rootFile, name));
        }
        electronClouds = loadElectronClouds(rootFile,
                                            cfg.electronCloudPrefix,
                                            cfg.nDeValues,
                                            cfg.nDriftTimes);
    } catch (const std::exception& ex) {
        cerr << ex.what() << endl;
        rootFile->Close();
        delete rootFile;
        return 1;
    }

    const size_t nVariantsPerEnergy =
          cfg.beamSpotHists.size() * cfg.thresholds.size()
        * cfg.stoppingFiles.size() * static_cast<size_t>(cfg.nDeValues);
    cout << "Scan size per energy: " << nVariantsPerEnergy
         << "  (" << cfg.beamSpotHists.size()   << " beamSpots x "
                 << cfg.thresholds.size()       << " thresholds x "
                 << cfg.stoppingFiles.size()    << " stoppingFiles x "
                 << cfg.nDeValues               << " De values)"
         << endl;

    // ------------------------------------------------------------------
    // Parallelize across energies.  Each worker thread:
    //   1. Clones the shared histograms into its own ThreadWorkspace.
    //   2. Atomically claims the next unprocessed energy index.
    //   3. Runs the full scan grid for that energy.
    //   4. Writes the result into a pre-sized slot (no ordering contention).
    // ------------------------------------------------------------------
    const size_t nEnergies = cfg.labEnergies.size();
    unsigned nThreads = (cfg.nThreads > 0)
                      ? static_cast<unsigned>(cfg.nThreads)
                      : std::thread::hardware_concurrency();
    if (nThreads == 0) nThreads = 1;                        // fallback
    if (nThreads > nEnergies)
        nThreads = static_cast<unsigned>(nEnergies);        // no idle threads

    cout << "Launching " << nThreads << " worker thread(s) for "
         << nEnergies << " energies." << endl;

    vector<EnergyResult> allResults(nEnergies);
    std::atomic<size_t>  nextJob{0};
    std::mutex           ioMutex;
    // Root seed for reproducibility-ish: each worker derives its TRandom3
    // seed from this plus its energy index.
    const unsigned long seedBase =
        static_cast<unsigned long>(std::time(nullptr));

    // Build ALL per-thread workspaces on the main thread, BEFORE launching
    // workers.  TH1::Clone() walks ROOT's global object/directory lists and
    // is not thread-safe; two threads cloning concurrently race and crash
    // inside TCollection::Clone.  Doing the cloning serially up-front
    // eliminates the race and costs only a few ms per workspace.
    //
    // There is an additional first-use cost to worry about: on the very
    // first TH1::Clone / TGraph construction, ROOT lazy-loads the relevant
    // class dictionaries through Cling.  If that first call happens inside
    // a worker thread it re-enters Cling from multiple threads and
    // segfaults deep inside TClingClassInfo.  The workspace-building loop
    // below guarantees the first TH1::Clone happens here on the main
    // thread; we also build a throwaway TGraph so that the TGraph
    // dictionary is loaded before any worker calls buildDepositGraph().
    {
        TGraph warmup;            // triggers TGraph dictionary load
        (void)warmup;
    }
    vector<ThreadWorkspace> workspaces;
    workspaces.reserve(nThreads);
    for (unsigned t = 0; t < nThreads; ++t) {
        workspaces.push_back(makeWorkspace(driftTimes, beamSpots,
                                           electronClouds,
                                           seedBase + 1009ul * t));
    }

    auto worker = [&](unsigned tid) {
        ThreadWorkspace& ws = workspaces[tid];
        while (true) {
            const size_t i = nextJob.fetch_add(1);
            if (i >= nEnergies) break;
            // Seed this job deterministically from (base, energy-index).
            ws.gen->SetSeed(seedBase + 7919ul * (i + 1));
            allResults[i] = runEnergyJob(cfg.labEnergies[i], cfg,
                                         ws, ioMutex);
        }
    };

    vector<std::thread> pool;
    pool.reserve(nThreads);
    if (nThreads <= 1) {
        worker(0);                                           // run inline
    } else {
        for (unsigned t = 0; t < nThreads; ++t)
            pool.emplace_back(worker, t);
        for (auto& th : pool) th.join();
    }

    writeSummary(cfg.summaryFile, cfg.nProtons, allResults);
    cout << "Wrote summary to " << cfg.summaryFile << endl;

    rootFile->Close();
    delete rootFile;
    return 0;
}

// Allow standalone compilation: `g++ SGMC.cpp ... -o SGMC`
#ifndef __CLING__
int main(int argc, char** argv) {
    const char* configFile = (argc > 1) ? argv[1] : "SGMC.config";
    return SGMC(configFile);
}
#endif
