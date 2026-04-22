// RecoverBeamSpot.cpp -------------------------------------------------------
//
// Fit a 2D Gaussian beam spot to the measured pad intensities of a five-pad
// detector arrangement: one central disc and four azimuthally-divided outer
// quadrants.  All runtime inputs are read from a plain-text config file;
// nothing is hardcoded for a specific experiment.
//
// Run as a ROOT macro:
//   root -l -b -q 'RecoverBeamSpot.cpp("RecoverBeamSpot.cfg")'
// or compile and run:
//   g++ -O2 $(root-config --cflags) RecoverBeamSpot.cpp $(root-config --libs) -o RecoverBeamSpot
//   ./RecoverBeamSpot RecoverBeamSpot.cfg
// ---------------------------------------------------------------------------

#include "TCanvas.h"
#include "TEllipse.h"
#include "TFitter.h"
#include "TH2D.h"
#include "TLine.h"
#include "TMath.h"
#include "TStyle.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::vector;

// ===========================================================================
// Configuration
// ===========================================================================
//
// Pad geometry:
//   The detector has one central circular pad of radius rCentral, surrounded
//   by an annular ring (rCentral -> rOuterPad) divided azimuthally into four
//   quadrant pads at theta = 0, 90, 180, 270 degrees.  For display purposes
//   two additional circles are drawn at rGuardInner and rGuardOuter with 45-
//   degree radial lines between them.
//
// Config keys (all optional; defaults in parentheses reproduce the original
// AstroBox-II geometry and the original integration grid):
//
//   Geometry:
//     rCentral       (14.2)  -- inner radius of the quadrant pads [mm]
//     rOuterPad      (27.6)  -- outer radius of the quadrant pads [mm]
//     rGuardInner    (40.0)  -- inner radius of the guard ring    [mm]
//     rGuardOuter    (50.0)  -- outer radius of the guard ring    [mm]
//
//   Numerical integration:
//     epsR           (1.0)   -- radial step [mm]
//     epsThetaDeg    (1.0)   -- azimuthal step [deg]
//
//   Measured pad intensities (REQUIRED unless `intensityFile` is given):
//     padsExperimental       -- whitespace-separated list of 5 counts:
//                               {central, quadrant0, quadrant1, quadrant2,
//                                quadrant3}.  The "quadrant0" pad covers
//                               0 < theta < 90 deg, quadrant1 covers 90 <
//                               theta < 180, etc.
//     intensityFile          -- alternative: read the 5 values from this
//                               plain-text file (one or more lines,
//                               whitespace-separated).  If both are given,
//                               intensityFile wins.
//
//   Minuit starting point and bounds:
//     xStart         (0.0)   -- initial beam-center x [mm]
//     yStart         (0.0)   -- initial beam-center y [mm]
//     RStart         (20.0)  -- initial beam 1/e^2 radius [mm]
//     xStep          (20.0)  -- Minuit initial step for x
//     yStep          (20.0)  -- Minuit initial step for y
//     RStep          (20.0)  -- Minuit initial step for R
//     xLow, xHigh    (-100, 100)  -- x bounds [mm]
//     yLow, yHigh    (-100, 100)  -- y bounds [mm]
//     RLow, RHigh    (10, 100)    -- R bounds [mm]
//
//   Output / display:
//     drawBeam       (0)     -- 1 to pop a TCanvas with the fitted beam
//                               spot and pad overlay, 0 for headless run.
//     verbose        (0)     -- 1 to print every function evaluation.
//
struct Config {
    // Geometry
    double rCentral     = 14.2;
    double rOuterPad    = 27.6;
    double rGuardInner  = 40.0;
    double rGuardOuter  = 50.0;

    // Integration
    double epsR         = 1.0;
    double epsThetaDeg  = 1.0;

    // Measured intensities (central, Q0, Q1, Q2, Q3)
    vector<double> padsExperimental;

    // Minuit starting point / bounds
    double xStart = 0.0, xStep = 20.0, xLow = -100.0, xHigh = 100.0;
    double yStart = 0.0, yStep = 20.0, yLow = -100.0, yHigh = 100.0;
    double RStart = 20.0, RStep = 20.0, RLow = 10.0,  RHigh = 100.0;

    // Output
    bool drawBeam = false;
    bool verbose  = false;
};

// ---- config parsing helpers (same idiom as the other project files) ------
static string trim(const string& s) {
    const char* ws = " \t\r\n";
    const auto a = s.find_first_not_of(ws);
    if (a == string::npos) return "";
    const auto b = s.find_last_not_of(ws);
    return s.substr(a, b - a + 1);
}

static bool parseConfigLine(const string& raw, string& key, string& value) {
    string line = raw;
    const auto hash = line.find('#');
    if (hash != string::npos) line.erase(hash);
    line = trim(line);
    if (line.empty()) return false;

    const auto eq = line.find('=');
    if (eq == string::npos) {
        cerr << "Warning: malformed config line (no '='): " << raw << endl;
        return false;
    }
    key   = trim(line.substr(0, eq));
    value = trim(line.substr(eq + 1));
    return !key.empty();
}

static vector<double> parseDoubleList(const string& value) {
    vector<double> out;
    std::istringstream iss(value);
    double v;
    while (iss >> v) out.push_back(v);
    return out;
}

static bool parseBool(const string& v) {
    const string s = trim(v);
    return !(s == "0" || s == "false" || s == "no" || s.empty());
}

// Pad intensities can live either inline in the config or in a separate
// whitespace-delimited text file (one run = 5 numbers).  Returns an empty
// vector if the file is unreadable.
static vector<double> readIntensityFile(const string& path) {
    ifstream in(path);
    if (!in.is_open()) {
        cerr << "Error: could not open intensity file: " << path << endl;
        return {};
    }
    vector<double> out;
    double v;
    while (in >> v) out.push_back(v);
    return out;
}

static Config loadConfig(const string& path) {
    Config cfg;
    string intensityFile;

    ifstream in(path);
    if (!in.is_open()) {
        cerr << "Error: could not open config file: " << path << endl;
        std::exit(1);
    }

    string raw;
    while (std::getline(in, raw)) {
        string key, value;
        if (!parseConfigLine(raw, key, value)) continue;

        if      (key == "rCentral")         cfg.rCentral    = std::stod(value);
        else if (key == "rOuterPad")        cfg.rOuterPad   = std::stod(value);
        else if (key == "rGuardInner")      cfg.rGuardInner = std::stod(value);
        else if (key == "rGuardOuter")      cfg.rGuardOuter = std::stod(value);
        else if (key == "epsR")             cfg.epsR        = std::stod(value);
        else if (key == "epsThetaDeg")      cfg.epsThetaDeg = std::stod(value);
        else if (key == "padsExperimental") cfg.padsExperimental = parseDoubleList(value);
        else if (key == "intensityFile")    intensityFile = value;
        else if (key == "xStart")           cfg.xStart = std::stod(value);
        else if (key == "yStart")           cfg.yStart = std::stod(value);
        else if (key == "RStart")           cfg.RStart = std::stod(value);
        else if (key == "xStep")            cfg.xStep  = std::stod(value);
        else if (key == "yStep")            cfg.yStep  = std::stod(value);
        else if (key == "RStep")            cfg.RStep  = std::stod(value);
        else if (key == "xLow")             cfg.xLow   = std::stod(value);
        else if (key == "xHigh")            cfg.xHigh  = std::stod(value);
        else if (key == "yLow")             cfg.yLow   = std::stod(value);
        else if (key == "yHigh")            cfg.yHigh  = std::stod(value);
        else if (key == "RLow")             cfg.RLow   = std::stod(value);
        else if (key == "RHigh")            cfg.RHigh  = std::stod(value);
        else if (key == "drawBeam")         cfg.drawBeam = parseBool(value);
        else if (key == "verbose")          cfg.verbose  = parseBool(value);
        else cerr << "Warning: unknown config key '" << key
                  << "' -- ignoring." << endl;
    }

    // intensityFile (if given) overrides inline padsExperimental
    if (!intensityFile.empty()) {
        vector<double> v = readIntensityFile(intensityFile);
        if (!v.empty()) cfg.padsExperimental = v;
    }

    if (cfg.padsExperimental.size() != 5) {
        cerr << "Error: padsExperimental must have exactly 5 values (got "
             << cfg.padsExperimental.size() << "). Either set "
             << "'padsExperimental = v0 v1 v2 v3 v4' in the config or "
             << "provide 'intensityFile = path'." << endl;
        std::exit(1);
    }
    return cfg;
}

static void printConfig(const Config& c) {
    cout << "---- RecoverBeamSpot configuration ----"                         << endl;
    std::printf("  Geometry [mm]: rCentral=%.2f, rOuterPad=%.2f, "
                "rGuardInner=%.1f, rGuardOuter=%.1f\n",
                c.rCentral, c.rOuterPad, c.rGuardInner, c.rGuardOuter);
    std::printf("  Integration : epsR=%.3f mm, epsTheta=%.3f deg\n",
                c.epsR, c.epsThetaDeg);
    std::printf("  Measured counts: %.1f (central) | %.1f %.1f %.1f %.1f "
                "(Q0..Q3)\n",
                c.padsExperimental[0], c.padsExperimental[1],
                c.padsExperimental[2], c.padsExperimental[3],
                c.padsExperimental[4]);
    std::printf("  Start: x=%.2f (step %.2f, [%.1f,%.1f])\n",
                c.xStart, c.xStep, c.xLow, c.xHigh);
    std::printf("         y=%.2f (step %.2f, [%.1f,%.1f])\n",
                c.yStart, c.yStep, c.yLow, c.yHigh);
    std::printf("         R=%.2f (step %.2f, [%.1f,%.1f])\n",
                c.RStart, c.RStep, c.RLow, c.RHigh);
    std::printf("  drawBeam=%d  verbose=%d\n",
                int(c.drawBeam), int(c.verbose));
    cout << "----------------------------------------" << endl;
}

// ===========================================================================
// Globals -- Minuit's FCN signature has no user-data pointer, so the fit
// parameters that aren't varied (geometry, measured intensities, verbosity)
// are held in this one place and consulted during each function evaluation.
// ===========================================================================
static const Config* g_cfg = nullptr;

// ===========================================================================
// Physics / geometry
// ===========================================================================

// Relative intensity of a 2D Gaussian beam at point (Vx, Vy).  The beam is
// centered at (beamX, beamY) with a "radius" R defined as the distance at
// which the intensity falls to 1/e^2 of the peak value.
static double intensity(double beamX, double beamY, double R,
                        double Vx, double Vy) {
    const double r2 = (beamX - Vx) * (beamX - Vx)
                    + (beamY - Vy) * (beamY - Vy);
    return std::exp(-2.0 * r2 / (R * R));
}

// Numerically integrate `intensity` over one pad (an annular sector between
// radii [rIn, rOut] and azimuthal angles [thLo, thHi]).  Returns the double
// integral of r * intensity(r, theta) dr dtheta, which is proportional to
// the charge collected on that pad for a uniform-efficiency response.
static double integratePad(double beamX, double beamY, double R,
                           double rIn, double rOut,
                           double thLo, double thHi,
                           double epsR, double epsTheta) {
    double total = 0.0;
    for (double r = rIn; r < rOut; r += epsR) {
        for (double th = thLo; th < thHi; th += epsTheta) {
            total += r * epsR * epsTheta
                   * intensity(beamX, beamY, R, r * std::cos(th),
                                                r * std::sin(th));
        }
    }
    return total;
}

// Calculate the four quadrant-pad intensities, normalized to the central
// pad intensity.  Quadrant k covers theta in [k*pi/2, (k+1)*pi/2].
static void calculatePadIntensities(const Config& c,
                                    double beamX, double beamY, double R,
                                    double outPads[4]) {
    const double epsTheta = c.epsThetaDeg * TMath::DegToRad();

    // Central pad (used for normalization).
    const double A = integratePad(beamX, beamY, R,
                                  0.0, c.rCentral,
                                  0.0, 2.0 * TMath::Pi(),
                                  c.epsR, epsTheta);

    for (int q = 0; q < 4; ++q) {
        const double thLo = q * TMath::Pi() / 2.0;
        const double thHi = (q + 1) * TMath::Pi() / 2.0;
        outPads[q] = integratePad(beamX, beamY, R,
                                  c.rCentral, c.rOuterPad,
                                  thLo, thHi,
                                  c.epsR, epsTheta);
        outPads[q] /= A;
    }
}

// ===========================================================================
// Minuit FCN
// ===========================================================================

// Chi-square between measured and calculated (quadrant / central) ratios.
// `par[0..2]` are the fit parameters (x, y, R).  Measured counts come from
// g_cfg->padsExperimental: index 0 is the central pad, indices 1..4 are
// the four quadrants.
static double chiSquare(const double* par) {
    double padsCalc[4];
    calculatePadIntensities(*g_cfg, par[0], par[1], par[2], padsCalc);

    // Normalized measurements + statistical errors on the ratios.
    const double A = g_cfg->padsExperimental[0];
    double chi2 = 0.0;
    for (int i = 0; i < 4; ++i) {
        const double B = g_cfg->padsExperimental[i + 1];
        const double padsMeas = B / A;
        const double padsErr  = std::sqrt(1.0 / A + 1.0 / B) * padsMeas;
        const double resid    = (padsCalc[i] - padsMeas) / padsErr;
        chi2 += resid * resid;
    }

    if (g_cfg->verbose) {
        std::printf("  eval: X0=%7.2f Y0=%7.2f R=%6.2f  chi2=%.4e\n",
                    par[0], par[1], par[2], chi2);
    }
    return chi2;
}

// Signature required by TFitter::SetFCN.
static void fcnForMinuit(int& /*nDim*/, double* /*gout*/,
                         double& result, double par[], int /*flg*/) {
    result = chiSquare(par);
}

// ===========================================================================
// Drawing (optional -- only invoked when cfg.drawBeam is true)
// ===========================================================================

// Owned canvas / primitives kept alive while ROOT renders them.
struct DrawingState {
    TCanvas*  canvas = nullptr;
    TH2D*     hist   = nullptr;
    TEllipse *eCentral=nullptr, *eOuterPad=nullptr, *eGuardOuter=nullptr,
             *eGuardInner=nullptr;
    TLine    *l1=nullptr, *l2=nullptr, *l3=nullptr, *l4=nullptr,
             *l5=nullptr, *l6=nullptr, *l7=nullptr, *l8=nullptr;
};

// Draw the pad outlines (two concentric rings for the fiducial area, plus a
// dashed guard-ring display circle and radial division lines).
static void drawPads(DrawingState& st, const Config& c) {
    st.eCentral = new TEllipse(0, 0, c.rCentral, c.rCentral);
    st.eCentral->SetFillStyle(0);
    st.eCentral->SetLineWidth(2);

    st.eGuardInner = new TEllipse(0, 0, c.rGuardInner, c.rGuardInner);
    st.eGuardInner->SetFillStyle(0);
    st.eGuardInner->SetLineWidth(2);

    st.eGuardOuter = new TEllipse(0, 0, c.rGuardOuter, c.rGuardOuter);
    st.eGuardOuter->SetFillStyle(0);
    st.eGuardOuter->SetLineWidth(2);

    st.eOuterPad = new TEllipse(0, 0, c.rOuterPad, c.rOuterPad);
    st.eOuterPad->SetFillStyle(0);
    st.eOuterPad->SetLineStyle(2);
    st.eOuterPad->SetLineColor(5);
    st.eOuterPad->SetLineWidth(3);

    // Radial divisions at 0, 90, 180, 270 deg between the central and
    // guard-outer circles.
    st.l1 = new TLine(0,  c.rCentral,  0,  c.rGuardOuter);
    st.l2 = new TLine( c.rCentral, 0,   c.rGuardOuter, 0);
    st.l3 = new TLine(0, -c.rCentral,  0, -c.rGuardOuter);
    st.l4 = new TLine(-c.rCentral, 0,  -c.rGuardOuter, 0);

    // Guard-ring 45-deg divisions (display only; not integrated).
    const double s = 1.0 / std::sqrt(2.0);
    st.l5 = new TLine( c.rGuardInner*s,  c.rGuardInner*s,
                       c.rGuardOuter*s,  c.rGuardOuter*s);
    st.l6 = new TLine( c.rGuardInner*s, -c.rGuardInner*s,
                       c.rGuardOuter*s, -c.rGuardOuter*s);
    st.l7 = new TLine(-c.rGuardInner*s,  c.rGuardInner*s,
                      -c.rGuardOuter*s,  c.rGuardOuter*s);
    st.l8 = new TLine(-c.rGuardInner*s, -c.rGuardInner*s,
                      -c.rGuardOuter*s, -c.rGuardOuter*s);

    TLine* lines[] = { st.l1, st.l2, st.l3, st.l4,
                       st.l5, st.l6, st.l7, st.l8 };
    for (TLine* l : lines) l->SetLineWidth(2);

    st.eCentral->Draw();
    st.eGuardInner->Draw();
    st.eGuardOuter->Draw();
    st.eOuterPad->Draw();
    for (TLine* l : lines) l->Draw();
}

// Render the best-fit beam intensity as a 2D histogram and overlay the pad
// geometry.  Also prints the calculated pad intensities and a "collimation"
// diagnostic (integrated flux inside the full pad area vs. within 2R of the
// beam center).
static DrawingState* drawBeam(const Config& c,
                              double beamX, double beamY, double R) {
    auto* st = new DrawingState();
    const int Nbins = 200;
    const double xmin = -100.0, xmax = 100.0;

    st->hist = new TH2D("beamSpot",
        Form("beam spot X_{0}=%.1f Y_{0}=%.1f R=%.1f mm", beamX, beamY, R),
        Nbins, xmin, xmax, Nbins, xmin, xmax);
    st->hist->SetXTitle("X (mm)");
    st->hist->SetYTitle("Y (mm)");

    for (int n = 1; n <= Nbins; ++n) {
        for (int m = 1; m <= Nbins; ++m) {
            const double x = st->hist->GetXaxis()->GetBinCenter(n);
            const double y = st->hist->GetYaxis()->GetBinCenter(m);
            st->hist->SetBinContent(n, m, intensity(beamX, beamY, R, x, y));
        }
    }

    st->canvas = new TCanvas("c1", "RecoverBeamSpot", 700, 500);
    st->hist->Draw("colz");
    drawPads(*st, c);

    // Diagnostics: pad intensities and a "collimation" number.
    double padsCalc[4];
    calculatePadIntensities(c, beamX, beamY, R, padsCalc);

    const double epsTheta = c.epsThetaDeg * TMath::DegToRad();
    const double padsArea = integratePad(beamX, beamY, R,
                                         0.0, c.rOuterPad,
                                         0.0, 2.0 * TMath::Pi(),
                                         c.epsR, epsTheta);
    const double beamArea = integratePad(beamX, beamY, R,
                                         0.0, 2.0 * R,
                                         0.0, 2.0 * TMath::Pi(),
                                         c.epsR, epsTheta);
    std::printf("\nCalculated ratios: %.3f %.3f %.3f %.3f  "
                "collimation = %.3f\n",
                padsCalc[0], padsCalc[1], padsCalc[2], padsCalc[3],
                (beamArea > 0 ? padsArea / beamArea : 0.0));
    return st;
}

// ===========================================================================
// Main entry point
// ===========================================================================
int RecoverBeamSpot(const char* configFile = "RecoverBeamSpot.cfg") {
    Config cfg = loadConfig(configFile);
    g_cfg = &cfg;
    printConfig(cfg);

    TFitter* minimizer = new TFitter(3);
    minimizer->SetFCN(fcnForMinuit);
    minimizer->SetParameter(0, "X", cfg.xStart, cfg.xStep, cfg.xLow, cfg.xHigh);
    minimizer->SetParameter(1, "Y", cfg.yStart, cfg.yStep, cfg.yLow, cfg.yHigh);
    minimizer->SetParameter(2, "R", cfg.RStart, cfg.RStep, cfg.RLow, cfg.RHigh);

    // SIMPLEX first to get near the minimum, then MIGRAD for the local polish.
    minimizer->ExecuteCommand("SIMPLEX", nullptr, 0);
    minimizer->ExecuteCommand("MIGRAD",  nullptr, 0);

    const double bestX = minimizer->GetParameter(0);
    const double bestY = minimizer->GetParameter(1);
    const double bestR = minimizer->GetParameter(2);

    const double par[3] = { bestX, bestY, bestR };
    const double chi2   = chiSquare(par);

    std::printf("\nBest fit: X0 = %.2f mm, Y0 = %.2f mm, R = %.2f mm  "
                "(chi^2 = %.3e)\n", bestX, bestY, bestR, chi2);

    if (cfg.drawBeam) {
        drawBeam(cfg, bestX, bestY, bestR);
    }

    // TFitter owns TMinuit under the hood; deleting it is safe now that
    // we've read out the results.
    delete minimizer;
    return 0;
}

// Allow standalone compilation.  CLING (root macro mode) defines __CLING__
// so main() is skipped when loaded as a macro.
#ifndef __CLING__
int main(int argc, char** argv) {
    const char* cfg = (argc > 1) ? argv[1] : "RecoverBeamSpot.cfg";
    return recoverBeamSpot(cfg);
}
#endif
