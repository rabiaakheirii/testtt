#ifndef BLASTWAVE_DRAW_PARAMS_H
#define BLASTWAVE_DRAW_PARAMS_H

#include "../core/Constants.h"
#include "../core/FitConfig.h"
#include "../core/ParamIO.h"

#include "../input/FormatOfEverything.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"

#include <cctype>
#include <iostream>
#include <string>
#include <vector>

namespace bw {

constexpr unsigned kParamSourceIndividual = 1u;
constexpr unsigned kParamSourceGlobal = 2u;

constexpr unsigned kAllParticlesMask = (1u << kMaxParticles) - 1u;
constexpr unsigned kPositiveHadronsMask = (1u << 0) | (1u << 2) | (1u << 4);

constexpr int kParamCanvasW = 1200;
constexpr int kParamCanvasH = 1000;

// PHENIX-style axis ranges (NpartDrawParams_syst.cc / publication plots)
constexpr double kTBetaXMin = 0.1;
constexpr double kTBetaXMax = 1.20;
constexpr double kTBetaYMin = 0.02;
constexpr double kTBetaYMax = 0.22;
constexpr double kTNpartLogXMin = 1.;
constexpr double kTNpartLogXMax = 10000.;
constexpr double kTNpartYMin = 0.05;
constexpr double kTNpartYMax = 0.25;
constexpr double kBetaNpartYMin = 0.;
constexpr double kBetaNpartYMax = 1.;

constexpr Color_t kGlobalFitColor = kRed + 1;
constexpr Style_t kGlobalFitMarker = 20;
constexpr double kGlobalFitMarkerSize = 3.0;

inline std::string ToLowerCopy(std::string s)
{
    for (char& ch : s)
    {
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    }
    return s;
}

inline const char* SystemDisplayName(const FitConfig& cfg)
{
    const std::string s = ToLowerCopy(cfg.system);
    if (s == "pal")
    {
        return "p+Al";
    }
    if (s == "heau")
    {
        return "He+Au";
    }
    if (s == "cuau")
    {
        return "Cu+Au";
    }
    if (s == "uu")
    {
        return "U+U";
    }
    if (s == "auau")
    {
        return "Au+Au";
    }
    return cfg.system.c_str();
}

struct ParamPlotOptions
{
    unsigned sources = kParamSourceIndividual | kParamSourceGlobal;
    unsigned particleMask = kAllParticlesMask;
    int globalCharge = 0;  // 0 = positive hadrons, 1 = negative
    bool drawPenaltyCurve = true;
};

inline TLegend* NewGlobalLegend()
{
    TLegend* leg = new TLegend(0.66, 0.88, 0.95, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.045);
    return leg;
}

inline TLegend* NewParticleLegend()
{
    TLegend* leg = new TLegend(0.66, 0.73, 0.95, 0.86);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetNColumns(2);
    leg->SetTextSize(0.042);
    return leg;
}

inline void DrawSystemLabels(const FitConfig& cfg, double x = 0.20, double ySystem = 0.88)
{
    const std::string label = std::string(SystemDisplayName(cfg))
                              + ", #sqrt{s_{NN}} = "
                              + std::to_string(int(cfg.sqrt_s_NN)) + " GeV";

    TLatex* sys = new TLatex(x, ySystem, label.c_str());
    sys->SetNDC();
    sys->SetTextFont(42);
    sys->SetTextSize(0.048);
    sys->SetTextAlign(12);
    sys->Draw();
}

inline Color_t ParticlePlotColor(int particle)
{
    static const Color_t kColors[kMaxParticles] = {
        kRed + 1, kRed - 7, kBlue + 1, kBlue - 7, kGreen + 2, kGreen - 2
    };
    return kColors[particle];
}

inline Style_t ParticleMarkerStyle(int particle)
{
    static const Style_t kStyles[kMaxParticles] = {20, 24, 21, 25, 22, 26};
    return kStyles[particle];
}

// mode examples:
//   "both" / ""           — individual + global, all particles
//   "individual"          — only per-particle fits
//   "global"              — only global fit
//   "both:pos"            — only pi+, K+, p
//   "individual:pi"       — pi+ and pi-
//   "global:neg"          — global (negative charge)
inline ParamPlotOptions ParseParamPlotOptions(const std::string& spec)
{
    ParamPlotOptions opt;
    std::string s = ToLowerCopy(spec);

    if (s.find("individual") != std::string::npos && s.find("global") == std::string::npos)
    {
        opt.sources = kParamSourceIndividual;
    }
    else if (s.find("global") != std::string::npos && s.find("individual") == std::string::npos)
    {
        opt.sources = kParamSourceGlobal;
    }

    if (s.find("pos") != std::string::npos)
    {
        opt.particleMask = kPositiveHadronsMask;
    }
    if (s.find("neg") != std::string::npos)
    {
        opt.particleMask = (1u << 1) | (1u << 3) | (1u << 5);
        opt.globalCharge = 1;
    }
    if (s.find("pi") != std::string::npos)
    {
        opt.particleMask = (1u << 0) | (1u << 1);
    }
    if (s.find(":k") != std::string::npos || s == "k")
    {
        opt.particleMask = (1u << 2) | (1u << 3);
    }
    if (s.find(":p") != std::string::npos || s == "p")
    {
        opt.particleMask = (1u << 4) | (1u << 5);
    }
    if (s.find("nopredict") != std::string::npos)
    {
        opt.drawPenaltyCurve = false;
    }

    return opt;
}

inline bool UseIndividual(const ParamPlotOptions& opt)
{
    return (opt.sources & kParamSourceIndividual) != 0u;
}

inline bool UseGlobal(const ParamPlotOptions& opt)
{
    return (opt.sources & kParamSourceGlobal) != 0u;
}

inline bool ParticleSelected(const ParamPlotOptions& opt, int particle)
{
    return particle >= 0 && particle < kMaxParticles
           && ((opt.particleMask >> particle) & 1u) != 0u;
}

inline void ApplyParamPlotStyle()
{
    gStyle->SetOptStat(0);
    gStyle->SetGridStyle(3);
    gStyle->SetGridColor(kGray);
}

inline TCanvas* NewParamCanvas(const std::string& name, const std::string& title)
{
    ApplyParamPlotStyle();
    TCanvas* canvas = new TCanvas(name.c_str(), title.c_str(), kParamCanvasW, kParamCanvasH);
    canvas->SetGrid();
    canvas->SetLeftMargin(0.12);
    canvas->SetRightMargin(0.03);
    canvas->SetTopMargin(0.03);
    canvas->SetBottomMargin(0.12);
    canvas->SetTickx();
    canvas->SetTicky();
    return canvas;
}

inline void FormatParamPad(double xmin, double xmax, double ymin, double ymax,
                           const char* xTitle, const char* yTitle,
                           int nDivX = -105, int nDivY = -205,
                           double offsetX = 1.0, double offsetY = 1.0)
{
    Format_Pad(xmin, xmax, ymin, ymax, xTitle, yTitle,
               offsetX, offsetY, 0.05, 0.05, "", nDivY, nDivX);
    gPad->SetTickx();
    gPad->SetTicky();
}

inline void StyleGraph(TGraph* gr, Color_t color, Style_t marker, double size = 2.3)
{
    gr->SetMarkerColor(color);
    gr->SetMarkerStyle(marker);
    gr->SetMarkerSize(size);
    gr->SetLineColor(color);
    gr->SetLineWidth(2);
}

inline void StyleGlobalGraph(TGraph* gr)
{
    StyleGraph(gr, kGlobalFitColor, kGlobalFitMarker, kGlobalFitMarkerSize);
}

inline void DrawPenaltyCurve(const FitConfig& cfg, double betaMin, double betaMax)
{
    if (!cfg.t_beta_penalty.enabled)
    {
        return;
    }
    const auto& pen = cfg.t_beta_penalty;
    TF1* f = new TF1("t_beta_penalty",
                     "[0] + x*[1] + x*x*[2]", betaMin, betaMax);
    f->SetParameters(pen.a0, pen.a1, pen.a2);
    f->SetLineColor(kMagenta + 1);
    f->SetLineStyle(2);
    f->SetLineWidth(2);
    f->Draw("L SAME");
}

inline std::string ParamPdfPath(const FitConfig& cfg, const std::string& filename)
{
    if (!filename.empty() && filename[0] == '/')
    {
        return filename;
    }
    return cfg.output_dir + "/" + filename;
}

inline void SaveCanvas(TCanvas* canvas, const std::string& path)
{
    canvas->SaveAs(path.c_str());
    std::cout << "Saved " << path << "\n";
}

inline void DrawTBetaPlot(const FitConfig& cfg,
                          const FitParams individual[kMaxParticles][kMaxCentralities],
                          const GlobalParams global[2][kMaxCentralities],
                          const ParamPlotOptions& opt,
                          const std::string& filename)
{
    EnsureDir(cfg.output_dir);
    const int nCent = cfg.NumCentralities();

    TCanvas* canvas = NewParamCanvas("c_t_beta", "T_{kin} vs #beta_{T}");
    FormatParamPad(kTBetaXMin, kTBetaXMax, kTBetaYMin, kTBetaYMax,
                   "#beta_{T} [GeV]", "T_{kin} [GeV]",
                   -105, -205, 1.0, 1.25);

    if (opt.drawPenaltyCurve)
    {
        DrawPenaltyCurve(cfg, kTBetaXMin, kTBetaXMax);
    }

    TLegend* legGlobal = nullptr;
    TLegend* legParts = nullptr;
    if (UseGlobal(opt))
    {
        legGlobal = NewGlobalLegend();
    }
    if (UseIndividual(opt))
    {
        legParts = NewParticleLegend();
    }

    bool first = true;

    if (UseGlobal(opt))
    {
        std::vector<double> bx(nCent), ty(nCent), ex(nCent), ey(nCent);
        int n = 0;
        for (int c = 0; c < nCent; ++c)
        {
            const GlobalParams& gp = global[opt.globalCharge][c];
            if (gp.beta <= 0. || gp.T <= 0.)
            {
                continue;
            }
            bx[n] = gp.beta;
            ty[n] = gp.T;
            ex[n] = 0.;
            ey[n] = 0.;
            ++n;
        }

        if (n > 0)
        {
            TGraphErrors* gr = new TGraphErrors(n, bx.data(), ty.data(), ex.data(), ey.data());
            StyleGlobalGraph(gr);
            if (first)
            {
                gr->Draw("P");
                first = false;
            }
            else
            {
                gr->Draw("P SAME");
            }
            legGlobal->AddEntry(gr, "Global fit", "p");
        }
    }

    if (UseIndividual(opt))
    {
        TGraphErrors* partGraphs[kMaxParticles] = {};
        for (int part = 0; part < kMaxParticles; ++part)
        {
            if (!ParticleSelected(opt, part))
            {
                continue;
            }

            std::vector<double> bx(nCent), ty(nCent), ex(nCent), ey(nCent);
            int n = 0;
            for (int c = 0; c < nCent; ++c)
            {
                const FitParams& p = individual[part][c];
                if (p.beta <= 0. || p.T <= 0.)
                {
                    continue;
                }
                bx[n] = p.beta;
                ty[n] = p.T;
                ex[n] = p.beta_err;
                ey[n] = p.T_err;
                ++n;
            }
            if (n == 0)
            {
                continue;
            }

            TGraphErrors* gr = new TGraphErrors(n, bx.data(), ty.data(), ex.data(), ey.data());
            StyleGraph(gr, ParticlePlotColor(part), ParticleMarkerStyle(part), 2.0);
            if (first)
            {
                gr->Draw("P E");
                first = false;
            }
            else
            {
                gr->Draw("P E SAME");
            }
            partGraphs[part] = gr;
        }

        for (int part = 0; part < kMaxParticles; ++part)
        {
            if (partGraphs[part])
            {
                legParts->AddEntry(partGraphs[part], kParticleTitles[part], "p");
            }
        }
    }

    if (legGlobal)
    {
        legGlobal->Draw();
    }
    if (legParts)
    {
        legParts->Draw();
    }
    DrawSystemLabels(cfg);

    SaveCanvas(canvas, ParamPdfPath(cfg, filename));
    delete canvas;
}

inline void DrawVsNpartPlot(const FitConfig& cfg,
                            const char* yTitle,
                            const std::string& filename,
                            const FitParams individual[kMaxParticles][kMaxCentralities],
                            const GlobalParams global[2][kMaxCentralities],
                            const ParamPlotOptions& opt,
                            double yMin, double yMax,
                            double (*yInd)(const FitParams&),
                            double (*eInd)(const FitParams&),
                            double (*yGlob)(const GlobalParams&))
{
    EnsureDir(cfg.output_dir);
    const int nCent = cfg.NumCentralities();

    TCanvas* canvas = NewParamCanvas("c_vs_npart", yTitle);
    canvas->SetLogx();
    FormatParamPad(kTNpartLogXMin, kTNpartLogXMax, yMin, yMax,
                   "N_{part}", yTitle);

    TLegend* legGlobal = nullptr;
    TLegend* legParts = nullptr;
    if (UseGlobal(opt))
    {
        legGlobal = NewGlobalLegend();
    }
    if (UseIndividual(opt))
    {
        legParts = NewParticleLegend();
    }

    bool first = true;

    if (UseGlobal(opt))
    {
        std::vector<double> x(nCent), y(nCent), ex(nCent), ey(nCent);
        int n = 0;
        for (int c = 0; c < nCent; ++c)
        {
            const double np = cfg.Npart(c);
            const double val = yGlob(global[opt.globalCharge][c]);
            if (np <= 0. || val <= 0.)
            {
                continue;
            }
            x[n] = np;
            y[n] = val;
            ex[n] = 0.;
            ey[n] = 0.;
            ++n;
        }
        if (n > 0)
        {
            TGraphErrors* gr = new TGraphErrors(n, x.data(), y.data(), ex.data(), ey.data());
            StyleGlobalGraph(gr);
            if (first)
            {
                gr->Draw("P");
                first = false;
            }
            else
            {
                gr->Draw("P SAME");
            }
            legGlobal->AddEntry(gr, "Global fit", "p");
        }
    }

    if (UseIndividual(opt))
    {
        TGraphErrors* partGraphs[kMaxParticles] = {};
        for (int part = 0; part < kMaxParticles; ++part)
        {
            if (!ParticleSelected(opt, part))
            {
                continue;
            }
            std::vector<double> x(nCent), y(nCent), ex(nCent), ey(nCent);
            int n = 0;
            for (int c = 0; c < nCent; ++c)
            {
                const double np = cfg.Npart(c);
                const FitParams& p = individual[part][c];
                const double val = yInd(p);
                if (np <= 0. || val <= 0.)
                {
                    continue;
                }
                x[n] = np;
                y[n] = val;
                ex[n] = 0.;
                ey[n] = eInd(p);
                ++n;
            }
            if (n == 0)
            {
                continue;
            }
            TGraphErrors* gr = new TGraphErrors(n, x.data(), y.data(), ex.data(), ey.data());
            StyleGraph(gr, ParticlePlotColor(part), ParticleMarkerStyle(part), 2.0);
            if (first)
            {
                gr->Draw("P E");
                first = false;
            }
            else
            {
                gr->Draw("P E SAME");
            }
            partGraphs[part] = gr;
        }

        for (int part = 0; part < kMaxParticles; ++part)
        {
            if (partGraphs[part])
            {
                legParts->AddEntry(partGraphs[part], kParticleTitles[part], "p");
            }
        }
    }

    if (legGlobal)
    {
        legGlobal->Draw();
    }
    if (legParts)
    {
        legParts->Draw();
    }
    DrawSystemLabels(cfg);
    SaveCanvas(canvas, ParamPdfPath(cfg, filename));
    delete canvas;
}

inline void DrawParamPlots(const FitConfig& cfg,
                           const FitParams individual[kMaxParticles][kMaxCentralities],
                           const GlobalParams global[2][kMaxCentralities],
                           const ParamPlotOptions& opt,
                           const std::string& tag = "")
{
    const std::string suffix = tag.empty() ? "" : ("_" + tag);

    DrawTBetaPlot(cfg, individual, global, opt, "BlastWave_T_beta" + suffix + ".pdf");

    DrawVsNpartPlot(cfg, "T_{kin} [GeV]", "BlastWave_T_Npart" + suffix + ".pdf",
                    individual, global, opt, kTNpartYMin, kTNpartYMax,
                    [](const FitParams& p) { return p.T; },
                    [](const FitParams& p) { return p.T_err; },
                    [](const GlobalParams& gp) { return gp.T; });

    DrawVsNpartPlot(cfg, "#beta_{T} [GeV]", "BlastWave_beta_Npart" + suffix + ".pdf",
                    individual, global, opt, kBetaNpartYMin, kBetaNpartYMax,
                    [](const FitParams& p) { return p.beta; },
                    [](const FitParams& p) { return p.beta_err; },
                    [](const GlobalParams& gp) { return gp.beta; });
}

inline void DrawParamPlotsFromFiles(const FitConfig& cfg,
                                    const ParamPlotOptions& opt,
                                    const std::string& tag = "")
{
    FitParams individual[kMaxParticles][kMaxCentralities]{};
    GlobalParams global[2][kMaxCentralities]{};

    if (UseIndividual(opt))
    {
        try
        {
            ReadIndividualParams(cfg, individual);
        }
        catch (const std::exception& e)
        {
            std::cout << "WARNING: " << e.what() << " — skip individual params\n";
            if (!UseGlobal(opt))
            {
                throw;
            }
        }
    }
    if (UseGlobal(opt))
    {
        try
        {
            ReadGlobalParams(cfg, global);
        }
        catch (const std::exception& e)
        {
            std::cout << "WARNING: " << e.what() << " — skip global params\n";
            if (!UseIndividual(opt))
            {
                throw;
            }
        }
    }

    DrawParamPlots(cfg, individual, global, opt, tag);
}

}  // namespace bw

#endif
