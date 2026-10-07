#ifndef BLASTWAVE_DRAW_SPECTRA_H
#define BLASTWAVE_DRAW_SPECTRA_H

#include "../core/BlastWaveModel.h"
#include "../core/Constants.h"
#include "../core/FitConfig.h"
#include "../core/ParamIO.h"
#include "../core/SpectrumData.h"
#include "../fit/GlobalFit.h"
#include "../fit/IndividualFit.h"

#include "../input/FormatOfEverything.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"

#include <functional>
#include <iostream>
#include <sstream>
#include <string>

namespace bw {

inline Color_t CentralityColor(int centrality)
{
    static const Color_t kColors[] = {
        kRed, kBlue, kGreen + 2, kBlack, kMagenta + 2, kBlue + 3, kCyan + 2
    };
    return kColors[centrality % 7];
}

inline TString XAxisTitle(const FitConfig& cfg)
{
    if (cfg.variable == "mt")
    {
        return "m_{T}-m_{0} [GeV/c^{2}]";
    }
    return "p_{T} [GeV/c]";
}

// Most central (c=0) on top: multiply by 10^(N-1-c).
inline int CentralityDrawExponent(int centrality, int nCentralities)
{
    return std::max(0, nCentralities - 1 - centrality);
}

inline double CentralityDrawScale(int centrality, int nCentralities)
{
    return TMath::Power(10.0, CentralityDrawExponent(centrality, nCentralities));
}

inline TString CentralityLegendLabel(const FitConfig& cfg, int centrality)
{
    const int exp = CentralityDrawExponent(centrality, cfg.NumCentralities());
    if (exp > 0)
    {
        return TString::Format("%s (#times10^{%d})",
                               cfg.centrality_labels[centrality].c_str(), exp);
    }
    return cfg.centrality_labels[centrality].c_str();
}

// TF1 is painted at SaveAs time; do not rely on temporary SetParameter scaling.
inline TGraph* MakeScaledFitGraph(TF1* f, double scale, int nPoints = 80)
{
    const double xlo = f->GetXmin();
    const double xhi = f->GetXmax();
    auto* g = new TGraph(nPoints + 1);
    for (int i = 0; i <= nPoints; ++i)
    {
        const double x = xlo + (xhi - xlo) * i / nPoints;
        g->SetPoint(i, x, f->Eval(x) * scale);
    }
    return g;
}

inline void FormatSpectraPadForConfig(const FitConfig& cfg,
                                      double ymin,
                                      double ymax,
                                      double texScale = 1.)
{
    const double ll = 0.01;
    const double rl = 2.49;
    const double pad_offset_x = 0.95;
    const double pad_offset_y = 0.95;
    const double pad_tsize = 0.09 * texScale;
    const double pad_lsize = 0.08 * texScale;
    const TString pad_title_y = "d^{2}N/(p_{T}dydp_{T})";

    Format_Pad(ll, rl, ymin, ymax,
               XAxisTitle(cfg), pad_title_y,
               pad_offset_x, pad_offset_y, pad_tsize, pad_lsize, "");
}

inline void SyncGlobalFunctions(GlobalFitter& fitter,
                                const GlobalParams global[2][kMaxCentralities],
                                int nCentralities)
{
    for (int part = 0; part < kMaxParticles; ++part)
    {
        for (int c = 0; c < nCentralities; ++c)
        {
            TF1* f = fitter.Function(part, c);
            if (!f)
            {
                continue;
            }
            const FitParams p = ParamsFromGlobal(part, c, global);
            f->SetParameters(p.C, p.T, p.beta, p.mass);
            f->SetLineColor(CentralityColor(c));
        }
    }
}

inline void SetupSpectraCanvas(TCanvas* canvas, int cols, int rows)
{
    canvas->Divide(cols, rows, 0, 0);
    const int nPads = cols * rows;
    for (int i = 0; i < nPads; ++i)
    {
        canvas->cd(i + 1);
        gPad->SetLogy();
        gPad->SetTickx(1);
        gPad->SetTicky(1);

        const bool leftCol = (i % cols == 0);
        const bool rightCol = ((i + 1) % cols == 0);
        const bool bottomRow = (i >= cols * (rows - 1));

        gPad->SetLeftMargin(leftCol ? 0.18 : 0.0);
        gPad->SetRightMargin(rightCol ? 0.0 : 0.0);
        gPad->SetTopMargin(0.0);
        gPad->SetBottomMargin(bottomRow ? 0.20 : 0.0);
    }
}

inline void DrawParticlePanel(const FitConfig& cfg,
                              const Dataset& data,
                              const std::function<TF1*(int particle, int centrality)>& getFunction,
                              int particle,
                              int padIndex,
                              TCanvas* canvas)
{
    canvas->cd(padIndex);
    const double texScale = (padIndex <= 3) ? 1.0 : 0.9;
    const double shiftX = (particle % 2 == 0) ? 0. : 0.08;
    const int nCent = cfg.NumCentralities();

    double ymin = 1e-4;
    double ymax = 1e5;
    cfg.PanelYLimits(particle, ymin, ymax);
    FormatSpectraPadForConfig(cfg, ymin, ymax, texScale);

    TLegend* leg = new TLegend(0.67 - shiftX, 0.68, 0.9 - shiftX, 0.97);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    //leg->SetNColumns(2);
    leg->SetTextSize(0.057 * texScale);

    TLatex* title = new TLatex();
    title->SetNDC();
    title->SetTextFont(42);
    title->SetTextSize(0.085 * texScale);
    title->SetLineWidth(2);
    title->DrawLatex(0.2, 0.86, kParticleTitles[particle]);

    bool first = true;
    for (int c = 0; c < nCent; ++c)
    {
        TGraphErrors* gr = data.Spectrum(particle, c);
        if (!gr || gr->GetN() <= 0)
        {
            continue;
        }

        const double scale = CentralityDrawScale(c, nCent);
        const Color_t color = CentralityColor(c);

        TGraphErrors* grDraw = static_cast<TGraphErrors*>(gr->Clone());
        if (grDraw->GetListOfFunctions())
        {
            grDraw->GetListOfFunctions()->Clear();
        }
        for (int i = 0; i < grDraw->GetN(); ++i)
        {
            double x = 0.;
            double y = 0.;
            grDraw->GetPoint(i, x, y);
            grDraw->SetPoint(i, x, y * scale);
            grDraw->SetPointError(i, grDraw->GetErrorX(i), grDraw->GetErrorY(i) * scale);
        }

        grDraw->SetMarkerColor(color);
        grDraw->SetMarkerStyle(8);
        grDraw->SetMarkerSize(0.7);
        grDraw->SetLineColor(color);
        if (first)
        {
            grDraw->Draw("P");
            first = false;
        }
        else
        {
            grDraw->Draw("P SAME");
        }

        TF1* f = getFunction(particle, c);
        if (f && data.HasSpectrum(particle, c))
        {
            TGraph* fitDraw = MakeScaledFitGraph(f, scale);
            fitDraw->SetLineColor(color);
            fitDraw->SetLineWidth(2);
            fitDraw->Draw("L SAME");
            leg->AddEntry(fitDraw, CentralityLegendLabel(cfg, c).Data(), "l");
        }
        else
        {
            leg->AddEntry(grDraw, CentralityLegendLabel(cfg, c).Data(), "p");
        }
    }

    leg->Draw();
}

inline void SaveIndividualFitPdf(const FitConfig& cfg,
                                 const Dataset& data,
                                 const std::function<TF1*(int particle, int centrality)>& getFunction,
                                 const std::string& filename)
{
    gStyle->SetOptStat(0);
    EnsureDir(cfg.output_dir);

    TCanvas* canvas = new TCanvas("c_bw_individual",
                                  ("BlastWave " + cfg.system).c_str(),
                                  1200, 1200);
    SetupSpectraCanvas(canvas, 2, 3);

    for (int part = 0; part < kMaxParticles; ++part)
    {
        DrawParticlePanel(cfg, data, getFunction, part, part + 1, canvas);
    }

    canvas->SaveAs(filename.c_str());
    std::cout << "Saved " << filename << "\n";
    delete canvas;
}

inline void SaveGlobalFitPdf(const FitConfig& cfg,
                             const Dataset& data,
                             GlobalFitter& fitter,
                             const GlobalParams global[2][kMaxCentralities],
                             const std::string& filename)
{
    gStyle->SetOptStat(0);
    EnsureDir(cfg.output_dir);
    SyncGlobalFunctions(fitter, global, cfg.NumCentralities());

    TCanvas* canvas = new TCanvas("c_bw_global",
                                  ("BlastWave Global " + cfg.system).c_str(),
                                  1200, 1200);
    SetupSpectraCanvas(canvas, 2, 3);

    for (int part = 0; part < kMaxParticles; ++part)
    {
        DrawParticlePanel(cfg, data,
                          [&](int p, int c) { return fitter.Function(p, c); },
                          part, part + 1, canvas);
    }

    canvas->SaveAs(filename.c_str());
    std::cout << "Saved " << filename << "\n";
    delete canvas;
}

inline void SaveIndividualFitPdfFromFile(const FitConfig& cfg,
                                           const Dataset& data,
                                           const std::string& filename)
{
    FitParams params[kMaxParticles][kMaxCentralities]{};
    ReadIndividualParams(cfg, params);

    TF1* integrand = new TF1("bw_draw_integrand", bwfitfunc, 0.01, 10., 5);
    MyIntegFunc integ(integrand);
    TF1* funcs[kMaxParticles][kMaxCentralities]{};

    const int nCent = cfg.NumCentralities();
    for (int part = 0; part < kMaxParticles; ++part)
    {
        const auto fitRange = cfg.FitRangeForParticle(part);
        for (int c = 0; c < nCent; ++c)
        {
            if (!data.HasSpectrum(part, c))
                continue;

            std::ostringstream name;
            name << "bw_draw_" << part << "_" << c;
            funcs[part][c] = new TF1(name.str().c_str(), integ,
                                     fitRange.xmin, fitRange.xmax, 4);
            const FitParams& p = params[part][c];
            funcs[part][c]->SetParameters(p.C, p.T, p.beta, kMasses[part]);
            funcs[part][c]->FixParameter(3, kMasses[part]);
        }
    }

    SaveIndividualFitPdf(cfg, data,
                         [&](int p, int c) { return funcs[p][c]; },
                         filename);

    for (int part = 0; part < kMaxParticles; ++part)
        for (int c = 0; c < kMaxCentralities; ++c)
            delete funcs[part][c];
    delete integrand;
}

inline std::string OutputPdfPath(const FitConfig& cfg, const std::string& filename)
{
    if (filename.empty())
    {
        return cfg.output_dir + "/BlastWave.pdf";
    }
    if (!filename.empty() && filename[0] == '/')
    {
        return filename;
    }
    return cfg.output_dir + "/" + filename;
}

inline void DrawHandNoFit(const FitConfig& cfg,
                          const Dataset& data,
                          const std::string& filename = "BlastWave_handNoFit.pdf")
{
    FitStrategy strategy{InitStrategy::HandNoFit, LimitsStrategy::Fix};
    IndividualFitter fitter(data, strategy);
    fitter.Run();
    SaveIndividualFitPdf(cfg, data,
                         [&](int p, int c) { return fitter.Function(p, c); },
                         OutputPdfPath(cfg, filename));
}

// Redraw individual fits from output/txtParams/BWparams.txt.
inline void DrawIndividualResults(const FitConfig& cfg,
                                  const Dataset& data,
                                  const std::string& filename = "")
{
    SaveIndividualFitPdfFromFile(cfg, data,
                                 OutputPdfPath(cfg, filename.empty() ? "BlastWave.pdf" : filename));
}

// Redraw global fit from output/txtParams/GlobalBWparams.txt.
inline void DrawGlobalResults(const FitConfig& cfg,
                              const Dataset& data,
                              const std::string& filename = "")
{
    GlobalParams global[2][kMaxCentralities]{};
    ReadGlobalParams(cfg, global);

    GlobalFitter fitter(data);
    fitter.InitFunctions();
    SaveGlobalFitPdf(cfg, data, fitter, global,
                     OutputPdfPath(cfg, filename.empty() ? "BlastWave_global.pdf" : filename));
}

}  // namespace bw

#endif
