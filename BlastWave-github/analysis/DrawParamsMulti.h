#ifndef BLASTWAVE_DRAW_PARAMS_MULTI_H
#define BLASTWAVE_DRAW_PARAMS_MULTI_H

#include "DrawParams.h"

#include <sstream>
#include <stdexcept>
#include <vector>

namespace bw {

enum class ParamDataSource
{
    Hand,
    Global,
    Individual,
    Both
};

struct SystemParamSeries
{
    FitConfig cfg;
    std::vector<double> npart;
    std::vector<double> T;
    std::vector<double> beta;
    std::vector<double> T_err;
    std::vector<double> beta_err;
};

struct SystemFitResults
{
    FitConfig cfg;
    FitParams individual[kMaxParticles][kMaxCentralities]{};
    GlobalParams global[2][kMaxCentralities]{};
    bool hasIndividual = false;
    bool hasGlobal = false;
};

inline std::string ToLowerTrim(std::string s)
{
    s = ToLowerCopy(s);
    const auto start = s.find_first_not_of(" \t");
    if (start == std::string::npos)
    {
        return "";
    }
    const auto end = s.find_last_not_of(" \t");
    return s.substr(start, end - start + 1);
}

inline std::string ConfigPathForSystem(const std::string& token)
{
    const std::string s = ToLowerTrim(token);
    if (s.empty())
    {
        throw std::runtime_error("Empty system name in list");
    }
    if (s.find('/') != std::string::npos || s.find(".json") != std::string::npos)
    {
        return token;
    }

    if (s == "pal" || s == "p+al" || s == "p_al")
    {
        return "config/phenix_pal.json";
    }
    if (s == "heau" || s == "he+au" || s == "he_au")
    {
        return "config/phenix_heau.json";
    }
    if (s == "cuau" || s == "cu+au" || s == "cu_au")
    {
        return "config/phenix_cuau.json";
    }
    if (s == "uu" || s == "u+u" || s == "u_u")
    {
        return "config/phenix_uu.json";
    }
    if (s == "auau" || s == "au+au" || s == "au_au")
    {
        return "config/phenix_auau.json";
    }

    return "config/phenix_" + s + ".json";
}

inline std::vector<std::string> SplitCommaList(const std::string& text)
{
    std::vector<std::string> out;
    std::stringstream ss(text);
    std::string item;
    while (std::getline(ss, item, ','))
    {
        item = ToLowerTrim(item);
        if (!item.empty())
        {
            out.push_back(item);
        }
    }
    if (out.empty())
    {
        throw std::runtime_error("System list is empty");
    }
    return out;
}

inline ParamDataSource ParseParamDataSource(const std::string& text)
{
    const std::string s = ToLowerTrim(text);
    if (s.empty() || s == "global" || s == "fit")
    {
        return ParamDataSource::Global;
    }
    if (s == "both")
    {
        return ParamDataSource::Both;
    }
    if (s == "hand" || s == "hand_params")
    {
        return ParamDataSource::Hand;
    }
    if (s == "individual" || s == "ind")
    {
        return ParamDataSource::Individual;
    }
    throw std::runtime_error("Unknown param source: " + text
                             + " (use global, both, individual, hand)");
}

inline std::string SourceTag(ParamDataSource source)
{
    switch (source)
    {
        case ParamDataSource::Hand: return "hand";
        case ParamDataSource::Global: return "global";
        case ParamDataSource::Individual: return "individual";
        case ParamDataSource::Both: return "both";
    }
    return "global";
}

inline Color_t SystemPlotColor(int index)
{
    static const Color_t kColors[] = {
        kAzure + 4, kOrange + 7, kSpring - 5, kViolet - 6, kGray + 2, kTeal + 2
    };
    return kColors[index % 6];
}

inline Style_t SystemMarkerStyle(int index)
{
    static const Style_t kStyles[] = {20, 21, 22, 23, 24, 25};
    return kStyles[index % 6];
}

inline Color_t GlobalSeriesColor(int systemIndex, int nSystems)
{
    return (nSystems == 1) ? kGlobalFitColor : SystemPlotColor(systemIndex);
}

inline SystemParamSeries LoadHandSeries(const FitConfig& cfg)
{
    std::cout << "WARNING: hand_params are init/reference values, not fit results. "
              << "Use source=global or source=both for fit output.\n";

    SystemParamSeries series;
    series.cfg = cfg;
    const int nCent = cfg.NumCentralities();
    series.npart.resize(nCent);
    series.T.resize(nCent);
    series.beta.resize(nCent);
    series.T_err.assign(nCent, 0.);
    series.beta_err.assign(nCent, 0.);

    for (int c = 0; c < nCent; ++c)
    {
        series.npart[c] = cfg.Npart(c);
        series.T[c] = (c < cfg.hand.T.size()) ? cfg.hand.T[c] : 0.;
        series.beta[c] = (c < cfg.hand.beta.size()) ? cfg.hand.beta[c] : 0.;
    }
    return series;
}

inline SystemParamSeries LoadGlobalSeries(const FitConfig& cfg, int charge)
{
    GlobalParams global[2][kMaxCentralities]{};
    ReadGlobalParams(cfg, global);

    SystemParamSeries series;
    series.cfg = cfg;
    const int nCent = cfg.NumCentralities();
    series.npart.resize(nCent);
    series.T.resize(nCent);
    series.beta.resize(nCent);
    series.T_err.assign(nCent, 0.);
    series.beta_err.assign(nCent, 0.);

    for (int c = 0; c < nCent; ++c)
    {
        series.npart[c] = cfg.Npart(c);
        series.T[c] = global[charge][c].T;
        series.beta[c] = global[charge][c].beta;
    }
    return series;
}

inline SystemParamSeries LoadIndividualSeries(const FitConfig& cfg, int particle)
{
    FitParams individual[kMaxParticles][kMaxCentralities]{};
    ReadIndividualParams(cfg, individual);

    SystemParamSeries series;
    series.cfg = cfg;
    const int nCent = cfg.NumCentralities();
    series.npart.resize(nCent);
    series.T.resize(nCent);
    series.beta.resize(nCent);
    series.T_err.resize(nCent);
    series.beta_err.resize(nCent);

    for (int c = 0; c < nCent; ++c)
    {
        series.npart[c] = cfg.Npart(c);
        series.T[c] = individual[particle][c].T;
        series.beta[c] = individual[particle][c].beta;
        series.T_err[c] = individual[particle][c].T_err;
        series.beta_err[c] = individual[particle][c].beta_err;
    }
    return series;
}

inline SystemFitResults LoadSystemFitResults(const FitConfig& cfg,
                                             bool needIndividual,
                                             bool needGlobal)
{
    SystemFitResults results;
    results.cfg = cfg;

    if (needIndividual)
    {
        try
        {
            ReadIndividualParams(cfg, results.individual);
            results.hasIndividual = true;
        }
        catch (const std::exception& e)
        {
            std::cout << "WARNING: " << e.what() << " — skip individual for "
                      << cfg.system << "\n";
        }
    }

    if (needGlobal)
    {
        try
        {
            ReadGlobalParams(cfg, results.global);
            results.hasGlobal = true;
        }
        catch (const std::exception& e)
        {
            throw std::runtime_error(std::string(e.what())
                                     + " Run global fit first for " + cfg.system);
        }
    }

    return results;
}

inline SystemParamSeries LoadSystemSeries(const FitConfig& cfg,
                                          ParamDataSource source,
                                          int charge,
                                          int particle)
{
    switch (source)
    {
        case ParamDataSource::Hand:
            return LoadHandSeries(cfg);
        case ParamDataSource::Global:
        case ParamDataSource::Both:
            return LoadGlobalSeries(cfg, charge);
        case ParamDataSource::Individual:
            return LoadIndividualSeries(cfg, particle);
    }
    return LoadGlobalSeries(cfg, charge);
}

inline std::vector<SystemParamSeries> LoadSystemSeriesList(
    const std::vector<std::string>& systemTokens,
    ParamDataSource source,
    int charge,
    int particle)
{
    std::vector<SystemParamSeries> out;
    out.reserve(systemTokens.size());
    for (const std::string& token : systemTokens)
    {
        const std::string path = ConfigPathForSystem(token);
        std::cout << "Loading " << path << "\n";
        const FitConfig cfg = FitConfig::Load(path);
        out.push_back(LoadSystemSeries(cfg, source, charge, particle));
    }
    return out;
}

inline std::vector<SystemFitResults> LoadSystemFitResultsList(
    const std::vector<std::string>& systemTokens,
    bool needIndividual,
    bool needGlobal)
{
    std::vector<SystemFitResults> out;
    out.reserve(systemTokens.size());
    for (const std::string& token : systemTokens)
    {
        const std::string path = ConfigPathForSystem(token);
        std::cout << "Loading fit results from " << path << "\n";
        const FitConfig cfg = FitConfig::Load(path);
        out.push_back(LoadSystemFitResults(cfg, needIndividual, needGlobal));
    }
    return out;
}

inline std::string SystemsTag(const std::vector<SystemParamSeries>& series)
{
    std::string tag;
    for (const auto& s : series)
    {
        if (!tag.empty())
        {
            tag += "_";
        }
        tag += s.cfg.system;
    }
    return tag;
}

inline std::string SystemsTagFromFit(const std::vector<SystemFitResults>& systems)
{
    std::string tag;
    for (const auto& s : systems)
    {
        if (!tag.empty())
        {
            tag += "_";
        }
        tag += s.cfg.system;
    }
    return tag;
}

inline TLegend* NewSystemsLegend()
{
    TLegend* leg = new TLegend(0.58, 0.70, 0.95, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.042);
    return leg;
}

inline void DrawMultiTBetaPlot(const std::vector<SystemParamSeries>& series,
                               const std::string& outputPath,
                               bool drawPenaltyCurve = true)
{
    EnsureDir(outputPath.substr(0, outputPath.find_last_of('/')));

    TCanvas* canvas = NewParamCanvas("c_t_beta_syst", "T_{kin} vs #beta_{T}");
    FormatParamPad(kTBetaXMin, kTBetaXMax, kTBetaYMin, kTBetaYMax,
                   "#beta_{T} [GeV]", "T_{kin} [GeV]",
                   -105, -205, 1.0, 1.25);

    if (drawPenaltyCurve && !series.empty())
    {
        DrawPenaltyCurve(series.front().cfg, kTBetaXMin, kTBetaXMax);
    }

    TLegend* leg = NewSystemsLegend();
    bool first = true;
    const auto nSystems = series.size();

    for (std::size_t i = 0; i < nSystems; ++i)
    {
        const auto& s = series[i];
        const auto nCent = s.npart.size();
        std::vector<double> bx(nCent), ty(nCent), ex(nCent), ey(nCent);
        int n = 0;
        for (std::size_t c = 0; c < nCent; ++c)
        {
            if (s.beta[c] <= 0. || s.T[c] <= 0.)
            {
                continue;
            }
            bx[n] = s.beta[c];
            ty[n] = s.T[c];
            ex[n] = s.beta_err[c];
            ey[n] = s.T_err[c];
            ++n;
        }
        if (n == 0)
        {
            continue;
        }

        TGraphErrors* gr = new TGraphErrors(n, bx.data(), ty.data(), ex.data(), ey.data());
        StyleGraph(gr,
                   GlobalSeriesColor(i, nSystems),
                   (nSystems == 1) ? kGlobalFitMarker : SystemMarkerStyle(i),
                   (nSystems == 1) ? kGlobalFitMarkerSize : 2.8);
        if (first)
        {
            gr->Draw("P");
            first = false;
        }
        else
        {
            gr->Draw("P SAME");
        }
        const std::string label = (nSystems == 1)
                                      ? "Global fit"
                                      : std::string(SystemDisplayName(s.cfg)) + " global";
        leg->AddEntry(gr, label.c_str(), "p");
    }

    leg->Draw();
    SaveCanvas(canvas, outputPath);
    delete canvas;
}

inline void DrawMultiTBetaPlotFromFit(const std::vector<SystemFitResults>& systems,
                                      const ParamPlotOptions& opt,
                                      const std::string& outputPath,
                                      bool drawPenaltyCurve = true)
{
    EnsureDir(outputPath.substr(0, outputPath.find_last_of('/')));

    TCanvas* canvas = NewParamCanvas("c_t_beta_syst", "T_{kin} vs #beta_{T}");
    FormatParamPad(kTBetaXMin, kTBetaXMax, kTBetaYMin, kTBetaYMax,
                   "#beta_{T} [GeV]", "T_{kin} [GeV]",
                   -105, -205, 1.0, 1.25);

    if (drawPenaltyCurve && !systems.empty())
    {
        DrawPenaltyCurve(systems.front().cfg, kTBetaXMin, kTBetaXMax);
    }

    TLegend* legGlobal = nullptr;
    TLegend* legParts = nullptr;
    if (UseGlobal(opt))
    {
        legGlobal = (systems.size() == 1) ? NewGlobalLegend() : NewSystemsLegend();
    }
    if (UseIndividual(opt))
    {
        legParts = NewParticleLegend();
    }

    bool first = true;
    const auto nSystems = systems.size();
    TGraphErrors* partGraphs[kMaxParticles] = {};

    for (std::size_t i = 0; i < nSystems; ++i)
    {
        const SystemFitResults& sys = systems[i];
        const FitConfig& cfg = sys.cfg;
        const int nCent = cfg.NumCentralities();

        if (UseGlobal(opt) && sys.hasGlobal)
        {
            std::vector<double> bx(nCent), ty(nCent), ex(nCent), ey(nCent);
            int n = 0;
            for (int c = 0; c < nCent; ++c)
            {
                const GlobalParams& gp = sys.global[opt.globalCharge][c];
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
                StyleGraph(gr,
                           GlobalSeriesColor(i, nSystems),
                           (nSystems == 1) ? kGlobalFitMarker : SystemMarkerStyle(i),
                           (nSystems == 1) ? kGlobalFitMarkerSize : 2.8);
                if (first)
                {
                    gr->Draw("P");
                    first = false;
                }
                else
                {
                    gr->Draw("P SAME");
                }
                const std::string label = (nSystems == 1)
                                              ? "Global fit"
                                              : std::string(SystemDisplayName(cfg)) + " global";
                legGlobal->AddEntry(gr, label.c_str(), "p");
            }
        }

        if (UseIndividual(opt) && sys.hasIndividual)
        {
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
                    const FitParams& p = sys.individual[part][c];
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
                if (nSystems == 1)
                {
                    partGraphs[part] = gr;
                }
            }
        }
    }

    if (legParts && nSystems == 1)
    {
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
    if (legParts && nSystems == 1)
    {
        legParts->Draw();
    }
    if (nSystems == 1)
    {
        DrawSystemLabels(systems.front().cfg);
    }

    SaveCanvas(canvas, outputPath);
    delete canvas;
}

inline void DrawMultiVsNpartPlot(const std::vector<SystemParamSeries>& series,
                                 const char* yTitle,
                                 const std::string& outputPath,
                                 double yMin,
                                 double yMax,
                                 bool useT)
{
    EnsureDir(outputPath.substr(0, outputPath.find_last_of('/')));

    TCanvas* canvas = NewParamCanvas("c_vs_npart_syst", yTitle);
    canvas->SetLogx();
    FormatParamPad(kTNpartLogXMin, kTNpartLogXMax, yMin, yMax,
                   "N_{part}", yTitle);

    TLegend* leg = NewSystemsLegend();
    bool first = true;
    const auto nSystems = series.size();

    for (std::size_t i = 0; i < nSystems; ++i)
    {
        const auto& s = series[i];
        const auto nCent = s.npart.size();
        std::vector<double> x(nCent), y(nCent), ex(nCent), ey(nCent);
        int n = 0;
        for (std::size_t c = 0; c < nCent; ++c)
        {
            const double np = s.npart[c];
            const double val = useT ? s.T[c] : s.beta[c];
            if (np <= 0. || val <= 0.)
            {
                continue;
            }
            x[n] = np;
            y[n] = val;
            ex[n] = 0.;
            ey[n] = useT ? s.T_err[c] : s.beta_err[c];
            ++n;
        }
        if (n == 0)
        {
            continue;
        }

        TGraphErrors* gr = new TGraphErrors(n, x.data(), y.data(), ex.data(), ey.data());
        StyleGraph(gr,
                   GlobalSeriesColor(i, nSystems),
                   (nSystems == 1) ? kGlobalFitMarker : SystemMarkerStyle(i),
                   (nSystems == 1) ? kGlobalFitMarkerSize : 2.8);
        if (first)
        {
            gr->Draw("P");
            first = false;
        }
        else
        {
            gr->Draw("P SAME");
        }
        const std::string label = (nSystems == 1)
                                      ? "Global fit"
                                      : std::string(SystemDisplayName(s.cfg)) + " global";
        leg->AddEntry(gr, label.c_str(), "p");
    }

    leg->Draw();
    SaveCanvas(canvas, outputPath);
    delete canvas;
}

inline void DrawMultiSystemParamPlots(const std::vector<std::string>& systemTokens,
                                      ParamDataSource source,
                                      const std::string& outputDir,
                                      const ParamPlotOptions& plotOpt,
                                      int charge = 0,
                                      int particle = 0,
                                      bool drawPenaltyCurve = true)
{
    EnsureDir(outputDir);
    const std::string sourceTag = SourceTag(source);

    if (source == ParamDataSource::Both)
    {
        ParamPlotOptions opt = plotOpt;
        opt.sources = kParamSourceIndividual | kParamSourceGlobal;
        opt.globalCharge = charge;

        const std::vector<SystemFitResults> systems =
            LoadSystemFitResultsList(systemTokens, true, true);
        const std::string tag = SystemsTagFromFit(systems);
        const std::vector<SystemParamSeries> globalSeries =
            LoadSystemSeriesList(systemTokens, ParamDataSource::Global, charge, particle);

        DrawMultiTBetaPlotFromFit(systems, opt,
                                  outputDir + "/BlastWave_T_beta_" + tag + "_" + sourceTag + ".pdf",
                                  drawPenaltyCurve);

        DrawMultiVsNpartPlot(globalSeries, "T_{kin} [GeV]",
                             outputDir + "/BlastWave_T_Npart_" + tag + "_" + sourceTag + ".pdf",
                             kTNpartYMin, kTNpartYMax, true);

        DrawMultiVsNpartPlot(globalSeries, "#beta_{T}",
                             outputDir + "/BlastWave_beta_Npart_" + tag + "_" + sourceTag + ".pdf",
                             kBetaNpartYMin, kBetaNpartYMax, false);
        return;
    }

    const std::vector<SystemParamSeries> series =
        LoadSystemSeriesList(systemTokens, source, charge, particle);
    const std::string tag = SystemsTag(series);

    DrawMultiTBetaPlot(series,
                       outputDir + "/BlastWave_T_beta_" + tag + "_" + sourceTag + ".pdf",
                       drawPenaltyCurve);

    DrawMultiVsNpartPlot(series, "T_{kin} [GeV]",
                         outputDir + "/BlastWave_T_Npart_" + tag + "_" + sourceTag + ".pdf",
                         kTNpartYMin, kTNpartYMax, true);

    DrawMultiVsNpartPlot(series, "#beta_{T}",
                         outputDir + "/BlastWave_beta_Npart_" + tag + "_" + sourceTag + ".pdf",
                         kBetaNpartYMin, kBetaNpartYMax, false);
}

}  // namespace bw

#endif
