// Overlay T(Npart), beta(Npart), T(beta) for several collision systems.
//
// Reads fit results from output/<system>/txtParams/ (same files as drawParams).
//
// Usage (from repository root):
//   root -l 'macros/draw_params_syst.C("UU")'                    // both = global + individual
//   root -l 'macros/draw_params_syst.C("CuAu,UU", "global")'     // global fit only
//   root -l 'macros/draw_params_syst.C("CuAu,UU", "both", "output/params_compare", "pos")'
//
// Arguments:
//   systems   — comma-separated names (CuAu, UU, pAl, HeAu, AuAu) or config paths
//   source    — both | global | individual | hand  (default: both)
//   outputDir — PDF output directory              (default: output/params_compare)
//   charge    — pos | neg for global fit          (default: pos)
//   particle  — pi | K | p for individual-only    (default: pi)
//
// source=hand uses hand_params from JSON (reference init values, not fit output).

#include "../analysis/DrawParamsMulti.h"

#include <iostream>
#include <string>

namespace {

int ParseCharge(const std::string& text)
{
    const std::string s = bw::ToLowerTrim(text);
    if (s.empty() || s == "pos" || s == "positive" || s == "+")
    {
        return 0;
    }
    if (s == "neg" || s == "negative" || s == "-")
    {
        return 1;
    }
    throw std::runtime_error("Unknown charge: " + text);
}

int ParseParticle(const std::string& text)
{
    const std::string s = bw::ToLowerTrim(text);
    if (s.empty() || s == "pi" || s == "pip" || s == "pi+")
    {
        return 0;
    }
    if (s == "pi-" || s == "pim")
    {
        return 1;
    }
    if (s == "k" || s == "kp" || s == "k+")
    {
        return 2;
    }
    if (s == "k-" || s == "km")
    {
        return 3;
    }
    if (s == "p" || s == "proton")
    {
        return 4;
    }
    if (s == "pbar" || s == "ap" || s == "anti-p")
    {
        return 5;
    }
    throw std::runtime_error("Unknown particle: " + text);
}

}  // namespace

void draw_params_syst(const char* systems = "CuAu,UU",
                      const char* source = "both",
                      const char* outputDir = "output/params_compare",
                      const char* charge = "pos",
                      const char* particle = "pi")
{
    using namespace bw;

    try
    {
        const std::vector<std::string> tokens = SplitCommaList(systems);
        const ParamDataSource src = ParseParamDataSource(source);
        const int ch = ParseCharge(charge);
        const int part = ParseParticle(particle);

        ParamPlotOptions opt = ParseParamPlotOptions("both:pos");
        if (src == ParamDataSource::Global)
        {
            opt.sources = kParamSourceGlobal;
            opt.globalCharge = ch;
        }
        else if (src == ParamDataSource::Individual)
        {
            opt.sources = kParamSourceIndividual;
            opt.particleMask = (1u << part);
        }
        else if (src == ParamDataSource::Both)
        {
            opt.sources = kParamSourceIndividual | kParamSourceGlobal;
            opt.globalCharge = ch;
        }

        std::cout << "Systems: " << systems << "\n";
        std::cout << "Source: " << source << " (fit files in output/<system>/txtParams/)\n";
        std::cout << "Output: " << outputDir << "\n";
        if (src == ParamDataSource::Global || src == ParamDataSource::Both)
        {
            std::cout << "Charge: " << charge << "\n";
        }
        if (src == ParamDataSource::Individual)
        {
            std::cout << "Particle: " << particle << " (" << kParticleNames[part] << ")\n";
        }

        DrawMultiSystemParamPlots(tokens, src, outputDir, opt, ch, part, true);
    }
    catch (const std::exception& e)
    {
        std::cerr << "draw_params_syst failed: " << e.what() << "\n";
    }

    gROOT->ProcessLine(".q");
}
