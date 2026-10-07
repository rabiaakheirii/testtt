// BlastWave fit entry point for PHENIX.
//
// Usage (from repository root):
//   root -l 'macros/run_fit.C("config/phenix_uu.json", "handNoFit")'
//   root -l 'macros/run_fit.C("config/phenix_uu.json", "handFit")'
//   root -l 'macros/run_fit.C("config/phenix_uu.json", "full")'
//   root -l 'macros/run_fit.C("config/phenix_uu.json", "individualGlobal")'
//   root -l 'macros/run_fit.C("config/phenix_uu.json", "draw")'
//   root -l 'macros/run_fit.C("config/phenix_uu.json", "drawParams")'
//   root -l 'macros/run_fit.C("config/phenix_uu.json", "drawParams:individual:pos")'
//   root -l 'macros/draw_params_syst.C("CuAu,UU", "global")'

#include "../analysis/DrawParams.h"
#include "../analysis/DrawSpectra.h"
#include "../core/FitConfig.h"
#include "../core/SpectrumData.h"
#include "../fit/FitPipeline.h"

#include <iostream>
#include <string>

void run_fit(const char* configPath = "config/phenix_uu.json",
             const char* mode = "load")
{
    using namespace bw;

    const FitConfig cfg = FitConfig::Load(configPath);
    std::cout << "Experiment: " << cfg.experiment << "  system: " << cfg.system
              << "  sqrt(s_NN)=" << cfg.sqrt_s_NN << " GeV\n";
    std::cout << "Spectra: " << cfg.spectra_path << "  variable=" << cfg.variable << "\n";
    std::cout << "Centralities: " << cfg.NumCentralities() << "\n";
    std::cout << "T-beta penalty: " << (cfg.t_beta_penalty.enabled ? "on" : "off") << "\n";
    std::cout << "PDF output dir: " << cfg.output_dir << "\n";

    Dataset data(cfg);

    for (int part = 0; part < kMaxParticles; ++part)
    {
        int nLoaded = 0;
        for (int c = 0; c < cfg.NumCentralities(); ++c)
        {
            if (data.HasSpectrum(part, c))
            {
                ++nLoaded;
            }
        }
        std::cout << "  particle " << kParticleNames[part] << ": " << nLoaded
                  << " spectra loaded\n";
    }

    const std::string m(mode);
    FitPipeline pipeline(data);

    if (m == "handNoFit")
    {
        pipeline.Run({PipelineStep::HandNoFit});
        DrawHandNoFit(cfg, data, "BlastWave_handNoFit.pdf");
    }
    else if (m == "handFit")
    {
        pipeline.Run({PipelineStep::HandFit});
        DrawIndividualResults(cfg, data, "BlastWave_handFit.pdf");
    }
    else if (m == "global")
    {
        pipeline.Run({PipelineStep::HandFit});
        pipeline.Run({PipelineStep::Global});
        DrawGlobalResults(cfg, data, "BlastWave_global.pdf");
        DrawParamPlotsFromFiles(cfg, ParseParamPlotOptions("both:pos"));
    }
    else if (m == "individualGlobal")
    {
        pipeline.Run({PipelineStep::Global, PipelineStep::IndividualFromGlobal});
        DrawIndividualResults(cfg, data, "BlastWave_individualGlobal.pdf");
        DrawParamPlotsFromFiles(cfg, ParseParamPlotOptions("both:pos"));
    }
    else if (m == "full")
    {
        pipeline.RunFull();
        DrawIndividualResults(cfg, data, "BlastWave.pdf");
        DrawGlobalResults(cfg, data, "BlastWave_global_refine.pdf");
        DrawParamPlotsFromFiles(cfg, ParseParamPlotOptions("both:pos"));
    }
    else if (m == "draw")
    {
        DrawIndividualResults(cfg, data);
    }
    else if (m == "drawGlobal")
    {
        DrawGlobalResults(cfg, data);
    }
    else if (m.rfind("drawParams", 0) == 0)
    {
        std::string subspec;
        if (m.size() > 10 && m[10] == ':')
        {
            subspec = m.substr(11);
        }
        DrawParamPlotsFromFiles(cfg, ParseParamPlotOptions(subspec));
    }
    else if (m != "load")
    {
        std::cout << "Unknown mode '" << mode << "'. "
                  << "Available: load, handNoFit, handFit, global, individualGlobal, full, "
                  << "draw, drawGlobal, drawParams[:individual|:global|:both][:pos|:neg|:pi|:k|:p][:nopredict]\n";
    }

    gROOT->ProcessLine(".q");
}
