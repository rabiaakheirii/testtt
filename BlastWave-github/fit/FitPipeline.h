#ifndef BLASTWAVE_FIT_PIPELINE_H
#define BLASTWAVE_FIT_PIPELINE_H

#include "GlobalFit.h"
#include "IndividualFit.h"
#include "../core/ParamIO.h"

#include <iostream>
#include <vector>

namespace bw {

enum class PipelineStep
{
    HandNoFit,            // step 1: hand params, curve only, no Minuit (case 3)
    HandFit,              // step 2: hand init + individual fit (case 2)
    Global,               // step 3: global fit pi+K+p
    IndividualFromGlobal, // step 4: individual from global params (case 0)
    GlobalRefine          // step 5: second global iteration
};

class FitPipeline
{
public:
    explicit FitPipeline(Dataset& data) : data_(data) {}

    void Run(const std::vector<PipelineStep>& steps)
    {
        for (auto step : steps)
        {
            RunStep(step);
        }
    }

    void RunFull()
    {
        Run({PipelineStep::HandFit,
             PipelineStep::Global,
             PipelineStep::IndividualFromGlobal,
             PipelineStep::GlobalRefine});
    }

    const FitParams (&IndividualResults() const)[kMaxParticles][kMaxCentralities]
    {
        return individual_;
    }

    const GlobalParams (&GlobalResults() const)[2][kMaxCentralities]
    {
        return global_;
    }

private:
    Dataset& data_;
    FitParams individual_[kMaxParticles][kMaxCentralities]{};
    GlobalParams global_[2][kMaxCentralities]{};

    void RunStep(PipelineStep step)
    {
        const auto& cfg = data_.Config();
        std::cout << "\n========== Pipeline step: " << int(step) << " ==========\n";

        switch (step)
        {
            case PipelineStep::HandNoFit:
            {
                std::cout << "handNoFit: set hand_params on BW functions (no Minuit)\n";
                FitStrategy strategy{InitStrategy::HandNoFit, LimitsStrategy::Fix};
                IndividualFitter fitter(data_, strategy);
                fitter.Run();
                fitter.CopyResults(individual_);
                break;
            }
            case PipelineStep::HandFit:
            {
                std::cout << "handFit: hand_params init + fit_limits scales\n";
                FitStrategy strategy{InitStrategy::HandFit, LimitsStrategy::Scaled};
                IndividualFitter fitter(data_, strategy);
                fitter.Run();
                fitter.CopyResults(individual_);
                WriteIndividualParams(cfg, individual_);
                break;
            }
            case PipelineStep::Global:
            {
                GlobalFitter fitter(data_, individual_);
                fitter.Run();
                fitter.CopyResults(global_);
                WriteGlobalParams(cfg, global_);
                break;
            }
            case PipelineStep::IndividualFromGlobal:
            {
                FitStrategy strategy{InitStrategy::Global, LimitsStrategy::Scaled};
                IndividualFitter fitter(data_, strategy);
                fitter.Run();
                fitter.CopyResults(individual_);
                WriteIndividualParams(cfg, individual_);
                break;
            }
            case PipelineStep::GlobalRefine:
            {
                GlobalFitter fitter(data_, individual_);
                fitter.Run();
                fitter.CopyResults(global_);
                WriteGlobalParams(cfg, global_);
                break;
            }
        }
    }
};

}  // namespace bw

#endif
