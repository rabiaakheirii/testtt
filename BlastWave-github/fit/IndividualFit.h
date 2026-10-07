#ifndef BLASTWAVE_INDIVIDUAL_FIT_H
#define BLASTWAVE_INDIVIDUAL_FIT_H

#include "../core/BlastWaveModel.h"
#include "../core/Constants.h"
#include "../core/FitConfig.h"
#include "../core/ParamIO.h"
#include "../core/SpectrumData.h"

#include "TF1.h"
#include "TGraphErrors.h"
#include "TVirtualFitter.h"

#include <sstream>
#include <string>

namespace bw {

enum class InitStrategy
{
    HandFit,  
    Global, 
    HandNoFit
};

enum class LimitsStrategy
{
    Fix,  
    Scaled // low * initParam, high * initParam
};

struct FitStrategy
{
    InitStrategy init;
    LimitsStrategy limits;
};

class IndividualFitter
{
public:
    IndividualFitter(const Dataset& data, FitStrategy fitStrategy = {})
        : data_(data), fitStrategy_(fitStrategy)
    {
        TVirtualFitter::SetDefaultFitter("Minuit");
        integrand_ = new TF1("bw_integrand", bwfitfunc, 0.01, 10., 5);
        integrand_->SetParNames("constant", "T", "beta", "mass", "pt");
        integ_ = new MyIntegFunc(integrand_);
    }

    ~IndividualFitter()
    {
        for (int p = 0; p < kMaxParticles; ++p)
            for (int c = 0; c < kMaxCentralities; ++c)
                delete fitfuncs_[p][c];
    }

    TF1* Function(int particle, int centrality) const
    {
        return fitfuncs_[particle][centrality];
    }

    void Run()
    {
        const int nCent = data_.Config().NumCentralities();
        for (int part = 0; part < kMaxParticles; ++part)
        {
            for (int c = 0; c < nCent; ++c)
            {
                if (!data_.HasSpectrum(part, c))
                    continue;

                FitOne(part, c);
            }
        }
    }

    void CopyResults(FitParams out[kMaxParticles][kMaxCentralities]) const
    {
        for (int p = 0; p < kMaxParticles; ++p)
        {
            for (int c = 0; c < kMaxCentralities; ++c)
            {
                out[p][c] = results_[p][c];
            }
        }
    }

private:
    const Dataset& data_;
    FitStrategy fitStrategy_;
    TF1* integrand_ = nullptr;
    MyIntegFunc* integ_ = nullptr;
    TF1* fitfuncs_[kMaxParticles][kMaxCentralities]{};
    FitParams results_[kMaxParticles][kMaxCentralities]{};
    GlobalParams global_[2][kMaxCentralities]{};

    TF1* MakeFitFunction(int particle, int centrality)
    {
        const auto fitRange = data_.Config().FitRangeForParticle(particle);
        std::ostringstream name;
        name << "bw_" << particle << "_" << centrality;
        return new TF1(name.str().c_str(), *integ_, fitRange.xmin, fitRange.xmax, 4);
    }


    FitParams HandParams(int particle, int centrality) const
    {
        const auto& cfg = data_.Config();

        FitParams p;
        p.mass = kMasses[particle];
        p.T = cfg.hand.T[centrality];
        p.beta = cfg.hand.beta[centrality];
        p.C = cfg.hand.C[particle][centrality];
        
        return p;
    }


    void FitOne(int particle, int centrality)
    {
        const auto& cfg = data_.Config();
        const auto fitRange = cfg.FitRangeForParticle(particle);
        auto* gr = data_.Spectrum(particle, centrality);
        
        if (!gr) return;
        
        delete fitfuncs_[particle][centrality];
        fitfuncs_[particle][centrality] = MakeFitFunction(particle, centrality);
        TF1* f = fitfuncs_[particle][centrality];

        // Init params
        FitParams init;
        switch (fitStrategy_.init)
        {
            case InitStrategy::Global:
                ReadGlobalParams(data_.Config(), global_);
                init = ParamsFromGlobal(particle, centrality, global_);
                break;

            default: // hand, noHand 
                init = HandParams(particle, centrality);
        }

        f->SetParameters(init.C, init.T, init.beta, init.mass);
        f->FixParameter(3, init.mass);
        
        // Apply limits
        const FitLimits& lim = data_.Config().fit_limits;
        switch (fitStrategy_.limits)
        {
            case LimitsStrategy::Scaled:
                f->SetParLimits(0, 0., init.C * lim.C_max_scale);
                f->SetParLimits(1, init.T * lim.T_min_scale, init.T * lim.T_max_scale);
                f->SetParLimits(2, init.beta * lim.beta_min_scale, init.beta * lim.beta_max_scale);
                break;

            case LimitsStrategy::Fix:
                f->SetParLimits(0, 0., init.C * lim.C_max_scale);
                f->SetParLimits(1, lim.T_min, lim.T_max);
                f->SetParLimits(2, lim.beta_min, lim.beta_max);
                break;
        }

        const bool runMinuit = fitStrategy_.init != InitStrategy::HandNoFit;

        if (runMinuit)
        {
            gr->Fit(f, "QR+SEX0", "", fitRange.xmin, fitRange.xmax);

            results_[particle][centrality].C = f->GetParameter(0);
            results_[particle][centrality].T = f->GetParameter(1);
            results_[particle][centrality].beta = f->GetParameter(2);
            results_[particle][centrality].mass = init.mass;
            results_[particle][centrality].C_err = f->GetParError(0);
            results_[particle][centrality].T_err = f->GetParError(1);
            results_[particle][centrality].beta_err = f->GetParError(2);
        }
        else
        {
            results_[particle][centrality] = init;
            results_[particle][centrality].C_err = 0.;
            results_[particle][centrality].T_err = 0.;
            results_[particle][centrality].beta_err = 0.;
        }
    }
};

}  // namespace bw

#endif
