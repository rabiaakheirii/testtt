#ifndef BLASTWAVE_GLOBAL_FIT_H
#define BLASTWAVE_GLOBAL_FIT_H

#include "../core/BlastWaveModel.h"
#include "../core/Constants.h"
#include "../core/FitConfig.h"
#include "../core/ParamIO.h"
#include "../core/SpectrumData.h"

#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "Fit/Fitter.h"
#include "Math/WrappedMultiTF1.h"
#include "TF1.h"
#include "TVirtualFitter.h"

#include <iostream>
#include <sstream>

namespace bw {

struct GlobalChi2Fcn
{
    GlobalChi2Fcn(ROOT::Math::IMultiGenFunction& f0,
                  ROOT::Math::IMultiGenFunction& f2,
                  ROOT::Math::IMultiGenFunction& f4,
                  const TBetaPenalty& penalty)
        : fChi2_0(&f0), fChi2_2(&f2), fChi2_4(&f4), penalty_(penalty)
    {}

    double operator()(const double* par) const
    {
        constexpr int kNFunc = 3;
        constexpr int kNPar = 4;
        double p[kNFunc][kNPar];

        for (int i = 0; i < kNFunc; ++i)
        {
            p[i][0] = par[2 + i];
            p[i][1] = par[0];
            p[i][2] = par[1];
            p[i][3] = kMasses[2 * i];
        }

        double chi2 = (*fChi2_0)(p[0]) + (*fChi2_2)(p[1]) + (*fChi2_4)(p[2]);

        if (penalty_.enabled)
        {
            const double T = par[0];
            const double beta = par[1];
            const double T_pred = PredictT(beta, penalty_.a0, penalty_.a1, penalty_.a2);
            const double diff = T - T_pred;
            chi2 += (diff * diff) / (penalty_.sigma * penalty_.sigma);
        }
        return chi2;
    }

    const ROOT::Math::IMultiGenFunction* fChi2_0;
    const ROOT::Math::IMultiGenFunction* fChi2_2;
    const ROOT::Math::IMultiGenFunction* fChi2_4;
    TBetaPenalty penalty_;
};

class GlobalFitter
{
public:
    explicit GlobalFitter(const Dataset& data,
                          const FitParams individual[kMaxParticles][kMaxCentralities] = nullptr)
        : data_(data), individual_(individual)
    {
        TVirtualFitter::SetDefaultFitter("Minuit");
        integrand_ = new TF1("bw_integrand_global", bwfitfunc, 0.01, 10., 5);
        integ_ = new MyIntegFunc(integrand_);
    }

    ~GlobalFitter()
    {
        for (int p = 0; p < kMaxParticles; ++p)
        {
            for (int c = 0; c < kMaxCentralities; ++c)
            {
                delete funcs_[p][c];
            }
        }
    }

    const GlobalParams& Result(int charge, int centrality) const
    {
        return results_[charge][centrality];
    }

    TF1* Function(int particle, int centrality) const
    {
        return funcs_[particle][centrality];
    }

    void Run()
    {
        const auto& cfg = data_.Config();
        for (int c = 0; c < cfg.NumCentralities(); ++c)
        {
            BuildFunctions(c);
            FitCentrality(c, 0);
            FitCentrality(c, 1);
        }
    }

    void InitFunctions()
    {
        const auto& cfg = data_.Config();
        for (int c = 0; c < cfg.NumCentralities(); ++c)
        {
            BuildFunctions(c);
        }
    }

    void CopyResults(GlobalParams out[2][kMaxCentralities]) const
    {
        for (int ch = 0; ch < 2; ++ch)
        {
            for (int c = 0; c < kMaxCentralities; ++c)
            {
                out[ch][c] = results_[ch][c];
            }
        }
    }

private:
    const Dataset& data_;
    const FitParams (*individual_)[kMaxCentralities] = nullptr;
    TF1* integrand_ = nullptr;
    MyIntegFunc* integ_ = nullptr;
    TF1* funcs_[kMaxParticles][kMaxCentralities]{};
    GlobalParams results_[2][kMaxCentralities]{};

    void BuildFunctions(int centrality)
    {
        for (int part = 0; part < kMaxParticles; ++part)
        {
            const auto fitRange = data_.Config().FitRangeForParticle(part);
            delete funcs_[part][centrality];
            std::ostringstream name;
            name << "bw_global_" << part << "_" << centrality;
            funcs_[part][centrality] =
                new TF1(name.str().c_str(), *integ_, fitRange.xmin, fitRange.xmax, 4);
            funcs_[part][centrality]->FixParameter(3, kMasses[part]);
        }
    }

    static bool HasIndividualInit(const FitParams& p)
    {
        return p.T > 0.01 && p.beta > 0.01;
    }

    void FitCentrality(int centrality, int charge)
    {
        const auto& cfg = data_.Config();

        const int pi = 0 + charge;
        const int k = 2 + charge;
        const int p = 4 + charge;

        ROOT::Math::WrappedMultiTF1 wf0(*funcs_[pi][centrality], 1);
        ROOT::Math::WrappedMultiTF1 wf2(*funcs_[k][centrality], 1);
        ROOT::Math::WrappedMultiTF1 wf4(*funcs_[p][centrality], 1);

        ROOT::Fit::DataOptions opt;
        ROOT::Fit::DataRange rangePi, rangeK, rangeP;

        const auto rPi = cfg.fit_range_pi;
        const auto rK = cfg.fit_range_K;
        const auto rP = cfg.fit_range_p;

        rangePi.SetRange(rPi.xmin, rPi.xmax);
        rangeK.SetRange(rK.xmin, rK.xmax);
        rangeP.SetRange(rP.xmin, rP.xmax);

        ROOT::Fit::BinData dataPi(opt, rangePi);
        ROOT::Fit::BinData dataK(opt, rangeK);
        ROOT::Fit::BinData dataP(opt, rangeP);
        ROOT::Fit::FillData(dataPi, data_.Spectrum(pi, centrality));
        ROOT::Fit::FillData(dataK, data_.Spectrum(k, centrality));
        ROOT::Fit::FillData(dataP, data_.Spectrum(p, centrality));

        ROOT::Fit::Chi2Function chi2_pi(dataPi, wf0);
        ROOT::Fit::Chi2Function chi2_k(dataK, wf2);
        ROOT::Fit::Chi2Function chi2_p(dataP, wf4);
        GlobalChi2Fcn globalChi2(chi2_pi, chi2_k, chi2_p, cfg.t_beta_penalty);

        double T0 = 0.108;
        double beta0 = 0.7;
        double Cpi = 10., Ck = 1., Cp = 0.1;
        if (centrality < cfg.hand.T.size())
        {
            T0 = cfg.hand.T[centrality];
        }
        if (centrality < cfg.hand.beta.size())
        {
            beta0 = cfg.hand.beta[centrality];
        }
        if (pi < cfg.hand.C.size() &&
            centrality < cfg.hand.C[pi].size())
        {
            Cpi = cfg.hand.C[pi][centrality];
        }
        if (k < cfg.hand.C.size() &&
            centrality < cfg.hand.C[k].size())
        {
            Ck = cfg.hand.C[k][centrality];
        }
        if (p < cfg.hand.C.size() &&
            centrality < cfg.hand.C[p].size())
        {
            Cp = cfg.hand.C[p][centrality];
        }

        if (individual_)
        {
            const FitParams& ip = individual_[pi][centrality];
            if (HasIndividualInit(ip))
            {
                T0 = ip.T;
                beta0 = ip.beta;
            }
            if (individual_[pi][centrality].C > 0.)
            {
                Cpi = individual_[pi][centrality].C;
            }
            if (individual_[k][centrality].C > 0.)
            {
                Ck = individual_[k][centrality].C;
            }
            if (individual_[p][centrality].C > 0.)
            {
                Cp = individual_[p][centrality].C;
            }
        }

        constexpr int kNpar = 5;
        double par0[kNpar] = {T0, beta0, Cpi, Ck, Cp};

        ROOT::Fit::Fitter fitter;
        fitter.Config().SetParamsSettings(kNpar, par0);
        fitter.Config().ParSettings(0).SetLimits(T0 * 0.9, T0 * 1.1);
        fitter.Config().ParSettings(1).SetLimits(beta0 * 0.9, beta0 * 1.1);
        fitter.Config().ParSettings(2).SetLimits(Cpi * 0.001, Cpi * 100.9);
        fitter.Config().ParSettings(3).SetLimits(Ck * 0.001, Ck * 100.9);
        fitter.Config().ParSettings(4).SetLimits(Cp * 0.001, Cp * 100.9);
        fitter.Config().MinimizerOptions().SetPrintLevel(0);
        fitter.Config().SetMinimizer("Minuit2", "Migrad");

        fitter.FitFCN(kNpar, globalChi2, nullptr,
                      dataPi.Size() + dataK.Size() + dataP.Size(), true);

        const ROOT::Fit::FitResult& result = fitter.Result();
        if (result.Status() != 0)
        {
            std::cout << "WARNING: global fit status " << result.Status()
                      << " centr=" << centrality << " charge=" << charge << "\n";
        }

        const double* fitPar = result.GetParams();
        results_[charge][centrality].T = fitPar[0];
        results_[charge][centrality].beta = fitPar[1];
        results_[charge][centrality].C_pi = fitPar[2];
        results_[charge][centrality].C_K = fitPar[3];
        results_[charge][centrality].C_p = fitPar[4];

        std::cout << "Global fit centr=" << centrality << " charge=" << charge
                  << "  T=" << fitPar[0] << " beta=" << fitPar[1]
                  << "  chi2/ndf=" << result.Chi2() / result.Ndf() << "\n";
    }
};

}  // namespace bw

#endif
