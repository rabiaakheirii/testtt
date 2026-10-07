#ifndef BLASTWAVE_MODEL_H
#define BLASTWAVE_MODEL_H

#include "TF1.h"
#include "TMath.h"

// Blast-wave integrand (E.J.Kim, flow_c12.C style).
// par[0]=C, par[1]=T, par[2]=beta, par[3]=mass, par[4]=mT (integration variable for MyIntegFunc).
inline double bwfitfunc(double* x, double* par)
{
    const double con = par[0];
    const double mass = par[3];
    const double mt = par[4] + mass;
    const double pt = std::sqrt(mt * mt - mass * mass);
    const double rho = TMath::ATanH(par[2]);

    return con * x[0] * mt
           * TMath::BesselI0(pt * TMath::SinH(rho) / par[1])
           * TMath::BesselK1(mt * TMath::CosH(rho) / par[1]);
}

struct MyIntegFunc
{
    explicit MyIntegFunc(TF1* f) : fFunc(f) {}

    TF1* fFunc;
    double param[5]{};

    double operator()(double* x, double* p)
    {
        constexpr double kRadius = 13.0;  // Rmax [fm]
        std::copy(p, p + 4, param);
        param[4] = *x;
        fFunc->SetParameters(param);
        return fFunc->Integral(0.0001, kRadius, 1.e-10);
    }
};

// T(beta) from PHENIX global-fit systematics (quadratic).
inline double PredictT(double beta, double a0, double a1, double a2)
{
    return a0 + beta * a1 + beta * beta * a2;
}

inline double GetT(double beta, double a0 = 0.129557, double a1 = 0.232774, double a2 = -0.339709)
{
    TF1 f("t_beta", "[0] + x * [1] + x * x * [2]", 0.25, 0.9);
    f.SetParameters(a0, a1, a2);
    return f.Eval(beta) - 0.01;
}

#endif
