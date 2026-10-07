#pragma once
#include "blastwave/params.hpp"
#include "particles/species.hpp"

namespace bwmu {

    // ------------------------------------------------------------------
    // Boltzmann-limit blast-wave spectrum
    //
    //   dN/(mT dmT) = (g/(2*pi)) * mT * exp(mu_i/T)
    //                 * integral_0^R r dr
    //                     I0( mT sinh(rho) / T )
    //                     K1( mT cosh(rho) / T )
    //
    // with rho(r) = atanh( beta_s * (r/R)^n )
    //
    // Returns dN/(mT dmT) in GeV^-2. The overall normalization (tau_f,
    // prefactors) is absorbed into the caller.
    // ------------------------------------------------------------------
    double dN_dmT_boltzmann(
        double mT,
        const Species& sp,
        const BWParams& bw,
        double mu_B,
        double mu_Q,
        double mu_S,
        double norm = 1.0);

    // ------------------------------------------------------------------
    // Quantum-statistics blast-wave spectrum via fugacity expansion.
    //
    //   f = 1 / (exp((u.p - mu)/T) -/+ 1)
    //
    // expands into a series over the fugacity z = exp(mu/T):
    //
    //   f = sum_{n>=1} (-/+ 1)^(n-1) z^n exp(-n u.p / T)
    //
    // so the spectrum becomes
    //
    //   dN/(mT dmT) = sum_{n=1}^{Nmax} (-/+ 1)^(n-1) z^n
    //                 * [ Boltzmann formula at T -> T/n ]
    //
    // The sign is '-' for bosons, '+' for fermions.
    // ------------------------------------------------------------------
    double dN_dmT_quantum(
        double mT,
        const Species& sp,
        const BWParams& bw,
        double mu_B,
        double mu_Q,
        double mu_S,
        double norm = 1.0,
        int    max_terms = 4);

} // namespace bwmu
