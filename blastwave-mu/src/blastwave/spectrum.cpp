#include "blastwave/spectrum.hpp"
#include "particles/chem_potential.hpp"
#include "utils/bessel.hpp"

#include <cmath>

namespace bwmu {

    // --------------------------------------------------------------
    // Boltzmann-limit shape (no chemical factor):
    //
    //   shape(mT, T, beta_s, R, n) = mT * integral_0^R r dr
    //                                  * I0(mT sinh(rho)/T)
    //                                  * K1(mT cosh(rho)/T)
    // --------------------------------------------------------------
    static double shape_boltzmann(double mT,
                                  double T,
                                  double beta_s,
                                  double R_fm,
                                  int    n_prof)
    {
        if (mT <= 0.0 || T <= 0.0 || R_fm <= 0.0) return 0.0;

        const int    N  = 200;
        const double dr = R_fm / N;
        double integral = 0.0;

        for (int i = 0; i <= N; ++i) {
            const double r    = i * dr;
            const double frac = r / R_fm;
            const double beta = beta_s * std::pow(frac, n_prof);
            if (beta >= 1.0) continue;

            const double rho   = std::atanh(beta);
            const double I0val = I0(mT * std::sinh(rho) / T);
            const double K1val = K1(mT * std::cosh(rho) / T);

            const double w = (i == 0 || i == N) ? 0.5 : 1.0;
            integral += w * r * I0val * K1val;
        }
        integral *= dr;

        return mT * integral;
    }

    // --------------------------------------------------------------
    // Boltzmann spectrum with chemical factor exp(mu_i / T)
    // --------------------------------------------------------------
    double dN_dmT_boltzmann(double mT,
                            const Species& sp,
                            const BWParams& bw,
                            double mu_B,
                            double mu_Q,
                            double mu_S,
                            double norm)
    {
        const double mu  = mu_of(sp, mu_B, mu_Q, mu_S);
        const double fac = std::exp(mu / bw.T);
        const double pre = sp.g / (2.0 * M_PI);

        return norm * pre * fac
             * shape_boltzmann(mT, bw.T, bw.beta_s, bw.R_fm, bw.n_prof);
    }

    // --------------------------------------------------------------
    // Quantum-statistics spectrum via fugacity expansion.
    //
    //   Bose   : f = sum_{n=1}^{inf} (+1)          z^n e^{-nE/T}
    //   Fermi  : f = sum_{n=1}^{inf} (-1)^{n-1}    z^n e^{-nE/T}
    //
    // with z = exp(mu/T). Each term is the Boltzmann formula evaluated
    // at T -> T/n, multiplied by z^n and the appropriate sign.
    // --------------------------------------------------------------
    double dN_dmT_quantum(double mT,
                          const Species& sp,
                          const BWParams& bw,
                          double mu_B,
                          double mu_Q,
                          double mu_S,
                          double norm,
                          int    max_terms)
    {
        const double mu  = mu_of(sp, mu_B, mu_Q, mu_S);
        const double pre = sp.g / (2.0 * M_PI);

        double total = 0.0;

        for (int n = 1; n <= max_terms; ++n) {
            // sign convention:
            //   bosons  : +1 for all n
            //   fermions: (-1)^(n-1)  ->  +1, -1, +1, -1, ...
            double sign;
            if (sp.boson) {
                sign = +1.0;
            } else {
                sign = (n % 2 == 1) ? +1.0 : -1.0;
            }

            const double Tn   = bw.T / n;
            const double zpow = std::exp(n * mu / bw.T);

            const double term = sign * zpow
                              * shape_boltzmann(mT, Tn,
                                                bw.beta_s, bw.R_fm,
                                                bw.n_prof);

            total += term;

            // early exit when the term is negligible compared to the sum
            if (n > 1 && std::abs(term) < 1e-14 * std::abs(total)) break;
        }

        return norm * pre * total;
    }

} // namespace bwmu
