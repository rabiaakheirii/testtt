#include "freezeout/blastwave.hpp"
#include "utils/bessel.hpp"
#include <cmath>

namespace qcd {

    double thermal_mT_spectrum(double mT, double T, double norm) {
        if (mT <= 0.0 || T <= 0.0) return 0.0;
        return norm * mT * K1(mT / T);
    }

    double blastwave_mT_spectrum(double mT, double T,
                                 double beta_s, double R, int n,
                                 double norm)
    {
        if (mT <= 0.0 || T <= 0.0 || R <= 0.0) return 0.0;

        const int    Nr = 200;
        const double dr = R / Nr;
        double integral = 0.0;

        for (int i = 0; i <= Nr; ++i) {
            const double r     = i * dr;
            const double frac  = (R > 0.0) ? r / R : 0.0;
            const double beta  = beta_s * std::pow(frac, n);
            if (beta >= 1.0) continue;
            const double rho   = std::atanh(beta);

            const double I0val = I0(mT * std::sinh(rho) / T);
            const double K1v = K1(mT * std::cosh(rho) / T);

            const double w  = (i == 0 || i == Nr) ? 0.5 : 1.0;
            integral += w * r * I0val * K1v;

        }
        integral *= dr;

        return norm * mT * integral;
    }

    double effective_temperature(double T, double beta) {
        return T * std::sqrt((1.0 + beta) / (1.0 - beta));
    }
}