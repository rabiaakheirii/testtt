#include "fit/tbeta_penalty.hpp"
#include <algorithm>

namespace bwmu {

    // Published T(beta_s) at sqrt(s_NN) = 9.2 GeV (NICA/MPD)
    //
    // From STAR BES PRC 96, 044904 (2017):
    //   beta_s = 0.55 -> T = 130 MeV
    //   beta_s = 0.60 -> T = 120 MeV
    //   beta_s = 0.65 -> T = 110 MeV
    //
    // Slope = -200 MeV per unit beta_s
    // Reference: T = 130 MeV at beta_s = 0.55
    double T_exp_of_beta(double beta_s) {
        constexpr double T_ref    = 0.130;   // 130 MeV at NICA
        constexpr double beta_ref = 0.55;
        constexpr double slope    = -0.200;

        double T = T_ref + slope * (beta_s - beta_ref);
        return std::max(T, 0.085);
    }

    double tbeta_penalty(double T, double beta_s, double lambda) {
        const double d = T - T_exp_of_beta(beta_s);
        return lambda * d * d;
    }

} // namespace bwmu
