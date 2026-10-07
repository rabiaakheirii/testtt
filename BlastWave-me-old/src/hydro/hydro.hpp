#pragma once
#include "hydro/eos.hpp"
#include <string>
#include <vector>

namespace qcd {

    struct HydroProfile {
        int    n_tau, n_r;
        double tau_min, tau_max;
        double r_min, r_max;
        std::vector<double> tau, r;
        std::vector<double> eps;   // [n_tau * n_r]
        std::vector<double> T;     // [n_tau * n_r]
        std::vector<double> ur;    // u^r
    };

    // 1+1D Bjorken hydro with ideal EoS
    //   tau step size dtau, r step dr
    //   initial energy density: eps0 * exp(-r^2 / (2 R^2))
    HydroProfile run_hydro_1p1d(const EoS& eos,
                                double eps0,
                                double R_fm,
                                double tau0, double tau_end, double dtau,
                                double r_max, double dr);

    void write_hydro_profile(const HydroProfile& h,
                             const std::string& filename);
}
