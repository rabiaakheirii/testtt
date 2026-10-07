#pragma once
#include <string>
#include <vector>

namespace qcd {

    struct Particle {
        std::string name;
        double mass;      // GeV
        double degeneracy;
        int baryon;
        int strangeness;
    };

    struct Config {
        // freeze-out
        double T_fo = 0.150;          // GeV
        double beta_s = 0.5;          // surface velocity
        double eta_max = 1.7;         // longitudinal flow
        int    n_profile = 2;         // beta_r(r) = beta_s * (r/R)^n
        double R = 5.0;               // fm (transverse radius)

        // percolation
        double string_density = 1.0;
        double percolation_threshold = 0.5;

        // hydro
        double eta_over_s = 0.08;
        double tau_0 = 1.0;           // fm/c
    };

}