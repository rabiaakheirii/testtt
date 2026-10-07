#pragma once
#include <vector>

namespace qcd {

    struct Cumulants {
        double C1;   // mean
        double C2;   // variance
        double C3;   // third central moment
        double C4;   // fourth cumulant
        double Ssigma;      // C3 / C2
        double kappa_sigma2; // C4 / C2
        int    n_events;
    };

    Cumulants compute_cumulants(const std::vector<int>& net_proton);
}
