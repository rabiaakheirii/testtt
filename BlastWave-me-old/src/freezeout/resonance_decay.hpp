#pragma once
#include <string>
#include <vector>

namespace qcd {

    struct Resonance {
        std::string name;
        double mass;      // GeV
        double width;     // GeV
        double pstar;     // pion momentum in resonance rest frame [GeV]
        double branching; // to pions
    };

    std::vector<Resonance> load_resonances(const std::string& filename);

    // Two-body decay contribution to the pion mT spectrum
    //  dn/(mT dmT)  at given mT, for a resonance with mass mR,
    //  daughter pion momentum pstar, temperature T.
    // Uses a simplified analytic approximation (paper Eq. 25-28).
    double resonance_decay_mT(double mT,
                              double mR,
                              double pstar,
                              double T);

    // Total pion mT spectrum = direct + all resonances
    double total_pion_mT(double mT,
                         double m_pi,
                         double T,
                         double norm,
                         const std::vector<Resonance>& resonances);
}