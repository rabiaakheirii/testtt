#include "freezeout/resonance_decay.hpp"
#include "freezeout/blastwave.hpp"
#include "utils/bessel.hpp"

#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>

namespace qcd {

    std::vector<Resonance> load_resonances(const std::string& filename) {
        std::vector<Resonance> out;
        std::ifstream in(filename);
        if (!in) {
            std::cerr << "[resonance] cannot open " << filename << "\n";
            return out;
        }
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            Resonance r;
            if (ss >> r.name >> r.mass >> r.width >> r.pstar >> r.branching)
                out.push_back(r);
        }
        std::cout << "[resonance] loaded " << out.size()
                  << " resonances from " << filename << "\n";
        return out;
    }

    // Simplified analytic two-body decay spectrum:
    //   dn/(mT dmT) ~ mT * exp(-mT / T_eff)   (paper Eq. 27-28)
    // with  T_eff = (p*/mR) * T.
    //
    // This captures the key feature: decay pions have a steeper
    // slope than direct pions, filling the low-mT region.
    double resonance_decay_mT(double mT,
                              double mR,
                              double pstar,
                              double T) {
        if (mT <= 0.0 || T <= 0.0 || mR <= 0.0) return 0.0;
        const double T_eff = (pstar / mR) * T;
        return mT * std::exp(-mT / T_eff);
    }

    double total_pion_mT(double mT,
                         double m_pi,
                         double T,
                         double norm,
                         const std::vector<Resonance>& resonances) {
        // direct pions
        double total = thermal_mT_spectrum(mT, T, norm);

        // resonance decays
        for (const auto& r : resonances) {
            total += norm * r.branching
                   * resonance_decay_mT(mT, r.mass, r.pstar, T);
        }
        return total;
    }
}