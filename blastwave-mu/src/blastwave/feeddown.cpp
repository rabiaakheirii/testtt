#include "blastwave/feeddown.hpp"
#include "blastwave/spectrum.hpp"
#include "particles/species.hpp"
#include <cmath>
#include <fstream>
#include <sstream>

namespace bwmu {

    namespace {
        std::vector<Resonance>& storage() {
            static std::vector<Resonance> r;
            return r;
        }
    }

    void set_resonances(const std::vector<Resonance>& r) {
        storage() = r;
    }
    const std::vector<Resonance>& get_resonances() {
        return storage();
    }

    std::vector<Resonance> load_resonances(const std::string& path) {
        std::vector<Resonance> out;
        std::ifstream in(path);
        if (!in) return out;
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            Resonance r;
            if (ss >> r.name >> r.mass >> r.g >> r.B >> r.Q >> r.S
                   >> r.pstar >> r.BR >> r.target) {
                out.push_back(r);
            }
        }
        return out;
    }

    double decay_contribution(double mT,
                              const std::string& target,
                              const BWParams& bw,
                              const FitParams& fp,
                              Statistics /*stats*/)
    {
        double total = 0.0;
        for (const auto& R : get_resonances()) {
            if (R.target != target) continue;

            // Parent chemical potential
            const double mu_R = R.B * fp.mu_B
                              + R.Q * fp.mu_Q
                              + R.S * fp.mu_S;

            // Two-temperature scheme:
            //   T_ch  = chemical freeze-out, sets the ABUNDANCE of the parent
            //   fp.T  = kinetic freeze-out, sets the m_T SHAPE
            constexpr double T_ch = 0.156;   // GeV, from lattice QCD
            constexpr double m_pi = 0.139570;

            // Boltzmann abundance ratio N_R / N_pi(direct):
            //   (m_R / m_pi)^(3/2) * exp(-(m_R - m_pi) / T_ch)
            //   * exp(mu_R / T_ch)
            // The spin degeneracy g_R is set on the parent Species below.
            // Boltzmann suppression factor, guarded against overflow
            const double exponent = -(R.mass - m_pi) / T_ch;
            const double suppression = (exponent > -50.0)
                                     ? std::exp(exponent)
                                     : 0.0;
            const double abundance =
                std::pow(R.mass / m_pi, 1.5)
              * suppression
              * std::exp(mu_R / T_ch);

            // Effective temperature for the decay pions (paper Eq. 28)
            const double T_eff = (R.pstar / R.mass) * fp.T;

            // Species object for the resonance; g = g_R
            Species parent;
            parent.name  = R.name;
            parent.mass  = R.mass;
            parent.g     = R.g;
            parent.B     = R.B;
            parent.Q     = R.Q;
            parent.S     = R.S;
            parent.boson = true;

            BWParams bw_decay = bw;
            bw_decay.T      = T_eff;
            bw_decay.beta_s = fp.beta_s;

            // BW shape of the parent (norm = 1; g_R enters via parent.g)
            const double shape = dN_dmT_boltzmann(mT, parent, bw_decay,
                                                   0.0, 0.0, 0.0, 1.0);

            total += R.BR * abundance * shape;
        }
        return total;
    }

} // namespace bwmu
