#include "particles/table.hpp"
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace bwmu {

    const std::vector<Species>& particle_table() {
        static const std::vector<Species> table = {
            //  name        mass       g    B   Q   S   boson
            { "pi+",       0.139570,  1,   0, +1,  0, true  },
            { "pi-",       0.139570,  1,   0, -1,  0, true  },
            { "K+",        0.493677,  1,   0, +1, +1, true  },
            { "K-",        0.493677,  1,   0, -1, -1, true  },
            { "p",         0.938272,  2,  +1, +1,  0, false },
            { "pbar",      0.938272,  2,  -1, -1,  0, false },
            { "Lambda",    1.115683,  2,  +1,  0, -1, false },
            { "Lambdabar", 1.115683,  2,  -1,  0, +1, false },
            { "Xi-",       1.321710,  2,  +1, -1, -2, false },
            { "Xibar+",    1.321710,  2,  -1, +1, +2, false },
            { "Omega-",    1.672450,  4,  +1, -1, -3, false },
            { "Omegabar+", 1.672450,  4,  -1, +1, +3, false },
            { "phi",       1.019461,  3,   0,  0,  0, true  },
            { "Jpsi",      3.096900,  3,   0,  0,  0, true  },
            // Light nuclei (for NICA/MPD at 9.2 GeV)
            { "d",         1.875612,  3,  +1, +1,  0, true  },   // deuteron
            { "t",         2.808921,  2,  +1, +1,  0, true  },   // triton
            { "He3",       2.808921,  2,  +1, +2,  0, true  },   // helium-3
            { "He4",       3.727379,  1,  +2, +2,  0, true  }    // helium-4 (alpha)
        };
        return table;
    }

    const Species& get_species(const std::string& name) {
        for (const auto& s : particle_table()) {
            if (s.name == name) return s;
        }
        throw std::runtime_error("get_species: unknown species '" + name + "'");
    }

    void print_particle_table() {
        std::cout << std::left
                  << std::setw(12) << "name"
                  << std::setw(10) << "mass[GeV]"
                  << std::setw(4)  << "g"
                  << std::setw(4)  << "B"
                  << std::setw(4)  << "Q"
                  << std::setw(4)  << "S"
                  << "stat\n";
        std::cout << std::string(45, '-') << "\n";
        for (const auto& s : particle_table()) {
            std::cout << std::left
                      << std::setw(12) << s.name
                      << std::setw(10) << std::fixed << std::setprecision(4) << s.mass
                      << std::setw(4)  << s.g
                      << std::setw(4)  << s.B
                      << std::setw(4)  << s.Q
                      << std::setw(4)  << s.S
                      << (s.boson ? "Bose" : "Fermi") << "\n";
        }
    }

} // namespace bwmu
