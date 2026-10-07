#pragma once
#include <string>

namespace bwmu {

    // One particle species with its quantum numbers.
    //
    // mass : rest mass in GeV
    // g    : spin-isospin degeneracy (2J+1, times isospin multiplicity)
    // B    : baryon number
    // Q    : electric charge (in units of e)
    // S    : strangeness
    // boson: true for mesons, false for baryons
    struct Species {
        std::string name;
        double mass;
        double g;
        int    B;
        int    Q;
        int    S;
        bool   boson;
    };

} // namespace bwmu
