#include "utils/math.hpp"
#include <cmath>
#include "core/constants.hpp"

namespace qcd {
    double transverse_mass(double pT, double m) {
        return std::sqrt(pT*pT + m*m);
    }

    double rapidity(double E, double pL) {
        return 0.5 * std::log((E + pL) / (E - pL));
    }

    double mT_from_pT(double pT, double m) {
        return transverse_mass(pT, m);
    }
}