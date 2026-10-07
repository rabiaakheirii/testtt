#pragma once
#include "particles/species.hpp"

namespace bwmu {

    // Chemical potential of a species:
    //   mu_i = B_i * mu_B + Q_i * mu_Q + S_i * mu_S
    // All in GeV.
    double mu_of(const Species& sp,
                 double mu_B,
                 double mu_Q,
                 double mu_S);

} // namespace bwmu
