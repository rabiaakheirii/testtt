#pragma once

namespace bwmu {

    // Blast-wave geometric and flow parameters.
    //
    //   T       : kinetic freeze-out temperature [GeV]
    //   beta_s  : surface transverse flow velocity [c]
    //   R_fm    : transverse radius [fm]
    //   n_prof  : profile exponent, beta_r(r) = beta_s * (r/R)^n
    //
    struct BWParams {
        double T      = 0.100;   // GeV
        double beta_s = 0.60;
        double R_fm   = 8.0;     // fm
        int    n_prof = 2;
    };

} // namespace bwmu
