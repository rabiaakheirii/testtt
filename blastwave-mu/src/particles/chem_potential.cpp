#include "particles/chem_potential.hpp"

namespace bwmu {

    double mu_of(const Species& sp,
                 double mu_B,
                 double mu_Q,
                 double mu_S)
    {
        return sp.B * mu_B
             + sp.Q * mu_Q
             + sp.S * mu_S;
    }

} // namespace bwmu
