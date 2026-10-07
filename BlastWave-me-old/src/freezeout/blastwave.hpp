#pragma once

namespace qcd {

    // Phase 1: pure thermal spectrum
    double thermal_mT_spectrum(double mT, double T, double norm);

    // Phase 4: full blast-wave spectrum (paper Eq. 7)
    //   dn/(mT dmT) = norm * mT * ∫_0^R r dr I0(mT sinh(rho)/T) K1(mT cosh(rho)/T)
    // with rho(r) = atanh( beta_s * (r/R)^n )
    double blastwave_mT_spectrum(double mT, double T,
                                 double beta_s, double R, int n,
                                 double norm);

    // Apparent slope temperature at high mT (paper Eq. 19)
    double effective_temperature(double T, double beta);
}