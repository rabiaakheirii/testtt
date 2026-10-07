#include "utils/bessel.hpp"
#include <gsl/gsl_sf_bessel.h>

namespace bwmu {

    double K0(double x) { return gsl_sf_bessel_K0(x); }
    double K1(double x) { return gsl_sf_bessel_K1(x); }
    double I0(double x) { return gsl_sf_bessel_I0(x); }
    double I1(double x) { return gsl_sf_bessel_I1(x); }

} // namespace bwmu
