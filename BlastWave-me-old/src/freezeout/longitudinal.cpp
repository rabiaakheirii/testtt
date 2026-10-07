#include "freezeout/longitudinal.hpp"
#include <cmath>

namespace qcd {

    double thermal_dn_dy(double y, double m, double T, double norm, double y0) {
        const double ch = std::cosh(y - y0);
        const double r  = m / T;
        return norm * T * T * T
             * (r * r + r * 2.0 / ch + 2.0 / (ch * ch))
             * std::exp(-r * ch);
    }

    double flowing_dn_dy(double y, double m, double T, double norm,
                         double eta_max, double y0) {
        const int    N  = 200;
        const double de = 2.0 * eta_max / N;
        double sum = 0.0;
        for (int i = 0; i <= N; ++i) {
            const double eta = -eta_max + i * de;
            const double w   = (i == 0 || i == N) ? 0.5 : 1.0;
            sum += w * thermal_dn_dy(y - eta, m, T, norm, y0);
        }
        return sum * de;
    }
}
