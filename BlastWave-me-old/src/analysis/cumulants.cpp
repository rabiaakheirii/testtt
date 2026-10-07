#include "analysis/cumulants.hpp"
#include <cmath>
#include <iostream>

namespace qcd {

    Cumulants compute_cumulants(const std::vector<int>& N) {
        Cumulants c{};
        c.n_events = static_cast<int>(N.size());
        if (N.empty()) return c;

        // C1 = mean
        double sum = 0.0;
        for (int x : N) sum += x;
        c.C1 = sum / N.size();

        // central moments
        double m2 = 0.0, m3 = 0.0, m4 = 0.0;
        for (int x : N) {
            const double d = x - c.C1;
            m2 += d * d;
            m3 += d * d * d;
            m4 += d * d * d * d;
        }
        m2 /= N.size();
        m3 /= N.size();
        m4 /= N.size();

        c.C2 = m2;
        c.C3 = m3;
        c.C4 = m4 - 3.0 * m2 * m2;

        c.Ssigma       = (c.C2 > 0.0) ? c.C3 / c.C2 : 0.0;
        c.kappa_sigma2 = (c.C2 > 0.0) ? c.C4 / c.C2 : 0.0;

        return c;
    }
}
