#include "utils/integrate.hpp"
#include <cmath>
#include <limits>

namespace bwmu {

    double trapezoid(const std::vector<double>& f,
                     const std::vector<double>& x)
    {
        const int n = static_cast<int>(f.size());
        if (n < 2 || static_cast<int>(x.size()) != n) return 0.0;

        double sum = 0.5 * (f[0] + f[n-1]) * (x[1] - x[0]);
        for (int i = 1; i < n - 1; ++i) {
            const double dxl = x[i] - x[i-1];
            const double dxr = x[i+1] - x[i];
            sum += 0.5 * f[i] * (dxl + dxr);
        }
        return sum;
    }

    double simpson(const std::vector<double>& f,
                   const std::vector<double>& x)
    {
        const int n = static_cast<int>(f.size());
        if (n < 3 || (n % 2) == 0) return std::numeric_limits<double>::quiet_NaN();
        if (static_cast<int>(x.size()) != n) return std::numeric_limits<double>::quiet_NaN();

        const double h = (x[n-1] - x[0]) / (n - 1);

        double sum = f[0] + f[n-1];
        for (int i = 1; i < n - 1; ++i) {
            sum += (i % 2 == 0 ? 2.0 : 4.0) * f[i];
        }
        return sum * h / 3.0;
    }

} // namespace bwmu
