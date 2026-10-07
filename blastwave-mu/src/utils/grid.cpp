#include "utils/grid.hpp"

namespace bwmu {

    Grid1D::Grid1D(double a, double b, int n_points)
        : x_min(a), x_max(b), n(n_points)
    {
        x.resize(n);
        const double d = (n > 1) ? (b - a) / (n - 1) : 0.0;
        for (int i = 0; i < n; ++i) x[i] = a + i * d;
    }

    double Grid1D::dx() const {
        return (n > 1) ? (x_max - x_min) / (n - 1) : 0.0;
    }

} // namespace bwmu
