#pragma once
#include <vector>

namespace bwmu {

    // Uniform 1D grid on [x_min, x_max] with n points
    struct Grid1D {
        double x_min = 0.0;
        double x_max = 0.0;
        int    n     = 0;
        std::vector<double> x;

        Grid1D() = default;
        Grid1D(double a, double b, int n_points);

        double dx() const;
        double operator[](int i) const { return x[i]; }
    };

} // namespace bwmu
