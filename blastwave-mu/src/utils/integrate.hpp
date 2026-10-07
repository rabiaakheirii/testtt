#pragma once
#include <vector>

namespace bwmu {

    // Trapezoidal rule on a uniform grid
    //   x.size() == f.size()
    double trapezoid(const std::vector<double>& f,
                     const std::vector<double>& x);

    // Simpson's 1/3 rule on a uniform grid
    //   x.size() == f.size() must be odd
    //   returns NAN if the input size is not suitable
    double simpson(const std::vector<double>& f,
                   const std::vector<double>& x);

} // namespace bwmu
