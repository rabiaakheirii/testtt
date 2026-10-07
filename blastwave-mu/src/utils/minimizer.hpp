#pragma once
#include <functional>
#include <vector>

namespace bwmu {

    using CostFn = std::function<double(const std::vector<double>&)>;

    struct MinimizeResult {
        std::vector<double> x;    // best parameter vector
        double f;                 // best cost
        int    iterations;
        int    evaluations;
        bool   converged;
    };

    // ------------------------------------------------------------------
    // Nelder-Mead simplex minimizer.
    //
    // Parameters
    //   f          : cost function, one evaluation per call
    //   x0         : initial guess, size n
    //   step       : initial simplex step per dimension, size n
    //   tol_f      : convergence tolerance on the function value
    //   tol_x      : convergence tolerance on the simplex size
    //   max_iter   : hard iteration cap
    //
    // Returns the best point found.
    //
    // Reference: Nelder & Mead (1965), Computer Journal 7, 308.
    // ------------------------------------------------------------------
    MinimizeResult nelder_mead(
        const CostFn& f,
        const std::vector<double>& x0,
        const std::vector<double>& step,
        double tol_f = 1e-8,
        double tol_x = 1e-8,
        int    max_iter = 5000);

} // namespace bwmu
