#include "utils/minimizer.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace bwmu {

    namespace {

        struct SimplexVertex {
            std::vector<double> x;
            double f;
        };

        double simplex_size(const std::vector<SimplexVertex>& s) {
            // max distance from the best vertex
            const auto& best = s.front().x;
            double d = 0.0;
            for (size_t i = 1; i < s.size(); ++i) {
                double d2 = 0.0;
                for (size_t k = 0; k < best.size(); ++k) {
                    const double dx = s[i].x[k] - best[k];
                    d2 += dx * dx;
                }
                d = std::max(d, std::sqrt(d2));
            }
            return d;
        }

    } // anonymous namespace

    MinimizeResult nelder_mead(const CostFn& f,
                               const std::vector<double>& x0,
                               const std::vector<double>& step,
                               double tol_f,
                               double tol_x,
                               int    max_iter)
    {
        const int n = static_cast<int>(x0.size());
        MinimizeResult res;
        res.evaluations = 0;
        res.iterations  = 0;
        res.converged   = false;

        auto eval = [&](const std::vector<double>& x) {
            ++res.evaluations;
            return f(x);
        };

        // Standard coefficients
        const double alpha = 1.0;   // reflection
        const double gamma = 2.0;   // expansion
        const double rho   = 0.5;   // contraction
        const double sigma = 0.5;   // shrink

        // Build initial simplex: n+1 vertices
        std::vector<SimplexVertex> simplex(n + 1);
        simplex[0].x = x0;
        simplex[0].f = eval(x0);
        for (int i = 0; i < n; ++i) {
            std::vector<double> x = x0;
            x[i] += (step[i] != 0.0) ? step[i] : 1e-3;
            simplex[i + 1].x = x;
            simplex[i + 1].f = eval(x);
        }

        for (int iter = 0; iter < max_iter; ++iter) {
            res.iterations = iter + 1;

            // Sort vertices by function value
            std::sort(simplex.begin(), simplex.end(),
                      [](const SimplexVertex& a, const SimplexVertex& b) {
                          return a.f < b.f;
                      });

            const double best_f = simplex.front().f;
            const double worst_f = simplex.back().f;
            const double second_worst_f = simplex[n - 1].f;

            // Convergence
            if (std::abs(worst_f - best_f) < tol_f &&
                simplex_size(simplex) < tol_x)
            {
                res.converged = true;
                break;
            }

            // Centroid of the best n vertices
            std::vector<double> centroid(n, 0.0);
            for (int i = 0; i < n; ++i) {
                for (int k = 0; k < n; ++k) {
                    centroid[k] += simplex[i].x[k];
                }
            }
            for (int k = 0; k < n; ++k) centroid[k] /= n;

            // Reflect
            std::vector<double> xr(n);
            for (int k = 0; k < n; ++k) {
                xr[k] = centroid[k] + alpha * (centroid[k] - simplex[n].x[k]);
            }
            const double fr = eval(xr);

            if (fr < second_worst_f && fr >= best_f) {
                simplex[n].x = xr;
                simplex[n].f = fr;
                continue;
            }

            // Expand
            if (fr < best_f) {
                std::vector<double> xe(n);
                for (int k = 0; k < n; ++k) {
                    xe[k] = centroid[k] + gamma * (xr[k] - centroid[k]);
                }
                const double fe = eval(xe);
                if (fe < fr) {
                    simplex[n].x = xe;
                    simplex[n].f = fe;
                } else {
                    simplex[n].x = xr;
                    simplex[n].f = fr;
                }
                continue;
            }

            // Contract
            std::vector<double> xc(n);
            if (fr < worst_f) {
                // outside contraction
                for (int k = 0; k < n; ++k) {
                    xc[k] = centroid[k] + rho * (xr[k] - centroid[k]);
                }
                const double fc = eval(xc);
                if (fc <= fr) {
                    simplex[n].x = xc;
                    simplex[n].f = fc;
                    continue;
                }
            } else {
                // inside contraction
                for (int k = 0; k < n; ++k) {
                    xc[k] = centroid[k] + rho * (simplex[n].x[k] - centroid[k]);
                }
                const double fc = eval(xc);
                if (fc < worst_f) {
                    simplex[n].x = xc;
                    simplex[n].f = fc;
                    continue;
                }
            }

            // Shrink
            for (int i = 1; i <= n; ++i) {
                for (int k = 0; k < n; ++k) {
                    simplex[i].x[k] =
                        simplex[0].x[k] + sigma * (simplex[i].x[k] - simplex[0].x[k]);
                }
                simplex[i].f = eval(simplex[i].x);
            }
        }

        // Final sort and return best
        std::sort(simplex.begin(), simplex.end(),
                  [](const SimplexVertex& a, const SimplexVertex& b) {
                      return a.f < b.f;
                  });

        res.x = simplex.front().x;
        res.f = simplex.front().f;
        return res;
    }

} // namespace bwmu
