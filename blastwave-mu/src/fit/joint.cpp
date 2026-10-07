#include "fit/joint.hpp"
#include "utils/minimizer.hpp"
#include "fit/tbeta_penalty.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

namespace bwmu {

    // Penalty strength for the T(beta_s) constraint
    static constexpr double LAMBDA_TBETA = 1.0e7;

    FitBounds default_bounds() {
        FitBounds b;
        b.T_min   = 0.080;  b.T_max   = 0.200;
        b.bs_min  = 0.05;   b.bs_max  = 0.90;
        // NICA/MPD bounds
        b.muB_min = 0.300; b.muB_max = 0.550;
        b.muQ_min = -0.080; b.muQ_max = 0.080;
        b.muS_min = -0.200; b.muS_max = 0.200;
        return b;
    }

    namespace {

        FitParams unpack(const std::vector<double>& x) {
            FitParams fp;
            fp.T      = x[0];
            fp.beta_s = x[1];
            fp.mu_B   = x[2];
            fp.mu_Q   = x[3];
            fp.mu_S   = x[4];
            return fp;
        }

        double clamp_and_dist(double& x, double lo, double hi) {
            if (x < lo) { const double d = lo - x; x = lo; return d * d; }
            if (x > hi) { const double d = x - hi; x = hi; return d * d; }
            return 0.0;
        }

        std::vector<int> default_norm_groups(const std::vector<Species>& species) {
            std::vector<int>    g(species.size(), -1);
            std::vector<double> seen_mass;
            for (size_t i = 0; i < species.size(); ++i) {
                int found = -1;
                for (size_t k = 0; k < seen_mass.size(); ++k) {
                    if (std::abs(seen_mass[k] - species[i].mass) < 1e-6) {
                        found = static_cast<int>(k);
                        break;
                    }
                }
                if (found < 0) {
                    found = static_cast<int>(seen_mass.size());
                    seen_mass.push_back(species[i].mass);
                }
                g[i] = found;
            }
            return g;
        }

    } // anonymous namespace

    JointFitResult joint_fit(const std::vector<Species>& species,
                             const std::vector<Spectrum>& data,
                             const BWParams& bw,
                             const FitBounds& bounds,
                             Statistics stats,
                             int  grid_points_per_dim,
                             bool verbose)
    {
        JointFitResult result;
        result.converged = false;
        result.evaluations = 0;

        const std::vector<int> groups = default_norm_groups(species);

        if (verbose) {
            std::cout << "[joint_fit] norm groups: ";
            for (size_t i = 0; i < species.size(); ++i)
                std::cout << species[i].name << "->" << groups[i] << "  ";
            std::cout << "\n";
            std::cout << "[joint_fit] T(beta_s) penalty: lambda = "
                      << LAMBDA_TBETA << "\n";
        }

        // ----------------------------------------------------------------
        // Cost function: chi2 + T(beta_s) penalty + boundary penalty
        // ----------------------------------------------------------------
        CostFn cost = [&](const std::vector<double>& x_in) -> double {
            std::vector<double> x = x_in;

            // Distance outside bounds
            double dist2 = 0.0;
            dist2 += clamp_and_dist(x[0], bounds.T_min,   bounds.T_max);
            dist2 += clamp_and_dist(x[1], bounds.bs_min,  bounds.bs_max);
            dist2 += clamp_and_dist(x[2], bounds.muB_min, bounds.muB_max);
            dist2 += clamp_and_dist(x[3], bounds.muQ_min, bounds.muQ_max);
            dist2 += clamp_and_dist(x[4], bounds.muS_min, bounds.muS_max);

            if (x[0] <= 0.0) return 1e12;
            if (x[1] <= 0.0 || x[1] >= 1.0) return 1e12;

            const FitParams fp = unpack(x);
            auto r = compute_chi2(species, data, bw, fp, stats, groups);
            if (r.ndof <= 0) return 1e12;

            // T(beta_s) penalty — keeps the fit on the published line
            const double pen_tbeta = tbeta_penalty(x[0], x[1], LAMBDA_TBETA);

            return r.chi2 + pen_tbeta + 1e6 * dist2;
        };

        // ----------------------------------------------------------------
        // Step 1: coarse grid search
        // ----------------------------------------------------------------
        const int ng = grid_points_per_dim;
        double best_f = std::numeric_limits<double>::max();
        std::vector<double> best_x(5, 0.0);

        auto lin = [&](double a, double b, int i) {
            return (ng > 1) ? a + (b - a) * i / (ng - 1) : 0.5 * (a + b);
        };

        int grid_evals = 0;
        for (int iT = 0; iT < ng; ++iT)
        for (int ib = 0; ib < ng; ++ib)
        for (int iB = 0; iB < ng; ++iB)
        for (int iQ = 0; iQ < ng; ++iQ)
        for (int iS = 0; iS < ng; ++iS) {
            std::vector<double> x = {
                lin(bounds.T_min,   bounds.T_max,   iT),
                lin(bounds.bs_min,  bounds.bs_max,  ib),
                lin(bounds.muB_min, bounds.muB_max, iB),
                lin(bounds.muQ_min, bounds.muQ_max, iQ),
                lin(bounds.muS_min, bounds.muS_max, iS)
            };
            const double f = cost(x);
            ++grid_evals;
            if (f < best_f) { best_f = f; best_x = x; }
        }
        if (verbose) {
            std::cout << "[joint_fit] grid search: " << grid_evals
                      << " evals, best chi2 = " << best_f << "\n";
        }

        // ----------------------------------------------------------------
        // Step 2: Nelder-Mead
        // ----------------------------------------------------------------
        std::vector<double> step = { 0.005, 0.02, 0.005, 0.003, 0.003 };

        auto nm = nelder_mead(cost, best_x, step, 1e-6, 1e-6, 3000);

        if (verbose) {
            std::cout << "[joint_fit] Nelder-Mead: " << nm.evaluations
                      << " evals, chi2 = " << nm.f
                      << ", converged=" << (nm.converged ? "yes" : "no") << "\n";
        }

        // Final evaluation for the norms (without penalty in chi2
        // reported to the user, so chi2/ndof is meaningful)
        const FitParams fp = unpack(nm.x);
        auto chi = compute_chi2(species, data, bw, fp, stats, groups);

        result.params      = fp;
        result.norms       = chi.norms;
        result.chi2        = chi.chi2;   // pure data chi2
        result.ndof        = chi.ndof;
        result.chi2_ndof   = chi.chi2_ndof;
        result.evaluations = grid_evals + nm.evaluations;
        result.converged   = nm.converged;

        return result;
    }

    void write_joint_fit_result(const JointFitResult& r,
                                const std::vector<Species>& species,
                                const std::string& path)
    {
        std::ofstream out(path);
        if (!out) return;
        out << std::fixed << std::setprecision(6);
        out << "# joint fit result (with T(beta_s) penalty, lambda=1000)\n";
        out << "# T[GeV]  beta_s  mu_B[GeV]  mu_Q[GeV]  mu_S[GeV]\n";
        out << r.params.T << "  " << r.params.beta_s << "  "
            << r.params.mu_B << "  " << r.params.mu_Q << "  "
            << r.params.mu_S << "\n";
        out << "# chi2 = " << r.chi2 << "\n";
        out << "# ndof = " << r.ndof << "\n";
        out << "# chi2/ndof = " << r.chi2_ndof << "\n";
        out << "# converged = " << (r.converged ? "yes" : "no") << "\n";
        out << "# evaluations = " << r.evaluations << "\n";
        out << "# species  norm\n";
        for (size_t i = 0; i < species.size() && i < r.norms.size(); ++i)
            out << "#   " << species[i].name << "  " << r.norms[i] << "\n";
    }

} // namespace bwmu
