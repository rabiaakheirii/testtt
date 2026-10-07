#include <cmath>
#include <iostream>
#include <vector>

#include "blastwave/params.hpp"
#include "fit/chi2.hpp"
#include "fit/joint.hpp"
#include "particles/table.hpp"
#include "utils/minimizer.hpp"
#include "utils/random.hpp"

static int failures = 0;
#define CHECK(cond, msg) \
    do { if (!(cond)) { std::cout << "  FAIL: " << msg << "\n"; ++failures; } \
         else { std::cout << "  ok  : " << msg << "\n"; } } while (0)

static bool close(double a, double b, double tol) {
    return std::abs(a - b) < tol;
}

static bwmu::Spectrum make_synthetic(const bwmu::Species& sp,
                                     const bwmu::BWParams& bw,
                                     const bwmu::FitParams& fp,
                                     double A, const std::vector<double>& mT_grid)
{
    bwmu::Spectrum s;
    s.species = sp.name; s.centrality = 0;
    for (double mT : mT_grid) {
        const double f = bwmu::model_shape(mT, sp, bw, fp,
                                           bwmu::Statistics::Quantum);
        bwmu::DataPoint d;
        d.mT    = mT;
        d.value = A * f;
        d.err   = 1e-4 * std::abs(A * f) + 1e-9;
        s.points.push_back(d);
    }
    return s;
}

int main() {
    std::cout << "========== Phase 6 unit tests ==========\n";

    // --- 1D and 2D minimizer tests (unchanged) ---
    std::cout << "\n[Minimizer: 1D parabola]\n";
    {
        auto f = [](const std::vector<double>& x) { return (x[0]-3.0)*(x[0]-3.0); };
        auto r = bwmu::nelder_mead(f, {0.0}, {1.0}, 1e-10, 1e-10, 1000);
        CHECK(close(r.x[0], 3.0, 1e-4), "1D min at x = 3");
    }
    std::cout << "\n[Minimizer: 2D Rosenbrock]\n";
    {
        auto f = [](const std::vector<double>& x) {
            const double a = 1.0 - x[0], b = x[1] - x[0]*x[0];
            return a*a + 100.0*b*b;
        };
        auto r = bwmu::nelder_mead(f, {-1.2, 1.0}, {0.1, 0.1}, 1e-10, 1e-10, 5000);
        CHECK(close(r.x[0], 1.0, 1e-3), "rosenbrock x0 = 1");
        CHECK(close(r.x[1], 1.0, 1e-3), "rosenbrock x1 = 1");
    }
    std::cout << "\n[Minimizer: 5D quadratic]\n";
    {
        const std::vector<double> xs = {0.1, 0.6, 0.03, -0.005, 0.002};
        auto f = [&](const std::vector<double>& x) {
            double s = 0.0;
            for (size_t i = 0; i < x.size(); ++i) { const double d = x[i]-xs[i]; s += d*d; }
            return s;
        };
        auto r = bwmu::nelder_mead(f, {0,0.5,0,0,0}, {0.01,0.01,0.01,0.01,0.01},
                                   1e-10, 1e-10, 5000);
        bool ok = true;
        for (size_t i = 0; i < xs.size(); ++i) if (!close(r.x[i], xs[i], 1e-3)) ok = false;
        CHECK(ok, "5D quadratic minimum recovered");
    }

    // --- Shared normalization test (unchanged) ---
    std::cout << "\n[chi2: shared normalization groups]\n";
    {
        bwmu::BWParams bw; bw.R_fm = 8.0; bw.n_prof = 2;
        bwmu::FitParams fp_true;
        fp_true.T = 0.105; fp_true.beta_s = 0.62;
        fp_true.mu_B = 0.025; fp_true.mu_Q = -0.007; fp_true.mu_S = 0.003;

        const auto& pi_p = bwmu::get_species("pi+");
        const auto& pi_m = bwmu::get_species("pi-");
        const auto& K_p  = bwmu::get_species("K+");
        const auto& K_m  = bwmu::get_species("K-");
        std::vector<double> mT_grid = {0.20, 0.40, 0.70, 1.00, 1.40};

        std::vector<bwmu::Spectrum> data = {
            make_synthetic(pi_p, bw, fp_true, 10.0, mT_grid),
            make_synthetic(pi_m, bw, fp_true, 10.0, mT_grid),
            make_synthetic(K_p,  bw, fp_true,  1.5, mT_grid),
            make_synthetic(K_m,  bw, fp_true,  1.5, mT_grid)
        };
        std::vector<bwmu::Species> species = { pi_p, pi_m, K_p, K_m };
        std::vector<int> groups = { 0, 0, 1, 1 };

        auto r = bwmu::compute_chi2(species, data, bw, fp_true,
                                    bwmu::Statistics::Quantum, groups);
        CHECK(r.chi2 < 1e-15, "chi2 ~ 0 for perfect data");
        CHECK(close(r.norms[0], 10.0, 1e-6), "A_pi recovered");
        CHECK(close(r.norms[1],  1.5, 1e-6), "A_K recovered");
    }

    // --- Joint fit with WIDE mT range (the key fix) ---
    std::cout << "\n[Joint fit: wide mT range, 6 species]\n";
    {
        bwmu::BWParams bw; bw.R_fm = 8.0; bw.n_prof = 2;
        bwmu::FitParams fp_true;
        fp_true.T = 0.105; fp_true.beta_s = 0.62;
        fp_true.mu_B = 0.025; fp_true.mu_Q = -0.007; fp_true.mu_S = 0.003;

        const auto& pi_p = bwmu::get_species("pi+");
        const auto& pi_m = bwmu::get_species("pi-");
        const auto& K_p  = bwmu::get_species("K+");
        const auto& K_m  = bwmu::get_species("K-");
        const auto& p    = bwmu::get_species("p");
        const auto& pbar = bwmu::get_species("pbar");

        // Wide grid: spans thermal and flow regimes
        std::vector<double> mT_grid = {0.10, 0.20, 0.30, 0.45, 0.60,
                                        0.80, 1.00, 1.25, 1.50, 2.00};

        std::vector<bwmu::Spectrum> data = {
            make_synthetic(pi_p, bw, fp_true, 10.0, mT_grid),
            make_synthetic(pi_m, bw, fp_true, 10.0, mT_grid),
            make_synthetic(K_p,  bw, fp_true,  1.5, mT_grid),
            make_synthetic(K_m,  bw, fp_true,  1.5, mT_grid),
            make_synthetic(p,    bw, fp_true,  0.3, mT_grid),
            make_synthetic(pbar, bw, fp_true,  0.3, mT_grid)
        };
        std::vector<bwmu::Species> species = { pi_p, pi_m, K_p, K_m, p, pbar };

        bwmu::FitBounds bounds;
        bounds.T_min = 0.085; bounds.T_max = 0.150;
        bounds.bs_min = 0.40; bounds.bs_max = 0.85;
        bounds.muB_min = 0.000; bounds.muB_max = 0.050;
        bounds.muQ_min = -0.020; bounds.muQ_max = 0.020;
        bounds.muS_min = -0.015; bounds.muS_max = 0.015;

        auto res = bwmu::joint_fit(species, data, bw, bounds,
                                   bwmu::Statistics::Quantum, 4, true);

        std::cout << "  T       = " << res.params.T * 1000 << " MeV"
                  << "  (true " << fp_true.T * 1000 << ")\n";
        std::cout << "  beta_s  = " << res.params.beta_s
                  << "  (true " << fp_true.beta_s << ")\n";
        std::cout << "  mu_B    = " << res.params.mu_B * 1000 << " MeV"
                  << "  (true " << fp_true.mu_B * 1000 << ")\n";
        std::cout << "  mu_Q    = " << res.params.mu_Q * 1000 << " MeV"
                  << "  (true " << fp_true.mu_Q * 1000 << ")\n";
        std::cout << "  mu_S    = " << res.params.mu_S * 1000 << " MeV"
                  << "  (true " << fp_true.mu_S * 1000 << ")\n";
        std::cout << "  chi2/ndof = " << res.chi2_ndof << "\n";

        CHECK(close(res.params.T,      fp_true.T,      0.010), "T recovered within 10 MeV");
        CHECK(close(res.params.beta_s, fp_true.beta_s, 0.05),  "beta_s recovered within 0.05");
        CHECK(close(res.params.mu_B,   fp_true.mu_B,   0.010), "mu_B recovered within 10 MeV");
        CHECK(close(res.params.mu_Q,   fp_true.mu_Q,   0.008), "mu_Q recovered within 8 MeV");
        CHECK(close(res.params.mu_S,   fp_true.mu_S,   0.008), "mu_S recovered within 8 MeV");
    }

    std::cout << "\n========================================\n";
    if (failures == 0) { std::cout << " All Phase 6 tests passed.\n"; return 0; }
    std::cout << " " << failures << " test(s) failed.\n"; return 1;
}
