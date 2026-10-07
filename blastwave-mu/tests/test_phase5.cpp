#include <cmath>
#include <iostream>
#include <vector>

#include "blastwave/params.hpp"
#include "blastwave/spectrum.hpp"
#include "fit/chi2.hpp"
#include "io/datapoint.hpp"
#include "particles/table.hpp"
#include "utils/random.hpp"

static int failures = 0;

#define CHECK(cond, msg)                                       \
    do {                                                       \
        if (!(cond)) {                                         \
            std::cout << "  FAIL: " << msg << "\n";            \
            ++failures;                                        \
        } else {                                               \
            std::cout << "  ok  : " << msg << "\n";            \
        }                                                      \
    } while (0)

static bool close(double a, double b, double rel = 1e-6) {
    if (a == 0.0 && b == 0.0) return true;
    const double denom = std::max(std::abs(a), std::abs(b));
    return std::abs(a - b) / denom < rel;
}

// Generate a synthetic spectrum from a known model with a known A_i
static bwmu::Spectrum make_synthetic(
    const bwmu::Species& sp,
    const bwmu::BWParams& bw,
    const bwmu::FitParams& fp,
    double A_true,
    const std::vector<double>& mT_grid,
    double noise_sigma,
    unsigned seed)
{
    bwmu::Spectrum s;
    s.species    = sp.name;
    s.centrality = 0;

    bwmu::RNG rng(seed);
    for (double mT : mT_grid) {
        const double f = bwmu::model_shape(mT, sp, bw, fp,
                                           bwmu::Statistics::Quantum);
        const double v = A_true * f;
        const double noise = (noise_sigma > 0)
            ? rng.gaussian(0.0, noise_sigma)
            : 0.0;
        bwmu::DataPoint d;
        d.mT    = mT;
        d.value = v + noise;
        d.err   = (noise_sigma > 0) ? noise_sigma : 0.01 * std::abs(v);
        if (d.err <= 0.0) d.err = 1e-6;
        s.points.push_back(d);
    }
    return s;
}

int main() {
    std::cout << "========== Phase 5 unit tests ==========\n";

    bwmu::BWParams bw;
    bw.R_fm   = 8.0;
    bw.n_prof = 2;
    // T and beta_s come from FitParams

    bwmu::FitParams fp;
    fp.T      = 0.100;
    fp.beta_s = 0.60;
    fp.mu_B   = 0.030;
    fp.mu_Q   = -0.005;
    fp.mu_S   = 0.002;

    const auto& pi_p  = bwmu::get_species("pi+");
    const auto& pi_m  = bwmu::get_species("pi-");
    const auto& K_p   = bwmu::get_species("K+");

    const std::vector<double> mT_grid = {
        0.15, 0.25, 0.35, 0.45, 0.55, 0.70, 0.90, 1.20
    };

    // --------------------------------------------------------------
    // Test 1: zero data -> chi2 = 0
    // --------------------------------------------------------------
    std::cout << "\n[Empty data]\n";
    {
        std::vector<bwmu::Species> sp = { pi_p };
        std::vector<bwmu::Spectrum> dt;   // empty
        auto res = bwmu::compute_chi2(sp, dt, bw, fp);
        CHECK(res.chi2 == 0.0, "chi2 = 0 when no data");
        CHECK(res.ndof <= 0, "ndof <= 0 when no data");
    }

    // --------------------------------------------------------------
    // Test 2: perfect data -> chi2 ~ 0
    // --------------------------------------------------------------
    std::cout << "\n[Perfect data]\n";
    {
        const double A_true = 42.0;
        auto spec = make_synthetic(pi_p, bw, fp, A_true, mT_grid, 0.0, 0);

        std::vector<bwmu::Species> sp = { pi_p };
        std::vector<bwmu::Spectrum> dt = { spec };

        auto res = bwmu::compute_chi2(sp, dt, bw, fp);
        std::cout << "  chi2 = " << res.chi2
                  << "  A = " << res.norms[0]
                  << "  A_true = " << A_true << "\n";
        CHECK(res.chi2 < 1e-20, "chi2 ~ 0 for perfect data");
        CHECK(close(res.norms[0], A_true, 1e-10),
              "recovered A_true");
        CHECK(res.ndof == static_cast<int>(mT_grid.size()) - 1,
              "ndof = n_points - n_species");
    }

    // --------------------------------------------------------------
    // Test 3: linearity — scaling A_true by k scales fitted A by k
    // --------------------------------------------------------------
    std::cout << "\n[Linearity in A]\n";
    {
        auto s1 = make_synthetic(pi_p, bw, fp,  1.0, mT_grid, 0.0, 0);
        auto s2 = make_synthetic(pi_p, bw, fp, 10.0, mT_grid, 0.0, 0);
        std::vector<bwmu::Species> sp = { pi_p };
        auto r1 = bwmu::compute_chi2(sp, { s1 }, bw, fp);
        auto r2 = bwmu::compute_chi2(sp, { s2 }, bw, fp);
        std::cout << "  A1 = " << r1.norms[0]
                  << "  A2 = " << r2.norms[0]
                  << "  ratio = " << r2.norms[0] / r1.norms[0] << "\n";
        CHECK(close(r2.norms[0] / r1.norms[0], 10.0, 1e-10),
              "A scales linearly");
    }

    // --------------------------------------------------------------
    // Test 4: multi-species — 2 species fit independently
    // --------------------------------------------------------------
    std::cout << "\n[Multi-species]\n";
    {
        auto s_pi = make_synthetic(pi_p, bw, fp,  5.0, mT_grid, 0.0, 0);
        auto s_K  = make_synthetic(K_p,  bw, fp,  0.5, mT_grid, 0.0, 0);
        std::vector<bwmu::Species> sp = { pi_p, K_p };
        std::vector<bwmu::Spectrum> dt = { s_pi, s_K };

        auto res = bwmu::compute_chi2(sp, dt, bw, fp);
        std::cout << "  A_pi = " << res.norms[0]
                  << "  A_K  = " << res.norms[1] << "\n";
        CHECK(close(res.norms[0], 5.0, 1e-10), "A_pi recovered");
        CHECK(close(res.norms[1], 0.5, 1e-10), "A_K recovered");
        CHECK(res.ndof == static_cast<int>(2 * mT_grid.size()) - 2,
              "ndof for 2 species");
    }

    // --------------------------------------------------------------
    // Test 5: chi2 grows when parameters are wrong
    // --------------------------------------------------------------
    std::cout << "\n[chi2 grows for wrong parameters]\n";
    {
        // Generate data with T = 0.100, beta_s = 0.60
        auto spec = make_synthetic(pi_p, bw, fp, 10.0, mT_grid, 0.0, 0);
        std::vector<bwmu::Species> sp = { pi_p };

        // Chi2 at true parameters
        auto r_true = bwmu::compute_chi2(sp, { spec }, bw, fp);

        // Chi2 at wrong T
        bwmu::FitParams fp_wrong = fp;
        fp_wrong.T = 0.150;
        auto r_wrongT = bwmu::compute_chi2(sp, { spec }, bw, fp_wrong);

        // Chi2 at wrong beta_s
        bwmu::FitParams fp_wrong_b = fp;
        fp_wrong_b.beta_s = 0.30;
        auto r_wrongB = bwmu::compute_chi2(sp, { spec }, bw, fp_wrong_b);

        std::cout << "  chi2(true)      = " << r_true.chi2 << "\n";
        std::cout << "  chi2(wrong T)   = " << r_wrongT.chi2 << "\n";
        std::cout << "  chi2(wrong bs)  = " << r_wrongB.chi2 << "\n";

        CHECK(r_true.chi2 < 1e-20, "chi2 at truth is ~0");
        CHECK(r_wrongT.chi2 > 1e3, "chi2 large for wrong T");
        CHECK(r_wrongB.chi2 > 1e3, "chi2 large for wrong beta_s");
    }

    // --------------------------------------------------------------
    // Test 6: quantum vs Boltzmann give different chi2
    // --------------------------------------------------------------
    std::cout << "\n[Quantum vs Boltzmann chi2]\n";
    {
        auto spec = make_synthetic(pi_p, bw, fp, 10.0, mT_grid, 0.0, 0);
        std::vector<bwmu::Species> sp = { pi_p };
        auto r_q = bwmu::compute_chi2(sp, { spec }, bw, fp,
                                      bwmu::Statistics::Quantum);
        auto r_b = bwmu::compute_chi2(sp, { spec }, bw, fp,
                                      bwmu::Statistics::Boltzmann);
        std::cout << "  chi2 quantum    = " << r_q.chi2 << "\n";
        std::cout << "  chi2 Boltzmann  = " << r_b.chi2 << "\n";
        std::cout << "  A quantum       = " << r_q.norms[0] << "\n";
        std::cout << "  A Boltzmann     = " << r_b.norms[0] << "\n";
        // The data was generated with QUANTUM statistics.
        //   - Fitting with quantum -> chi2 ~ 0 (shape matches)
        //   - Fitting with Boltzmann -> chi2 >> 0 (shape mismatch)
        // This is the physical signal we are trying to measure: the
        // Bose/Fermi correction is not negligible for pions.
        CHECK(r_q.chi2 < 1e-20, "quantum fit of quantum data: chi2 ~ 0");
        CHECK(r_b.chi2 > 100.0,
              "Boltzmann fit of quantum data: chi2 >> 0 (Bose enhancement)");
        CHECK(!close(r_q.norms[0], r_b.norms[0], 1e-3),
              "A_quantum != A_boltzmann (different shapes)");
        // The Boltzmann fit needs a larger normalisation to compensate
        // for the missing Bose enhancement at low mT.
        CHECK(r_b.norms[0] > r_q.norms[0],
              "A_boltzmann > A_quantum (compensates for missing Bose peak)");
    }

    // --------------------------------------------------------------
    // Test 7: chi2 with realistic noise
    // --------------------------------------------------------------
    std::cout << "\n[Noisy data]\n";
    {
        const double sigma = 0.05;   // 5% noise
        auto spec = make_synthetic(pi_p, bw, fp, 10.0, mT_grid, sigma, 7);
        std::vector<bwmu::Species> sp = { pi_p };
        auto res = bwmu::compute_chi2(sp, { spec }, bw, fp);
        std::cout << "  chi2/ndof = " << res.chi2_ndof << "\n";
        // For Gaussian noise, chi2/ndof ~ 1 on average
        CHECK(res.chi2_ndof < 5.0, "chi2/ndof reasonable with 5% noise");
    }

    // --------------------------------------------------------------
    // Test 8: species with no data -> skipped, norm = 0
    // --------------------------------------------------------------
    std::cout << "\n[Species without data]\n";
    {
        auto spec = make_synthetic(pi_p, bw, fp, 10.0, mT_grid, 0.0, 0);
        std::vector<bwmu::Species> sp = { pi_p, K_p };  // K_p has no data
        std::vector<bwmu::Spectrum> dt = { spec };
        auto res = bwmu::compute_chi2(sp, dt, bw, fp);
        std::cout << "  A_pi = " << res.norms[0]
                  << "  A_K  = " << res.norms[1] << "\n";
        CHECK(res.norms[1] == 0.0, "missing species has norm 0");
    }

    std::cout << "\n========================================\n";
    if (failures == 0) {
        std::cout << " All Phase 5 tests passed.\n";
        return 0;
    } else {
        std::cout << " " << failures << " test(s) failed.\n";
        return 1;
    }
}
