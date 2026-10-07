#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "fit/chi2_surface.hpp"
#include "blastwave/params.hpp"
#include "fit/chi2.hpp"
#include "particles/table.hpp"

using namespace bwmu;

static int failures = 0;
#define CHECK(cond, msg) \
    do { if (!(cond)) { std::cout << "  FAIL: " << msg << "\n"; ++failures; } \
         else { std::cout << "  ok  : " << msg << "\n"; } } while (0)

static Spectrum make_synth(const Species& sp, const BWParams& bw,
                           const FitParams& fp, double A,
                           const std::vector<double>& mTs)
{
    Spectrum s; s.species = sp.name; s.centrality = 0;
    for (double mT : mTs) {
        const double f = model_shape(mT, sp, bw, fp, Statistics::Quantum);
        DataPoint d;
        d.mT    = mT;
        d.value = A * f;
        d.err   = 1e-4 * std::abs(A * f) + 1e-9;
        s.points.push_back(d);
    }
    return s;
}

int main() {
    std::cout << "========== Phase 10 unit tests ==========\n\n";

    BWParams bw; bw.R_fm = 12.0; bw.n_prof = 2;
    FitParams truth;
    truth.T = 0.105; truth.beta_s = 0.62;
    truth.mu_B = 0.025; truth.mu_Q = -0.007; truth.mu_S = 0.003;

    const auto& pi_p = get_species("pi+");
    const auto& pi_m = get_species("pi-");
    const auto& K_p  = get_species("K+");
    const auto& K_m  = get_species("K-");

    std::vector<double> mTs = {0.10, 0.20, 0.30, 0.45, 0.60,
                               0.80, 1.00, 1.25, 1.50, 2.00};

    std::vector<Spectrum> data = {
        make_synth(pi_p, bw, truth, 10.0, mTs),
        make_synth(pi_m, bw, truth, 10.0, mTs),
        make_synth(K_p,  bw, truth,  1.5, mTs),
        make_synth(K_m,  bw, truth,  1.5, mTs)
    };
    std::vector<Species> species = { pi_p, pi_m, K_p, K_m };

    std::cout << "[Test 1: grid size]\n";
    auto surface = compute_chi2_surface(species, data, bw, truth,
                                        0.080, 0.150, 8,
                                        0.40,  0.80, 11,
                                        Statistics::Quantum);
    CHECK(surface.size() == 88, "8 x 11 = 88 grid points");

    std::cout << "\n[Test 2: chi2 values]\n";
    bool all_ok = true;
    for (const auto& p : surface)
        if (!std::isfinite(p.chi2) || p.chi2 < 0.0) { all_ok = false; break; }
    CHECK(all_ok, "all chi2 finite and non-negative");

    std::cout << "\n[Test 3: minimum near truth]\n";
    double best_chi2 = 1e300, best_T = 0.0, best_bs = 0.0;
    for (const auto& p : surface)
        if (p.chi2 < best_chi2) { best_chi2 = p.chi2; best_T = p.T; best_bs = p.beta_s; }
    std::cout << "  grid minimum: T=" << best_T*1000 << " MeV, beta_s=" << best_bs << "\n";
    CHECK(std::abs(best_T - truth.T) < 0.010, "T min within 10 MeV");
    CHECK(std::abs(best_bs - truth.beta_s) < 0.030, "beta_s min within 0.03");

    std::cout << "\n[Test 4: chi2 rises away from minimum]\n";
    double dmin = 1e9; Chi2SurfacePoint p_near;
    for (const auto& p : surface) {
        const double d = std::hypot(p.T - truth.T, p.beta_s - truth.beta_s);
        if (d < dmin) { dmin = d; p_near = p; }
    }
    Chi2SurfacePoint p_far = surface.front();
    for (const auto& p : surface) if (p.T < p_far.T) p_far = p;
    CHECK(p_far.chi2 > p_near.chi2, "chi2 at edge > chi2 at minimum");

    std::cout << "\n[Test 5: file output]\n";
    const std::string out = "/tmp/test_chi2_surface.dat";
    scan_chi2_surface(species, data, bw, truth,
                      0.080, 0.150, 8,
                      0.40,  0.80, 11,
                      Statistics::Quantum, out);
    std::ifstream f(out);
    CHECK(f.good(), "output file created");
    f.close();

    std::cout << "\n[Test 6: minimum inside the box]\n";
    CHECK(best_T > 0.080 && best_T < 0.150, "T minimum strictly inside");
    CHECK(best_bs > 0.40 && best_bs < 0.80, "beta_s minimum strictly inside");

    std::cout << "\n========================================\n";
    if (failures == 0) { std::cout << " All Phase 10 tests passed.\n"; return 0; }
    std::cout << " " << failures << " test(s) failed.\n"; return 1;
}
