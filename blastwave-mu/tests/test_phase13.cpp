#include <cmath>
#include <iostream>
#include <vector>

#include "analysis/dndy.hpp"
#include "blastwave/params.hpp"
#include "particles/table.hpp"

using namespace bwmu;
static int failures = 0;
#define CHECK(cond, msg) \
    do { if (!(cond)) { std::cout << "  FAIL: " << msg << "\n"; ++failures; } \
         else { std::cout << "  ok  : " << msg << "\n"; } } while (0)

int main() {
    std::cout << "========== Phase 13 unit tests ==========\n\n";

    BWParams bw; bw.R_fm = 12.0; bw.n_prof = 2;
    FitParams fp;
    fp.T = 0.150; fp.beta_s = 0.62;
    fp.mu_B = 0.025; fp.mu_Q = -0.005; fp.mu_S = 0.003;

    const double eta_max = 2.0;

    std::cout << "[Test 1: dN/dy returns positive values]\n";
    const auto& pi_p = get_species("pi+");
    const double v0 = dN_dy(0.0, pi_p, bw, fp, eta_max, 1.0);
    const double v1 = dN_dy(1.0, pi_p, bw, fp, eta_max, 1.0);
    const double v2 = dN_dy(2.0, pi_p, bw, fp, eta_max, 1.0);
    std::cout << "  dN/dy(y=0) = " << v0 << "\n";
    std::cout << "  dN/dy(y=1) = " << v1 << "\n";
    std::cout << "  dN/dy(y=2) = " << v2 << "\n";
    CHECK(v0 > 0 && v1 > 0 && v2 > 0, "all positive");
    CHECK(v0 >= v1 && v1 >= v2, "monotonically decreasing");

    std::cout << "\n[Test 2: symmetry — dN/dy is even in y]\n";
    const double vp = dN_dy(1.5, pi_p, bw, fp, eta_max, 1.0);
    const double vm = dN_dy(-1.5, pi_p, bw, fp, eta_max, 1.0);
    std::cout << "  dN/dy(+1.5) = " << vp << "\n";
    std::cout << "  dN/dy(-1.5) = " << vm << "\n";
    CHECK(std::abs(vp - vm) / vp < 1e-6, "even in y");

    std::cout << "\n[Test 3: at y=0, larger eta_max gives larger value]\n";
    const double v_small = dN_dy(0.0, pi_p, bw, fp, 1.0, 1.0);
    const double v_large = dN_dy(0.0, pi_p, bw, fp, 3.0, 1.0);
    std::cout << "  eta_max=1.0: " << v_small << "\n";
    std::cout << "  eta_max=3.0: " << v_large << "\n";
    CHECK(v_large > v_small, "larger eta_max -> larger y=0 yield");

    std::cout << "\n[Test 4: particle/antiparticle ratio correct]\n";
    const auto& pi_m = get_species("pi-");
    const double r_pi_p = dN_dy(0.0, pi_p, bw, fp, eta_max, 1.0);
    const double r_pi_m = dN_dy(0.0, pi_m, bw, fp, eta_max, 1.0);
    const double expected = std::exp(2.0 * fp.mu_Q / fp.T);
    std::cout << "  pi+/pi- = " << r_pi_p / r_pi_m
              << "  (expected " << expected << ")\n";
    CHECK(std::abs(r_pi_p / r_pi_m - expected) / expected < 1e-6,
          "pi+/pi- = exp(2 mu_Q / T)");

    std::cout << "\n[Test 5: curve generation]\n";
    auto curve = compute_dndy_curve(pi_p, bw, fp, eta_max,
                                    -3.0, 3.0, 31, 1.0);
    CHECK(curve.size() == 31, "31 points");
    CHECK(curve.front().first == -3.0, "starts at -3");
    CHECK(curve.back().first  ==  3.0, "ends at 3");

    std::cout << "\n========================================\n";
    if (failures == 0) { std::cout << " All Phase 13 tests passed.\n"; return 0; }
    std::cout << " " << failures << " test(s) failed.\n"; return 1;
}
