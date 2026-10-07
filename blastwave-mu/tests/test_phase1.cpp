#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "utils/bessel.hpp"
#include "utils/integrate.hpp"
#include "utils/grid.hpp"
#include "utils/random.hpp"

static int failures = 0;

#define CHECK(cond, msg)                                        \
    do {                                                        \
        if (!(cond)) {                                          \
            std::cout << "  FAIL: " << msg << "\n";             \
            ++failures;                                         \
        } else {                                                \
            std::cout << "  ok  : " << msg << "\n";             \
        }                                                       \
    } while (0)

static bool close(double a, double b, double tol) {
    return std::abs(a - b) < tol;
}

int main() {
    std::cout << "========== Phase 1 unit tests ==========\n";

    // --- Bessel ---
    std::cout << "\n[Bessel functions]\n";
    CHECK(close(bwmu::K1(1.0),  0.601907, 1e-5), "K1(1.0) = 0.601907");
    CHECK(close(bwmu::K0(1.0),  0.421024, 1e-5), "K0(1.0) = 0.421024");
    CHECK(close(bwmu::I0(1.0),  1.266066, 1e-5), "I0(1.0) = 1.266066");
    CHECK(close(bwmu::I1(1.0),  0.565159, 1e-5), "I1(1.0) = 0.565159");

    // --- Trapezoid / Simpson ---
    std::cout << "\n[Integration]\n";
    // integral of sin(x) from 0 to pi = 2
    const int n = 1001;
    bwmu::Grid1D g(0.0, M_PI, n);
    std::vector<double> f(n);
    for (int i = 0; i < n; ++i) f[i] = std::sin(g[i]);

    const double trap = bwmu::trapezoid(f, g.x);
    const double simp = bwmu::simpson  (f, g.x);
    std::cout << "  trapezoid = " << trap << "\n";
    std::cout << "  simpson   = " << simp << "\n";
    CHECK(close(trap, 2.0, 1e-4), "trapezoid of sin on [0,pi] = 2");
    CHECK(close(simp, 2.0, 1e-8), "simpson   of sin on [0,pi] = 2");

    // --- Grid ---
    std::cout << "\n[Grid]\n";
    bwmu::Grid1D g2(0.0, 1.0, 11);
    CHECK(g2.n == 11,        "grid has 11 points");
    CHECK(close(g2[0],  0.0, 1e-12), "grid[0]  = 0");
    CHECK(close(g2[10], 1.0, 1e-12), "grid[10] = 1");
    CHECK(close(g2.dx(), 0.1, 1e-12), "grid spacing = 0.1");

    // --- RNG ---
    std::cout << "\n[RNG]\n";
    bwmu::RNG rng(42);
    const double u1 = rng.uniform();
    const double u2 = rng.uniform();
    CHECK(u1 >= 0.0 && u1 < 1.0, "uniform() in [0,1)");
    CHECK(u2 >= 0.0 && u2 < 1.0, "uniform() in [0,1)");
    CHECK(u1 != u2,              "consecutive draws differ");

    // Reproducibility: same seed -> same sequence
    bwmu::RNG rA(7), rB(7);
    bool same = true;
    for (int i = 0; i < 100; ++i)
        if (std::abs(rA.uniform() - rB.uniform()) > 1e-15) same = false;
    CHECK(same, "same seed reproduces sequence");

    // Gaussian mean and stddev over many samples
    bwmu::RNG rG(123);
    const int N = 200000;
    double s = 0.0, s2 = 0.0;
    for (int i = 0; i < N; ++i) {
        const double x = rG.gaussian(1.5, 2.0);
        s  += x;
        s2 += x * x;
    }
    const double mean = s / N;
    const double var  = s2 / N - mean * mean;
    std::cout << "  gaussian mean = " << mean << "\n";
    std::cout << "  gaussian std  = " << std::sqrt(var) << "\n";
    CHECK(close(mean,             1.5, 0.02), "gaussian mean ~ 1.5");
    CHECK(close(std::sqrt(var),   2.0, 0.02), "gaussian std  ~ 2.0");

    std::cout << "\n========================================\n";
    if (failures == 0) {
        std::cout << " All Phase 1 tests passed.\n";
        return 0;
    } else {
        std::cout << " " << failures << " test(s) failed.\n";
        return 1;
    }
}
