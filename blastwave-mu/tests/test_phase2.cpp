#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "particles/table.hpp"
#include "particles/chem_potential.hpp"

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

static bool close(double a, double b, double tol = 1e-12) {
    return std::abs(a - b) < tol;
}

int main() {
    std::cout << "========== Phase 2 unit tests ==========\n";

    std::cout << "\n[Table completeness]\n";
    CHECK(bwmu::particle_table().size() == 14, "table has 14 species");

    std::cout << "\n[Quantum numbers]\n";
    CHECK(bwmu::get_species("pi+").Q == +1, "pi+ has Q = +1");
    CHECK(bwmu::get_species("pi+").S ==  0, "pi+ has S = 0");
    CHECK(bwmu::get_species("pi+").B ==  0, "pi+ has B = 0");

    CHECK(bwmu::get_species("K-").S == -1,  "K- has S = -1");
    CHECK(bwmu::get_species("K-").Q == -1,  "K- has Q = -1");

    CHECK(bwmu::get_species("p").B == +1,   "p has B = +1");
    CHECK(bwmu::get_species("pbar").B == -1,"pbar has B = -1");

    CHECK(bwmu::get_species("Lambda").S == -1, "Lambda has S = -1");
    CHECK(bwmu::get_species("Xi-").S == -2,    "Xi- has S = -2");
    CHECK(bwmu::get_species("Omega-").S == -3, "Omega- has S = -3");

    std::cout << "\n[Statistics]\n";
    CHECK(bwmu::get_species("pi+").boson == true,  "pi+ is a boson");
    CHECK(bwmu::get_species("K+").boson  == true,  "K+ is a boson");
    CHECK(bwmu::get_species("p").boson   == false, "p is a fermion");
    CHECK(bwmu::get_species("Lambda").boson == false, "Lambda is a fermion");

    std::cout << "\n[Masses]\n";
    CHECK(close(bwmu::get_species("pi+").mass, 0.139570, 1e-6), "pi mass");
    CHECK(close(bwmu::get_species("K+").mass,  0.493677, 1e-6), "K mass");
    CHECK(close(bwmu::get_species("p").mass,   0.938272, 1e-6), "p mass");

    std::cout << "\n[Chemical potentials]\n";
    const auto& pi_minus = bwmu::get_species("pi-");
    const auto& K_minus  = bwmu::get_species("K-");
    const auto& pbar     = bwmu::get_species("pbar");

    // mu(pi-) = -mu_Q
    CHECK(close(bwmu::mu_of(pi_minus, 0.0, -0.005, 0.002), +0.005),
          "mu(pi-) = -mu_Q = +5 MeV");

    // mu(K-) = -mu_Q - mu_S
    CHECK(close(bwmu::mu_of(K_minus, 0.0, -0.005, 0.002), +0.003),
          "mu(K-) = -mu_Q - mu_S = +3 MeV");

    // mu(pbar) = -mu_B - mu_Q
    CHECK(close(bwmu::mu_of(pbar, 0.030, -0.005, 0.0), -0.025),
          "mu(pbar) = -mu_B - mu_Q = -25 MeV");

    // All-zero case
    CHECK(close(bwmu::mu_of(bwmu::get_species("pi+"), 0.0, 0.0, 0.0), 0.0),
          "mu(pi+) = 0 when all potentials vanish");

    std::cout << "\n[Lookup]\n";
    bool threw = false;
    try { bwmu::get_species("not-a-particle"); }
    catch (const std::runtime_error&) { threw = true; }
    CHECK(threw, "unknown species throws");

    std::cout << "\n[Table print]\n";
    bwmu::print_particle_table();

    std::cout << "\n========================================\n";
    if (failures == 0) {
        std::cout << " All Phase 2 tests passed.\n";
        return 0;
    } else {
        std::cout << " " << failures << " test(s) failed.\n";
        return 1;
    }
}
