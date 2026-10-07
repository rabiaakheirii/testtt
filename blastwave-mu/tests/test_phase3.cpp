#include <cmath>
#include <iostream>
#include <vector>

#include "blastwave/params.hpp"
#include "blastwave/spectrum.hpp"
#include "particles/table.hpp"

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

int main() {
    std::cout << "========== Phase 3 unit tests ==========\n";

    bwmu::BWParams bw;
    bw.T       = 0.100;   // 100 MeV
    bw.beta_s  = 0.60;
    bw.R_fm    = 8.0;
    bw.n_prof  = 2;

    const auto& pi_p  = bwmu::get_species("pi+");
    const auto& pi_m  = bwmu::get_species("pi-");
    const auto& K_p   = bwmu::get_species("K+");
    const auto& K_m   = bwmu::get_species("K-");
    const auto& p     = bwmu::get_species("p");
    const auto& pbar  = bwmu::get_species("pbar");

    std::cout << "\n[BW parameters]\n";
    std::cout << "  T      = " << bw.T * 1000 << " MeV\n";
    std::cout << "  beta_s = " << bw.beta_s << "\n";
    std::cout << "  R      = " << bw.R_fm << " fm\n";
    std::cout << "  n      = " << bw.n_prof << "\n";

    // --------------------------------------------------------------
    // Test 1: at mu_B = mu_Q = mu_S = 0, Boltzmann spectrum must be
    // positive and monotonically decreasing in mT
    // --------------------------------------------------------------
    std::cout << "\n[Boltzmann: monotonically decreasing]\n";
    const double mT1 = 0.2;
    const double mT2 = 0.5;
    const double mT3 = 1.0;
    const double s1 = bwmu::dN_dmT_boltzmann(mT1, pi_p, bw, 0, 0, 0);
    const double s2 = bwmu::dN_dmT_boltzmann(mT2, pi_p, bw, 0, 0, 0);
    const double s3 = bwmu::dN_dmT_boltzmann(mT3, pi_p, bw, 0, 0, 0);
    std::cout << "  pi+ mT=0.2 : " << s1 << "\n";
    std::cout << "  pi+ mT=0.5 : " << s2 << "\n";
    std::cout << "  pi+ mT=1.0 : " << s3 << "\n";
    CHECK(s1 > 0.0, "spectrum positive at mT=0.2");
    CHECK(s1 > s2 && s2 > s3, "spectrum decreasing in mT");

    // --------------------------------------------------------------
    // Test 2: quantum and Boltzmann agree at mu = 0 for protons
    //         (proton: m >> T, so fugacity expansion converges in 1-2 terms)
    // --------------------------------------------------------------
    std::cout << "\n[Quantum vs Boltzmann: proton, mu=0]\n";
    for (double mT : {0.95, 1.2, 1.5}) {
        const double boltz = bwmu::dN_dmT_boltzmann(mT, p, bw, 0, 0, 0);
        const double quant = bwmu::dN_dmT_quantum  (mT, p, bw, 0, 0, 0);
        const double ratio = quant / boltz;
        std::cout << "  mT=" << mT
                  << "  boltz=" << boltz
                  << "  quant=" << quant
                  << "  ratio=" << ratio << "\n";
        CHECK(std::abs(ratio - 1.0) < 0.02,
              "quantum ~ boltzmann for proton at mu=0");
    }

    // --------------------------------------------------------------
    // Test 3: pi+/pi- ratio at Boltzmann level must equal
    //         exp(2 mu_Q / T)
    // --------------------------------------------------------------
    std::cout << "\n[pi+/pi- ratio at Boltzmann level]\n";
    const double mu_Q = -0.005;   // -5 MeV
    const double mT   = 0.3;
    const double pi_p_s = bwmu::dN_dmT_boltzmann(mT, pi_p, bw, 0, mu_Q, 0);
    const double pi_m_s = bwmu::dN_dmT_boltzmann(mT, pi_m, bw, 0, mu_Q, 0);
    const double ratio_model = pi_p_s / pi_m_s;
    const double ratio_expect = std::exp(2.0 * mu_Q / bw.T);
    std::cout << "  ratio model  = " << ratio_model  << "\n";
    std::cout << "  ratio expect = " << ratio_expect << "\n";
    CHECK(close(ratio_model, ratio_expect, 1e-10),
          "pi+/pi- = exp(2 mu_Q / T)");

    // --------------------------------------------------------------
    // Test 4: K+/K- ratio = exp(2 (mu_Q + mu_S) / T)
    // --------------------------------------------------------------
    std::cout << "\n[K+/K- ratio at Boltzmann level]\n";
    const double mu_S = 0.002;    // +2 MeV
    const double K_p_s = bwmu::dN_dmT_boltzmann(mT, K_p, bw, 0, mu_Q, mu_S);
    const double K_m_s = bwmu::dN_dmT_boltzmann(mT, K_m, bw, 0, mu_Q, mu_S);
    const double ratioK_model  = K_p_s / K_m_s;
    const double ratioK_expect = std::exp(2.0 * (mu_Q + mu_S) / bw.T);
    std::cout << "  ratio model  = " << ratioK_model  << "\n";
    std::cout << "  ratio expect = " << ratioK_expect << "\n";
    CHECK(close(ratioK_model, ratioK_expect, 1e-10),
          "K+/K- = exp(2(mu_Q+mu_S)/T)");

    // --------------------------------------------------------------
    // Test 5: pbar/p ratio = exp(-2 (mu_B + mu_Q) / T)
    // --------------------------------------------------------------
    std::cout << "\n[pbar/p ratio at Boltzmann level]\n";
    const double mu_B = 0.030;    // +30 MeV
    const double pbar_s = bwmu::dN_dmT_boltzmann(mT, pbar, bw, mu_B, mu_Q, 0);
    const double p_s    = bwmu::dN_dmT_boltzmann(mT, p,    bw, mu_B, mu_Q, 0);
    const double ratioP_model  = pbar_s / p_s;
    const double ratioP_expect = std::exp(-2.0 * (mu_B + mu_Q) / bw.T);
    std::cout << "  ratio model  = " << ratioP_model  << "\n";
    std::cout << "  ratio expect = " << ratioP_expect << "\n";
    CHECK(close(ratioP_model, ratioP_expect, 1e-10),
          "pbar/p = exp(-2(mu_B+mu_Q)/T)");

    // --------------------------------------------------------------
    // Test 6: at mu = 0, Boltzmann spectrum is invariant under
    //         particle <-> antiparticle for the same mass
    // --------------------------------------------------------------
    std::cout << "\n[Particle = antiparticle at mu=0]\n";
    CHECK(close(bwmu::dN_dmT_boltzmann(mT, pi_p, bw, 0, 0, 0),
                bwmu::dN_dmT_boltzmann(mT, pi_m, bw, 0, 0, 0), 1e-12),
          "pi+ = pi- at mu=0");
    CHECK(close(bwmu::dN_dmT_boltzmann(mT, K_p, bw, 0, 0, 0),
                bwmu::dN_dmT_boltzmann(mT, K_m, bw, 0, 0, 0), 1e-12),
          "K+ = K- at mu=0");

    // --------------------------------------------------------------
    // Test 7: quantum > Boltzmann for pions (Bose enhancement)
    // --------------------------------------------------------------
    std::cout << "\n[Bose enhancement for pions]\n";
    for (double mT : {0.15, 0.25, 0.50}) {
        const double boltz = bwmu::dN_dmT_boltzmann(mT, pi_p, bw, 0, 0, 0);
        const double quant = bwmu::dN_dmT_quantum  (mT, pi_p, bw, 0, 0, 0);
        const double ratio = quant / boltz;
        std::cout << "  mT=" << mT
                  << "  quant/boltz = " << ratio << "\n";
        CHECK(ratio > 1.0, "quantum > boltzmann for pions");
        if (mT > 0.5) {
            CHECK(ratio < 1.05, "enhancement < 5% at high mT");
        }
    }

    // --------------------------------------------------------------
    // Test 8: quantum spectra are positive and monotonically decreasing
    //         for all species considered
    // --------------------------------------------------------------
    std::cout << "\n[Quantum monotonic decrease]\n";
    struct Row { const char* name; const bwmu::Species* sp; };
    std::vector<Row> rows = {
        {"pi+", &pi_p}, {"K+", &K_p}, {"p", &p}
    };
    for (const auto& row : rows) {
        const double q1 = bwmu::dN_dmT_quantum(0.3, *row.sp, bw, 0.0, 0.0, 0.0);
        const double q2 = bwmu::dN_dmT_quantum(0.8, *row.sp, bw, 0.0, 0.0, 0.0);
        const double q3 = bwmu::dN_dmT_quantum(1.5, *row.sp, bw, 0.0, 0.0, 0.0);
        std::cout << "  " << row.name
                  << "  mT=0.3 : " << q1
                  << "  mT=0.8 : " << q2
                  << "  mT=1.5 : " << q3 << "\n";
        CHECK(q1 > 0 && q1 > q2 && q2 > q3,
              std::string(row.name) + " quantum decreasing");
    }

    std::cout << "\n========================================\n";
    if (failures == 0) {
        std::cout << " All Phase 3 tests passed.\n";
        return 0;
    } else {
        std::cout << " " << failures << " test(s) failed.\n";
        return 1;
    }
}
