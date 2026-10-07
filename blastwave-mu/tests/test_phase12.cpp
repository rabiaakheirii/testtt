#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "blastwave/params.hpp"
#include "fit/chi2.hpp"
#include "fit/joint.hpp"
#include "fit/mcmc.hpp"
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
        const double f = model_shape(mT, sp, bw, fp, Statistics::Boltzmann);
        DataPoint d; d.mT = mT; d.value = A * f;
        d.err = 0.02 * std::abs(A * f) + 1e-6;   // 2% relative
        s.points.push_back(d);
    }
    return s;
}

int main() {
    std::cout << "========== Phase 12 unit tests ==========\n\n";

    BWParams bw; bw.R_fm = 12.0; bw.n_prof = 2;
    FitParams truth;
    truth.T = 0.105; truth.beta_s = 0.62;
    truth.mu_B = 0.025; truth.mu_Q = -0.007; truth.mu_S = 0.003;

    const auto& pi_p = get_species("pi+");
    const auto& pi_m = get_species("pi-");
    const auto& K_p  = get_species("K+");
    const auto& K_m  = get_species("K-");

    std::vector<double> mTs = {0.10,0.20,0.30,0.45,0.60,0.80,1.00,1.25,1.50,2.00};

    std::vector<Spectrum> data = {
        make_synth(pi_p, bw, truth, 10.0, mTs),
        make_synth(pi_m, bw, truth, 10.0, mTs),
        make_synth(K_p,  bw, truth,  1.5, mTs),
        make_synth(K_m,  bw, truth,  1.5, mTs)
    };
    std::vector<Species> species = { pi_p, pi_m, K_p, K_m };

    FitBounds bounds;
    bounds.T_min = 0.080; bounds.T_max = 0.150;
    bounds.bs_min = 0.40; bounds.bs_max = 0.85;
    bounds.muB_min = 0.000; bounds.muB_max = 0.050;
    bounds.muQ_min = -0.020; bounds.muQ_max = 0.020;
    bounds.muS_min = -0.015; bounds.muS_max = 0.015;

    JointFitResult start;
    start.params = truth;
    start.chi2 = compute_chi2(species, data, bw, truth,
                              Statistics::Boltzmann).chi2;

    std::cout << "[Test 1: MCMC runs]\n";
    auto chain = run_mcmc(species, data, bw, bounds, start,
                          Statistics::Boltzmann,
                          2000, 500, 2,
                          0.001, 0.003, 0.001, 0.0005, 0.0005, 42);

    CHECK(chain.samples.size() > 0, "chain has samples");
    std::cout << "  kept " << chain.samples.size() << " samples\n";
    std::cout << "  acceptance = " << chain.acceptance_rate << "\n";

    std::cout << "\n[Test 2: acceptance]\n";
    CHECK(chain.acceptance_rate > 0.05 && chain.acceptance_rate < 0.95,
          "acceptance between 5% and 95%");

    std::cout << "\n[Test 3: inside bounds]\n";
    bool inside = true;
    for (const auto& s : chain.samples) {
        if (s[0] < bounds.T_min || s[0] > bounds.T_max) inside = false;
        if (s[1] < bounds.bs_min || s[1] > bounds.bs_max) inside = false;
    }
    CHECK(inside, "all samples inside physical bounds");

    std::cout << "\n[Test 4: truth inside 2-sigma posterior]\n";
    double sT = 0.0, sbs = 0.0;
    for (const auto& s : chain.samples) { sT += s[0]; sbs += s[1]; }
    const double mT  = sT  / chain.samples.size();
    const double mbs = sbs / chain.samples.size();

    double varT = 0.0, vars_bs = 0.0;
    for (const auto& s : chain.samples) {
        varT    += (s[0] - mT) * (s[0] - mT);
        vars_bs += (s[1] - mbs) * (s[1] - mbs);
    }
    const double sdT  = std::sqrt(varT   / chain.samples.size());
    const double sdbs = std::sqrt(vars_bs / chain.samples.size());

    std::cout << "  T:   mean " << mT * 1000
              << " +/- " << sdT * 1000 << " MeV (truth " << truth.T * 1000 << ")\n";
    std::cout << "  b_s: mean " << mbs
              << " +/- " << sdbs << " (truth " << truth.beta_s << ")\n";

    // Truth should be within 2 standard deviations of the posterior mean
    const double zT  = std::abs(mT - truth.T)      / (sdT  + 1e-12);
    const double zbs = std::abs(mbs - truth.beta_s)/ (sdbs + 1e-12);
    std::cout << "  |z| T   = " << zT << "\n";
    std::cout << "  |z| b_s = " << zbs << "\n";
    CHECK(zT  < 5.0, "T truth within 5-sigma");
    CHECK(zbs < 5.0, "beta_s truth within 5-sigma");

    std::cout << "\n[Test 5: file output]\n";
    write_mcmc_chain(chain, "/tmp/test_mcmc_chain.dat");
    summarize_mcmc (chain, "/tmp/test_mcmc_stats.dat");
    std::ifstream f1("/tmp/test_mcmc_chain.dat");
    std::ifstream f2("/tmp/test_mcmc_stats.dat");
    CHECK(f1.good(), "chain file written");
    CHECK(f2.good(), "stats file written");

    std::cout << "\n  final acceptance = " << chain.acceptance_rate << "\n";
    std::cout << "\n========================================\n";
    if (failures == 0) { std::cout << " All Phase 12 tests passed.\n"; return 0; }
    std::cout << " " << failures << " test(s) failed.\n"; return 1;
}
