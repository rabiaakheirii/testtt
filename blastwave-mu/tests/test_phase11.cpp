#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "fit/ptbin_scan.hpp"
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
                           const std::vector<double>& pTs)
{
    Spectrum s; s.species = sp.name; s.centrality = 0;
    for (double pT : pTs) {
        const double mT = std::sqrt(pT*pT + sp.mass*sp.mass);
        const double mT_m0 = mT - sp.mass;
        const double f = model_shape(mT_m0, sp, bw, fp, Statistics::Quantum);
        DataPoint d;
        d.mT    = mT_m0;
        d.value = A * f;
        d.err   = 1e-4 * std::abs(A * f) + 1e-9;
        s.points.push_back(d);
    }
    return s;
}

int main() {
    std::cout << "========== Phase 11 unit tests ==========\n\n";

    BWParams bw; bw.R_fm = 12.0; bw.n_prof = 2;
    FitParams truth;
    truth.T = 0.105; truth.beta_s = 0.62;
    truth.mu_B = 0.025; truth.mu_Q = -0.007; truth.mu_S = 0.003;

    const auto& pi_p = get_species("pi+");
    const auto& pi_m = get_species("pi-");
    const auto& K_p  = get_species("K+");
    const auto& K_m  = get_species("K-");
    const auto& p    = get_species("p");

    // Wide pT grids
    std::vector<double> pTs_pi, pTs_K, pTs_p;
    for (double x = 0.20; x <= 2.00; x += 0.10) pTs_pi.push_back(x);
    for (double x = 0.30; x <= 1.80; x += 0.10) pTs_K.push_back(x);
    for (double x = 0.50; x <= 2.50; x += 0.10) pTs_p.push_back(x);

    std::vector<Spectrum> data = {
        make_synth(pi_p, bw, truth, 10.0, pTs_pi),
        make_synth(pi_m, bw, truth, 10.0, pTs_pi),
        make_synth(K_p,  bw, truth,  1.5, pTs_K),
        make_synth(K_m,  bw, truth,  1.5, pTs_K),
        make_synth(p,    bw, truth,  0.3, pTs_p)
    };
    std::vector<Species> species = { pi_p, pi_m, K_p, K_m, p };

    FitBounds bounds;
    bounds.T_min = 0.080; bounds.T_max = 0.150;
    bounds.bs_min = 0.40; bounds.bs_max = 0.85;
    bounds.muB_min = 0.000; bounds.muB_max = 0.050;
    bounds.muQ_min = -0.020; bounds.muQ_max = 0.020;
    bounds.muS_min = -0.015; bounds.muS_max = 0.015;

    std::vector<std::pair<double,double>> windows = {
        {0.20, 0.80},
        {0.30, 1.00},
        {0.40, 1.20}
    };

    std::cout << "[Test 1: scan returns correct number of windows]\n";
    auto results = scan_pT_windows(species, data, bw, bounds,
                                   windows, Statistics::Quantum, false);
    CHECK(results.size() == windows.size(), "3 windows processed");

    std::cout << "\n[Test 2: all windows converge]\n";
    bool all_converged = true;
    for (const auto& r : results) {
        if (r.chi2_ndof <= 0.0 || r.params.T <= 0.0) {
            all_converged = false;
        }
    }
    CHECK(all_converged, "all windows have valid fit");

    std::cout << "\n[Test 3: parameters within physical range]\n";
    bool all_physical = true;
    for (const auto& r : results) {
        if (r.params.T < 0.05 || r.params.T > 0.30) all_physical = false;
        if (r.params.beta_s < 0.0 || r.params.beta_s > 1.0) all_physical = false;
    }
    CHECK(all_physical, "all T in [50, 300] MeV, beta_s in [0, 1]");

    std::cout << "\n[Test 4: output file written]\n";
    const std::string out = "/tmp/test_ptbin.dat";
    write_ptbin_scan(results, out);
    std::ifstream f(out);
    CHECK(f.good(), "output file created");
    f.close();

    std::cout << "\n[Test 5: T stability across pT windows]\n";
    // Print T values
    for (const auto& r : results) {
        std::cout << "  pT[" << r.pT_min << "," << r.pT_max
                  << "]: T = " << r.params.T * 1000 << " MeV\n";
    }
    double T_min = 1e9, T_max = -1e9;
    for (const auto& r : results) {
        T_min = std::min(T_min, r.params.T);
        T_max = std::max(T_max, r.params.T);
    }
    const double T_spread_MeV = (T_max - T_min) * 1000;
    std::cout << "  T spread: " << T_spread_MeV << " MeV\n";
    CHECK(T_spread_MeV < 50.0, "T stable within 50 MeV across pT windows");

    std::cout << "\n========================================\n";
    if (failures == 0) { std::cout << " All Phase 11 tests passed.\n"; return 0; }
    std::cout << " " << failures << " test(s) failed.\n"; return 1;
}
