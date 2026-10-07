#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "blastwave/feeddown.hpp"
#include "blastwave/params.hpp"
#include "fit/chi2.hpp"
#include "particles/table.hpp"

using namespace bwmu;

int main() {
    std::cout << "========== Phase 8 unit tests ==========\n\n";

    auto res = load_resonances("data/resonances/resonances_QS.csv");
    std::cout << "[test] loaded " << res.size() << " resonances\n";

    int failures = 0;
    if (res.size() < 20) {
        std::cout << "FAIL: too few resonances loaded (need 20+)\n";
        ++failures;
    } else {
        std::cout << "ok  : resonance table loaded\n";
    }

    // Parameters for test
    BWParams bw; bw.R_fm = 12.0; bw.n_prof = 2;
    FitParams fp;
    fp.T = 0.100; fp.beta_s = 0.60;
    fp.mu_B = 0.025; fp.mu_Q = -0.007; fp.mu_S = 0.003;

    const auto& pi_p = get_species("pi+");
    const std::vector<double> mTs = {0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 0.80, 1.20, 1.80};

    // ---- Direct-only ----
    set_resonances({});
    std::vector<double> direct;
    for (double mT : mTs)
        direct.push_back(model_shape(mT, pi_p, bw, fp, Statistics::Quantum));

    // ---- With feed-down ----
    set_resonances(res);
    std::vector<double> total;
    for (double mT : mTs)
        total.push_back(model_shape(mT, pi_p, bw, fp, Statistics::Quantum));

    std::cout << "\nmT[GeV]   direct       +decay       ratio\n";
    std::cout << "---------------------------------------------------\n";
    for (size_t i = 0; i < mTs.size(); ++i) {
        const double r = direct[i] > 0 ? total[i] / direct[i] : 0.0;
        std::cout << std::fixed << std::setprecision(2)
                  << std::setw(6) << mTs[i] << "   "
                  << std::scientific << std::setprecision(3)
                  << std::setw(11) << direct[i] << "  "
                  << std::setw(11) << total[i]  << "  "
                  << std::fixed << std::setprecision(3) << r << "\n";
    }

    // Physics check: at low mT, feed-down enhances pions by >30%
    const double ratio_low = total[1] / direct[1];
    if (ratio_low > 1.3) {
        std::cout << "ok  : low-mT enhancement > 30%\n";
    } else {
        std::cout << "FAIL: low-mT enhancement too small (" << ratio_low << ")\n";
        ++failures;
    }

    // Physics check: at high mT, feed-down contributes < 50%
    const double ratio_high = total[8] / direct[8];
    if (ratio_high < 1.5) {
        std::cout << "ok  : high-mT enhancement small\n";
    } else {
        std::cout << "FAIL: high-mT enhancement too large (" << ratio_high << ")\n";
        ++failures;
    }

    // Pion/antipion ratio should now include resonances with their own charges
    const auto& pi_m = get_species("pi-");
    const double pi_p_tot = model_shape(0.15, pi_p, bw, fp, Statistics::Quantum);
    const double pi_m_tot = model_shape(0.15, pi_m, bw, fp, Statistics::Quantum);
    const double r = pi_p_tot / pi_m_tot;
    std::cout << "\npi+/pi- at mT=0.15 with feed-down: " << r << "\n";
    std::cout << "(expected near exp(2 mu_Q/T) = "
              << std::exp(2*fp.mu_Q/fp.T) << ")\n";
    if (r > 0.7 && r < 1.3) {
        std::cout << "ok  : pi+/pi- ratio reasonable\n";
    } else {
        std::cout << "FAIL: pi+/pi- ratio out of range\n";
        ++failures;
    }

    std::cout << "\n========================================\n";
    if (failures == 0) {
        std::cout << " All Phase 8 tests passed.\n";
        return 0;
    }
    std::cout << " " << failures << " test(s) failed.\n";
    return 1;
}
