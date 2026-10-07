#include "analysis/critical_point.hpp"
#include "utils/random.hpp"
#include <cmath>
#include <fstream>
#include <iostream>

namespace qcd {

    ScanPoint simulate_event_sample(double sqrt_s,
                                    int n_events,
                                    unsigned seed)
    {
        ScanPoint sp;
        sp.sqrt_s = sqrt_s;

        // phenomenological mu_B(sqrt_s) in GeV
        sp.mu_B = 1.31 / (1.0 + 0.27 * sqrt_s);

        // Distance from critical point (mu_CEP ~ 0.4 GeV at sqrt_s ~ 7 GeV)
        const double mu_CEP = 0.40;
        const double dmu    = (sp.mu_B - mu_CEP) / 0.15;
        sp.cs2 = 0.05 + 0.28 * std::tanh(dmu * dmu) + 0.02 * std::exp(-dmu*dmu);

        // Build toy net-proton distribution:
        //  - Gaussian with mean ~ mu_B * A (A = 20)
        //  - width sigma = sqrt(C2)  where C2 scales like 1/cs2 near CEP
        const double mean  = 20.0 * sp.mu_B;
        const double sigma = 4.0 + 8.0 * std::exp(-dmu*dmu);

        RNG rng(seed);
        std::vector<int> N(n_events);
        for (int i = 0; i < n_events; ++i) {
            // Box-Muller
            const double u1 = std::max(rng.uniform(), 1e-12);
            const double u2 = rng.uniform();
            const double z  = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
            N[i] = static_cast<int>(std::round(mean + sigma * z));
        }

        sp.cum = compute_cumulants(N);
        return sp;
    }

    void scan_beam_energy(double sqrt_s_min, double sqrt_s_max, int n_steps,
                          int n_events_per_point,
                          unsigned seed,
                          const std::string& output_file)
    {
        std::ofstream out(output_file);
        out << "# sqrt_s[GeV]  mu_B[GeV]  cs2  C1  C2  C3  C4  Ssigma  kappa_sigma2\n";

        const double dlog = (std::log(sqrt_s_max) - std::log(sqrt_s_min)) / (n_steps - 1);

        for (int i = 0; i < n_steps; ++i) {
            const double sqrt_s = std::exp(std::log(sqrt_s_min) + i * dlog);
            auto sp = simulate_event_sample(sqrt_s, n_events_per_point,
                                            seed + i * 1000);
            out << sp.sqrt_s << "  " << sp.mu_B << "  " << sp.cs2 << "  "
                << sp.cum.C1 << "  " << sp.cum.C2 << "  "
                << sp.cum.C3 << "  " << sp.cum.C4 << "  "
                << sp.cum.Ssigma << "  " << sp.cum.kappa_sigma2 << "\n";

            std::cout << "[cumulants] sqrt_s=" << sp.sqrt_s
                      << "  mu_B=" << sp.mu_B
                      << "  cs2="   << sp.cs2
                      << "  Ssigma=" << sp.cum.Ssigma
                      << "  kappa_sigma2=" << sp.cum.kappa_sigma2 << "\n";
        }
        std::cout << "[cumulants] wrote " << output_file << "\n";
    }
}
