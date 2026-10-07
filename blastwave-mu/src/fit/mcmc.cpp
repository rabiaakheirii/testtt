#include "fit/mcmc.hpp"
#include "utils/random.hpp"
#include "fit/tbeta_penalty.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

namespace bwmu {

    namespace {

        std::vector<int> build_norm_groups(const std::vector<Species>& species) {
            std::vector<int>    g(species.size(), -1);
            std::vector<double> seen_mass;
            for (size_t i = 0; i < species.size(); ++i) {
                int found = -1;
                for (size_t k = 0; k < seen_mass.size(); ++k) {
                    if (std::abs(seen_mass[k] - species[i].mass) < 1e-6) {
                        found = static_cast<int>(k); break;
                    }
                }
                if (found < 0) {
                    found = static_cast<int>(seen_mass.size());
                    seen_mass.push_back(species[i].mass);
                }
                g[i] = found;
            }
            return g;
        }

        bool inside(const std::vector<double>& x, const FitBounds& b) {
            return x[0] >= b.T_min   && x[0] <= b.T_max
                && x[1] >= b.bs_min  && x[1] <= b.bs_max
                && x[2] >= b.muB_min && x[2] <= b.muB_max
                && x[3] >= b.muQ_min && x[3] <= b.muQ_max
                && x[4] >= b.muS_min && x[4] <= b.muS_max;
        }

        FitParams unpack(const std::vector<double>& x) {
            FitParams fp;
            fp.T = x[0]; fp.beta_s = x[1];
            fp.mu_B = x[2]; fp.mu_Q = x[3]; fp.mu_S = x[4];
            return fp;
        }

        std::vector<double> pack(const FitParams& fp) {
            return { fp.T, fp.beta_s, fp.mu_B, fp.mu_Q, fp.mu_S };
        }
    }

    McmcChain run_mcmc(const std::vector<Species>& species,
                       const std::vector<Spectrum>& data,
                       const BWParams& bw,
                       const FitBounds& bounds,
                       const JointFitResult& start,
                       Statistics stats,
                       int n_steps, int burn_in, int thinning,
                       double step_T, double step_bs,
                       double step_muB, double step_muQ, double step_muS,
                       unsigned seed)
    {
        McmcChain chain;
        chain.param_names = {"T", "beta_s", "mu_B", "mu_Q", "mu_S"};
        chain.n_total = n_steps;

        const auto groups = build_norm_groups(species);

        constexpr double LAMBDA_TBETA = 1000.0;

        auto log_post = [&](const std::vector<double>& x) -> double {
            if (!inside(x, bounds)) return -std::numeric_limits<double>::infinity();
            FitParams fp = unpack(x);
            auto r = compute_chi2(species, data, bw, fp, stats, groups);
            if (r.ndof <= 0) return -std::numeric_limits<double>::infinity();
            const double penalty = tbeta_penalty(x[0], x[1], LAMBDA_TBETA);
            return -0.5 * (r.chi2 + penalty);
        };

        RNG rng(seed);
        std::vector<double> x = pack(start.params);
        double lp = log_post(x);

        if (!std::isfinite(lp)) {
            std::cerr << "[mcmc] starting point rejected\n";
            return chain;
        }

        const std::vector<double> step = { step_T, step_bs, step_muB, step_muQ, step_muS };
        int accepted = 0;

        const int n_keep = std::max(0, (n_steps - burn_in) / thinning);
        chain.samples.reserve(n_keep);
        chain.log_post.reserve(n_keep);

        for (int i = 0; i < n_steps; ++i) {
            std::vector<double> x_new(5);
            for (int k = 0; k < 5; ++k)
                x_new[k] = x[k] + step[k] * rng.gaussian();

            const double lp_new = log_post(x_new);
            double log_r = lp_new - lp;
            bool accept = false;
            if (std::isfinite(log_r)) {
                if (log_r >= 0.0) accept = true;
                else accept = (std::log(rng.uniform()) < log_r);
            }

            if (accept) { x = x_new; lp = lp_new; ++accepted; }

            if (i >= burn_in && (i - burn_in) % thinning == 0) {
                chain.samples.push_back(x);
                chain.log_post.push_back(lp);
            }

            if (i > 0 && i % 2000 == 0) {
                std::cout << "  [mcmc] " << i << " / " << n_steps
                          << "  accept = " << std::fixed << std::setprecision(2)
                          << (100.0 * accepted / (i + 1)) << "%\n";
            }
        }

        chain.n_accepted = accepted;
        chain.acceptance_rate = double(accepted) / n_steps;

        std::cout << "[mcmc] " << n_steps << " steps, "
                  << "acceptance = " << std::fixed << std::setprecision(3)
                  << chain.acceptance_rate
                  << ", kept " << chain.samples.size() << " samples\n";
        return chain;
    }

    void write_mcmc_chain(const McmcChain& chain, const std::string& out_file) {
        std::ofstream out(out_file);
        if (!out) return;
        out << std::fixed << std::setprecision(8);
        out << "# T[GeV]  beta_s  mu_B[GeV]  mu_Q[GeV]  mu_S[GeV]  log_post\n";
        for (size_t i = 0; i < chain.samples.size(); ++i) {
            for (double v : chain.samples[i]) out << v << "  ";
            out << chain.log_post[i] << "\n";
        }
        std::cout << "[mcmc] wrote " << out_file << "\n";
    }

    void summarize_mcmc(const McmcChain& chain, const std::string& out_file) {
        std::ofstream out(out_file);
        if (!out) return;
        out << "# param  mean  std  median  q16  q84\n";
        out << std::fixed << std::setprecision(6);

        const int n = static_cast<int>(chain.samples.size());
        if (n == 0) { out << "# empty chain\n"; return; }

        for (int k = 0; k < 5; ++k) {
            std::vector<double> vals(n);
            for (int i = 0; i < n; ++i) vals[i] = chain.samples[i][k];

            double sum = 0.0;
            for (double v : vals) sum += v;
            const double mean = sum / n;

            double var = 0.0;
            for (double v : vals) var += (v - mean) * (v - mean);
            const double sd = std::sqrt(var / (n - 1));

            std::sort(vals.begin(), vals.end());
            const double median = vals[n / 2];
            const double q16    = vals[static_cast<int>(0.16 * n)];
            const double q84    = vals[static_cast<int>(0.84 * n)];

            out << chain.param_names[k] << "  " << mean << "  " << sd << "  "
                << median << "  " << q16 << "  " << q84 << "\n";
        }
        std::cout << "[mcmc] wrote " << out_file << "\n";
    }

} // namespace bwmu
