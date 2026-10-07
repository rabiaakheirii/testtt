#include "fit/chi2_surface.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace bwmu {

    namespace {
        std::vector<int> build_norm_groups(const std::vector<Species>& species) {
            std::vector<int>    g(species.size(), -1);
            std::vector<double> seen_mass;
            for (size_t i = 0; i < species.size(); ++i) {
                int found = -1;
                for (size_t k = 0; k < seen_mass.size(); ++k) {
                    if (std::abs(seen_mass[k] - species[i].mass) < 1e-6) {
                        found = static_cast<int>(k);
                        break;
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
    }

    std::vector<Chi2SurfacePoint> compute_chi2_surface(
        const std::vector<Species>& species,
        const std::vector<Spectrum>& data,
        const BWParams& bw,
        const FitParams& best,
        double T_min, double T_max, int n_T,
        double b_min, double b_max, int n_b,
        Statistics stats)
    {
        std::vector<Chi2SurfacePoint> out;
        out.reserve(n_T * n_b);
        const auto groups = build_norm_groups(species);

        for (int i = 0; i < n_T; ++i) {
            const double T = (n_T > 1)
                ? T_min + (T_max - T_min) * i / (n_T - 1)
                : 0.5 * (T_min + T_max);
            for (int j = 0; j < n_b; ++j) {
                const double bs = (n_b > 1)
                    ? b_min + (b_max - b_min) * j / (n_b - 1)
                    : 0.5 * (b_min + b_max);
                FitParams fp = best;
                fp.T      = T;
                fp.beta_s = bs;
                auto r = compute_chi2(species, data, bw, fp, stats, groups);
                Chi2SurfacePoint p;
                p.T         = T;
                p.beta_s    = bs;
                p.chi2      = r.chi2;
                p.chi2_ndof = r.chi2_ndof;
                out.push_back(p);
            }
        }
        return out;
    }

    void scan_chi2_surface(
        const std::vector<Species>& species,
        const std::vector<Spectrum>& data,
        const BWParams& bw,
        const FitParams& best,
        double T_min, double T_max, int n_T,
        double b_min, double b_max, int n_b,
        Statistics stats,
        const std::string& out_file)
    {
        auto points = compute_chi2_surface(species, data, bw, best,
                                           T_min, T_max, n_T,
                                           b_min, b_max, n_b,
                                           stats);
        std::ofstream out(out_file);
        if (!out) {
            std::cerr << "[chi2_surface] cannot write " << out_file << "\n";
            return;
        }
        out << std::fixed << std::setprecision(6);
        out << "# T[GeV]  beta_s  chi2  chi2_ndof\n";
        for (const auto& p : points) {
            out << p.T         << "  "
                << p.beta_s    << "  "
                << p.chi2      << "  "
                << p.chi2_ndof << "\n";
        }
        std::cout << "[chi2_surface] wrote " << out_file
                  << "  (" << points.size() << " points)\n";
    }

} // namespace bwmu
