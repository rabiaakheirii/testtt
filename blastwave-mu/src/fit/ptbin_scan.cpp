#include "fit/ptbin_scan.hpp"
#include "particles/table.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace bwmu {

    namespace {

        // Filter a spectrum by pT range, converting (mT - m0) to pT.
        Spectrum filter_by_pT(const Spectrum& s,
                              const Species& sp,
                              double pT_min, double pT_max)
        {
            Spectrum out;
            out.species    = s.species;
            out.centrality = s.centrality;
            for (const auto& d : s.points) {
                const double mT = d.mT + sp.mass;
                const double pT2 = mT * mT - sp.mass * sp.mass;
                if (pT2 < 0.0) continue;
                const double pT = std::sqrt(pT2);
                if (pT >= pT_min && pT <= pT_max) {
                    out.points.push_back(d);
                }
            }
            return out;
        }

    } // anonymous namespace

    std::vector<PtBinResult> scan_pT_windows(
        const std::vector<Species>& species,
        const std::vector<Spectrum>& data,
        const BWParams& bw,
        const FitBounds& bounds,
        const std::vector<std::pair<double,double>>& windows,
        Statistics stats,
        bool verbose)
    {
        std::vector<PtBinResult> results;

        for (size_t w = 0; w < windows.size(); ++w) {
            const double pT_min = windows[w].first;
            const double pT_max = windows[w].second;

            // Filter every species's data to this pT window
            std::vector<Spectrum> filtered;
            int n_points_total = 0;
            for (size_t i = 0; i < species.size(); ++i) {
                const Spectrum* src = nullptr;
                for (const auto& s : data) {
                    if (s.species == species[i].name) { src = &s; break; }
                }
                if (src == nullptr) continue;

                Spectrum f = filter_by_pT(*src, species[i], pT_min, pT_max);
                if (f.points.size() < 2) continue;   // need at least 2 pts

                filtered.push_back(f);
                n_points_total += static_cast<int>(f.points.size());
            }

            if (verbose) {
                std::cout << "\n[ptbin " << w << "] pT = ["
                          << pT_min << ", " << pT_max << "]  species="
                          << filtered.size()
                          << "  points=" << n_points_total << "\n";
            }

            PtBinResult res;
            res.window_index = static_cast<int>(w);
            res.pT_min       = pT_min;
            res.pT_max       = pT_max;
            res.n_points     = n_points_total;

            if (filtered.size() < 3) {
                // Too few species with data
                res.params      = FitParams{};
                res.chi2        = 0.0;
                res.chi2_ndof   = 0.0;
                if (verbose) std::cout << "  skipping: too few species\n";
                results.push_back(res);
                continue;
            }

            // Build species list matching the filtered data
            std::vector<Species> species_used;
            for (const auto& f : filtered) {
                species_used.push_back(get_species(f.species));
            }

            auto fit = joint_fit(species_used, filtered, bw, bounds,
                                 stats, 3, verbose);

            res.params    = fit.params;
            res.chi2      = fit.chi2;
            res.chi2_ndof = fit.chi2_ndof;
            results.push_back(res);
        }
        return results;
    }

    void write_ptbin_scan(const std::vector<PtBinResult>& results,
                          const std::string& out_file)
    {
        std::ofstream out(out_file);
        if (!out) {
            std::cerr << "[ptbin_scan] cannot write " << out_file << "\n";
            return;
        }
        out << std::fixed << std::setprecision(6);
        out << "# win  pT_min  pT_max  T[GeV]  beta_s  mu_B[GeV]  mu_Q[GeV]  mu_S[GeV]  chi2_ndof  n_points\n";
        for (const auto& r : results) {
            out << r.window_index << "  "
                << r.pT_min       << "  "
                << r.pT_max       << "  "
                << r.params.T     << "  "
                << r.params.beta_s<< "  "
                << r.params.mu_B  << "  "
                << r.params.mu_Q  << "  "
                << r.params.mu_S  << "  "
                << r.chi2_ndof    << "  "
                << r.n_points     << "\n";
        }
        std::cout << "[ptbin_scan] wrote " << out_file
                  << "  (" << results.size() << " windows)\n";
    }

} // namespace bwmu
