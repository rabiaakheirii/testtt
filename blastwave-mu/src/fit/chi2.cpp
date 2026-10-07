#include "fit/chi2.hpp"
#include "blastwave/spectrum.hpp"
#include "blastwave/feeddown.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace bwmu {

    double model_shape(double mT,
                       const Species& sp,
                       const BWParams& bw,
                       const FitParams& fp,
                       Statistics stats)
    {
        // T and beta_s come from fp, not bw (see Phase 6 notes).
        BWParams bw_local = bw;
        bw_local.T      = fp.T;
        bw_local.beta_s = fp.beta_s;

        double direct;
        if (stats == Statistics::Boltzmann) {
            direct = dN_dmT_boltzmann(mT, sp, bw_local,
                                      fp.mu_B, fp.mu_Q, fp.mu_S,
                                      1.0);
        } else {
            direct = dN_dmT_quantum(mT, sp, bw_local,
                                    fp.mu_B, fp.mu_Q, fp.mu_S,
                                    1.0, 4);
        }

        // If no resonances are loaded, return the direct contribution only.
        // This preserves the behaviour of Phases 5-7 tests.
        if (get_resonances().empty()) return direct;

        // Add feed-down for pions and kaons
        if (sp.name == "pi+" || sp.name == "pi-" ||
            sp.name == "K+"  || sp.name == "K-") {
            return direct + decay_contribution(mT, sp.name, bw_local, fp, stats);
        }
        return direct;
    }

    Chi2Result compute_chi2(
        const std::vector<Species>& species,
        const std::vector<Spectrum>& data,
        const BWParams& bw,
        const FitParams& fp,
        Statistics stats,
        const std::vector<int>& norm_group)
    {
        Chi2Result res;

        // Default: one group per species
        std::vector<int> groups = norm_group;
        if (groups.empty()) {
            groups.resize(species.size());
            for (size_t i = 0; i < species.size(); ++i)
                groups[i] = static_cast<int>(i);
        }

        const int n_groups = *std::max_element(groups.begin(), groups.end()) + 1;
        res.norms.assign(n_groups, 0.0);
        res.n_pts.assign(n_groups, 0);

        BWParams bw_local = bw;
        bw_local.T      = fp.T;
        bw_local.beta_s = fp.beta_s;

        // ------------------------------------------------------------
        // Step 1: accumulate linear-least-squares sums per group
        // ------------------------------------------------------------
        std::vector<double> num(n_groups, 0.0), den(n_groups, 0.0);

        for (size_t i = 0; i < species.size(); ++i) {
            const Spectrum* spec = nullptr;
            for (const auto& s : data) {
                if (s.species == species[i].name) { spec = &s; break; }
            }
            if (spec == nullptr) continue;

            const int g = groups[i];
            for (const auto& d : spec->points) {
                if (d.err <= 0.0) continue;
                const double f = model_shape(d.mT, species[i],
                                             bw_local, fp, stats);
                const double w = 1.0 / (d.err * d.err);
                num[g] += d.value * f * w;
                den[g] += f * f * w;
                res.n_pts[g] += 1;
            }
        }

        for (int g = 0; g < n_groups; ++g) {
            res.norms[g] = (den[g] > 0.0) ? num[g] / den[g] : 0.0;
        }

        // ------------------------------------------------------------
        // Step 2: accumulate chi2
        // ------------------------------------------------------------
        double total_chi2 = 0.0;
        int    total_pts  = 0;

        for (size_t i = 0; i < species.size(); ++i) {
            const Spectrum* spec = nullptr;
            for (const auto& s : data) {
                if (s.species == species[i].name) { spec = &s; break; }
            }
            if (spec == nullptr) continue;

            const int g = groups[i];
            const double A = res.norms[g];

            for (const auto& d : spec->points) {
                if (d.err <= 0.0) continue;
                const double f = model_shape(d.mT, species[i],
                                             bw_local, fp, stats);
                const double r = (d.value - A * f) / d.err;
                total_chi2 += r * r;
                ++total_pts;
            }
        }

        res.chi2 = total_chi2;
        res.ndof = total_pts - n_groups;
        res.chi2_ndof = (res.ndof > 0) ? total_chi2 / res.ndof : 0.0;

        return res;
    }

} // namespace bwmu
