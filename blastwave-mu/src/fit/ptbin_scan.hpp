#pragma once
#include "blastwave/params.hpp"
#include "fit/chi2.hpp"
#include "fit/joint.hpp"
#include "io/datapoint.hpp"
#include "particles/species.hpp"
#include <string>
#include <utility>
#include <vector>

namespace bwmu {

    struct PtBinResult {
        int    window_index;
        double pT_min, pT_max;
        FitParams params;
        double chi2, chi2_ndof;
        int    n_points;    // total data points used in this window
    };

    // ------------------------------------------------------------------
    // Fit all species simultaneously in each pT window.
    //
    // For every window (pT_min, pT_max):
    //   - filter each species's points to that pT range
    //   - run joint_fit
    //   - record (T, beta_s, mu_B, mu_Q, mu_S, chi2/ndof)
    // ------------------------------------------------------------------
    std::vector<PtBinResult> scan_pT_windows(
        const std::vector<Species>& species,
        const std::vector<Spectrum>& data,
        const BWParams& bw,
        const FitBounds& bounds,
        const std::vector<std::pair<double,double>>& windows,
        Statistics stats,
        bool verbose = true);

    void write_ptbin_scan(
        const std::vector<PtBinResult>& results,
        const std::string& out_file);

} // namespace bwmu
