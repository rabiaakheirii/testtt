#pragma once
#include "blastwave/params.hpp"
#include "fit/chi2.hpp"
#include "io/datapoint.hpp"
#include "particles/species.hpp"
#include <string>
#include <vector>

namespace bwmu {

    struct Chi2SurfacePoint {
        double T;          // GeV
        double beta_s;
        double chi2;
        double chi2_ndof;
    };

    std::vector<Chi2SurfacePoint> compute_chi2_surface(
        const std::vector<Species>& species,
        const std::vector<Spectrum>& data,
        const BWParams& bw,
        const FitParams& best,
        double T_min, double T_max, int n_T,
        double b_min, double b_max, int n_b,
        Statistics stats = Statistics::Quantum);

    void scan_chi2_surface(
        const std::vector<Species>& species,
        const std::vector<Spectrum>& data,
        const BWParams& bw,
        const FitParams& best,
        double T_min, double T_max, int n_T,
        double b_min, double b_max, int n_b,
        Statistics stats,
        const std::string& out_file);

} // namespace bwmu
