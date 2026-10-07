#pragma once
#include "blastwave/params.hpp"
#include "fit/chi2.hpp"
#include "io/datapoint.hpp"
#include "particles/species.hpp"

#include <string>
#include <vector>

namespace bwmu {

    // Full result of a joint fit
    struct JointFitResult {
        FitParams params;               // T, beta_s, mu_B, mu_Q, mu_S
        std::vector<double> norms;      // per-species normalisations
        double chi2;
        int    ndof;
        double chi2_ndof;
        int    evaluations;
        bool   converged;
    };

    // Parameter bounds for the grid search and Nelder-Mead
    struct FitBounds {
        double T_min,  T_max;
        double bs_min, bs_max;
        double muB_min, muB_max;
        double muQ_min, muQ_max;
        double muS_min, muS_max;
    };

    // Default bounds for PHENIX/ALICE physics
    FitBounds default_bounds();

    // ------------------------------------------------------------------
    // Joint fit over (T, beta_s, mu_B, mu_Q, mu_S).
    //
    // Algorithm:
    //   1. Coarse 5D grid search to bracket the minimum
    //   2. Nelder-Mead refinement from the best grid point
    //
    // The per-species normalisations are solved analytically inside
    // compute_chi2 (see Phase 5).
    // ------------------------------------------------------------------
    JointFitResult joint_fit(
        const std::vector<Species>& species,
        const std::vector<Spectrum>& data,
        const BWParams& bw,
        const FitBounds& bounds = default_bounds(),
        Statistics stats = Statistics::Quantum,
        int  grid_points_per_dim = 5,
        bool verbose = true);

    // Write the result to a text file
    void write_joint_fit_result(const JointFitResult& r,
                                const std::vector<Species>& species,
                                const std::string& path);

} // namespace bwmu
