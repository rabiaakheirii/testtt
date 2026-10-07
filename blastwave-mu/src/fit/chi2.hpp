#pragma once
#include "blastwave/params.hpp"
#include "io/datapoint.hpp"
#include "particles/species.hpp"

#include <string>
#include <vector>

namespace bwmu {

    struct FitParams {
        double T      = 0.100;
        double beta_s = 0.60;
        double mu_B   = 0.000;
        double mu_Q   = 0.000;
        double mu_S   = 0.000;
    };

    struct Chi2Result {
        double chi2      = 0.0;
        int    ndof      = 0;
        double chi2_ndof = 0.0;
        std::vector<double> norms;    // one entry per normalization group
        std::vector<int>    n_pts;    // one entry per normalization group
    };

    enum class Statistics {
        Boltzmann,
        Quantum
    };

    // ------------------------------------------------------------------
    // Compute chi2.
    //
    // norm_group[i] : index of the normalization group of species i.
    //   Species in the same group share a single normalization A.
    //   Example: {pi+, pi-, K+, K-, p, pbar} with
    //            norm_group = {0, 0, 1, 1, 2, 2}
    //   gives one A for the pion pair, one for the kaon pair, one for
    //   the proton pair. This is the physically correct treatment:
    //   volume, tau_f and degeneracy are identical for particle and
    //   antiparticle, so A must be shared.
    //
    // If `norm_group` is empty, each species is given its own group
    // (backward-compatible with Phase 5 behaviour).
    // ------------------------------------------------------------------
    Chi2Result compute_chi2(
        const std::vector<Species>& species,
        const std::vector<Spectrum>& data,
        const BWParams& bw,
        const FitParams& fp,
        Statistics stats = Statistics::Quantum,
        const std::vector<int>& norm_group = {});

    double model_shape(
        double mT,
        const Species& sp,
        const BWParams& bw,
        const FitParams& fp,
        Statistics stats);

} // namespace bwmu
