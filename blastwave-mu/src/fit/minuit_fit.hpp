#pragma once
#include "blastwave/params.hpp"
#include "fit/chi2.hpp"
#include "fit/joint.hpp"
#include "io/datapoint.hpp"
#include "particles/species.hpp"
#include <vector>

namespace bwmu {

    // Joint fit using ROOT's MINUIT2.
    // Bounds handled natively by MINUIT. T(beta_s) penalty added.
    JointFitResult joint_fit_minuit(
        const std::vector<Species>& species,
        const std::vector<Spectrum>& data,
        const BWParams& bw,
        const FitBounds& bounds,
        Statistics stats = Statistics::Quantum,
        double lambda_tbeta = 1.0e4,
        bool verbose = true);

} // namespace bwmu
