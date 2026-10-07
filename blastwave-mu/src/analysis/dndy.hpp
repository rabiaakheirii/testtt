#pragma once
#include "blastwave/params.hpp"
#include "fit/chi2.hpp"
#include "particles/species.hpp"
#include <string>
#include <utility>
#include <vector>

namespace bwmu {

    // ------------------------------------------------------------------
    // Rapidity distribution dN/dy from the 4D Cooper-Frye integral:
    //
    //   dN/dy = (g τ_f / (2π cosh y)) · e^{μ/T}
    //         × ∫ m_T dm_T ∫ r dr I₀(p_T sinh ρ / T)
    //             × ∫_{-η_max}^{+η_max} dη
    //                 cosh(y - η) · exp(-m_T cosh ρ cosh(y - η) / T)
    //
    // The double integral over (m_T, r) is done by trapezoidal rule.
    // The η integral is done by trapezoidal rule over [-η_max, +η_max].
    //
    // The overall normalization `norm` is a single scaling factor.
    // Returns dN/dy in arbitrary units (same as `norm`).
    // ------------------------------------------------------------------
    double dN_dy(double y,
                 const Species& sp,
                 const BWParams& bw,
                 const FitParams& fp,
                 double eta_max,
                 double norm = 1.0);

    // Compute dN/dy on a grid of y values.
    // Returns vector of (y, dN/dy) pairs.
    std::vector<std::pair<double, double>>
    compute_dndy_curve(const Species& sp,
                       const BWParams& bw,
                       const FitParams& fp,
                       double eta_max,
                       double y_min, double y_max, int n_y,
                       double norm = 1.0);

    // Write dN/dy curves for a list of species to a single file
    void write_dndy(const std::vector<Species>& species,
                    const BWParams& bw,
                    const std::vector<FitParams>& params,
                    const std::vector<double>& norms,
                    double eta_max,
                    double y_min, double y_max, int n_y,
                    const std::string& out_file);

} // namespace bwmu
