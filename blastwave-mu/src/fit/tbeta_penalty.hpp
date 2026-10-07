#pragma once

namespace bwmu {

    // Empirical T(beta_s) relation from published PHENIX
    // PRC 69, 034909 (2004). Linear interpolation between
    // published points.
    //
    // Returns T in GeV.
    double T_exp_of_beta(double beta_s);

    // Penalty = lambda * (T - T_exp(beta_s))^2
    // Adds to chi2 to keep the fit on the published T(beta) curve.
    double tbeta_penalty(double T, double beta_s, double lambda);

} // namespace bwmu
