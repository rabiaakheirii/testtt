#pragma once
#include "fit/chi2.hpp"
#include "fit/joint.hpp"
#include "blastwave/params.hpp"
#include "io/datapoint.hpp"
#include "particles/species.hpp"
#include <string>
#include <vector>

namespace bwmu {

    struct McmcChain {
        std::vector<std::vector<double>> samples;   // [n_keep][5]
        std::vector<double>              log_post;
        int    n_total         = 0;
        int    n_accepted      = 0;
        double acceptance_rate = 0.0;
        std::vector<std::string> param_names;
    };

    McmcChain run_mcmc(
        const std::vector<Species>& species,
        const std::vector<Spectrum>& data,
        const BWParams& bw,
        const FitBounds& bounds,
        const JointFitResult& start,
        Statistics stats,
        int n_steps  = 20000,
        int burn_in  = 4000,
        int thinning = 5,
        double step_T   = 0.003,
        double step_bs  = 0.010,
        double step_muB = 0.002,
        double step_muQ = 0.001,
        double step_muS = 0.001,
        unsigned seed   = 12345);

    void write_mcmc_chain(const McmcChain& chain, const std::string& out_file);
    void summarize_mcmc  (const McmcChain& chain, const std::string& out_file);

} // namespace bwmu
