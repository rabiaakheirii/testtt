#pragma once
#include "blastwave/params.hpp"
#include "fit/chi2.hpp"
#include <string>
#include <vector>

namespace bwmu {

    struct Resonance {
        std::string name;
        double mass;      // GeV
        double g;         // spin degeneracy (2J+1)
        int    B, Q, S;
        double pstar;     // pion momentum in resonance rest frame [GeV]
        double BR;        // branching ratio into the target
        std::string target;  // "pi+", "pi-", "K+", "K-"
    };

    // Read the CSV, ignoring lines that start with '#'
    std::vector<Resonance> load_resonances(const std::string& path);

    // Global storage — set once at program start
    void set_resonances(const std::vector<Resonance>& r);
    const std::vector<Resonance>& get_resonances();

    // Sum of decays contributing to the given daughter species
    double decay_contribution(double mT,
                              const std::string& target,
                              const BWParams& bw,
                              const FitParams& fp,
                              Statistics stats);

} // namespace bwmu
