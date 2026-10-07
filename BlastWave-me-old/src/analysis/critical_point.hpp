#pragma once
#include "analysis/cumulants.hpp"
#include <string>
#include <vector>

namespace qcd {

    struct ScanPoint {
        double sqrt_s;     // collision energy [GeV]
        double mu_B;       // baryon chemical potential [GeV]
        double cs2;        // effective sound speed^2
        Cumulants cum;
    };

    // Generate toy events for given sqrt_s, compute cumulants
    ScanPoint simulate_event_sample(double sqrt_s,
                                    int    n_events,
                                    unsigned seed);

    // Scan over sqrt_s and write to file
    void scan_beam_energy(double sqrt_s_min, double sqrt_s_max, int n_steps,
                          int n_events_per_point,
                          unsigned seed,
                          const std::string& output_file);
}
