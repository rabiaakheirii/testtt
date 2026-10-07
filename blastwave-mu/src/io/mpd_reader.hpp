#pragma once
#include "io/datapoint.hpp"
#include "particles/species.hpp"
#include <string>
#include <vector>

namespace bwmu {

    // ------------------------------------------------------------------
    // Reader for MPD/NICA simulation output.
    //
    // Expected format (text or ROOT-exported TXT):
    //   # centrality bin: 0-20%
    //   # particle: pi+
    //   pT      dN/(2pi pT dpT dy)    stat_err    sys_err
    //   0.15    1.23e+02              4.5e+00     1.2e+00
    //   0.25    8.90e+01              3.2e+00     8.9e-01
    //   ...
    //   # particle: pi-
    //   ...
    //
    // Three centrality bins: [0-20%], [20-40%], [40-80%].
    // ------------------------------------------------------------------
    std::vector<Spectrum> read_mpd(const std::string& path,
                                   int centrality_index,   // 0, 1, 2
                                   int n_centrality = 3);

    Spectrum read_mpd_species(const std::string& path,
                              int centrality_index,
                              const std::string& species,
                              int n_centrality = 3);

} // namespace bwmu
