#pragma once
#include "io/datapoint.hpp"
#include <string>
#include <vector>

namespace bwmu {

    // ------------------------------------------------------------------
    // Reader for ALICE-style text files.
    //
    // Expected format (one row per data point):
    //     pT  value  err
    // Lines starting with '#' are comments. Blank lines are skipped.
    //
    // The caller supplies the species name and mass, so the reader
    // does not need to know which particle it is.
    // ------------------------------------------------------------------
    Spectrum read_alice(const std::string& path,
                        int centrality,
                        const std::string& species);

} // namespace bwmu
