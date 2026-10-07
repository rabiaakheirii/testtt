#pragma once
#include "io/datapoint.hpp"
#include <string>
#include <vector>

namespace bwmu {

    // ------------------------------------------------------------------
    // Reader for PHENIX Au+Au files in the 6-block format.
    //
    // File structure:
    //   <n_points>            (integer on its own line)
    //   pT  val0  err0  val1  err1  ... val11 err11
    //   ...
    //   <blank line>
    //   <n_points>
    //   ...
    //
    // Block index -> species:
    //   0 -> pi+    1 -> pi-    2 -> K+    3 -> K-    4 -> p    5 -> pbar
    //
    // Each block has `n_centrality` columns (12 by default).
    // We keep only one centrality bin per call.
    //
    // Returns one Spectrum per species, for the requested centrality.
    // ------------------------------------------------------------------
    std::vector<Spectrum> read_phenix(const std::string& path,
                                      int centrality,
                                      int n_centrality = 12);

    // ------------------------------------------------------------------
    // Same, but returns only one species by name.
    // ------------------------------------------------------------------
    Spectrum read_phenix_species(const std::string& path,
                                 int centrality,
                                 const std::string& species,
                                 int n_centrality = 12);

} // namespace bwmu
