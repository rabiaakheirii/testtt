#pragma once
#include "particles/species.hpp"
#include <string>
#include <vector>

namespace bwmu {

    // Full PDG-based species table for the analysis.
    // Order: pi+, pi-, K+, K-, p, pbar, Lambda, Lambdabar,
    //        Xi-, Xibar+, Omega-, Omegabar+, phi, J/psi
    const std::vector<Species>& particle_table();

    // Lookup by name (throws std::runtime_error if not found)
    const Species& get_species(const std::string& name);

    // Print the table to stdout
    void print_particle_table();

} // namespace bwmu
