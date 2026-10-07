#pragma once
#include "io/datapoint.hpp"
#include <string>
#include <vector>

namespace bwmu {

    // Write one spectrum to a CSV file
    void write_spectrum_csv(const std::string& path,
                            const Spectrum& s);

    // Write multiple spectra to one CSV with a species column
    void write_spectra_csv(const std::string& path,
                           const std::vector<Spectrum>& spectra);

} // namespace bwmu
