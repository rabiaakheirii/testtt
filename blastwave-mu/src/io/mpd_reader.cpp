#include "io/mpd_reader.hpp"
#include "particles/table.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace bwmu {

    std::vector<Spectrum> read_mpd(const std::string& path,
                                   int centrality_index,
                                   int /*n_centrality*/)
    {
        std::ifstream in(path);
        if (!in) throw std::runtime_error("read_mpd: cannot open " + path);

        std::vector<Spectrum> out;
        Spectrum current;
        bool in_species = false;
        int  current_cent = -1;

        std::string line;
        while (std::getline(in, line)) {
            if (line.empty()) continue;

            // ---- centrality header ----
            if (line.find("# centrality") != std::string::npos) {
                // Parse the index from "centrality: N" if present
                // Otherwise increment
                ++current_cent;
                continue;
            }

            // ---- particle header ----
            if (line[0] == '#') {
                if (line.find("# particle:") != std::string::npos) {
                    if (in_species) out.push_back(current);
                    current = Spectrum{};
                    current.centrality = centrality_index;

                    // Extract species name
                    const auto pos = line.find(':');
                    std::string name = line.substr(pos + 1);
                    // Trim
                    while (!name.empty() && name.front() == ' ') name.erase(0, 1);
                    while (!name.empty() && name.back()  == ' ') name.pop_back();
                    current.species = name;
                    in_species = true;
                }
                continue;
            }

            if (!in_species) continue;

            // ---- data line ----
            std::istringstream ss(line);
            double pT, val, err, sys_err;
            if (!(ss >> pT >> val >> err)) continue;
            if (!(ss >> sys_err)) sys_err = 0.0;

            // Only keep the requested centrality
            if (current_cent != centrality_index) continue;

            // Convert pT -> mT - m0
            const Species& sp = get_species(current.species);
            const double mT = std::sqrt(pT*pT + sp.mass*sp.mass);

            DataPoint d;
            d.mT    = mT - sp.mass;
            d.value = val;
            // Combine statistical and systematic in quadrature
            d.err   = std::sqrt(err*err + sys_err*sys_err);
            if (d.err <= 0.0 || d.value <= 0.0) continue;

            current.points.push_back(d);
        }
        if (in_species) out.push_back(current);

        // Silent (caller will print summary if needed)
        return out;
    }

    Spectrum read_mpd_species(const std::string& path,
                              int centrality_index,
                              const std::string& species,
                              int n_centrality)
    {
        auto all = read_mpd(path, centrality_index, n_centrality);
        for (auto& s : all) if (s.species == species) return s;
        throw std::runtime_error("read_mpd_species: not found: " + species);
    }

} // namespace bwmu
