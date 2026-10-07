#include "io/alice_reader.hpp"
#include "particles/table.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace bwmu {

    Spectrum read_alice(const std::string& path,
                        int centrality,
                        const std::string& species)
    {
        const Species& sp = get_species(species);

        std::ifstream in(path);
        if (!in) {
            throw std::runtime_error(
                "read_alice: cannot open '" + path + "'");
        }

        Spectrum s;
        s.centrality = centrality;
        s.species    = species;

        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            double pT, value, err;
            if (!(ss >> pT >> value >> err)) continue;
            if (value <= 0.0 || err <= 0.0) continue;

            const double mT = std::sqrt(pT*pT + sp.mass*sp.mass);

            DataPoint d;
            d.mT    = mT - sp.mass;
            d.value = value;
            d.err   = err;
            s.points.push_back(d);
        }

        std::cout << "[alice_reader] " << path << "  species=" << species
                  << "  points=" << s.points.size() << "\n";
        return s;
    }

} // namespace bwmu
