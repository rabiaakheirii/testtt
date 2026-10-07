#include "io/csv_writer.hpp"
#include <fstream>
#include <iostream>

namespace bwmu {

    void write_spectrum_csv(const std::string& path,
                            const Spectrum& s)
    {
        std::ofstream out(path);
        if (!out) {
            std::cerr << "[csv_writer] cannot write " << path << "\n";
            return;
        }
        out << "# species=" << s.species
            << "  centrality=" << s.centrality << "\n";
        out << "# mT-m0[GeV]  dN/(mT dmT)[GeV^-2]  err\n";
        for (const auto& d : s.points) {
            out << d.mT << "  " << d.value << "  " << d.err << "\n";
        }
        std::cout << "[csv_writer] wrote " << path
                  << "  (" << s.points.size() << " points)\n";
    }

    void write_spectra_csv(const std::string& path,
                           const std::vector<Spectrum>& spectra)
    {
        std::ofstream out(path);
        if (!out) {
            std::cerr << "[csv_writer] cannot write " << path << "\n";
            return;
        }
        out << "# species  centrality  mT-m0[GeV]  dN/(mT dmT)[GeV^-2]  err\n";
        for (const auto& s : spectra) {
            for (const auto& d : s.points) {
                out << s.species << "  "
                    << s.centrality << "  "
                    << d.mT << "  " << d.value << "  " << d.err << "\n";
            }
        }
        std::cout << "[csv_writer] wrote " << path << "\n";
    }

} // namespace bwmu
