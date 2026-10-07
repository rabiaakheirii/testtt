#include "io/phenix_reader.hpp"
#include "particles/table.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace bwmu {

    namespace {

        const std::vector<std::string> BLOCK_SPECIES = {
            "pi+", "pi-", "K+", "K-", "p", "pbar"
        };

        struct RawBlock {
            int n_header = 0;
            std::vector<double> pT;
            std::vector<std::vector<std::pair<double,double>>> cols;
        };

        // Return true if `s` (trimmed) is a single positive integer
        bool is_int_line(const std::string& s, int& value) {
            std::string t;
            for (char c : s) if (!std::isspace((unsigned char)c)) t += c;
            if (t.empty()) return false;
            for (char c : t) if (!std::isdigit((unsigned char)c)) return false;
            value = std::atoi(t.c_str());
            return true;
        }

        std::vector<RawBlock> parse_blocks(const std::string& path) {
            std::ifstream in(path);
            if (!in) {
                throw std::runtime_error(
                    "read_phenix: cannot open '" + path + "'");
            }

            std::vector<RawBlock> blocks;
            RawBlock cur;
            bool in_block = false;

            std::string line;
            while (std::getline(in, line)) {
                if (!line.empty() && line.back() == '\r') line.pop_back();
                if (line.empty()) continue;

                // Skip comment lines starting with '#'
                if (line[0] == '#') continue;

                int n_int = 0;
                if (is_int_line(line, n_int)) {
                    // New block: flush the previous one
                    if (in_block) {
                        blocks.push_back(cur);
                        cur = RawBlock{};
                    }
                    cur.n_header = n_int;
                    in_block = true;
                    continue;
                }

                if (!in_block) continue;

                std::istringstream row(line);
                double pT;
                if (!(row >> pT)) continue;

                std::vector<std::pair<double,double>> pairs;
                double v, e;
                while (row >> v >> e) {
                    pairs.emplace_back(v, e);
                }
                if (pairs.empty()) continue;

                cur.pT.push_back(pT);
                cur.cols.push_back(std::move(pairs));
            }
            if (in_block) blocks.push_back(cur);

            return blocks;
        }

    } // anonymous namespace

    std::vector<Spectrum> read_phenix(const std::string& path,
                                      int centrality,
                                      int n_centrality)
    {
        (void)n_centrality;   // unused but kept for API compatibility

        auto blocks = parse_blocks(path);

        std::cout << "[phenix_reader] parsed " << blocks.size() << " blocks\n";
        for (size_t i = 0; i < blocks.size(); ++i) {
            const auto& b = blocks[i];
            const int nc = b.cols.empty() ? 0
                        : static_cast<int>(b.cols[0].size());
            std::cout << "   block " << i << ": header=" << b.n_header
                      << "  rows=" << b.pT.size()
                      << "  columns=" << nc << "\n";
        }

        std::vector<Spectrum> out;
        for (size_t b = 0; b < blocks.size() &&
                           b < BLOCK_SPECIES.size(); ++b)
        {
            const std::string& sp_name = BLOCK_SPECIES[b];
            const Species& sp = get_species(sp_name);

            Spectrum s;
            s.centrality = centrality;
            s.species    = sp_name;

            for (size_t i = 0; i < blocks[b].pT.size(); ++i) {
                if (centrality < 0 ||
                    centrality >= static_cast<int>(blocks[b].cols[i].size()))
                    continue;

                const auto& p = blocks[b].cols[i][centrality];
                const double value = p.first;
                const double err   = std::abs(p.second);
                if (value <= 0.0 || err <= 0.0) continue;

                const double pT = blocks[b].pT[i];
                const double mT = std::sqrt(pT*pT + sp.mass*sp.mass);

                DataPoint d;
                d.mT    = mT - sp.mass;
                d.value = value;
                d.err   = err;
                s.points.push_back(d);
            }

            out.push_back(std::move(s));
        }
        return out;
    }

    Spectrum read_phenix_species(const std::string& path,
                                 int centrality,
                                 const std::string& species,
                                 int n_centrality)
    {
        auto all = read_phenix(path, centrality, n_centrality);
        for (auto& s : all) {
            if (s.species == species) return s;
        }
        throw std::runtime_error(
            "read_phenix_species: species '" + species + "' not found");
    }

} // namespace bwmu
