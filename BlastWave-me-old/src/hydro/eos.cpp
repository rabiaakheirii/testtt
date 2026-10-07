#include "hydro/eos.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

namespace qcd {

    EoS::EoS(const std::string& filename) {
        std::ifstream in(filename);
        if (!in) { std::cerr << "[eos] cannot open " << filename << "\n"; return; }
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            EosPoint p;
            if (ss >> p.eps >> p.P >> p.T >> p.cs2) table_.push_back(p);
        }
        std::cout << "[eos] loaded " << table_.size()
                  << " points from " << filename << "\n";
    }

    double EoS::interpolate(double eps, int col) const {
        if (table_.empty()) return 0.0;
        if (eps <= table_.front().eps) {
            switch (col) {
                case 0: return table_.front().P;
                case 1: return table_.front().T;
                case 2: return table_.front().cs2;
            }
        }
        if (eps >= table_.back().eps) {
            switch (col) {
                case 0: return table_.back().P;
                case 1: return table_.back().T;
                case 2: return table_.back().cs2;
            }
        }
        for (size_t i = 0; i + 1 < table_.size(); ++i) {
            const auto& a = table_[i];
            const auto& b = table_[i + 1];
            if (eps >= a.eps && eps <= b.eps) {
                const double t = (eps - a.eps) / (b.eps - a.eps);
                switch (col) {
                    case 0: return a.P   + t * (b.P   - a.P);
                    case 1: return a.T   + t * (b.T   - a.T);
                    case 2: return a.cs2 + t * (b.cs2 - a.cs2);
                }
            }
        }
        return 0.0;
    }

    double EoS::pressure(double eps)    const { return interpolate(eps, 0); }
    double EoS::temperature(double eps) const { return interpolate(eps, 1); }
    double EoS::cs2(double eps)         const { return interpolate(eps, 2); }

    double EoS::eps_from_T(double T) const {
        for (size_t i = 0; i + 1 < table_.size(); ++i) {
            const auto& a = table_[i];
            const auto& b = table_[i + 1];
            if (T >= a.T && T <= b.T) {
                const double t = (T - a.T) / (b.T - a.T);
                return a.eps + t * (b.eps - a.eps);
            }
        }
        return table_.empty() ? 0.0 : table_.back().eps;
    }
}
