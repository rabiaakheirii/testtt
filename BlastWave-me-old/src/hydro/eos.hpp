#pragma once
#include <string>
#include <vector>

namespace qcd {

    struct EosPoint {
        double eps;   // energy density [GeV/fm^3]
        double P;     // pressure      [GeV/fm^3]
        double T;     // temperature   [GeV]
        double cs2;   // sound speed^2
    };

    class EoS {
    public:
        explicit EoS(const std::string& filename);

        double pressure(double eps)  const;
        double temperature(double eps) const;
        double cs2(double eps)       const;
        double eps_from_T(double T)  const;

        const std::vector<EosPoint>& table() const { return table_; }

    private:
        std::vector<EosPoint> table_;
        double interpolate(double eps, int col) const;
    };
}
