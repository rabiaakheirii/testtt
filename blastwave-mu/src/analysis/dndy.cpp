#include "analysis/dndy.hpp"
#include "particles/chem_potential.hpp"
#include "utils/bessel.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace bwmu {

    double dN_dy(double y,
                 const Species& sp,
                 const BWParams& bw,
                 const FitParams& fp,
                 double eta_max,
                 double norm)
    {
        // Integration grids
        const int    n_mT  = 60;
        const int    n_r   = 60;
        const int    n_eta = 60;

        const double mT_min = sp.mass + 1e-4;   // start slightly above rest mass
        const double mT_max = 4.0;              // GeV
        const double dmT    = (mT_max - mT_min) / n_mT;

        const double r_min = 0.0;
        const double r_max = bw.R_fm;
        const double dr    = (r_max - r_min) / n_r;

        const double eta_min = -eta_max;
        const double deta    = 2.0 * eta_max / n_eta;

        const double mu    = mu_of(sp, fp.mu_B, fp.mu_Q, fp.mu_S);
        const double chem  = std::exp(mu / fp.T);

        double total = 0.0;

        for (int iT = 0; iT <= n_mT; ++iT) {
            const double mT = mT_min + iT * dmT;
            const double pT2 = mT*mT - sp.mass*sp.mass;
            if (pT2 < 0.0) continue;
            const double pT = std::sqrt(pT2);
            const double wT = (iT == 0 || iT == n_mT) ? 0.5 : 1.0;

            for (int ir = 0; ir <= n_r; ++ir) {
                const double r = r_min + ir * dr;
                const double frac = (r_max > 0.0) ? r / r_max : 0.0;
                const double beta = bw.beta_s * std::pow(frac, bw.n_prof);
                if (beta >= 1.0) continue;
                const double rho = std::atanh(beta);
                const double wr = (ir == 0 || ir == n_r) ? 0.5 : 1.0;

                const double a   = mT * std::cosh(rho) / fp.T;
                const double b   = pT * std::sinh(rho) / fp.T;
                const double I0b = I0(b);

                // Integrate over eta
                double inner = 0.0;
                for (int ie = 0; ie <= n_eta; ++ie) {
                    const double eta = eta_min + ie * deta;
                    const double ch  = std::cosh(y - eta);
                    const double we  = (ie == 0 || ie == n_eta) ? 0.5 : 1.0;
                    inner += we * ch * std::exp(-a * ch);
                }
                inner *= deta * I0b;

                total += wT * mT * wr * r * inner;
            }
        }
        total *= dmT * dr;

        // Prefactors: g / (2 pi cosh y)
        // (tau_f is absorbed into `norm`)
        const double pre = sp.g / (2.0 * M_PI * std::cosh(y));

        return norm * pre * chem * total;
    }

    std::vector<std::pair<double, double>>
    compute_dndy_curve(const Species& sp,
                       const BWParams& bw,
                       const FitParams& fp,
                       double eta_max,
                       double y_min, double y_max, int n_y,
                       double norm)
    {
        std::vector<std::pair<double, double>> out;
        out.reserve(n_y);
        for (int i = 0; i < n_y; ++i) {
            const double y = (n_y > 1)
                ? y_min + (y_max - y_min) * i / (n_y - 1)
                : 0.5 * (y_min + y_max);
            out.emplace_back(y, dN_dy(y, sp, bw, fp, eta_max, norm));
        }
        return out;
    }

    void write_dndy(const std::vector<Species>& species,
                    const BWParams& bw,
                    const std::vector<FitParams>& params,
                    const std::vector<double>& norms,
                    double eta_max,
                    double y_min, double y_max, int n_y,
                    const std::string& out_file)
    {
        std::ofstream out(out_file);
        if (!out) {
            std::cerr << "[dndy] cannot write " << out_file << "\n";
            return;
        }
        out << std::fixed << std::setprecision(6);
        out << "# species  y  dN/dy\n";

        for (size_t i = 0; i < species.size(); ++i) {
            const auto& sp = species[i];
            const auto& fp = (i < params.size()) ? params[i] : FitParams{};
            const double n = (i < norms.size()) ? norms[i] : 1.0;

            std::cout << "[dndy] computing " << sp.name << " ...\n";
            auto curve = compute_dndy_curve(sp, bw, fp, eta_max,
                                            y_min, y_max, n_y, n);
            for (const auto& p : curve) {
                out << sp.name << "  " << p.first << "  " << p.second << "\n";
            }
            out << "\n";
        }
        std::cout << "[dndy] wrote " << out_file << "\n";
    }

} // namespace bwmu
