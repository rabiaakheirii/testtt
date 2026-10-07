#include "hydro/hydro.hpp"
#include <cmath>
#include <fstream>
#include <iostream>

namespace qcd {

    HydroProfile run_hydro_1p1d(const EoS& eos,
                                double eps0,
                                double R_fm,
                                double tau0, double tau_end, double dtau,
                                double r_max, double dr)
    {
        HydroProfile H;
        H.n_tau = static_cast<int>((tau_end - tau0) / dtau) + 1;
        H.n_r   = static_cast<int>(r_max / dr) + 1;
        H.tau_min = tau0; H.tau_max = tau_end;
        H.r_min   = 0.0;  H.r_max = r_max;

        H.tau.resize(H.n_tau);
        H.r.resize(H.n_r);
        for (int i = 0; i < H.n_tau; ++i) H.tau[i] = tau0 + i * dtau;
        for (int j = 0; j < H.n_r;   ++j) H.r[j]   = j * dr;

        H.eps.assign(H.n_tau * H.n_r, 0.0);
        H.T.assign  (H.n_tau * H.n_r, 0.0);
        H.ur.assign (H.n_tau * H.n_r, 0.0);

        // initial energy density: Gaussian in r
        for (int j = 0; j < H.n_r; ++j) {
            const double r = H.r[j];
            const double e = eps0 * std::exp(-r*r / (2.0 * R_fm * R_fm));
            H.eps[0 * H.n_r + j] = e;
            H.T  [0 * H.n_r + j] = eos.temperature(e);
        }

        // 1+1D Bjorken: pure longitudinal expansion, no transverse flow
        //   d(eps)/d(tau) = -(eps + P) / tau
        for (int i = 0; i + 1 < H.n_tau; ++i) {
            for (int j = 0; j < H.n_r; ++j) {
                const double e = H.eps[i * H.n_r + j];
                const double P = eos.pressure(e);
                const double dedtau = -(e + P) / H.tau[i];
                const double e_next = e + dtau * dedtau;
                H.eps[(i+1) * H.n_r + j] = std::max(e_next, 0.0);
                H.T  [(i+1) * H.n_r + j] = eos.temperature(H.eps[(i+1) * H.n_r + j]);
                H.ur [(i+1) * H.n_r + j] = 0.0;   // no transverse flow in Bjorken
            }
        }

        std::cout << "[hydro] grid = " << H.n_tau << " x " << H.n_r
                  << "  tau = " << tau0 << " -> " << tau_end
                  << " fm/c   eps0 = " << eps0 << " GeV/fm^3\n";
        return H;
    }

    void write_hydro_profile(const HydroProfile& H,
                             const std::string& filename)
    {
        std::ofstream out(filename);
        out << "# tau[fm/c]  r[fm]  eps[GeV/fm^3]  T[GeV]  u^r\n";
        for (int i = 0; i < H.n_tau; ++i) {
            for (int j = 0; j < H.n_r; ++j) {
                const int k = i * H.n_r + j;
                out << H.tau[i] << "  " << H.r[j] << "  "
                    << H.eps[k] << "  " << H.T[k] << "  " << H.ur[k] << "\n";
            }
            out << "\n";
        }
        std::cout << "[hydro] wrote " << filename << "\n";
    }
}
