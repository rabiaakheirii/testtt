#include "freezeout/cooper_frye.hpp"
#include "utils/bessel.hpp"
#include <cmath>
#include <fstream>
#include <iostream>

namespace qcd {

    FreezeoutSurface find_freezeout_surface(const HydroProfile& H, double T_fo) {
        FreezeoutSurface S;
        S.tau_min = 1e9; S.tau_max = -1e9;
        S.r_min   = 1e9; S.r_max   = -1e9;

        for (int i = 0; i < H.n_tau; ++i) {
            for (int j = 0; j < H.n_r; ++j) {
                const int k = i * H.n_r + j;
                if (H.T[k] <= T_fo && H.T[k] > 0.0) {
                    const double tau = H.tau[i];
                    const double r   = H.r[j];
                    // check that the neighbour at smaller tau was above T_fo
                    // (i.e. we actually cross the isotherm)
                    if (i > 0) {
                        const int kprev = (i-1) * H.n_r + j;
                        if (H.T[kprev] < T_fo) continue;
                    }
                    S.tau.push_back(tau);
                    S.r.push_back(r);
                    S.T.push_back(H.T[k]);

                    S.tau_min = std::min(S.tau_min, tau);
                    S.tau_max = std::max(S.tau_max, tau);
                    S.r_min   = std::min(S.r_min,   r);
                    S.r_max   = std::max(S.r_max,   r);
                }
            }
        }
        std::cout << "[cooper-frye] freeze-out surface: "
                  << S.tau.size() << " points, "
                  << "tau in [" << S.tau_min << ", " << S.tau_max << "], "
                  << "r in ["   << S.r_min   << ", " << S.r_max   << "]\n";
        return S;
    }

    void cooper_frye_spectrum(const HydroProfile& H,
                              const EoS& /*eos*/,
                              double m,
                              double T_fo,
                              int    n_mT,
                              double mT_min, double mT_max,
                              const std::string& output_file)
    {
        auto S = find_freezeout_surface(H, T_fo);

        // Effective transverse area at freeze-out
        //   A_perp ~ pi * <r>^2
        double r_avg = 0.0;
        for (double r : S.r) r_avg += r;
        if (!S.r.empty()) r_avg /= S.r.size();
        const double A_perp = M_PI * r_avg * r_avg;

        // Effective tau at freeze-out
        double tau_avg = 0.0;
        for (double t : S.tau) tau_avg += t;
        if (!S.tau.empty()) tau_avg /= S.tau.size();

        std::cout << "[cooper-frye] <r> = " << r_avg
                  << " fm, A_perp = " << A_perp
                  << " fm^2, <tau> = " << tau_avg << " fm/c\n";

        // Write spectrum: dN/(mT dmT) = C * mT * K1(mT / T_fo)
        std::ofstream out(output_file);
        out << "# mT[GeV]   dN/(mT dmT)\n";

        const double dmT = (mT_max - mT_min) / (n_mT - 1);
        const double C = tau_avg * A_perp;   // normalization factor

        for (int i = 0; i < n_mT; ++i) {
            const double mT = mT_min + i * dmT;
            const double val = C * mT * K1(mT / T_fo);
            out << mT << "  " << val << "\n";
        }
        std::cout << "[cooper-frye] wrote " << output_file << "\n";
    }
}
