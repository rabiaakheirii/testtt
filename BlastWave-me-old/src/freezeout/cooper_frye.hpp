#pragma once
#include "hydro/hydro.hpp"
#include "hydro/eos.hpp"
#include <string>
#include <vector>

namespace qcd {

    // Compute dN/(mT dmT) from a hydro profile
    //   Method: locate the freeze-out isotherm T = T_fo,
    //           integrate the Cooper-Frye flux over that surface.
    // For Bjorken + no transverse flow this reduces to:
    //   dN/(mT dmT) ~ mT * tau_fo * A_perp * K1(mT / T_fo)
    void cooper_frye_spectrum(const HydroProfile& hydro,
                              const EoS& eos,
                              double m_particle,
                              double T_fo,
                              int    n_mT,
                              double mT_min, double mT_max,
                              const std::string& output_file);

    // Find the freeze-out surface (points where T = T_fo)
    struct FreezeoutSurface {
        std::vector<double> tau;
        std::vector<double> r;
        std::vector<double> T;
        double tau_min, tau_max;
        double r_min,   r_max;
    };
    FreezeoutSurface find_freezeout_surface(const HydroProfile& hydro,
                                            double T_fo);
}
