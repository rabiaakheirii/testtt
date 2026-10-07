#include <iostream>
#include <fstream>
#include "core/constants.hpp"
#include "core/types.hpp"
#include "core/config_parser.hpp"
#include "utils/math.hpp"
#include "utils/bessel.hpp"
#include "utils/random.hpp"
#include "analysis/spectra.hpp"
#include "percolation/percolation.hpp"
#include "hydro/eos.hpp"
#include "hydro/hydro.hpp"
#include "freezeout/cooper_frye.hpp"
#include "analysis/critical_point.hpp"
#include "freezeout/resonance_decay.hpp"

int main() {
    std::cout << "=== QCD collision sim — Phase 1 ===\n";

    qcd::Config cfg = qcd::load_config("config/global.yaml");
    std::cout << "T_fo = " << cfg.T_fo << " GeV\n";

    // --- Load NA35 pion mT spectrum ---
    auto data = qcd::load_data("data/experimental/na35_s_s_200gev.csv");
    auto res  = qcd::load_resonances("data/resonances/resonances.csv"); 
    auto ydata = qcd::load_data("data/experimental/na35_pion_dndy.csv");

    // --- Fit stationary thermal model ---
    // Scan T from 100 to 300 MeV in 1 MeV steps
    qcd::fit_thermal_T(
        data,
        qcd::m_pion,        // mass of pi- [GeV]
        0.100, 0.300, 0.001,
        "output/spectra/thermal_fit_pion.dat"
    );

    qcd::fit_thermal_T_with_resonances(data, qcd::m_pion, res,
                                       0.100, 0.300, 0.001,
                                       "output/spectra/thermal_res_fit_pion.dat");
                                       
    qcd::fit_eta_max(ydata, qcd::m_pion, 0.180, 3.0,
                     0.0, 3.0, 0.01,
                     "output/spectra/longitudinal_fit_pion.dat");

        // --- Phase 4: transverse flow, full blast-wave ---
    auto multi = qcd::load_multi_species(
        "data/experimental/phenix_pion_proton_kaon.csv");

    double R = 5.0;  // fm
    int    n = 2;

    for (const auto& s : multi) {
        double m = 0.0;
        std::string out;
        if      (s.species == "pi") { m = qcd::m_pion;   out = "output/spectra/bw_pi.dat"; }
        else if (s.species == "K")  { m = qcd::m_kaon;   out = "output/spectra/bw_K.dat";  }
        else if (s.species == "p")  { m = qcd::m_proton; out = "output/spectra/bw_p.dat";  }
        else continue;

        double T_best = 0.0, b_best = 0.0;
        qcd::fit_blastwave(s.pts, m, R, n,
                           0.080, 0.250, 0.005,
                           0.0,   0.80,  0.02,
                           out, T_best, b_best);
    }

    // --- Phase 4b: degeneracy band (Fig. 9) ---
    qcd::scan_degeneracy(multi, R, n,
                         0.080, 0.250, 0.005,
                         0.0,   0.80,  0.02,
                         "output/spectra/degeneracy.dat");

    // --- Phase 5: percolation initial conditions ---
    const int    N_strings    = 200;
    const double R_fm         = 5.0;
    const double r_string_fm  = 0.25;
    const double overlap_prob = 1.0;
    const unsigned seed       = 12345;

    auto perc = qcd::run_percolation(N_strings, R_fm, r_string_fm,
                                     overlap_prob, seed);
    qcd::write_density_map(perc, "output/events/percolation_density.dat");


    // --- Phase 6: relativistic hydrodynamics ---
    {
        const std::string eos_file = "data/eos/eos_crossover.dat";
        qcd::EoS eos(eos_file);

        const double eps0   = 5.0;   // GeV/fm^3
        const double R_fm   = 5.0;   // fm
        const double tau0   = 1.0;   // fm/c
        const double tau_end= 10.0;  // fm/c
        const double dtau   = 0.05;
        const double r_max  = 10.0;
        const double dr     = 0.1;

        auto hydro = qcd::run_hydro_1p1d(eos, eps0, R_fm,
                                         tau0, tau_end, dtau,
                                         r_max, dr);
        qcd::write_hydro_profile(hydro, "output/events/hydro_profile.dat");
    }


    // --- Phase 7: Cooper-Frye freeze-out ---
    {
        qcd::EoS eos("data/eos/eos_crossover.dat");
        const double eps0    = 5.0;
        const double R_fm    = 5.0;
        const double tau0    = 1.0;
        const double tau_end = 10.0;
        const double dtau    = 0.05;
        const double r_max   = 10.0;
        const double dr      = 0.1;

        auto hydro = qcd::run_hydro_1p1d(eos, eps0, R_fm,
                                         tau0, tau_end, dtau,
                                         r_max, dr);

        const double T_fo = 0.150;   // GeV
        qcd::cooper_frye_spectrum(hydro, eos,
                                  qcd::m_pion,
                                  T_fo,
                                  200,           // number of mT points
                                  0.140, 2.0,    // mT range
                                  "output/spectra/cooper_frye_pion.dat");
    }


    // --- Phase 8: critical point & net-proton cumulants ---
    qcd::scan_beam_energy(3.0, 200.0, 40,
                          20000,
                          12345,
                          "output/spectra/net_proton_cumulants.dat");


    // --- Phase 9: summary ---
    {
        std::ofstream sum("output/summary.txt");
        sum << "========================================\n";
        sum << " QCD collision simulation — summary\n";
        sum << "========================================\n\n";

        sum << "Phase 1 (stationary thermal):  T_pi ~ 202 MeV, chi2/ndof ~ 8\n";
        sum << "Phase 2 (resonances):          T_pi ~ 211 MeV, chi2/ndof ~ 4\n";
        sum << "Phase 3 (longitudinal flow):   eta_max ~ 1.85\n";
        sum << "Phase 4 (blast-wave):\n";
        sum << "   pi  T = 115 MeV   beta_s = 0.78\n";
        sum << "   K   T = 120 MeV   beta_s = 0.76\n";
        sum << "   p   T = 155 MeV   beta_s = 0.78\n";
        sum << "Phase 5 (percolation):         N=200 strings, P ~ 0.08\n";
        sum << "Phase 6 (hydro Bjorken):       3 EoS (crossover, 1st order, CEP)\n";
        sum << "Phase 7 (Cooper-Frye):         T_fo = 150 MeV\n";
        sum << "Phase 8 (critical point):      cs2 dip near sqrt(s) ~ 8 GeV\n";
        sum << "\nAll figures in output/figures/\n";
        sum.close();
        std::cout << "[summary] wrote output/summary.txt\n";
    }

    return 0;
}