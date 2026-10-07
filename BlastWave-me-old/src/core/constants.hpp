#pragma once

namespace qcd {
    // Physical constants (natural units: hbar = c = kB = 1)
    constexpr double PI = 3.14159265358979323846;
    constexpr double GeV = 1.0;          // energy unit
    constexpr double MeV = 1e-3 * GeV;

    // Particle masses (in GeV)
    constexpr double m_pion   = 0.13957;
    constexpr double m_kaon   = 0.49367;
    constexpr double m_proton = 0.93827;
    constexpr double m_lambda = 1.11568;
    constexpr double m_rho    = 0.77526;
    constexpr double m_Delta  = 1.232;
}