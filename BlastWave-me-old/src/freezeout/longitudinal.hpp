#pragma once

namespace qcd {
    // thermal rapidity distribution centered at y0
    double thermal_dn_dy(double y, double m, double T, double norm, double y0);

    // convolution with longitudinal flow, centered at y0
    double flowing_dn_dy(double y, double m, double T, double norm,
                         double eta_max, double y0);
}
