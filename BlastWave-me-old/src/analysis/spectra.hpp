#pragma once
#include <string>
#include <vector>
#include "freezeout/resonance_decay.hpp"   // for Resonance

namespace qcd {

    struct DataPoint {
        double x;
        double value;
        double err;
    };

    std::vector<DataPoint> load_data(const std::string& filename);

    void fit_thermal_T(const std::vector<DataPoint>& data,
                       double m,
                       double T_min, double T_max, double T_step,
                       const std::string& output_file);

    void fit_thermal_T_with_resonances(const std::vector<DataPoint>& data,
                                       double m,
                                       const std::vector<Resonance>& resonances,
                                       double T_min, double T_max, double T_step,
                                       const std::string& output_file);
    
    void fit_eta_max(const std::vector<DataPoint>& data,
                     double m, double T, double y0,
                     double eta_min, double eta_max_max, double eta_step,
                     const std::string& output_file);
        struct SpeciesData {
        std::string species;          // "pi", "K", "p"
        std::vector<DataPoint> pts;
    };

    // Load multi-species CSV
    std::vector<SpeciesData> load_multi_species(const std::string& filename);

    // Phase 4: fit (T, beta_s) on one species with full blast-wave
    void fit_blastwave(const std::vector<DataPoint>& data,
                       double m, double R, int n,
                       double T_min, double T_max, double T_step,
                       double b_min, double b_max, double b_step,
                       const std::string& output_file,
                       double& T_best, double& b_best);

    // Phase 4b: degeneracy scan — for each T, find best beta_s
    void scan_degeneracy(const std::vector<SpeciesData>& all,
                         double R, int n,
                         double T_min, double T_max, double T_step,
                         double b_min, double b_max, double b_step,
                         const std::string& output_file);
} // namespace qcd