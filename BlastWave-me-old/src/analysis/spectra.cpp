#include "analysis/spectra.hpp"
#include "freezeout/blastwave.hpp"
#include "freezeout/resonance_decay.hpp"
#include "utils/bessel.hpp"
#include "freezeout/longitudinal.hpp"
#include "freezeout/blastwave.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>
#include <algorithm>

namespace qcd {

    // -------- data loader --------
    std::vector<DataPoint> load_data(const std::string& filename) {
        std::vector<DataPoint> out;
        std::ifstream in(filename);
        if (!in) {
            std::cerr << "[spectra] cannot open " << filename << "\n";
            return out;
        }
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            DataPoint d;
            if (ss >> d.x >> d.value >> d.err) out.push_back(d);
        }
        std::cout << "[spectra] loaded " << out.size()
                  << " points from " << filename << "\n";
        return out;
    }

    // -------- Phase 1: pure thermal fit --------
    void fit_thermal_T(const std::vector<DataPoint>& data,
                       double m,
                       double T_min, double T_max, double T_step,
                       const std::string& output_file)
    {
        if (data.empty()) return;

        double best_T = 0.0, best_A = 0.0;
        double best_chi2 = std::numeric_limits<double>::max();

        for (double T = T_min; T <= T_max; T += T_step) {
            double num = 0.0, den = 0.0;
            for (const auto& d : data) {
                const double mT  = d.x + m;
                const double mdl = mT * K1(mT / T);
                const double w   = 1.0 / (d.err * d.err);
                num += d.value * mdl * w;
                den += mdl * mdl * w;
            }
            if (den <= 0.0) continue;
            const double A = num / den;

            double chi2 = 0.0;
            for (const auto& d : data) {
                const double mT  = d.x + m;
                const double mdl = A * mT * K1(mT / T);
                const double r   = (d.value - mdl) / d.err;
                chi2 += r * r;
            }
            if (chi2 < best_chi2) {
                best_chi2 = chi2;
                best_T    = T;
                best_A    = A;
            }
        }

        const int ndof = static_cast<int>(data.size()) - 1;

        std::cout << "----------------------------------------\n";
        std::cout << "Stationary thermal fit\n";
        std::cout << "  best T   = " << best_T * 1000.0 << " MeV\n";
        std::cout << "  norm A   = " << best_A << "\n";
        std::cout << "  chi2     = " << best_chi2 << "\n";
        std::cout << "  chi2/ndof= " << best_chi2 / ndof << "\n";
        std::cout << "----------------------------------------\n";

        std::ofstream out(output_file);
        out << "# mT-m0[GeV]  data  err  model(T=" << best_T << ")\n";
        for (const auto& d : data) {
            const double mT  = d.x + m;
            const double mdl = best_A * mT * K1(mT / best_T);
            out << d.x << "  " << d.value << "  " << d.err
                << "  " << mdl << "\n";
        }
        std::cout << "[spectra] wrote " << output_file << "\n";
    }

    // -------- Phase 2: thermal + resonance decays --------
    void fit_thermal_T_with_resonances(const std::vector<DataPoint>& data,
                                       double m,
                                       const std::vector<Resonance>& resonances,
                                       double T_min, double T_max, double T_step,
                                       const std::string& output_file)
    {
        if (data.empty()) return;

        double best_T = 0.0, best_A = 0.0;
        double best_chi2 = std::numeric_limits<double>::max();

        for (double T = T_min; T <= T_max; T += T_step) {
            double num = 0.0, den = 0.0;
            for (const auto& d : data) {
                const double mT  = d.x + m;
                const double mdl = total_pion_mT(mT, m, T, 1.0, resonances);
                const double w   = 1.0 / (d.err * d.err);
                num += d.value * mdl * w;
                den += mdl * mdl * w;
            }
            if (den <= 0.0) continue;
            const double A = num / den;

            double chi2 = 0.0;
            for (const auto& d : data) {
                const double mT  = d.x + m;
                const double mdl = A * total_pion_mT(mT, m, T, 1.0, resonances);
                const double r   = (d.value - mdl) / d.err;
                chi2 += r * r;
            }
            if (chi2 < best_chi2) {
                best_chi2 = chi2;
                best_T    = T;
                best_A    = A;
            }
        }

        const int ndof = static_cast<int>(data.size()) - 1;

        std::cout << "----------------------------------------\n";
        std::cout << "Thermal + resonance fit\n";
        std::cout << "  best T   = " << best_T * 1000.0 << " MeV\n";
        std::cout << "  norm A   = " << best_A << "\n";
        std::cout << "  chi2     = " << best_chi2 << "\n";
        std::cout << "  chi2/ndof= " << best_chi2 / ndof << "\n";
        std::cout << "----------------------------------------\n";

        std::ofstream out(output_file);
        out << "# mT-m0  data  err  model  direct  decays\n";
        for (const auto& d : data) {
            const double mT     = d.x + m;
            const double direct = best_A * thermal_mT_spectrum(mT, best_T, 1.0);
            const double total  = best_A * total_pion_mT(mT, m, best_T, 1.0, resonances);
            const double decays = total - direct;
            out << d.x << "  " << d.value << "  " << d.err
                << "  " << total << "  " << direct << "  " << decays << "\n";
        }
        std::cout << "[spectra] wrote " << output_file << "\n";
    }
    void fit_eta_max(const std::vector<DataPoint>& data,
                     double m, double T, double y0,
                     double eta_min, double eta_max_max, double eta_step,
                     const std::string& output_file)
    {
        if (data.empty()) return;

        double best_eta = 0.0, best_A = 0.0;
        double best_chi2 = std::numeric_limits<double>::max();

        for (double eta = eta_min; eta <= eta_max_max; eta += eta_step) {
            double num = 0.0, den = 0.0;
            for (const auto& d : data) {
                const double mdl = flowing_dn_dy(d.x, m, T, 1.0, eta, y0);
                const double w   = 1.0 / (d.err * d.err);
                num += d.value * mdl * w;
                den += mdl * mdl * w;
            }
            if (den <= 0.0) continue;
            const double A = num / den;

            double chi2 = 0.0;
            for (const auto& d : data) {
                const double mdl = A * flowing_dn_dy(d.x, m, T, 1.0, eta, y0);
                const double r   = (d.value - mdl) / d.err;
                chi2 += r * r;
            }
            if (chi2 < best_chi2) {
                best_chi2 = chi2;
                best_eta  = eta;
                best_A    = A;
            }
        }

        const int ndof = static_cast<int>(data.size()) - 1;

        std::cout << "----------------------------------------\n";
        std::cout << "Longitudinal flow fit\n";
        std::cout << "  T        = " << T * 1000.0 << " MeV (fixed)\n";
        std::cout << "  eta_max  = " << best_eta << "\n";
        std::cout << "  norm A   = " << best_A << "\n";
        std::cout << "  chi2/ndof= " << best_chi2 / ndof << "\n";
        std::cout << "----------------------------------------\n";

        std::ofstream out(output_file);
        out << "# y  data  err  flow(eta=" << best_eta
            << ")  thermal\n";
        for (const auto& d : data) {
            const double flow    = best_A * flowing_dn_dy(d.x, m, T, 1.0, best_eta, y0);
            const double thermal = best_A * thermal_dn_dy(d.x, m, T, 1.0, y0);
            out << d.x << "  " << d.value << "  " << d.err
                << "  " << flow << "  " << thermal << "\n";
        }
        std::cout << "[spectra] wrote " << output_file << "\n";
    }
    

    // ---------- multi-species loader ----------
    std::vector<SpeciesData> load_multi_species(const std::string& filename) {
        std::vector<SpeciesData> out;
        std::ifstream in(filename);
        if (!in) { std::cerr << "cannot open " << filename << "\n"; return out; }
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            std::string sp; DataPoint d;
            if (!(ss >> sp >> d.x >> d.value >> d.err)) continue;

            auto it = std::find_if(out.begin(), out.end(),
                [&](const SpeciesData& s){ return s.species == sp; });
            if (it == out.end()) {
                SpeciesData s; s.species = sp; s.pts.push_back(d);
                out.push_back(s);
            } else {
                it->pts.push_back(d);
            }
        }
        std::cout << "[spectra] multi-species: " << out.size() << " species\n";
        return out;
    }

    // ---------- Phase 4: blast-wave fit for one species ----------
    void fit_blastwave(const std::vector<DataPoint>& data,
                       double m, double R, int n,
                       double T_min, double T_max, double T_step,
                       double b_min, double b_max, double b_step,
                       const std::string& output_file,
                       double& T_best, double& b_best)
    {
        if (data.empty()) return;
        double best_chi2 = std::numeric_limits<double>::max();
        double best_A = 0.0;
        T_best = 0.0; b_best = 0.0;

        for (double T = T_min; T <= T_max; T += T_step) {
            for (double b = b_min; b <= b_max; b += b_step) {
                double num = 0.0, den = 0.0;
                for (const auto& d : data) {
                    const double mT  = d.x + m;
                    const double mdl = blastwave_mT_spectrum(mT, T, b, R, n, 1.0);
                    const double w   = 1.0 / (d.err * d.err);
                    num += d.value * mdl * w;
                    den += mdl * mdl * w;
                }
                if (den <= 0.0) continue;
                const double A = num / den;

                double chi2 = 0.0;
                for (const auto& d : data) {
                    const double mT  = d.x + m;
                    const double mdl = A * blastwave_mT_spectrum(mT, T, b, R, n, 1.0);
                    const double r   = (d.value - mdl) / d.err;
                    chi2 += r * r;
                }
                if (chi2 < best_chi2) {
                    best_chi2 = chi2;
                    T_best = T; b_best = b; best_A = A;
                }
            }
        }

        const int ndof = static_cast<int>(data.size()) - 2;

        std::cout << "----------------------------------------\n";
        std::cout << "Blast-wave fit (m=" << m << " GeV)\n";
        std::cout << "  T       = " << T_best * 1000.0 << " MeV\n";
        std::cout << "  beta_s  = " << b_best << "\n";
        std::cout << "  chi2/ndof = " << best_chi2 / ndof << "\n";
        std::cout << "----------------------------------------\n";

        std::ofstream out(output_file);
        out << "# mT-m0  data  err  model\n";
        for (const auto& d : data) {
            const double mT  = d.x + m;
            const double mdl = best_A * blastwave_mT_spectrum(mT, T_best, b_best, R, n, 1.0);
            out << d.x << "  " << d.value << "  " << d.err << "  " << mdl << "\n";
        }
        std::cout << "[spectra] wrote " << output_file << "\n";
    }

    // ---------- Phase 4b: degeneracy scan ----------
    void scan_degeneracy(const std::vector<SpeciesData>& all,
                         double R, int n,
                         double T_min, double T_max, double T_step,
                         double b_min, double b_max, double b_step,
                         const std::string& output_file)
    {
        std::ofstream out(output_file);
        out << "# species  T[GeV]  beta_s  chi2_min\n";

        for (const auto& s : all) {
            double m = 0.0;
            if (s.species == "pi") m = 0.13957;
            else if (s.species == "K")  m = 0.49367;
            else if (s.species == "p")  m = 0.93827;
            else continue;

            for (double T = T_min; T <= T_max; T += T_step) {
                double best_chi2 = std::numeric_limits<double>::max();
                double best_b = 0.0;
                for (double b = b_min; b <= b_max; b += b_step) {
                    double num = 0.0, den = 0.0;
                    for (const auto& d : s.pts) {
                        const double mT  = d.x + m;
                        const double mdl = blastwave_mT_spectrum(mT, T, b, R, n, 1.0);
                        const double w   = 1.0 / (d.err * d.err);
                        num += d.value * mdl * w;
                        den += mdl * mdl * w;
                    }
                    if (den <= 0.0) continue;
                    const double A = num / den;

                    double chi2 = 0.0;
                    for (const auto& d : s.pts) {
                        const double mT  = d.x + m;
                        const double mdl = A * blastwave_mT_spectrum(mT, T, b, R, n, 1.0);
                        const double r   = (d.value - mdl) / d.err;
                        chi2 += r * r;
                    }
                    if (chi2 < best_chi2) { best_chi2 = chi2; best_b = b; }
                }
                out << s.species << "  " << T << "  " << best_b
                    << "  " << best_chi2 << "\n";
            }
        }
        std::cout << "[spectra] wrote degeneracy scan " << output_file << "\n";
    }
} // namespace qcd