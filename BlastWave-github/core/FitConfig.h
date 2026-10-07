#ifndef BLASTWAVE_FIT_CONFIG_H
#define BLASTWAVE_FIT_CONFIG_H

#include "Constants.h"
#include "../external/nlohmann/json.hpp"

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace bw {

struct FitRange
{
    double xmin = 0.2;
    double xmax = 1.0;
};

struct TBetaPenalty
{
    bool enabled = false;
    double a0 = 0.129557;
    double a1 = 0.232774;
    double a2 = -0.339709;
    double sigma = 0.01;
};

struct HandParams
{
    std::vector<double> T;
    std::vector<double> beta;
    // C[particle][centrality]
    std::vector<std::vector<double>> C;
};

// PDF spectra panels: log Y limits per particle (pi+, pi-, K+, K-, p, pbar).
struct DrawSettings
{
    std::vector<double> y_min;
    std::vector<double> y_max;
};

// Окна Minuit для individual-фита (см. fit_limits в JSON).
struct FitLimits
{
    double T_min_scale = 0.5;
    double T_max_scale = 1.5;
    double beta_min_scale = 0.5;
    double beta_max_scale = 1.5;
    double C_max_scale = 100000.;
    double T_min = 0.08;
    double T_max = 0.22;
    double beta_min = 0.3;
    double beta_max = 0.9;
};

struct FitConfig
{
    std::string experiment;
    std::string system;
    double sqrt_s_NN = 200.0;
    std::string variable = "mt";          // "mt" or "pt"
    std::string spectra_format = "phenix_txt";
    std::string spectra_path;
    std::vector<std::string> centrality_labels;
    std::vector<double> npart;  // N_part per centrality bin (same order as centrality_labels)
    // Per species group: pi, K, p (applied to +/- pairs).
    FitRange fit_range_pi;
    FitRange fit_range_K;
    FitRange fit_range_p;
    HandParams hand;
    FitLimits fit_limits;
    TBetaPenalty t_beta_penalty;
    DrawSettings draw;
    std::string output_dir = "output";

    FitRange FitRangeForParticle(int particle) const
    {
        switch (SpeciesGroup(particle))
        {
            case 0: return fit_range_pi;
            case 1: return fit_range_K;
            default: return fit_range_p;
        }
    }

    int NumCentralities() const
    {
        return centrality_labels.size();
    }

    double Npart(int centrality) const
    {
        if (centrality >= 0 && centrality < npart.size())
        {
            return npart[centrality];
        }
        return 0.;
    }

    void PanelYLimits(int particle, double& ymin, double& ymax) const
    {
        static const double kDefaultMin[kMaxParticles] = {
            0.003, 0.003, 0.003, 0.003, 0.03, 0.003
        };
        static const double kDefaultMax[kMaxParticles] = {
            1e5, 1e5, 2e5, 2e5, 5e6, 5e5
        };

        if (draw.y_min.size() == kMaxParticles && draw.y_max.size() == kMaxParticles)
        {
            ymin = draw.y_min[particle];
            ymax = draw.y_max[particle];
            return;
        }

        ymin = kDefaultMin[particle];
        ymax = kDefaultMax[particle];
    }

    static FitConfig Load(const std::string& path);
};

inline FitRange FitRangeFromJson(const nlohmann::json& fitRange, const char* key)
{
    const auto& a = fitRange.at(key);
    return {a[0].get<double>(), a[1].get<double>()};
}

inline FitConfig FitConfig::Load(const std::string& path)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open config: " + path);

    nlohmann::json j;
    in >> j;

    FitConfig cfg;
    cfg.experiment = j.at("experiment").get<std::string>();
    cfg.system = j.at("system").get<std::string>();
    cfg.sqrt_s_NN = j.value("sqrt_s_NN", 200.0);
    cfg.variable = j.value("variable", "mt");
    cfg.output_dir = j.value("output_dir", "output");

    const auto& spectra = j.at("spectra");
    cfg.spectra_format = spectra.value("format", "phenix_txt");
    cfg.spectra_path = spectra.at("path").get<std::string>();

    cfg.centrality_labels = j.at("centrality").get<std::vector<std::string>>();
    cfg.npart = j.at("npart").get<std::vector<double>>();

    const auto& fitRange = j.at("fit_range");
    cfg.fit_range_pi = FitRangeFromJson(fitRange, "pi");
    cfg.fit_range_K = FitRangeFromJson(fitRange, "K");
    cfg.fit_range_p = FitRangeFromJson(fitRange, "p");

    const auto& hand = j.at("hand_params");
    cfg.hand.T = hand.at("T").get<std::vector<double>>();
    cfg.hand.beta = hand.at("beta").get<std::vector<double>>();
    cfg.hand.C = hand.at("C").get<std::vector<std::vector<double>>>();

    const auto& penalty = j.value("t_beta_penalty", nlohmann::json::object());
    cfg.t_beta_penalty.enabled = penalty.value("enabled", false);
    cfg.t_beta_penalty.a0 = penalty.value("a0", 0.129557);
    cfg.t_beta_penalty.a1 = penalty.value("a1", 0.232774);
    cfg.t_beta_penalty.a2 = penalty.value("a2", -0.339709);
    cfg.t_beta_penalty.sigma = penalty.value("sigma", 0.01);

    if (j.contains("draw"))
    {
        const auto& drawSec = j.at("draw");
        cfg.draw.y_min = drawSec.value("y_min", std::vector<double>{});
        cfg.draw.y_max = drawSec.value("y_max", std::vector<double>{});
    }

    if (j.contains("fit_limits"))
    {
        const auto& lim = j.at("fit_limits");
        cfg.fit_limits.T_min_scale = lim.value("T_min_scale", cfg.fit_limits.T_min_scale);
        cfg.fit_limits.T_max_scale = lim.value("T_max_scale", cfg.fit_limits.T_max_scale);
        cfg.fit_limits.beta_min_scale = lim.value("beta_min_scale", cfg.fit_limits.beta_min_scale);
        cfg.fit_limits.beta_max_scale = lim.value("beta_max_scale", cfg.fit_limits.beta_max_scale);
        cfg.fit_limits.C_max_scale = lim.value("C_max_scale", cfg.fit_limits.C_max_scale);
        cfg.fit_limits.T_min = lim.value("T_min", cfg.fit_limits.T_min);
        cfg.fit_limits.T_max = lim.value("T_max", cfg.fit_limits.T_max);
        cfg.fit_limits.beta_min = lim.value("beta_min", cfg.fit_limits.beta_min);
        cfg.fit_limits.beta_max = lim.value("beta_max", cfg.fit_limits.beta_max);
    }

    return cfg;
}

}  // namespace bw

#endif
