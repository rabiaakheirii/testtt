#ifndef BLASTWAVE_PARAM_IO_H
#define BLASTWAVE_PARAM_IO_H

#include "Constants.h"
#include "FitConfig.h"

#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>

/*=========================================================
ParamIO.h - чтение и запись параметров фитов в txt файлы
=========================================================*/
namespace bw {

// Результат individual-фита одной частицы в одной центральности.
struct FitParams
{
    double C = 0.;
    double T = 0.;
    double beta = 0.;
    double mass = 0.;
    double C_err = 0.;
    double T_err = 0.;
    double beta_err = 0.;
};

// Результат global-фита: общие T и beta, свои C для pi / K / p.
struct GlobalParams
{
    double T = 0.;
    double beta = 0.;
    double C_pi = 0.;
    double C_K = 0.;
    double C_p = 0.;
};

// output/CuAu/txtParams/BWparams.txt
inline std::string IndividualParamsFile(const FitConfig& cfg)
{
    return cfg.output_dir + "/txtParams/BWparams.txt";
}

// output/CuAu/txtParams/GlobalBWparams.txt
inline std::string GlobalParamsFile(const FitConfig& cfg)
{
    return cfg.output_dir + "/txtParams/GlobalBWparams.txt";
}

// Записать individual-параметры после handFit / individualFromGlobal.
// Строка: particle  centrality  C  T  T_err  beta  beta_err
inline void WriteIndividualParams(const FitConfig& cfg,
                                  const FitParams results[kMaxParticles][kMaxCentralities])
{
    mkdir((cfg.output_dir + "/txtParams").c_str(), 0755);

    std::ofstream out(IndividualParamsFile(cfg));
    for (int part = 0; part < kMaxParticles; ++part)
    {
        for (int c = 0; c < cfg.NumCentralities(); ++c)
        {
            const FitParams& p = results[part][c];
            out << part << "  " << c << "  " << p.C << "  "
                << p.T << "   " << p.T_err << "   "
                << p.beta << "   " << p.beta_err << "\n";
        }
    }
    std::cout << "Wrote " << IndividualParamsFile(cfg) << "\n";
}

// Прочитать individual-параметры (для drawParams / draw).
inline void ReadIndividualParams(const FitConfig& cfg,
                                 FitParams results[kMaxParticles][kMaxCentralities])
{
    std::ifstream in(IndividualParamsFile(cfg));
    int part = 0;
    int c = 0;
    while (in >> part >> c >> results[part][c].C
               >> results[part][c].T >> results[part][c].T_err
               >> results[part][c].beta >> results[part][c].beta_err)
    {
        results[part][c].mass = kMasses[part];
    }
}

// Записать global-параметры для обоих зарядов (0 = pi+/K+/p, 1 = pi-/K-/pbar).
// Строка: charge  centrality  T  beta  C_pi  C_K  C_p
inline void WriteGlobalParams(const FitConfig& cfg,
                              const GlobalParams global[2][kMaxCentralities])
{
    mkdir((cfg.output_dir + "/txtParams").c_str(), 0755);

    std::ofstream out(GlobalParamsFile(cfg));
    for (int charge = 0; charge < 2; ++charge)
    {
        for (int c = 0; c < cfg.NumCentralities(); ++c)
        {
            const GlobalParams& gp = global[charge][c];
            out << charge << "  " << c << "  "
                << gp.T << "  " << gp.beta << "  "
                << gp.C_pi << "   " << gp.C_K << "   " << gp.C_p << "\n";
        }
    }
    std::cout << "Wrote " << GlobalParamsFile(cfg) << "\n";
}

// Прочитать global-параметры.
inline void ReadGlobalParams(const FitConfig& cfg,
                             GlobalParams global[2][kMaxCentralities])
{
    std::ifstream in(GlobalParamsFile(cfg));
    int charge = 0;
    int c = 0;
    while (in >> charge >> c >> global[charge][c].T >> global[charge][c].beta
               >> global[charge][c].C_pi >> global[charge][c].C_K >> global[charge][c].C_p)
    {
    }
}

// Собрать FitParams одной частицы из global-результата (для отрисовки кривых).
inline FitParams ParamsFromGlobal(int particle, int centrality,
                                  const GlobalParams global[2][kMaxCentralities])
{
    const int charge = particle % 2;
    const GlobalParams& gp = global[charge][centrality];

    FitParams p;
    p.T = gp.T;
    p.beta = gp.beta;
    p.mass = kMasses[particle];

    const int species = SpeciesGroup(particle);
    if (species == 0)
    {
        p.C = gp.C_pi;
    }
    else if (species == 1)
    {
        p.C = gp.C_K;
    }
    else
    {
        p.C = gp.C_p;
    }
    return p;
}

// Для DrawParams / DrawSpectra — создать папку output/<system>.
inline void EnsureDir(const std::string& path)
{
    mkdir(path.c_str(), 0755);
}

}  // namespace bw

#endif
