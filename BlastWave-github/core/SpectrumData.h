#ifndef BLASTWAVE_SPECTRUM_DATA_H
#define BLASTWAVE_SPECTRUM_DATA_H

#include "Constants.h"
#include "FitConfig.h"

#include "TGraphErrors.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace bw {

// Загружает спектры PHENIX из txt-файлов в TGraphErrors[частица][центральность].
//
// phenix_txt:       Spectra_particle_<i>_<i>.txt в cfg.spectra_path
// phenix_auau_txt:  один файл spectra.txt (все 6 частиц подряд)
class Dataset
{
public:
    explicit Dataset(const FitConfig& cfg) : cfg_(cfg)
    {
        if (cfg_.spectra_format == "phenix_auau_txt")
        {
            LoadAuAu();
        }
        else
        {
            for (int part = 0; part < kMaxParticles; ++part)
            {
                LoadParticle(part);
            }
        }
    }

    const FitConfig& Config() const { return cfg_; }

    TGraphErrors* Spectrum(int particle, int centrality) const
    {
        return spectra_[particle][centrality];
    }

    bool HasSpectrum(int particle, int centrality) const
    {
        return spectra_[particle][centrality] != nullptr;
    }

private:
    FitConfig cfg_;
    TGraphErrors* spectra_[kMaxParticles][kMaxCentralities]{};

    TGraphErrors* MakeGraph(int nPoints, const double* x, const double* y, const double* ey) const
    {
        std::vector<double> ex(nPoints, 0.05);
        return new TGraphErrors(nPoints, x, y, ex.data(),ey);
    }

    const double* XAxis(int particle, const std::vector<double>& pT,
                        const std::vector<double>& mT) const
    {
        return (cfg_.variable == "pt") ? pT.data() : mT.data();
    }

    // Формат файла: nPoints, затем для каждой центральности nPoints строк:
    // pT  yield  stat_err  syst_err
    void LoadParticle(int particle)
    {
        const std::string path = cfg_.spectra_path + "/Spectra_particle_"
                                 + std::to_string(particle) + "_"
                                 + std::to_string(particle) + ".txt";

        std::ifstream in(path);
        int nPoints = 0;
        in >> nPoints;

        std::vector<double> pT(nPoints), mT(nPoints), y(nPoints), ey(nPoints);

        for (int c = 0; c < cfg_.NumCentralities(); ++c)
        {
            for (int i = 0; i < nPoints; ++i)
            {
                double stat = 0.;
                double syst = 0.;
                in >> pT[i] >> y[i] >> stat >> syst;
                ey[i] = stat;
                mT[i] = MtMinusM(particle, pT[i]);
            }
            spectra_[particle][c] = 
                MakeGraph(nPoints, XAxis(particle, pT, mT), y.data(), ey.data());
        }

        std::cout << "Loaded " << path << "\n";
    }

    // Формат: для каждой частицы — nPoints, затем nPoints строк:
    // pT  y_c0  ey_c0  y_c1  ey_c1  ...
    void LoadAuAu()
    {
        std::ifstream in(cfg_.spectra_path);
        const int nCent = cfg_.NumCentralities();

        for (int particle = 0; particle < kMaxParticles; ++particle)
        {
            int nPoints = 0;
            in >> nPoints;

            std::vector<double> pT(nPoints), mT(nPoints);
            std::vector<std::vector<double>> y(nCent, std::vector<double>(nPoints));
            std::vector<std::vector<double>> ey(nCent, std::vector<double>(nPoints));

            for (int i = 0; i < nPoints; ++i)
            {
                in >> pT[i];
                mT[i] = MtMinusM(particle, pT[i]);
                for (int c = 0; c < nCent; ++c)
                {
                    in >> y[c][i] >> ey[c][i];
                    if (particle == 2 || particle == 3)
                    {
                        y[c][i] *= 10.;  // K+/K- в AuAu-файле в других единицах
                    }
                }
            }

            for (int c = 0; c < nCent; ++c)
            {
                spectra_[particle][c] = 
                    MakeGraph(nPoints, XAxis(particle, pT, mT), y[c].data(), ey[c].data());
            }
        }

        std::cout << "Loaded AuAu spectra: " << cfg_.spectra_path << "\n";
    }
};

}  // namespace bw

#endif
