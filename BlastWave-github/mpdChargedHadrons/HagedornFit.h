#include "def.h"
#include "WriteReadFiles.h"

double HagedornFunction(double *x, double *par) 
{
    // https://arxiv.org/pdf/2112.03187
    // x[0] - это pT
    double pT = x[0];
    double A = par[0]; // Нормировочный коэффициент
    double n = par[1]; // Параметр n
    double T = par[2]; // Параметр T
    double betaT = par[3];
    double m = par[4]; // масса
    double gammaT = 1 / sqrt(1 - betaT * betaT);

    double mT = sqrt(pT * pT + m * m);

    // Функция Хагедорна
    return A * pow(1.0 + gammaT * (mT - pT * betaT) / (n * T), -n);
    // return A * pow(1.0 + mT / (n * T), -n);
}       

void FitHagedorn( void )
{
    double outParams[N_PARTS][N_CENTR][5];
    double outParamsErr[N_PARTS][N_CENTR][5];

    double xmin[] = {0.4, 0.4, 0.3, 0.3, 0.4, 0.4};
    double xmax[] = {1.5, 1.2, 1.5, 1.5, 1.5, 1.5};

    TF1 *hagedorn[N_PARTS][N_CENTR];
    
    for (int part: PARTS)
    {
        for (int centr: CENTR)
        {
            hagedorn[part][centr] = new TF1("hagedorn", HagedornFunction, xmin[part], xmax[part], 5);
            hagedorn[part][centr]->SetParNames("A", "n", "T", "betaT", "m");
            hagedorn[part][centr]->FixParameter(4, masses[part]);

            hagedorn[part][centr]->SetParameters(1000, 10, 0.1, 0.5, masses[part]);
            hagedorn[part][centr]->SetParLimits(0, 0, 10000);
            hagedorn[part][centr]->SetParLimits(1, 0, 15);
            hagedorn[part][centr]->SetParLimits(2, 0.08, 0.12);
            hagedorn[part][centr]->SetParLimits(3, 0.3, 0.8);
 
            grSpectra[part][centr]->Fit(hagedorn[part][centr],  "QR+", "", xmin[part], xmax[part]);
            
            double *params = hagedorn[part][centr]->GetParameters();
            const double *paramsErr = hagedorn[part][centr]->GetParErrors();
            
            std::copy(params, params + 5, outParams[part][centr]);
            std::copy(paramsErr, paramsErr + 5, outParamsErr[part][centr]);

            cout << part << " " << centr << " " 
                 << params[0] << "  " << params[1] << "  " << params[2] << "  " << params[3] << endl;
        }
    }
    WriteParams(outParams, outParamsErr, "output/txtParams/HagedornParams.txt");
}