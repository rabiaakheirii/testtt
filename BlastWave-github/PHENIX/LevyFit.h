#include "def.h"
#include "WriteReadFiles.h"

// levy080 = new TF1("levy",  "([0]*[3]*([1]-1)*([1]-2))*(pow((sqrt(x*x+0.13957*0.13957)+[2])/([2]+0.13957),-[1]))/(2*3.14159)/([2]+0.13957*([1]-1))/([2]+0.13957)",
//                               0.5, 2.5);

double LevyFunction(double *x, double *par) 
{
    // https://www.researchgate.net/publication/367678541_Nonextensive_Statistics_in_High_Energy_Collisions
    // x[0] - это pT
    double pT = x[0];
    double A = par[0]; // Нормировочный коэффициент
    double n = par[1]; // Параметр n
    double T = par[2]; // Параметр T
    double m = par[3]; // масса
    double mT = sqrt(pT * pT + m * m);
    double nT = n * T;
    
    // Функция Леви
    return A * (n - 1) * (n - 2) * pT / (nT * (nT + m * (n - 2))) \
             * pow(1 + (mT - m)/ nT, -n);
    
    // Другая формула для Леви
    // return A * pow(1.0 + (pT * pT) / (n * T), -n);
}                              

void FitLevy( void )
{
    double outParams[N_PARTS][N_CENTR][4];
    double outParamsErr[N_PARTS][N_CENTR][4];
    
    double xmin[] = {0.2, 0.2, 0.3, 0.3, 0.4, 0.4};
    double xmax[] = {1.2, 1.2, 1.5, 1.5, 1.5, 1.5};

    TF1 *levy[N_PARTS][N_CENTR];
    
    for (int part: PARTS)
    {
        for (int centr: CENTR)
        {
            levy[part][centr] = new TF1("levy", LevyFunction, xmin[part], xmax[part], 4);
            levy[part][centr]->SetParNames("A", "n", "T", "m");
            levy[part][centr]->FixParameter(3, masses[part]);

            if (part <= 1) // pi
            {
                levy[part][centr]->SetParameters(100, 4, 0.1, masses[part]);
                levy[part][centr]->SetParLimits(0, 0, 1500);
                levy[part][centr]->SetParLimits(1, 0, 15);
                levy[part][centr]->SetParLimits(2, 0.08, 0.15);
            }
            if (part == 2 || part == 3) //K
            {
                levy[part][centr]->SetParameters(10, 11, 0.1, masses[part]);
                levy[part][centr]->SetParLimits(0, 0, 100);
                levy[part][centr]->SetParLimits(1, 0, 30);
                levy[part][centr]->SetParLimits(2, 0.08, 0.15);
            }
            if (part >= 4) //p
            {
                levy[part][centr]->SetParameters(10, 15, 0.1, masses[part]);
                levy[part][centr]->SetParLimits(0, 0, 100);
                levy[part][centr]->SetParLimits(1, 0, 30);
                levy[part][centr]->SetParLimits(2, 0.08, 0.2);
            }
 
            grSpectra[part][centr]->Fit(levy[part][centr],  "QR+", "", xmin[part], xmax[part]);
            
            double *params = levy[part][centr]->GetParameters();
            const double *paramsErr = levy[part][centr]->GetParErrors();
        
            std::copy(params, params + 4, outParams[part][centr]);
            std::copy(paramsErr, paramsErr + 4, outParamsErr[part][centr]);

            cout << part << " " << centr << " " 
                 << params[0] << "  " << params[1] << "  " << params[2] << endl;
        }
    }

    WriteParams(outParams, outParamsErr, "output/txtParams/LevyParams.txt");
}
