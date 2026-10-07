#include "def.h"
#include "WriteReadFiles.h"

using namespace std;

void GetContourPlots( int part, int centr )
{
    for (int s = 1; s < N_SIGMA; s++)
    {
        gMinuit->SetErrorDef(s * 8); //note 4 and not 2!
        contour[part][centr][s] = (TGraph*)gMinuit->Contour(40, 2, 1);
        if (contour[part][centr][s]) 
        {
            contour[part][centr][s]->SetLineColor(partColors[part] + centr);
            contour[part][centr][s]->SetFillStyle(0);
        }
        else 
            cout << "ERROR " << part << " " <<centr << " " << s << endl;
        
    }
}

class BlastWaveFit {
public:

    bool isContour = false;
    bool isDraw = true;
    double outParams[N_PARTS][N_CENTR][4];
    double outParamsErr[N_PARTS][N_CENTR][4];
    double paramsSystematics[N_PARTS][N_CENTR][4];
    double lLimitMult = 0.5, rLimitMult = 1.5; // for parLimits in case 4 (Systematic)
    double lLimitMultPi = 0.5, rLimitMultPi = 1.; // for parLimits in case 4 (Systematic Pi meson)

    void Fit( int systN, int initParamsType = 0 )
    {    
        cout << " ========== FIT(BWFit.h) ========== " << endl;

        TVirtualFitter::SetDefaultFitter("Minuit");  
        
        auto funcx = new TF1("funcx", bwfitfunc, 0.01, 10, 5);
        funcx->SetParameters(2,1);
        funcx->SetParNames("constant", "T", "beta", "mass", "pt");
        MyIntegFunc *integ = new MyIntegFunc(funcx);

        for (int part: PARTS)
        {  
            for (int centr: CENTR)
            {
                // cout << "PART: " << part << "   CENTR: " << centr << endl;
                string ifuncxName = "MyIntegFunc_" + to_string(part) + "_" + to_string(centr);
                ifuncx[part][centr] = new TF1("ifuncx", integ, xmin[part], xmax[part], 4, ifuncxName.c_str());

                switch(initParamsType)
                {
                    case 0: { /* DEFAULT */
                        // ================== version1 Params from Global fit ============================
                        double parResults[5];
                        TString systNamesT[5] = {"AuAu", "pAl", "HeAu", "CuAu", "UU"};
                        ReadGlobalParams(paramsGlobal, "output/txtParams/GlobalBWparams_" + systNamesT[systN] + ".txt");
                        getGlobalParams(part, centr, parResults);
                        if (parResults[0] == 0)
                            continue;
                            
                        parResults[2] = (parResults[2] > 0.6) ? 0.6 : parResults[2];
                        ifuncx[part][centr]->SetParameters(parResults);
                        for (int par = 0; par < 3; par++)
                        {
                            ifuncx[part][centr]->SetParLimits(par, parResults[par] * 0.5, parResults[par] * 1.5);
                            if (part <= 1 && par == 2) 
                                ifuncx[part][centr]->SetParLimits(par, parResults[par] * 0.5, parResults[par]);
                        }

                        ifuncx[part][centr]->FixParameter(3, masses[part]);
                        grSpectra[part][centr]->Fit(ifuncx[part][centr],  "QR+", "", xmin[part], xmax[part]);
                        
                        break;
                    }
                        // ====================================================================================

                    case 1: {
                        // ================== version2 Params from individual fit results =====================
                        double parResults[4];
                        ReadParams(part, centr, parResults);
                        if (parResults[0] == 0)
                            continue;
                            
                        ifuncx[part][centr]->SetParameters(parResults);
                        for (int par = 0; par < 3; par++)
                        {
                            ifuncx[part][centr]->SetParLimits(par, parResults[par] * 0.6, parResults[par] * 1.5);
                        }

                        ifuncx[part][centr]->FixParameter(3, masses[part]);
                        grSpectra[part][centr]->Fit(ifuncx[part][centr],  "QR+", "", xmin[part], xmax[part]);
                        break;
                    }
                        // ====================================================================================

                    case 2: {
                        // ================= version 2 Params with limits =================================
                        double customParams[4] = {con[part], 0.09, 0.75, masses[part]};
                        ifuncx[part][centr]->SetParameters(customParams);
                        ifuncx[part][centr]->SetParLimits(0, conmin[part], conmax[part]);
                        ifuncx[part][centr]->SetParLimits(1, 0.08, 0.14);	
                        ifuncx[part][centr]->SetParLimits(2, 0.4, 0.8);	
                        ifuncx[part][centr]->FixParameter(3, masses[part]);	//	mass
                        grSpectra[part][centr]->Fit(ifuncx[part][centr],  "QR+", "", xmin[part], xmax[part]);
                        break;
                    }   
                        // ====================================================================================

                    case 3: {
                        // ================= version 3 hand Params without Fit =============================+==
                        double handParams[4] = {handConst[part][centr], handT[centr], handBeta[centr], masses[part]};
                        ifuncx[part][centr]->SetParameters(handParams);
                        break;
                    }
                        // ====================================================================================
                    
                    case 4: {
                        // ================= version 4 Params For Systematics =================================    
                        paramsSystematics[part][centr][2] = (paramsSystematics[part][centr][2] > 0.95) ? 0.95 : paramsSystematics[part][centr][2];
                        ifuncx[part][centr]->SetParameters(paramsSystematics[part][centr]);
                        for (int par = 0; par < 3; par++)
                        {
                            if (paramsSystematics[part][centr][par] * rLimitMult > 0.95) rLimitMult = 0.95 / paramsSystematics[part][centr][par];
                            ifuncx[part][centr]->SetParLimits(par, paramsSystematics[part][centr][par] * lLimitMult, paramsSystematics[part][centr][par] * rLimitMult);
                            cout << paramsSystematics[part][centr][par] * lLimitMult << "   " << paramsSystematics[part][centr][par] * rLimitMult << endl;
                            // if (part <= 1 && par == 2) {
                            //     double rl = paramsSystematics[part][centr][par] * rLimitMultPi;

                            //     if (rl > 0.95) rl = 0.95;
                            //     ifuncx[part][centr]->SetParLimits(par, paramsSystematics[part][centr][par] * lLimitMultPi, rl);
                            // }
                                
                        }

                        ifuncx[part][centr]->FixParameter(3, masses[part]);
                        grSpectra[part][centr]->Fit(ifuncx[part][centr],  "QR+", "", xmin[part], xmax[part]);
                        break;
                    }
                        // ====================================================================================
                }
                
                ifuncx[part][centr]->SetLineColor(centrColors[centr]);
                
                double *params = ifuncx[part][centr]->GetParameters();
                const double *paramsErr = ifuncx[part][centr]->GetParErrors();
                std::copy(params, params + 4, outParams[part][centr]);
                std::copy(paramsErr, paramsErr + 4, outParamsErr[part][centr]);

            // +++++++++ Metrics ++++++++++++++++++++++++++++++++++++
                #if 0
                // NormalizeErrors on chi2/NDF
                double chi2 = ifuncx[part][centr]->GetChisquare();
                double ndf = ifuncx[part][centr]->GetNDF();
                double chi2Ndf = chi2 / ndf;
                for (int i = 0; i < 3; i++ )
                {
                    outParamsErr[part][centr][i] *= sqrt(chi2Ndf);
                }

                int N = 0, fitN = 0;
                double *x, *y, d = 0;
                x = grSpectra[part][centr]->GetX();
                y = grSpectra[part][centr]->GetY();
                N = grSpectra[part][centr]->GetN();

                for (int i = 0; i < N; i++)
                {
                    if (x[i] >= xmin[part] && x[i] <= xmax[part])
                    {
                        d += pow((y[i] - ifuncx[part][centr]->Eval(x[i])) / y[i], 2);
                        fitN++;
                    }
                }
                d = sqrt(d) / fitN;

                cout << part << " " << centr << "  " << d << " " << chi2Ndf << endl;
            #endif
            }
        }
    }
};



