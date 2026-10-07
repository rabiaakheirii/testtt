#include "def.h"
#include "WriteReadFiles.h"

using namespace std;

void GetContourPlots( int part, int centr )
{
    int addColor[7] = {-9, -7, -4, 1, 2, 3, -1};
    Color_t partColorsContour[6] = {kRed, kRed, kBlue, kBlue, kGreen, kGreen};
    for (int s = 1; s < N_SIGMA; s++)
    {
        gMinuit->SetErrorDef(s * 256); //note 4 and not 2!
        contour[part][centr][s] = (TGraph*)gMinuit->Contour(16, 2, 1);
        if (contour[part][centr][s]) 
        {
            contour[part][centr][s]->SetLineColor(partColorsContour[part] + addColor[centr]);
            contour[part][centr][s]->SetLineWidth(3);
            contour[part][centr][s]->SetFillStyle(3001);
            contour[part][centr][s]->SetFillColorAlpha(partColorsContour[part]+ addColor[centr], 0.7);
        }
        else 
            cout << "ERROR " << part << " " <<centr << " " << s << endl;
        
    }
}

class BlastWaveFit {
public:

    bool isContour = true;
    bool isDraw = true;
    double outParams[N_PARTS][N_CENTR][4];
    double outParamsErr[N_PARTS][N_CENTR][4];
    double paramsSystematics[N_PARTS][N_CENTR][4];
    double lLimitMult = 0.5, rLimitMult = 1.5; // for parLimits in case 4 (Systematic)
    double lLimitMultPi = 0.5, rLimitMultPi = 1.; // for parLimits in case 4 (Systematic Pi meson)

    double Iexp[N_PARTS][N_CENTR];
    double Ibw[N_PARTS][N_CENTR];
    double IbwInt[N_PARTS][N_CENTR];
    double Idiff[N_PARTS][N_CENTR]; // Integrals and their difference

    void Fit( int initParamsType = 0 )
    {    
        // ++++++ Read data +++++++++++++++++++++++++++++++++++++

        string inputFileName = "postprocess_mpdpid10";
        // SetSpectra(inputFileName, "mt");
        SetSpectraMalaev("mt");

        // +++++++++ Fit +++++++++++++++++++++++++++++++++++++++
        TF1 * f = new TF1("f","[0] + x * [1] + x * x * [2]", 0.25, 0.9);
        f->SetParameters(0.129557, 0.232774, -0.339709);

        TVirtualFitter::SetDefaultFitter("Minuit");  
        // TMinuit* minuit = new TMinuit(5); 
        auto funcx = new TF1("funcx", bwfitfunc, 0.01, 10, 5);
        funcx->SetParameters(2,1);
        funcx->SetParNames("constant", "T", "beta", "mass", "pt");
        MyIntegFunc *integ = new MyIntegFunc(funcx);

        for (int part: PARTS)
        {  
            for (int centr: CENTR_MALAEV)
            {
                
                cout << "PART: " << part << "   CENTR: " << centr << endl;
                string ifuncxName = "MyIntegFunc_" + to_string(part) + "_" + to_string(centr);

                // cout << xmin[part] << "  " << xmax[part] << endl;
                ifuncx[part][centr] = new TF1("ifuncx", integ, xmin[part], xmax[part], 4, ifuncxName.c_str());
                ifuncx[part][centr]->SetLineColor(centrColors[centr]);  

                switch(initParamsType)
                {
                    case 0: { /* DEFAULT */
                        // ================== version1 Params from Global fit ============================
                        cout << "Case 0" << endl;
                        double parResults[5];
                        ReadGlobalParams(paramsGlobal);
                        getGlobalParams(part, centr, parResults);
                        if (parResults[0] == 0)
                            continue;
                            
                        parResults[2] = (parResults[2] > 0.6) ? 0.6 : parResults[2];
                        ifuncx[part][centr]->SetParameters(parResults);
                        
                        for (int par = 0; par < 3; par++)
                        {
                            double parMax = (par == 1) ?  parResults[par] :  parResults[par] * 1.5;
                            ifuncx[part][centr]->SetParLimits(par, parResults[par] * 0.5,parMax);
                            // if (part <= 1 && par == 2) 
                            //     ifuncx[part][centr]->SetParLimits(par, parResults[par] * 0.5, parResults[par]);

                            cout << par << "  " << parResults[par] * 0.5 << "  " << parResults[par] << endl;
                        }


                        ifuncx[part][centr]->FixParameter(3, masses[part]);

                        TFitResultPtr fitResult = grSpectra[part][centr]->Fit(ifuncx[part][centr], "QR+SEX0", "", xmin[part], xmax[part]);

                        // cout << "FIT RESULTS  " << ifuncx[part][centr]->GetParameter(1) << " " << ifuncx[part][centr]->GetParameter(2) << endl;
                        // if(fitResult->Status() != 0 || fitResult->Chi2() < 1e-6) {
                        //     cout << "BAD FIT: part " << part << " centr " << centr 
                        //         << " status: " << fitResult->Status() << endl;
                        //     continue;
                        // }

                        // TMatrixDSym corr = fitResult->GetCorrelationMatrix();
                        // if(corr.GetNrows() >= 3) { // Добавляем проверку размерности
                        //     double T_beta_corr = corr(1,2);
                        //     if(fabs(T_beta_corr) > 0.8) {
                        //         outParamsErr[part][centr][1] *= 1.2;
                        //         outParamsErr[part][centr][2] *= 1.2;
                        //     }
                        // }

                        // // Проверяем валидность результата
                        // if (fitResult->IsValid()) {

                        //     double chi2 = fitResult->Chi2();
                        //     int ndf = fitResult->Ndf();
                        //     double chi2_ndf = (ndf > 0) ? chi2 / ndf : -1;

                        //     std::cout << part << " "
                        //               << centr << " "
                        //               << "Chi2/NDF = " << chi2_ndf 
                        //               << " (Chi2 = " << chi2 
                        //               << ", NDF = " << ndf << ")\n" 
                        //               << std::endl;
                        // }
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
                        grSpectra[part][centr]->Fit(ifuncx[part][centr],  "QR+SEX0", "", xmin[part], xmax[part]);
                        break;
                    }
                        // ====================================================================================

                    case 2: {
                        // ================= version 2 Params with limits =================================
                        double customParams[4] = {handConst[part][centr], 0.1, 0.75, masses[part]};
                        ifuncx[part][centr]->SetParameters(customParams);
                        ifuncx[part][centr]->SetParLimits(0, 0, handConst[part][centr] * 10);
                        ifuncx[part][centr]->SetParLimits(1, 0.09, 0.15);	
                        ifuncx[part][centr]->SetParLimits(2, 0.4, 0.8);	
                        ifuncx[part][centr]->FixParameter(3, masses[part]);	//	mass

                        cout << "BW 2" << endl;
                        TFitResultPtr fitResult = 
                            grSpectra[part][centr]->Fit(ifuncx[part][centr], "QR+SEX0", "", 0.1, xmax[part]);

                        if(fitResult->Status() != 0 || fitResult->Chi2() < 1e-6) {
                            cout << "BAD FIT: part " << part << " centr " << centr 
                                << " status: " << fitResult->Status() << endl;
                            continue;
                        }
                        break;
                    }   
                        // ====================================================================================

                    case 3: {
                        // ================= version 3 hand Params without Fit =============================+==
                        isContour = false;
                        double handTeval = f->Eval(handBeta[centr]);
                        double handParams[4] = {handConst[part][centr], handTeval, handBeta[centr], masses[part]};
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
                
                if (isContour) GetContourPlots(part,  centr);    
                
        // +++++++++ Metrics ++++++++++++++++++++++++++++++++++++

                double *params = ifuncx[part][centr]->GetParameters();
                const double *paramsErr = ifuncx[part][centr]->GetParErrors();
                std::copy(params, params + 4, outParams[part][centr]);
                std::copy(paramsErr, paramsErr + 4, outParamsErr[part][centr]);

                cout << "out " << params[1] << "  " << outParams[part][centr][1]<<endl;
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
            }
        }
    }

    void GetIntegralDiff( void )
    {
        for (int part: PARTS)
        {
            for (int centr: CENTR_MALAEV)
            {
                Iexp[part][centr] = 0;
                for (int i = 0; i < grSpectra[part][centr]->GetN(); i++)
                {
                    double x = grSpectra[part][centr]->GetPointX(i);
                    if (x < xmax[part] && x > xmin[part])
                    {
                        Iexp[part][centr] += grSpectra[part][centr]->GetPointY(i);
                        Ibw[part][centr] += ifuncx[part][centr]->Eval(x);
                    }
                }
                
                IbwInt[part][centr] = ifuncx[part][centr]->Integral(xmin[part], xmax[part]);

                Idiff[part][centr] = abs(Iexp[part][centr] - Ibw[part][centr]);

                // cout << "Part: " << part << "  CENTR: " << centr 
                //      << "  Iexp: " << Iexp[part][centr] << "  Ibw: " << Ibw[part][centr] << " Idiff: " << Idiff[part][centr] << endl;
            }
        }

        // Write to file
        ofstream txtFile;
        txtFile.open("output/txtParams/DiffIntegrals.txt");

        for (int part: PARTS)
        {
            for (int centr: CENTR)
            {
                txtFile << part << "  " << centr << "  " 
                        << Iexp[part][centr] << "  " 
                        << Ibw[part][centr] << "   " 
                        << Idiff[part][centr] << endl;
            }
        }
        
        txtFile.close();
    }
};



