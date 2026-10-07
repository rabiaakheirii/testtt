#include "input/headers/def.h"
#include "input/headers/WriteReadFiles.h"
#include "input/headers/BlastWaveFit.h"

using namespace std;

void setParamsForSys ( int systematicType, double parResults[4], BlastWaveFit *bwFit )
{
    if ( systematicType == 0 ) return;

    switch(systematicType)
    {
        case 1: 
            parResults[1] *= 0.9;
            break; // Добавить
        case 2: 
            parResults[1] *= 1.1;
            break; // Добавить
        case 3:
            parResults[2] *= 0.9;
            break; // Добавить
        case 4:
            parResults[2] *= 1.1;            
            break; // Добавить
        case 5:
            parResults[0] *= 0.9;
            break; // Добавить
        case 6:
            parResults[0] *= 1.1;
            break; // Добавить
        case 7:
            bwFit->lLimitMult *= 0.9;
            bwFit->rLimitMult *= 0.9;
            break; // Добавить
        case 8:
            bwFit->lLimitMult *= 1.1;
            bwFit->rLimitMult *= 1.1;
            break; // Добавить
        case 9:
            bwFit->lLimitMultPi *= 1.1;        
            bwFit->rLimitMultPi *= 0.9;
            break; // Добавить
        case 10:
            bwFit->lLimitMultPi *= 0.9;
            bwFit->rLimitMultPi *= 1.1;
            break; // Добавить
    }
}


void BlastWaveSystematic( void )
{
    bool isDraw = true;

    double systErr[N_PARTS][N_CENTR][4];
    
    // Добавляем объявления systErr_low и systErr_high
    double systErr_low[N_PARTS][N_CENTR][4] = {0};
    double systErr_high[N_PARTS][N_CENTR][4] = {0};

    for (int part: PARTS) {
        for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
            int centr = CENTR_SYST[systN][j];
            
                for (int p = 0; p < 4; p++) {
                    systErr[part][centr][p] = 0;
                }
        }
    }

                
    // ++++++ Read data +++++++++++++++++++++++++++++++++++++

    if (systN == 0) 
        ReadFromFileAuAu(); // Для системы AuAu
    else                    // Для других систем
        for (int part: PARTS) ReadFromFile(part, systN);

    // +++++++++ Fit +++++++++++++++++++++++++++++++++++++++

    BlastWaveFit *bwFitRef = new BlastWaveFit();
    bwFitRef->Fit(0);
    WriteParams(systN, bwFitRef->outParams, bwFitRef->outParamsErr, true, "output/parameters/SYST_" + systNamesT[systN] + ".txt");

    
    for (int systematicType = 1; systematicType <= 10; systematicType++)
    {
        if(systematicType >=7 && systN != 0) continue;
        cout << "Starting systematicType: " << systematicType << endl; // Отладка
        BlastWaveFit *bwFit = new BlastWaveFit();
        cout << "  [Debug] BlastWaveFit object created" << endl;
        
        for (int part: PARTS)
        {
            cout << "  Processing part: " << part << endl;
            for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
                int centr = CENTR_SYST[systN][j];
                cout << "    Processing centr: " << centr << " (j=" << j << "/" << N_CENTR_SYST[systN] << ")" << endl;

                double parResultsRef[4];
                ReadParams(part, centr, parResultsRef, "output/parameters/SYST_" + systNamesT[systN] + ".txt"); 
                cout << "      T=" << parResultsRef[1] << ", beta=" << parResultsRef[2] << endl;

                ReadGlobalParams(systN, paramsGlobal, "output/parameters/ALL_GlobalBWparams_" + systNamesT[systN] + ".txt");
                // ReadGlobalParams(systN, paramsGlobal, "output/parameters/ALL_FinalBWparams_" + systNamesT[systN] + ".txt");
                getGlobalParams(part, centr, bwFit->paramsSystematics[part][centr]);
                setParamsForSys(systematicType, bwFit->paramsSystematics[part][centr], bwFit);
            }
        }

        cout << "  [Debug] Before Fit()" << endl;
        if (systematicType == -1) 
            bwFit->Fit(0);
        else
            bwFit->Fit(4);
        cout << "  [Debug] After Fit()" << endl;
    
        for (int part: PARTS)
        {
            for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
                int centr = CENTR_SYST[systN][j];

                double parResultsRef[4];
                ReadParams(part, centr, parResultsRef, "output/parameters/SYST_" + systNamesT[systN] + ".txt"); 

                for (int p = 0; p < 3; p++) {
                    // systErr[part][centr][p] +=  pow(bwFit->outParams[part][centr][p] / parResultsRef[p] - 1, 2);

                    double ref_value = parResultsRef[p];
                    double syst_value = bwFit->outParams[part][centr][p];

                    // Рассчитываем асимметричные ошибки
                    if(systematicType%2 == 1) { // Lower variations
                        double delta = std::max(ref_value - syst_value, 0.0);
                        systErr_low[part][centr][p] += delta*delta;
                    } else { // Upper variations
                        double delta = std::max(syst_value - ref_value, 0.0);
                        systErr_high[part][centr][p] += delta*delta;
                    }

                    const double CL_scale = 1.0; // Для 1σ
                    systErr[part][centr][p] = 0.5*(sqrt(systErr_low[part][centr][p]))/CL_scale 
                                            + 0.5*(sqrt(systErr_high[part][centr][p]))/CL_scale;
                }

                cout << "Syst " << systematicType 
                     << " T: " << bwFit->outParams[part][centr][1] 
                     << " beta: " << bwFit->outParams[part][centr][2] 
                     << " chi2: " << ifuncx[part][centr]->GetChisquare() << endl;
                    

                // cout << "PART: " << part << "   CENTR: " << centr << "  "
                //      << ifuncx[part][centr]->GetChisquare() / ifuncx[part][centr]->GetNDF()<< " " 
                //      << ifuncx[part][centr]->GetProb() << endl;

                cout << part << "  " << centr << "  "
                     << parResultsRef[1] << "  " << parResultsRef[2] << "  | "
                     << bwFit->outParams[part][centr][1] << "  " << sqrt(systErr[part][centr][1]) << "  |  "
                     << bwFit->outParams[part][centr][2] << "  " << sqrt(systErr[part][centr][2]) << endl;
            }
        }

        delete bwFit; // Освободить память
        cout << "  [Debug] bwFit deleted" << endl;
    }
    
    for (int part: PARTS)
    {
        for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
            int centr = CENTR_SYST[systN][j];

            for (int p = 0; p < 3; p++)
                // systErr[part][centr][p] = pow(systErr[part][centr][p], 0.5) / 2.;
                systErr[part][centr][p] = 0.5*(sqrt(systErr_low[part][centr][p]) + sqrt(systErr_high[part][centr][p]));

            cout << "PART: " << part << "   CENTR: " << centr
                 << "   systErr T: " << systErr[part][centr][1]
                 << "   systErr beta: " << systErr[part][centr][2] << endl;
        }
    }

    WriteParamsSyst(systN, bwFitRef->outParams, bwFitRef->outParamsErr, systErr, "output/parameters/SYSTErr_" + systNamesT[systN] + ".txt");
   
    if (!isDraw)
        return;

    // ++++++ Draw spectra +++++++++++++++++++++++++++++++++++++

    TCanvas *c2 = new TCanvas("c2", "c2", 29, 30, 1200, 1200);
    Format_Canvas(c2, 2, 3, 0);

    for (int i: PARTS)
    {
        c2->SetLogy();
        c2->cd(i + 1);  
        
        double shiftX = (i % 2 == 0) ? 0 : 0.1;
        double texScale = (i < 3) ? 1 : 0.9;
        TLegend *legend = new TLegend(0.5 - shiftX, 0.7, 0.95 - shiftX, 0.9); //1 column
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetNColumns(2);
        legend->SetTextSize(0.072 * texScale);

        TLatex *titleTex = new TLatex(0.4, 500, partTitles[i].c_str());
        titleTex->SetTextFont(42);
        titleTex->SetTextSize(0.09);
        titleTex->SetLineWidth(2 * texScale);

        FormatSpectraPad(texScale);

        for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
            int centr = CENTR_SYST[systN][j];

            if (!ifuncx[i][centr]) continue;

            grSpectra[i][centr]->SetMarkerColor(centrColors[centr]);
            grSpectra[i][centr]->SetMarkerSize(1);
            grSpectra[i][centr]->SetMarkerStyle(8);
            grSpectra[i][centr]->Draw("P SAME");      
            ifuncx[i][centr]->Draw("SAME");
            legend->AddEntry(grSpectra[i][centr], centrTitles[centr].c_str(), "p");        
        }
        legend->Draw();
        titleTex->Draw();    
    }

    // c2->SaveAs("output/BlastWaveSyst.png");
    delete c2;
}

