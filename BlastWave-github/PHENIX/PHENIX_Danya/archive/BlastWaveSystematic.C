#include "def.h"
#include "WriteReadFiles.h"
#include "BlastWaveFit.h"

using namespace std;

void setParamsForSys ( int systematicType, double parResults[4], BlastWaveFit *bwFit )
{
    if (systematicType == 0 ) return;

    switch(systematicType)
    {
        case 1: 
            parResults[1] *= 0.8;
        case 2: 
            parResults[1] *= 1.2;
        case 3:
            parResults[2] *= 0.8;
        case 4:
            parResults[2] *= 1.2;            
        case 5:
            parResults[0] *= 0.1;
        case 6:
            parResults[0] *= 10;
        case 7:
            bwFit->lLimitMult *= 0.8;
            bwFit->rLimitMult *= 0.8;
        case 8:
            bwFit->lLimitMult *= 1.2;
            bwFit->rLimitMult *= 1.2;
        case 9:
            bwFit->lLimitMultPi *= 1.2;        
            bwFit->rLimitMultPi *= 0.8;
        case 10:
            bwFit->lLimitMultPi *= 1.2;
            bwFit->rLimitMultPi *= 1.2;
        // case 8:
        //     xmin[part] *= 1.1;
        //     xmax[part] *= 1.1;
    }
}


void BlastWaveSystematic( void )
{
    bool isDraw = true;

    double systErr[N_PARTS][N_CENTR][4]; 

    for (int part: PARTS)
        for (int centr: CENTR) 
            for (int p = 0; p < 4; p++)
                systErr[part][centr][p] = 0;
                
    // ++++++ Read data +++++++++++++++++++++++++++++++++++++

    string inputFileName = "postprocess_mpdpid10";
    SetSpectra(inputFileName, "mt");

    // +++++++++ Fit +++++++++++++++++++++++++++++++++++++++

    BlastWaveFit *bwFitRef = new BlastWaveFit();
    bwFitRef->Fit(0);
    WriteParams(bwFitRef->outParams, bwFitRef->outParamsErr);

    
    for (int systematicType: {0, -1})
    {
        BlastWaveFit *bwFit = new BlastWaveFit();

        for (int part: PARTS)
        {
            for (int centr: CENTR) 
            {
                double parResultsRef[4];
                ReadParams(part, centr, parResultsRef); 

                ReadGlobalParams(paramsGlobal);
                getGlobalParams(part, centr, bwFit->paramsSystematics[part][centr]);
                setParamsForSys(systematicType, bwFit->paramsSystematics[part][centr], bwFit);
            }
        }

        if (systematicType == -1) 
            bwFit->Fit(2);
        else
            bwFit->Fit(4);
    
        for (int part: PARTS)
        {
            for (int centr: CENTR) 
            {
                double parResultsRef[4];
                ReadParams(part, centr, parResultsRef); 

                for (int p = 0; p < 3; p++)
                    systErr[part][centr][p] +=  pow(bwFit->outParams[part][centr][p] / parResultsRef[p] - 1, 2);

                // cout << "PART: " << part << "   CENTR: " << centr << "  "
                //      << ifuncx[part][centr]->GetChisquare() / ifuncx[part][centr]->GetNDF()<< " " 
                //      << ifuncx[part][centr]->GetProb() << endl;

                cout << part << "  " << centr << "  "
                     << parResultsRef[1] << "  " << parResultsRef[2] << "  | "
                     << bwFit->outParams[part][centr][1] << "  " << sqrt(systErr[part][centr][1]) << "  |  "
                     << bwFit->outParams[part][centr][2] << "  " << sqrt(systErr[part][centr][2]) << endl;
            }
        }
    }
    
    for (int part: PARTS)
    {
        for (int centr: CENTR) 
        {
            for (int p = 0; p < 3; p++)
                systErr[part][centr][p] = pow(systErr[part][centr][p], 0.5) / 2.;

            cout << "PART: " << part << "   CENTR: " << centr
                 << "   systErr T: " << systErr[part][centr][1]
                 << "   systErr beta: " << systErr[part][centr][2] << endl;
        }
    }

    WriteParamsSyst(bwFitRef->outParams, bwFitRef->outParamsErr, systErr);
   
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
        for (int centr: CENTR)
        {
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

    c2->SaveAs("output/BlastWaveSyst.pdf");
    delete c2;
}

