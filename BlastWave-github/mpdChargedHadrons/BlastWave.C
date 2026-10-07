#include "def.h"
#include "WriteReadFiles.h"
#include "BlastWaveFit.h"

using namespace std;

void BlastWave( void )
{
    bool isContour = true;
    bool isDraw = true;

    // ++++++ Read data +++++++++++++++++++++++++++++++++++++

    SetSpectraMalaev("mt");

    // +++++++++ Fit +++++++++++++++++++++++++++++++++++++++
    BlastWaveFit *bwFit = new BlastWaveFit();
    bwFit->isContour = isContour;
    bwFit->Fit(0);
    bwFit->GetIntegralDiff();
    
    WriteParams(bwFit->outParams, bwFit->outParamsErr);
   
    if (!isDraw)
        return;

    // ++++++ Draw spectra +++++++++++++++++++++++++++++++++++++

    TCanvas *c2 = new TCanvas("c2", "c2", 29, 30, 2000, 1200);
    Format_Canvas(c2, 3, 2, 0);
    double shiftXarr[6] = {0, 0.1, 0.18, 0, 0.1, 0.18};
      
    int padN = 0;
    for (int i: {0, 2, 4, 1, 3, 5})
    {
        padN++;
        c2->SetLogy();
        c2->cd(padN);  
      
        double shiftX = shiftXarr[padN - 1];
        
        // double shiftX = (i % 2 == 0) ? 0 : 0.1;
        double texScale = (padN <= 3) ? 1 : 0.9;
        cout << "texscale " << texScale << endl;
        TLegend *legend = new TLegend(0.45 - shiftX, 0.6, 0.95 - shiftX, 0.9); //1 column
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetNColumns(2);
        legend->SetTextSize(0.065 * texScale);

        TLatex *titleTex = new TLatex(0.4, 500, partTitles[i].c_str());
        titleTex->SetTextFont(42);
        titleTex->SetTextSize(0.09);
        titleTex->SetLineWidth(2 * texScale);

        FormatSpectraPad(texScale);
        for (int centr: CENTR_MALAEV)
        {
            
            // if (!ifuncx[i][centr]) continue;
            grSpectra[i][centr]->SetMarkerColor(centrColors[centr]);
            grSpectra[i][centr]->SetMarkerSize(0.5);
            grSpectra[i][centr]->SetMarkerStyle(8);
            grSpectra[i][centr]->SetLineColor(centrColors[centr]);
            grSpectra[i][centr]->Draw("P SAME");      
            // if (ifuncx[i][centr]) ifuncx[i][centr]->Draw("SAME");
            legend->AddEntry(grSpectra[i][centr], centrTitlesMalaev[centr].c_str(), "l");   
        }
        legend->Draw();
        titleTex->Draw();    
    }

    c2->SaveAs("output/BlastWave_Malaev.pdf");
    delete c2;

    if (!bwFit->isContour) {
        gROOT->ProcessLine(".q");
        return;
    }

    //++++++++ Draw Contour plots ++++++++++++++++++++++++++++++

    TCanvas *c3 = new TCanvas("c3", "c3", 29, 30, 1200, 1200);
    c3->Divide(1,1);
    c3->cd(1);
    
    TPad *c3_1= (TPad*) c3->GetListOfPrimitives()->FindObject("c3_1");
    c3_1->SetLeftMargin(0.1154148);
    c3_1->SetRightMargin(0.02);
    c3_1->SetTopMargin(0.02);
    c3_1->SetBottomMargin(0.09);
    c3_1->SetGrid();

    double ll = 0.001, rl = 1.5, pad_min = 0.01, pad_max = 0.180, 
        pad_offset_x = 0.85, pad_offset_y = 1.2, 
        pad_tsize = 0.05, pad_lsize=0.04;
    TString pad_title_y = "T";
    TString pad_title_x = "#beta";
    Format_Pad(ll, rl, pad_min, pad_max, pad_title_x, pad_title_y, pad_offset_x, pad_offset_y, pad_tsize, pad_lsize, "", 8);        

    TLegend *legendContour = new TLegend(0.2, 0.1, 0.85, 0.479);
    legendContour->SetBorderSize(0);
    legendContour->SetFillStyle(0);
    legendContour->SetNColumns(3);
    legendContour->SetTextSize(0.03);

    for (int centr: CENTR_MALAEV)
    {
        for (int part: {0,2,4})
        {
            string legendText = partTitles[part] + ", " + centrTitlesMalaev[centr];
            legendContour->AddEntry(contour[part][centr][1], legendText.c_str(), "l"); 

            for (int nsigma = 1; nsigma < N_SIGMA; nsigma++)
            {
                if (contour[part][centr][nsigma])
                    contour[part][centr][nsigma]->Draw("lf"); //"CONT3"
                else 
                    cout << "ERROR " << part << " " <<centr << " " << nsigma << endl;
            }
        }
    }
    legendContour->Draw();

    c3->SaveAs("output/BlastWave_contour.pdf");

    gROOT->ProcessLine(".q");
}

