#include "def.h"
#include "Cumulative.h"
#include "HagedornFit.h"
#include "LevyFit.h"

using namespace std;

void spectra( int systN = 0 )
{
    string fitFunc = "levy";

    // ++++++ Read data +++++++++++++++++++++++++++++++++++++
    string inputFileName = "postprocess_mpdpid10";
    // string inputFileName = "postprocess_test-XeW";
    SetSpectra(inputFileName, "pt");

    if (fitFunc == "levy")
         FitLevy();
    else
        FitHagedorn();

    // ++++++ Draw spectra +++++++++++++++++++++++++++++++++++++

    TCanvas *c2 = new TCanvas("c2", "c2", 29, 30, 1200, 1200);
    Format_Canvas(c2, 2, 3, 0);

    for (int i: PARTS)
    {
        c2->SetLogy();
        c2->cd(i + 1);  
        
        double shiftX = (i % 2 == 0) ? 0 : 0.1;
        double texScale = (i < 3) ? 1 : 0.9;

        TLegend *legend = new TLegend(0.55 - shiftX, 0.7, 0.98 - shiftX, 0.9); //1 column
        legend->SetNColumns(2);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.07 * texScale);

        TLatex *titleTex = new TLatex(0.6, 500, partTitles[i].c_str());
        titleTex->SetTextFont(42);
        titleTex->SetTextSize(0.08);
        titleTex->SetLineWidth(2 * texScale);

        FormatSpectraPad(texScale);
        for (int centr: CENTR)
        {
            if (!grSpectra[i][centr]) continue;
            grSpectra[i][centr]->Draw("SAME");      
            legend->AddEntry(grSpectra[i][centr], centrTitles[centr].c_str(), "l");        
        }
        legend->Draw();
        titleTex->Draw();    
        
       // DrawCumulativeBorder(i, PAD_MIN, PAD_MAX);
    }

    c2->SaveAs(("output/spectra_" + fitFunc + "_" + inputFileName + ".pdf").c_str());
    delete c2;
}

