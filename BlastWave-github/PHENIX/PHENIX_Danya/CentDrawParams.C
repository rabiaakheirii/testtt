#include "input/headers/def.h"
#include "input/headers/WriteReadFiles.h"


double avgT, avgTerr;
double avgUt, avgUtErr;


void DrawParam(string paramName = "T", bool isSyst = true)
{
    TGraphErrors *gr[N_PARTS], *grSys[N_PARTS];
    double xerr[N_CENTR], xerrSys[N_CENTR];

    for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
        int centr = CENTR_SYST[systN][j];
        xerr[j] = 0., xerrSys[j] = 1;
    }

    for (int part: PARTS)
    {
        if (paramName == "T")
            gr[part] = new TGraphErrors(N_CENTR, centrX, Tpar[part], xerr, Tpar_err[part]);
        else if (paramName == "beta")
            gr[part] = new TGraphErrors(N_CENTR, centrX, utPar[part], xerr, utPar_err[part]);
        
        gr[part]->SetMarkerStyle(8);
        gr[part]->SetMarkerSize(2);
        gr[part]->SetMarkerColor(partColors[part]);

        if (isSyst)
        {
            if (paramName == "T")
                grSys[part] = new TGraphErrors(N_CENTR, centrX, Tpar[part], xerrSys, Tpar_sys[part]);
            else if (paramName == "beta")
                grSys[part] = new TGraphErrors(N_CENTR, centrX, utPar[part], xerrSys, utPar_sys[part]);
            
            grSys[part]->SetLineColorAlpha(partColors[part], 0.6);
            grSys[part]->SetFillStyle(0);
            grSys[part]->SetFillColorAlpha(partColors[part], 0.5);
            grSys[part]->SetLineWidth(2);
            grSys[part]->SetMarkerColorAlpha(partColors[part], 0.6);
        }
    }

    TCanvas *c2 = new TCanvas("c2", "c2", 30, 30, 1200, 1000);
    c2->cd();
    c2->SetGrid();
    // c2->SetLogx();

    double ll = 0, rl = 100., pad_min = 0., pad_max = (paramName == "T") ? 0.3 : 1., 
        pad_offset_x = 1., pad_offset_y = 1., 
        pad_tsize = 0.05, pad_lsize=0.05;
    TString pad_title_y = (paramName == "T") ? "T [GeV]" : "#beta";
    TString pad_title_x = "centr. [%]";
    Format_Pad(ll, rl, pad_min, pad_max, pad_title_x, pad_title_y, pad_offset_x, pad_offset_y, pad_tsize, pad_lsize, "", 8);        
    
    TLegend *legend = new TLegend(0.2, 0.15, 0.4, 0.30);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetNColumns(2);
    legend->SetTextSize(0.05);

    if (paramName == "T")
    {
        TLine *lineT = new TLine(ll, avgT, rl, avgT); 
        lineT->SetLineColor(kBlack);
        lineT->SetLineWidth(2);
        lineT->SetLineStyle(9);
        lineT->Draw("SAME");
    }
    else if (paramName == "beta")
    {
        TLine *lineT = new TLine(ll, avgUt, rl, avgUt); 
        lineT->SetLineColor(kBlack);
        lineT->SetLineWidth(2);
        lineT->SetLineStyle(9);
        lineT->Draw("SAME");
    }

    for (int part: PARTS)
    {
       gr[part]->Draw("P SAME");
       if (isSyst) grSys[part]->Draw("P2");

       legend->AddEntry(gr[part], partTitles[part].c_str(), "P");
    }
    legend->Draw();
    c2->SaveAs(("output/pics/BWparFinal_" + paramName + ".png").c_str());
}


void CentDrawParams ( void )
{
    ReadParam(1, Tpar, Tpar_err, Tpar_sys);
    ReadParam(2, utPar, utPar_err, utPar_sys);

    int count = 0;
    for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
        int centr = CENTR_SYST[systN][j];
        for (int part = 0; part < 6; part++)
        {
            avgT += Tpar[part][centr];
            avgTerr += pow(Tpar_sys[part][centr], 2);
    
            avgUt += utPar[part][centr];
            avgUtErr += pow(utPar_sys[part][centr], 2);
    
            count++;
        }
    }
    avgT = avgT / double(count);
    avgTerr = sqrt(avgTerr) / double(count);
    
    avgUt = avgUt / double(count);
    avgUtErr = sqrt(avgUtErr) / double(count);
    
    cout << "T = " << avgT << " ± " << avgTerr << endl;
    cout << "beta = " << avgUt << " ± " << avgUtErr << endl;


    DrawParam("T");
    DrawParam("beta");

    gROOT->ProcessLine(".q");
}