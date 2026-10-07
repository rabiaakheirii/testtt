#include "TGraph.h"
#include "input/FormatOfEverything.h"
TGraph *grTbeta;
int centrN[5] = {4, 5, 5, 4, 11};
TString systName[5] = {"pAl", "HeAu", "CuAu", "UU", "AuAu"};

void ReadPHENIXData( int systN )
{
    TString filename = "input/PHENIX_Tbeta/ALL_GlobalBWparams_" + systName[systN] + ".txt";
    ifstream f;
    f.open(filename);

    double tmp;
    int tmpInt;
    for (int charge: {0, 1})
    {
        for (int centr = 0; centr < centrN[systN]; centr++)
        {
            double T, beta, C;
            f >> tmpInt >> tmpInt >> T >> beta >> tmp >> tmp >> tmp;
            int Npoint = grTbeta ? grTbeta->GetN() : 1;
            cout << Npoint << "  " << T << "  " << beta << endl;
            grTbeta->SetPoint(Npoint, beta,T);
        }
    }
    f.close();
}

void TbetaDependency ( void )
{
    grTbeta = new TGraph(1);
    for (int syst = 0; syst < 5; syst++)
        ReadPHENIXData(syst);
    
    TF1 * f = new TF1("f","[0] + x * [1] + x * x * [2]", 0.25, 0.9);
    f->SetParameters(0.2, 1, -1);
    grTbeta->Fit(f, "R");
    
    double par[3];
    f->GetParameters(par);
    cout << par[0] << ", " << par[1] << ", " << par[2] << endl;
    grTbeta->SetMarkerStyle(8);
    grTbeta->SetMarkerSize(2);
    grTbeta->SetMarkerColor(kBlue);


    TCanvas *c2 = new TCanvas("c2", "c2", 29, 30, 1200, 1000);
    c2->cd();
    c2->SetGrid();
    // c2->SetLogx();
    double ll = 0, rl = 1., pad_min = 0., pad_max = 0.2, 
        pad_offset_x = 1., pad_offset_y = 1., 
        pad_tsize = 0.05, pad_lsize=0.05;
    TString pad_title_y = "T [GeV]";
    TString pad_title_x = "#beta_{T}";
    
    Format_Pad(ll, rl, pad_min, pad_max, pad_title_x, pad_title_y, pad_offset_x, pad_offset_y, pad_tsize, pad_lsize, "", 8);        

    TLatex *titleTex = new TLatex(0.2, 0.05, "PHENIX");
    titleTex->SetTextFont(42);
    titleTex->SetTextSize(0.09);
    titleTex->SetLineWidth(2);
    titleTex->Draw();
    grTbeta->Draw("P SAME");

    c2->SaveAs("output/TbetaDependency.pdf");
}