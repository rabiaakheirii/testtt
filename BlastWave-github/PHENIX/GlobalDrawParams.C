#include "def.h"
#include "WriteReadFiles.h"

TGraphErrors *gr[2][5]; // [charge][system]
TString systNamesT[5] = {"AuAu", "pAl", "HeAu", "CuAu", "UU"};
int SYSTS[] = {0, 1, 2, 3, 4};

void DrawParam(TString paramName = "T")
{
    TCanvas *c2 = new TCanvas("c2", "c2", 29, 30, 1200, 1000);
    c2->cd();
    c2->SetGrid();
    c2->SetLogx();
    double ll = 10, rl = 500., pad_min = 0., pad_max = (paramName == "T") ? 0.3 : 1., 
        pad_offset_x = 1., pad_offset_y = 1., 
        pad_tsize = 0.05, pad_lsize=0.05;
    TString pad_title_y = "T [GeV]";
    TString pad_title_x = "p_{T} [GeV/c]";
    Format_Pad(ll, rl, pad_min, pad_max, pad_title_x, pad_title_y, pad_offset_x, pad_offset_y, pad_tsize, pad_lsize, "", 8);        
    
    TLegend *legend = new TLegend(0.2, 0.7, 0.6, 0.85);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetNColumns(2);
    legend->SetTextSize(0.05);

    for (int systN: SYSTS)
    {
        for (int charge: {0, 1})
        {
            gr[charge][systN]->Draw("P SAME");
            legend->AddEntry(gr[charge][systN], systNamesT[systN], "P");
        }
    }
    
    legend->Draw();
    cout << "test" << endl;
    TString name = "output/BWparamGlobal_" + paramName + ".pdf";
    c2->SaveAs(name);
}


void SetGraphs( int systN, TString paramName )
{
    cout << systNamesT[systN] << endl;
    TString filename = "output/GlobalBWparams_" + systNamesT[systN] + ".txt";
    ReadGlobalParams(paramsGlobal, "output/txtParams/GlobalBWparams_" + systNamesT[systN] + ".txt");
    cout << filename << "  " <<  N_CENTR_SYST[systN] << endl;
    
    for (int charge: {0, 1})
    {
        for (int centr = 0; centr < N_CENTR_SYST[systN]; centr++)
        {
            // cout << charge << " " << centr << " " << Tpar[charge][centr] << endl;

            Tpar[charge][centr] = paramsGlobal[charge][centr][0];
            utPar[charge][centr] = paramsGlobal[charge][centr][1];

        }
    }
    cout <<"DEFINE GRAPHS " << N_CENTR_SYST[systN] << endl;
    
    double xerr[MAX_CENTR];
    for (int i: CENTR_SYST[systN]) xerr[i] = 0.;

    for (int charge: {0, 1})
    {
        if (paramName == "T")
            gr[charge][systN] = new TGraphErrors(N_CENTR_SYST[systN], Npart[systN], Tpar[charge], xerr, xerr); // Tpar_err[part]);
        else if (paramName == "ut")
            gr[charge][systN] = new TGraphErrors(N_CENTR_SYST[systN], Npart[systN], utPar[charge], xerr, xerr); // utPar_err[part]);
        
        gr[charge][systN]->SetMarkerStyle(8);
        gr[charge][systN]->SetMarkerSize(2);
        gr[charge][systN]->SetMarkerColor(systColors[systN]);
    }

}

void GlobalDrawParams ( void )
{
    
    for (int systN: SYSTS)
    {
        cout << systN << " " << systNamesT[systN] << endl;
        SetGraphs(systN, "ut");
    }
    
    DrawParam("ut");
}