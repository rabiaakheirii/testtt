#include "TF1.h"
#include "TLegend.h"
#include "input/FormatOfEverything.h"
using namespace std;

const double PI = 3.1415;
const int N_CENTR = 3;

string particles[6] = {"pip", "pim", "kp", "km", "p", "ap"};
Color_t centrColors[5] = {kMagenta + 1, kBlue, kGreen + 2};
string centrTitles[5] = {"0-20%", "20-40%", "40-60%", "60-80%"};
TH1D *hSpectra[6][5];
TH1D *hRatio[N_CENTR];


void ChemicalPotentialStrange( void )
{
    // ++++++ Read data +++++++++++++++++++++++++++++++++++++

    // TFile *f = new TFile("input/nuclei_ptspectra.root");
    TFile *f = new TFile("input/postprocess_mpdpid10.root");
    TDirectory *fd;

    for (int i = 0; i < 6; i++)
    {
        fd = (TDirectory*)f->Get(particles[i].c_str());
        fd->cd();
        for (int centr = 0; centr < N_CENTR; centr++)
        {
            string name = "h__pt_" + particles[i] +"_centrality" + to_string(centr) + "_mc_y-0.5_0.5";
            cout << name << endl;
            hSpectra[i][centr] = (TH1D *)fd->Get(name.c_str());    
        }
    }

    TF1 *fitF = new TF1("fitF", "[0]", 0, 3);
    fitF->SetParameter(0, 0.1);
    double ratioP[N_CENTR];

    for (int centr = 0; centr < N_CENTR; centr++)
    {
        hRatio[centr] = (TH1D *)hSpectra[3][centr]->Clone("RatioP");
        hRatio[centr]->Divide(hSpectra[2][centr]);

        fitF->SetLineColor(centrColors[centr]);
        hRatio[centr]->Fit("fitF");
        
        ratioP[centr] = fitF->GetParameter(0);
        cout << "centr  " << centr << " " << ratioP[centr] << endl;
    }

    // ++++++ Draw data +++++++++++++++++++++++++++++++++++++

    TCanvas *c2 = new TCanvas("c2", "c2", 29, 30, 1200, 1000);
    c2->cd();
    c2->SetGrid();
    double ll = 0.00001, rl = 1.9, pad_min = 0., pad_max = 1, 
        pad_offset_x = 1., pad_offset_y = 1., 
        pad_tsize = 0.05, pad_lsize=0.05;
    TString pad_title_y = "K^{-}/K^{+}";
    TString pad_title_x = "p_{T} [GeV/c]";
    Format_Pad(ll, rl, pad_min, pad_max, pad_title_x, pad_title_y, pad_offset_x, pad_offset_y, pad_tsize, pad_lsize, "", 8);        
    
    TLegend *legend = new TLegend(0.2, 0.65, 0.5, 0.85);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.04);

    for (int centr = 0; centr < N_CENTR; centr++)
    {
        hRatio[centr]->SetLineColor(centrColors[centr]);
        hRatio[centr]->Draw("SAME");
        legend->AddEntry(hRatio[centr], centrTitles[centr].c_str(), "l");

    }
    legend->Draw();
    c2->SaveAs("output/RatioK.pdf");

    double Tch = 177;
    double Tk = 110;
    double RatioValue = fitF->GetParameter(0);
    double muS_ch = -Tch / 2 * log(RatioValue);
    double muS_k = -Tk / 2 * log(RatioValue);
    cout << "muS chemical: " << muS_ch << "   muS kinetic: " << muS_k << endl;
}