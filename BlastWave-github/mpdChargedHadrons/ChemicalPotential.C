#include "TF1.h"
#include "TLegend.h"
#include "input/FormatOfEverything.h"
using namespace std;

const double PI = 3.1415;
const int N_CENTR = 3;


int markerStyle[6] = {8, 21, 22, 23, 33};
string particles[6] = {"pip", "pim", "kp", "km", "p", "ap"};
Color_t centrColors[5] = {kMagenta + 1, kBlue, kGreen + 2};
string centrTitles[5] = {"0-20%", "20-40%", "40-60%", "60-80%"};
TH1D *hSpectra[6][5];
TH1D *hRatioP[N_CENTR];

double ratPHENIX = 0.73, ratMPD;

double tempFunc( double x )  
{ 
    double B = pow(220, 4);
    double Sq = sqrt(340 * PI * PI * B + 55 * pow(x, 4));
    if (Sq - 15 * x * x > 0)
        return  (1 / PI) * sqrt(3/34.) * sqrt( Sq - 15 * x * x); 
    else 
        return 0;
}

void DrawLineOnPhaseDiagram(double parValue, Color_t color, string title, TLegend *legend)
{
    TF1 *expTmuFunc = new TF1("ratioF", "-2. * x / log([0])", 0, 500);
    expTmuFunc->SetParameter(0, parValue);
    expTmuFunc->SetLineColor(color);
    expTmuFunc->SetLineWidth(3);
    expTmuFunc->Draw("SAME");
    legend->AddEntry(expTmuFunc, title.c_str(), "l");
}

void DrawKineticFreezeOut()
{
    TF1 *expTmuFunc = new TF1("ratioF", "-2. * x / log([0])", 0, 500);
    double T0[2] = {118, 108.652}, mu[2];
    expTmuFunc->SetParameter(0, ratPHENIX);
    mu[0] = expTmuFunc->GetX(T0[0]); // PHENIX
    expTmuFunc->SetParameter(0, ratMPD);
    mu[1] = expTmuFunc->GetX(T0[1]); // MPD

    cout << " \n Kinetic freeze-out" << endl;
    cout << "PHENIX  " << mu[0] << "  " << T0[0] << endl;
    cout << "MPD     " << mu[1] << "  " << T0[1] << endl;

    TGraph *gr = new TGraph(2, mu, T0);
    gr->SetMarkerStyle(8);
    gr->SetMarkerSize(3);
    gr->SetMarkerColor(kBlue);
    gr->Draw("P SAME");
}

void DrawChemicalFreezeOut()
{
    TF1 *expTmuFunc = new TF1("ratioF", "-2. * x / log([0])", 0, 500);
    double mu[2] = {24.87, 247.}, T[2];
    expTmuFunc->SetParameter(0, ratPHENIX);
    T[0] = expTmuFunc->Eval(mu[0]); // PHENIX
    expTmuFunc->SetParameter(0, ratMPD);
    T[1] = expTmuFunc->Eval(mu[1]); // MPD

    cout << "\n Chemical freeze-out" << endl;
    cout << "PHENIX  " << mu[0] << "  " << T[0] << endl;
    cout << "MPD     " << mu[1] << "  " << T[1] << endl;

    TGraph *gr = new TGraph(2, mu, T);
    gr->SetMarkerStyle(8);
    gr->SetMarkerSize(3);
    gr->SetMarkerColor(kRed);
    gr->Draw("P SAME");
}
void DrawPhaseDiagram( double ratioP[N_CENTR] )
{
    TCanvas *c3 = new TCanvas("c3", "c3", 29, 30, 1200, 1200);
    c3->cd();
    c3->SetGrid();
    double ll = 0, rl = 500, pad_min = 0, pad_max = 199, 
        pad_offset_x = 1., pad_offset_y = 1., 
        pad_tsize = 0.05, pad_lsize=0.05;
    TString pad_title_y = "T";
    TString pad_title_x = "#mu";
    Format_Pad(ll, rl, pad_min, pad_max, pad_title_x, pad_title_y, pad_offset_x, pad_offset_y, pad_tsize, pad_lsize, "", 4);        

    auto f = new TF1("tempFunc", "tempFunc(x)",0, 500);
    f->SetLineWidth(3);
    f->SetMarkerSize(3);
    f->Draw("SAME");

    TLegend *legend = new TLegend(0.4, 0.15, 0.6, 0.35);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.04);

    for (int centr = 0; centr < N_CENTR; centr++)
    {
        DrawLineOnPhaseDiagram(ratioP[centr], centrColors[centr], centrTitles[centr], legend);
    }
    DrawLineOnPhaseDiagram(ratPHENIX, kBlack, "PHENIX", legend);
    DrawKineticFreezeOut();
    DrawChemicalFreezeOut();

    TString title = "#frac{#sqrt{3/34}}{#pi} #sqrt{#sqrt{340#pi^{2}(220)^{4} + 55#mu^{4}}-15#mu^{2}}";
    TLatex *titleTex = new TLatex(50, 170, title);
    titleTex->SetTextFont(42);
    titleTex->SetTextSize(0.05);
    titleTex->SetLineWidth(2);
    titleTex->Draw();

    legend->Draw();

    c3->SaveAs("output/PhaseDiagram.pdf");
}

double CalculateIntRatio( int centr )
{
    int Nbins = 33; // hRatioP[centr]->GetNbinsX();
    double Int = 0;
    for (int bin = 1; bin < Nbins; bin++)
    {
        Int += hRatioP[centr]->GetBinContent(bin);
    }
    return Int / (Nbins - 1);
}

void ChemicalPotential( void )
{
    // ++++++ Read data +++++++++++++++++++++++++++++++++++++

    // TFile *f = new TFile("input/nuclei_ptspectra.root");
    TFile *f = new TFile("input/postprocess_mpdpid10.root");
    TDirectory *fd;

    for (int i = 4; i < 6; i++)
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
        hRatioP[centr] = (TH1D *)hSpectra[5][centr]->Clone("RatioP");
        hRatioP[centr]->Divide(hSpectra[4][centr]);

        ratioP[centr] = CalculateIntRatio(centr);
        // fitF->SetLineColor(centrColors[centr]);
        // hRatioP[centr]->Fit("fitF");
        // 
        // ratioP[centr] = fitF->GetParameter(0);
        cout << "centr  " << centr << " " << ratioP[centr] << endl;
    }

    // ++++++ Draw data +++++++++++++++++++++++++++++++++++++

    TCanvas *c2 = new TCanvas("c2", "c2", 29, 30, 1200, 1000);
    //c2->Divide(1,1);
    
    Format_Canvas(c2, 1, 1, 1);
    c2->cd();
    c2->SetGrid();
    double ll = 0.00001, rl = 1.9, pad_min = 0., pad_max = 0.08, 
        pad_offset_x = 1., pad_offset_y = 1., 
        pad_tsize = 0.05, pad_lsize=0.05;
    TString pad_title_y = "#bar{p}/p";
    TString pad_title_x = "p_{T} [GeV/c]";
    Format_Pad(ll, rl, pad_min, pad_max, pad_title_x, pad_title_y, pad_offset_x, pad_offset_y, pad_tsize, pad_lsize, "", 8);        
    
    /*c2->Divide(1,1);
    c2->cd(1);
    TPad pad1("pad1", "pad1", 0, 0., 1, 1);
    pad1.SetLogy(true);
    pad1.Draw();
    pad1.cd();

    pad1.SetLeftMargin(0.13);
    pad1.SetRightMargin(0.02);
    pad1.SetTopMargin(0.02);
    pad1.SetBottomMargin(0.02);
*/
    TLegend *legend = new TLegend(0.2, 0.65, 0.5, 0.85);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.04);

    for (int centr = 0; centr < N_CENTR; centr++)
    {
        hRatioP[centr]->SetMarkerColor(centrColors[centr]);
        hRatioP[centr]->SetMarkerSize(2);
        hRatioP[centr]->SetMarkerStyle(markerStyle[centr]);
        hRatioP[centr]->SetLineColor(centrColors[centr]);
        hRatioP[centr]->Draw("SAME");
        legend->AddEntry(hRatioP[centr], centrTitles[centr].c_str(), "pl");

    }
    legend->Draw();
    c2->SaveAs("output/RatioP.pdf");
    ratMPD = ratioP[0];
    DrawPhaseDiagram(ratioP);
}