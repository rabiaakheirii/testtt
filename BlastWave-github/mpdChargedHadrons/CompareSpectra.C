#include "def.h"
#include "WriteReadFiles.h"
#include "BlastWaveFit.h"

using namespace std;

TH1D *hSpectraM[6][12], *hSpectraK[6][12], *hdiff[6][12];

void SetHistSpectraMalaev(string type = "pt")
{
    int centrBins[] = {0, 10, 20, 30, 40, 50, 60, 90};
    int N_CENTR_M = 8;
    string inputFileName; 
    
    string partMalaev[6] = {"Pion_pl", "Pion_mn", "Kaon_pl", "Kaon_mn", "Proton_pl", "Proton_mn"};

    for (int charge = 0; charge < 2; charge++)
    {
        inputFileName = (charge == 0) ? "Malaev_Spectra_positive_BiBi_9_2_GeV" : "Malaev_Spectra_negative_BiBi_9_2_GeV";
        TFile *f = new TFile(("input/" + inputFileName + ".root").c_str());

        for (int i = charge; i < 6; i += 2)
        {
            for (int centr = 0; centr < N_CENTR_M - 1; centr++)
            {
                string name = partMalaev[i] +"_rec_" + to_string(centrBins[centr]) + "_" + to_string(centrBins[centr + 1]);
                hSpectraM[i][centr] = (TH1D *)f->Get(name.c_str());    

                for (int bin = 1; bin < hSpectraM[i][centr]->GetNbinsX(); bin++)
                {
                    double content = hSpectraM[i][centr]->GetBinContent(bin);
                    double ptval = hSpectraM[i][centr]->GetBinCenter(bin);
                    double width = hSpectraM[i][centr]->GetBinWidth(bin);
                    hSpectraM[i][centr]->SetBinContent(bin, content * ptval * width);
                }
                hSpectraM[i][centr]->Rebin(2);
                for (int bin = 1; bin < hSpectraM[i][centr]->GetNbinsX(); bin++)
                {
                    double content = hSpectraM[i][centr]->GetBinContent(bin);
                    double ptval = hSpectraM[i][centr]->GetBinCenter(bin);
                    double width = hSpectraM[i][centr]->GetBinWidth(bin);
                    hSpectraM[i][centr]->SetBinContent(bin, content /(ptval * width));
                }
                if (!hSpectraM[i][centr]) continue;

                // hSpectraM[i][centr]->SetMarkerColor(centrColors[centr]);
                // hSpectraM[i][centr]->SetMarkerStyle(8); 
                // hSpectraM[i][centr]->SetMarkerSize(1);
            }
        }
        
    }
}


void SetHistSpectraKireev(string inputFileName = "postprocess_mpdpid10", string type = "pt")
{
    TFile *f1 = new TFile(("input/" + inputFileName + ".root").c_str());
    TDirectory *fd;

    for (int i = 0; i < 6; i++)
    {
        cout << i << endl;
        fd = (TDirectory*)f1->Get(particles[i].c_str());
        fd->cd();
        for (int centr = 0; centr < N_CENTR - 1; centr++)
        {
            string name = "h__pt_" + particles[i] +"_centrality" + to_string(centr) + "_mc_y-0.5_0.5";
            hSpectraK[i][centr] = (TH1D *)fd->Get(name.c_str());    
            if (!hSpectraK[i][centr]) continue;
            
            // hSpectraK[i][centr]->SetMarkerColor(centrColors[centr]);
            // hSpectraK[i][centr]->SetMarkerStyle(21);
            // hSpectraK[i][centr]->SetMarkerSize(1);
        }
    }
}

void CompareSpectra( void )
{
    SetHistSpectraMalaev();
    SetHistSpectraKireev();
    
    int part = 0;
    TCanvas *c2 = new TCanvas("c2", "c2", 1000, 1000);
    Format_Canvas(c2, 1, 1, 0);
    c2->cd(1);

    TPad pad1("pad1", "pad1", 0, 0.3, 1, 1);
    pad1.SetLogy(true);
    pad1.Draw();
    pad1.cd();

    pad1.SetLeftMargin(0.13);
    pad1.SetRightMargin(0.02);
    pad1.SetTopMargin(0.02);
    pad1.SetBottomMargin(0.2);

    double ll = 0.01, rl = 2.49, pad_min = 0.00009, pad_max = 29999, 
            pad_offset_x = 1., pad_offset_y = 1., 
            pad_tsize = 0.06, pad_lsize=0.06;
    TString pad_title_y = "d^{2}N/(p_{T}dydp_{T})";
    TString pad_title_x = "p_{T} [GeV/c]";
    Format_Pad(ll, rl, pad_min, pad_max, pad_title_x, pad_title_y, pad_offset_x, pad_offset_y, pad_tsize, pad_lsize, "");        

    for (int centr: CENTR)
    {
        if (hSpectraM[part][centr]) 
        {
            hSpectraM[part][centr]->SetMarkerColor(kBlue);
            hSpectraM[part][centr]->SetMarkerStyle(8); 
            hSpectraM[part][centr]->SetMarkerSize(1);
            hSpectraM[part][centr]->Draw("SAME");
        }
            
              
        if (hSpectraK[part][centr])
        {
            hSpectraK[part][centr]->SetMarkerColor(kRed);
            hSpectraK[part][centr]->SetMarkerStyle(21); 
            hSpectraK[part][centr]->SetMarkerSize(1);
            hSpectraK[part][centr]->Draw("SAME");
        }
    }

    TLegend *legend = new TLegend(0.5, 0.7, 0.95, 0.9); //1 column
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.05);
    legend->AddEntry(hSpectraM[0][0], "Malaev", "p");
    legend->AddEntry(hSpectraK[0][0], "Kireev", "p");
    legend->Draw();

    c2->cd();
    TPad pad2("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2.Draw();
    pad2.cd();

    pad_min = 0, pad_max = 2, 
    pad_title_y = "Kireev / Malaev";
    gStyle->SetOptFit(0000);
    Format_Pad(ll, rl, pad_min, pad_max, pad_title_x, pad_title_y, 0.5, 0.5, 0.1, 0.1, "");        

    for (int centr: CENTR)
    {
        hdiff[part][centr] = (TH1D *)hSpectraM[part][centr]->Clone("hdiff");
        hdiff[part][centr]->Divide(hSpectraK[part][centr]);
        hdiff[part][centr]->SetMarkerSize(1);
        hdiff[part][centr]->SetMarkerColor(centrColors[centr]);
        hdiff[part][centr]->SetMarkerStyle(34);
        hdiff[part][centr]->Draw("SAME");
    }

    c2->SaveAs("output/CompareSpectra.pdf");
}