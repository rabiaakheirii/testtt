#ifndef __DEF_H_
#define __DEF_H_

#include "input/FormatOfEverything.h"

#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include <TH1.h>
#include <TFile.h>
#include "TLegend.h"
#include "BlastWave.h"
#include "TLine.h"


const int MAX_CENTR = 10;
const int MAX_PARTS = 6;
const int N_CENTR = 7;
const int N_PARTS = 6;
const int N_SIGMA = 3;
const int PARTS[] = {0, 1, 2, 3, 4, 5}; // pi^+, pi-, K+, K-, p, antiproton
const int PARTS_POS[] = {0, 2, 4};
const int PARTS_NEG[] = {1, 3, 5};
const int PARTS_ALL[] = {0, 1, 2, 3, 4, 5};
const int CENTR[] = {0, 1, 2, 3, 4, 5};
const int CENTR_MALAEV[] = {0, 1, 2, 3, 4, 5, 6};

TH1D *hSpectra[6][12];
TGraphErrors *grSpectra[6][12];

string particles[6] = {"pip", "pim", "kp", "km", "p", "ap"};
string partTitles[6] = { "#pi^{+}","#pi^{#minus}","K^{+}","K^{#minus}","p", "#bar{p}"};
double masses[6] = {0.13957061, 0.13957061, 0.493667, 0.493667, 0.938272, 0.938272};
Color_t centrColors[11] = {kRed, kBlue, kGreen + 2, kBlack, kMagenta, kBlue+2, kTeal -1, kBlack, kBlack, kBlack, kBlack};
Color_t partColors[6] = {kMagenta + 1, kMagenta + 1, kBlue, kBlue, kGreen + 2, kGreen + 2};
string centrTitlesMalaev[] = {"0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-90%"};
string centrTitles[10] = {"0-10%", "10-20%", "20-30%", "30-40%", "40-60%", "60-80%"};
double centrX[10] = {5, 15, 25, 35, 45, 70};
// string centrTitles[10] = {"0-20%", "20-40%", "40-60%", "60-80%", "40-50%", "50-60%", "60-70%", "70-80%"};


// =============== Для  BlastWave ======================
// Kireev
// double xmin[] = {0.5, 0.5, 0.12, 0.4, 0.2, 0.12};
// double xmax[] = {1., 1., 1, 1, 1, 1};

// // Malaev по pt
// double xmin[] = {0.4, 0.4, 0.5, 0.5, 0.5, 0.5};
// double xmax[] = {0.8, 0.8, 1, 1, 1.2, 1.2};

// Malaev по pt
double xmin[] = {0.2, 0.2, 0.2, 0.2, 0.2, 0.2};
double xmax[] = {0.8, 0.8, 1, 1, 1.2, 1.2};

// double xmin[] = {0.001, 0.001, 0.001, 0.001, 0.001, 0.001};
// double xmax[] = {2, 2, 2, 2, 2, 2};

TGraph *contour[MAX_PARTS][N_CENTR][N_SIGMA];
TF1 *ifuncx[MAX_PARTS][N_CENTR], *ifuncxGlobal[MAX_PARTS][N_CENTR];

double paramsGlobal[2][N_CENTR][5]; // [2] - charge, 5 - количество параметров 1) T 2) ut 3) const pi 4) const K 5) const p
// по частицам
double con[]    = {10, 10, 1, 1, 0.1, 0.001};  
double conmin[] = {0, 0, 0, 0, 0, 0};
double conmax[] = {10000, 10000, 50000, 500000, 500000, 500000};

//For GlobalFit
double conGlobal[]    = {100, 100, 120, 60, 0.01, 0.0001};  
double conminGlobal[] = {0, 0, 0, 0, 0, 0};
double conmaxGlobal[] = {5000, 5000, 5000, 5000, 5000, 50};

// Try to find constants by hand
// double handT[] = {0.09, 0.095, 0.1, 0.1, 0.105, 0.105};
// double handBeta[] = {0.75, 0.7, 0.6, 0.55, 0.5, 0.5};

// From AuAu CH_PHENIX
double handT[] = {0.11, 0.116, 123.0, 0.132, 0.142, 0.153,};
double handBeta[] = {0.75, 0.7, 0.6, 0.51, 0.45, 0.4, 0.4};

// Kireev
// double handConst[MAX_PARTS][MAX_CENTR] = 
//     {
//         {1000, 600, 500, 300, 200, 100},
//         {1000, 600, 500, 300, 200, 100},
//         {700, 600, 520, 420, 200, 60},
//         {600, 500, 420, 320, 100, 30},
//         {30000, 15000, 10000, 5000, 4000, 3000},
//         {1000, 800, 650, 500, 200, 80},
//     };

// Malaev
double handConst[MAX_PARTS][MAX_CENTR] = 
    {
        {250, 100, 53, 35, 20, 10, 1},
        {250, 100, 53, 35, 20, 10, 1},
        {350, 120, 53, 35, 20, 10, 1},
        {350, 120, 53, 35, 20, 10, 1},
        {30000, 6000, 1000, 500, 200, 100, 10},
        {70, 60, 52, 42, 20, 6, 1},
    };

void getGlobalParams( int part, int centr, double parResults[4] )
{
    int charge = part % 2; 

    parResults[0] = paramsGlobal[charge][centr][2 + part / 2];
    parResults[1] = paramsGlobal[charge][centr][0]; 
    parResults[2] = paramsGlobal[charge][centr][1]; 
    parResults[3] = masses[part];     
}

//=======================================================
//Kireev
void SetSpectra(string inputFileName = "postprocess_mpdpid10", string type = "pt")
{
    TFile *f = new TFile(("input/" + inputFileName + ".root").c_str());
    TDirectory *fd;

    for (int i = 0; i < 6; i++)
    {
        fd = (TDirectory*)f->Get(particles[i].c_str());
        fd->cd();
        for (int centr = 0; centr < N_CENTR; centr++)
        {
            string name = "h__pt_" + particles[i] +"_centrality" + to_string(centr) + "_mc_y-0.5_0.5";
            hSpectra[i][centr] = (TH1D *)fd->Get(name.c_str());    
            
            if (!hSpectra[i][centr]) continue;
            const int N_BINS = hSpectra[i][centr]->GetNbinsX();
            double mT[1000], pT[1000], sp[1000], sp_err[1000], xerr[1000];
            for (int bin = 1; bin < N_BINS; bin++)
            {
                sp[bin - 1] = hSpectra[i][centr]->GetBinContent(bin);
                sp_err[bin - 1] = hSpectra[i][centr]->GetBinError(bin);
                xerr[bin - 1] = 0.;
                pT[bin - 1] = hSpectra[i][centr]->GetBinCenter(bin);
                mT[bin - 1] = sqrt(pT[bin - 1] * pT[bin - 1] + masses[i] * masses[i]) - masses[i];
                // cout << i << "  " << pT[bin - 1] << "  " << sp_err[bin - 1] << endl;
            }
    
            if (type == "mt") 
                grSpectra[i][centr] = new TGraphErrors(N_BINS - 1, mT, sp, xerr, sp_err);
            else if (type == "pt")
                grSpectra[i][centr] = new TGraphErrors(hSpectra[i][centr]);

            grSpectra[i][centr]->SetLineColor(centrColors[centr]);
        }
    }
}

void SetSpectraMalaev(string type = "pt")
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
                // cout << i << " " << centr << "  " <<  N_CENTR_M - 1 << endl;
                string name = partMalaev[i] +"_rec_" + to_string(centrBins[centr]) + "_" + to_string(centrBins[centr + 1]);
                hSpectra[i][centr] = (TH1D *)f->Get(name.c_str());    
                
                if (!hSpectra[i][centr]) continue;
                const int N_BINS = hSpectra[i][centr]->GetNbinsX();
                double mT[1000], pT[1000], sp[1000], sp_err[1000], xerr[1000];
                for (int bin = 1; bin < N_BINS; bin++)
                {
                    sp[bin - 1] = hSpectra[i][centr]->GetBinContent(bin);

                    // Increase errors for countour chi2 plots
                    double IncreaseErr = 1;
                    // double IncreaseErr = i < 2 ? 8 : 3;
                    // if (i < 2 && centr == 0) IncreaseErr = 10;
                    // if (i > 3 && centr == 0) IncreaseErr = 5;
                    // if ((i == 2 || i == 3) && centr > 4) IncreaseErr = 5;
                    sp_err[bin - 1] = IncreaseErr * hSpectra[i][centr]->GetBinError(bin);

                    xerr[bin - 1] = 0.;
                    pT[bin - 1] = hSpectra[i][centr]->GetBinCenter(bin);
                    mT[bin - 1] = sqrt(pT[bin - 1] * pT[bin - 1] + masses[i] * masses[i]) - masses[i];
                    // cout << i << "  " << pT[bin - 1] << "  " << sp_err[bin - 1] << endl;
                }
        
                if (type == "mt") 
                    grSpectra[i][centr] = new TGraphErrors(N_BINS - 1, mT, sp, xerr, sp_err);
                else if (type == "pt")
                    grSpectra[i][centr] = new TGraphErrors(hSpectra[i][centr]);

                grSpectra[i][centr]->SetLineColor(centrColors[centr]);
            }
        }
        delete f;
    }
    
}

#endif /* __DEF_H_ */