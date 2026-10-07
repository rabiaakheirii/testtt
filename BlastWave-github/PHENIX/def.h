#ifndef __DEF_H_
#define __DEF_H_

#include "BlastWave.h"   
#include "input/FormatOfEverything.h"

#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include <TH1.h>
#include <TFile.h>
#include "TLegend.h"
#include "TMinuit.h"
#include "TVirtualFitter.h"

const int MAX_CENTR = 20;
const int MAX_PARTS = 6;
const int N_CENTR = 4;
const int N_PARTS = 6;
const int N_SIGMA = 2;
const int PARTS[] = {0, 1, 2, 3, 4, 5};
const int PARTS_POS[] = {0, 2, 4};
const int PARTS_NEG[] = {1, 3, 5};
const int PARTS_ALL[] = {0, 1, 2, 3, 4, 5};

const int CENTR[] = {0, 1, 2, 3}; //, 5, 6, 7, 8 ,9, 10, 11};
const int N_CENTR_SYST[5] = {12, 4, 5, 5, 4};
const int CENTR_SYST[5][20] = {
        {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
        {0, 1, 2, 3},
        {0, 1, 2, 3, 4},
        {0, 1, 2, 3, 4},
        {0, 1, 2, 3},
};

TGraphErrors *grSpectra[6][MAX_CENTR];
TGraph *contour[N_PARTS][N_CENTR][N_SIGMA];
TF1 *ifuncx[N_PARTS][N_CENTR], *ifuncxGlobal[MAX_PARTS][N_CENTR]; 
int Ncentr[5] = {12, 4, 5, 5, 4};

double constPar[N_PARTS][N_CENTR],
       Tpar[N_PARTS][N_CENTR], Tpar_err[N_PARTS][N_CENTR], 
       utPar[N_PARTS][N_CENTR], utPar_err[N_PARTS][N_CENTR];
double paramsGlobal[2][N_CENTR][5]; // [2] - charge, 5 - количество параметров 1) T 2) ut 3) const pi 4) const K 5) const p
double paramsGlobalAllParts[N_CENTR][5];

double Npart[5][MAX_CENTR] =
        { //systN, centr
                {109.1, 351.4, 299.0, 253.9, 215.3, 166.6, 114.2, 74.4 , 45.5 , 25.7 , 13.4 ,6.3 }, //AuAu
                {3.1, 4.35, 3.3, 2.7, 0},         //pAl
                {11.34, 21.84, 15.38, 9.51, 4.87},//HeAu
                {57.0, 154.8, 80.4, 34.9, 7.5},   //CuAu
                {330, 159.0, 61.6, -17.8},       //UU
        };
        
string systNames[5] = {"AuAu", "pAl", "HeAu", "CuAu", "UU"};
// TString systNamesT[5] = {"AuAu", "pAl", "HeAu", "CuAu", "UU"};
string particles[6] = {"pip", "pim", "kp", "km", "p", "ap"};
string partTitles[6] = { "#pi^{+}","#pi^{#minus}","K^{+}","K^{#minus}","p", "#bar{p}"};
double masses[6] = {0.13957061, 0.13957061, 0.493667, 0.493667, 0.938272, 0.938272};
Color_t centrColors[11] = {kRed, kBlue, kGreen + 2, kBlack, kMagenta+2, kBlue + 3, kCyan + 2, kRed + 3, kBlack, kBlack, kBlack};
Color_t partColors[6] = {kRed, kRed, kBlue, kBlue, kGreen + 2, kGreen + 2};
Color_t systColors[6] = {kBlack, kBlue, kGreen + 2, kRed + 2, kMagenta};
string centrTitles[5] = {"0-20%", "20-40%", "40-60%", "60-80%"};
string centrTitlesAuAu[MAX_CENTR] = {"MB", "0–5%", "5–10%", "10–15%", "15–20%", "20–30%", "30–40%", "40–50%", "50–60%", "60–70%", "70–80%", "80–92%"};
double TAuAu[MAX_CENTR] = {132, 107.8, 109.8, 113.3, 116.5, 123, 132, 142, 153, 163, 168, 179};


TString NpartStr[5][MAX_CENTR] =
        {
                {"3.1", "4.35", "3.3", "2.7"},
                {"11.34", "21.84", "15.38", "9.51", "4.87"},
                {"57.0", "154.8", "80.4", "34.9", "7.5"},
                {"330", "159.0", "61.6", "17.8"},
                { "109.1", "351.4", "325.2", "299.0", "253.9", "234.6", "215.3", "166.6", "114.2", "74.4 ", "45.5 ", "25.7 ", "19.5 ", "14.5 ", "13.4 ", "9.5 ", "6.3 ", "14.5 "}
        };
TString centrStr[5][MAX_CENTR] =
        {
                {"0-72%", "0-20%", "20-40%", "40-72%"},
                {"0-88%", "0-20%", "20-40%", "40-60%", "60-88%"},
                {"0-80%", "0-20%", "20-40%", "40-60%", "60-80%"},
                {"0-80%", "0-20%", "20-40%", "40-60%", "60-80%"},
                {"MB", "0-10", "5-10", "10-15", "10-20", "15-20", "20-30", "30-40", "40-50", "50-60", "60-70", "60-80", "60-92", "70-80", "70-92", "80-92", "60-92"}
        };


// ======= Fit paramteres ============
// по центральностям
double xmin[] = {0.4, 0.2, 0.12, 0.4, 0.2, 0.12};
double xmax[] = {1, 1., 1, 1, 1, 1};

double T[] = {0.132,  0.1078, 0.1098, 0.1133, 0.1165, 0.123, 0.132, 0.142, 0.153, 0.163, 0,168, 0.179};
double Tmin[] = {0.06, 0.06, 0.06, 0.06, 0.06, 0.06};
double Tmax[] = {0.22, 0.22, 0.22, 0.22, 0.22, 0.22};

double betaAuAu[] = {0.71, 0.773, 0.769, 0.763, 0.754, 0.738, 0.71, 0.67, 0.614, 0.555, 0.497, 0.399};
double beta_[MAX_CENTR] = {0.673, 0.769, 0.614, 0.497, 0.399};
double betamin[] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
double betamax[] = {0.9, 0.9, 0.8, 0.9, 0.9, 0.9};

// по частицам
double con[]    = {10, 10, 1, 1, 0.1, 0.1};  
double conmin[] = {0, 0, 0, 0, 0, 0};
double conmax[] = {500, 500, 50, 50, 5000, 5000};

//For GlobalFit
double conGlobal[]    = {100, 100, 120, 60, 0.01, 0.01};  
double conminGlobal[] = {0, 0, 0, 0, 0, 0};
double conmaxGlobal[] = {5000, 5000, 5000, 5000, 5000, 5000};

// Try to find constants by hand
double handT[] = {0.09, 0.095, 0.1, 0.1, 0.105, 0.105};
double handBeta[] = {0.75, 0.7, 0.6, 0.55, 0.5, 0.5};
double handConst[MAX_PARTS][MAX_CENTR] = 
    {
        {1000, 600, 500, 300, 200, 100},
        {1000, 600, 500, 300, 200, 100},
        {700, 600, 520, 420, 200, 60},
        {600, 500, 420, 320, 100, 30},
        {30000, 15000, 10000, 5000, 4000, 3000},
        {1000, 800, 650, 500, 200, 80},
    };


double GetMt( int part, double pT )
{
    return sqrt(pT * pT  + masses[part] * masses[part]) - masses[part];
}

void getGlobalParams( int part, int centr, double parResults[4] )
{
    int charge = part % 2;    
    parResults[0] = paramsGlobal[charge][centr][2 + part / 2];
    parResults[1] = paramsGlobal[charge][centr][0]; 
    parResults[2] = paramsGlobal[charge][centr][1]; 
    parResults[3] = masses[part];     
}

#endif /* __DEF_H_ */