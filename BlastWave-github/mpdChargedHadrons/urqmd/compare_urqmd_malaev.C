#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

using namespace std;

// ================= CONFIG =================
const double yMin = -0.5;
const double yMax =  0.5;
const double dy   = yMax - yMin;

const int nPtBins = 50;
const double ptMin = 0.0;
const double ptMax = 2.0;

// центральности (%)
vector<pair<int,int>> centBins = {
    {0,10},{10,20},{20,30},{30,40},{40,50},{50,60},{60,90}
};

// ==========================================

int GetParticleType(double m, int chg)
{
    // массы (GeV)
    const double mpi = 0.13957;
    const double mK  = 0.493677;
    const double mp  = 0.938272;

    // допуски (важно! UrQMD не даёт идеальные значения)
    const double dPi = 0.03;
    const double dK  = 0.05;
    const double dP  = 0.08;

    // нейтральные сразу убираем
    if (chg == 0) return -1;

    if (fabs(m - mpi) < dPi) return 1; // pion
    if (fabs(m - mK ) < dK ) return 2; // kaon
    if (fabs(m - mp ) < dP ) return 0; // proton

    return -1;
}

// charge: +1 / -1
int GetCharge(int chg) {
    return chg;
}

// центральность через b
int GetCentralityBin(double b, double bmax) {
    double c = (b*b)/(bmax*bmax) * 100.0;

    for (int i=0;i<centBins.size();i++) {
        if (c >= centBins[i].first && c < centBins[i].second)
            return i;
    }
    return -1;
}

// ==========================================

void compare_urqmd_malaev(
    const char* f14file = "/home/dasha/urqmd-3.4/output_XeW/urqmdXeW.f14",
    const char* malaev_neg = "../input/Malaev_Spectra_negative_XeW_2_5_GeV.root",
    const char* malaev_pos = "../input/Malaev_Spectra_positive_XeW_2_5_GeV.root")
{
    gStyle->SetOptStat(0);

    // ===== histograms =====
    TH1D* h[3][7][2]; // particle, centrality, charge

    const char* partName[3] = {"Proton","Pion","Kaon"};
    const char* chName[2]   = {"mn","pl"};

    for (int p=0;p<3;p++) {
        for (int c=0;c<centBins.size();c++) {
            for (int ch=0;ch<2;ch++) {
                h[p][c][ch] = new TH1D(
                    Form("h_%s_%d_%s",partName[p],c,chName[ch]),
                    "",
                    nPtBins, ptMin, ptMax
                );
            }
        }
    }

    // ===== read f14 =====
    ifstream in(f14file);
    string line;

    double b = 0;
    double bmax = 15.0;

    int Nevt[7] = {0};

    while (getline(in,line)) {

        // impact parameter
        if (line.find("impact_parameter_real") != string::npos) {
            sscanf(line.c_str(),
                   "impact_parameter_real/min/max(fm): %lf %*lf %lf",
                   &b, &bmax);
        }

        // начало события
        if (line.find("event#") != string::npos) {

            int cent = GetCentralityBin(b,bmax);
            if (cent<0) continue;

            Nevt[cent]++;

            // ищем pvec
            while (getline(in,line)) {
                if (line.find("pvec") != string::npos) break;
            }

            // пропускаем 2 строки
            getline(in,line);
            getline(in,line);

            // читаем частицы
            while (getline(in,line)) {

                if (line.empty()) break;
                if (!isdigit(line[0]) && line[0] != ' ') break;

                double r0,rx,ry,rz;
                double E,px,py,pz,m;
                int ityp,i3,chg;

                sscanf(line.c_str(),
                    "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d",
                    &r0,&rx,&ry,&rz,
                    &E,&px,&py,&pz,&m,
                    &ityp,&i3,&chg);

                int ptype = GetParticleType(m, chg);
                if (ptype < 0) continue;

                int charge = (chg<0) ? 0 : 1;

                double pt = sqrt(px*px + py*py);
                double p  = sqrt(px*px + py*py + pz*pz);
                double y  = 0.5 * log((E+pz)/(E-pz));

                if (y < yMin || y > yMax) continue;

                h[ptype][cent][charge]->Fill(pt);
            }
        }
    }

    // ===== нормировка =====
    for (int p=0;p<3;p++) {
        for (int c=0;c<centBins.size();c++) {
            for (int ch=0;ch<2;ch++) {

                if (Nevt[c]==0) continue;

                double scale = 1.0 / (Nevt[c] * dy);

                for (int b=1;b<=nPtBins;b++) {
                    double content = h[p][c][ch]->GetBinContent(b);
                    double width   = h[p][c][ch]->GetBinWidth(b);
                    h[p][c][ch]->SetBinContent(b, content/(width)*scale);
                }
            }
        }
    }

    // ===== open Malaev =====
    TFile *f_neg = new TFile(malaev_neg);
    TFile *f_pos = new TFile(malaev_pos);

    int colors[7] = {kRed,kBlue,kGreen+2,kMagenta,kCyan+2,kOrange,kBlack};

    // ===== plotting =====
    for (int p=0;p<3;p++) {

        for (int ch=0;ch<2;ch++) {

            TCanvas *c1 = new TCanvas(
                Form("c_%s_%s",partName[p],chName[ch]),
                "",800,600);

            TLegend *leg = new TLegend(0.6,0.6,0.88,0.88);

            bool first = true;

            for (int c=0;c<centBins.size();c++) {

                h[p][c][ch]->SetLineColor(colors[c]);
                h[p][c][ch]->SetLineWidth(1);
                h[p][c][ch]->SetMarkerColor(colors[c]);
                h[p][c][ch]->SetMarkerStyle(20);

                if (first) {
                    h[p][c][ch]->SetTitle(
                        Form("%s (%s);p_{T} (GeV/c);1/N d^{2}N/dp_{T}dy",
                             partName[p],chName[ch]));
                    h[p][c][ch]->Draw("lE");
                    first = false;
                } else {
                    h[p][c][ch]->Draw("lE SAME");
                }

                leg->AddEntry(h[p][c][ch],
                    Form("UrQMD %d-%d%%",
                         centBins[c].first,centBins[c].second),
                    "lp");

                // ===== Malaev =====
                TString hname = Form("%s_%s_rec_%d_%d",
                    partName[p], chName[ch],
                    centBins[c].first,
                    centBins[c].second);

                cout << hname << endl;
                TH1D* hM = (TH1D*)(
                    (ch==0 ? f_neg : f_pos)->Get(hname)
                );

                if (hM) {
                    hM->SetLineStyle(2);
                    hM->SetLineColor(colors[c]);
                    hM->Draw("HIST SAME");

                    cout << "found "<< endl;
                    leg->AddEntry(hM,
                        Form("Malaev %d-%d%%",
                             centBins[c].first,
                             centBins[c].second),
                        "l");
                }
            }

            leg->Draw();
            c1->SetLogy();

            c1->SaveAs(
                Form("%s_%s_comparison.png",
                     partName[p], chName[ch]));
        }
    }

    cout << "DONE" << endl;
}