#include "TTree.h"
#include "TFile.h"


void urqmd_Tmu( int start = 0, int end = 0 )
{
    int t;
    float e, T, muB;

    TFile *f = new TFile("tree_hyper.root");
    TTree *tree = (TTree *)f->Get("tree");

    tree->SetBranchAddress("t", &t);
    tree->SetBranchAddress("e", &e);
    tree->SetBranchAddress("T", &T);
    tree->SetBranchAddress("muB", &muB);

    int nentries = (int)tree->GetEntries();
    cout << nentries << endl;
    
    if (end == 0) end = nentries;

    TH2F *h = new TH2F("h", "T vs #mu_{B}", 100, 0, 0.3, 100, 0, 0.3);

    for (int i = start; i < end; i++)
    {
        tree->GetEntry(i);
        h->Fill(muB, T);
        // cout << i << " " << muB << " " << T << endl;
    }

    TCanvas *c = new TCanvas("c", "c", 29, 30, 1200, 900);
    c->Divide(1,1);
    gStyle->SetOptStat(0000);

    h->GetXaxis()->SetTitle("#mu_{B} [GeV]");
    h->GetYaxis()->SetTitle("T [GeV]");
    h->Draw("");
    
    char outname[50];
    sprintf(outname, "output/urqmd_TmuB_hyper_%i_%i.pdf", start, end);
    c->SaveAs(outname);

}