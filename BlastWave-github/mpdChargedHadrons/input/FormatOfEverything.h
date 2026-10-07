//
// Created by vlad on 28.09.2019.
//

#ifndef FORMATOFEVERYTHING_H
#define FORMATOFEVERYTHING_H

#include "TGraphErrors.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include <iostream>
#include <fstream>
#include <TLatex.h>
#include <TStyle.h>
#include <TMarker.h>
#include <TPolyLine.h>
#include <TASImage.h>

void Format_Graph(TGraph *gr, int mark_style, Float_t mark_size, Color_t mark_col, int line_style, Float_t line_wd, Color_t line_col, Float_t alpha);
void Format_Latex(TLatex *lat, int font, Float_t size, Float_t line_wd);
void Format_typeC(TPolyLine *pl, Float_t line_wd, Color_t line_fill_col, Float_t alpha);
void Format_Marker(TMarker *mk, Float_t mark_size, Color_t mark_col, Float_t alpha);
void Format_Canvas(TCanvas *c2, int divide_x, int divide_y, float space);
void Format_Pad(double_t left, double_t right, double_t min, double_t max, const char * title_x, const char * title_y, double_t offset_x, double_t offset_y, double_t Tsize, double_t Lsize, const char * title, int NdivisionsY, int NdivisionsX );


void Format_Graph(TGraph *gr, int mark_style=0, Float_t mark_size=0, Color_t mark_col=0, int line_style=0, Float_t line_wd=0, Color_t line_col=0, Float_t alpha=0) {
  gr->SetMarkerStyle(mark_style);
  gr->SetMarkerSize(mark_size);
  gr->SetMarkerColorAlpha(mark_col, alpha);
  gr->SetLineStyle(line_style);
  gr->SetLineWidth(line_wd);
  gr->SetLineColorAlpha(line_col, alpha);
}

void Format_Latex(TLatex *lat, int font=0, Float_t size=0, Float_t line_wd=0) {
  lat->SetTextFont(font);
  lat->SetTextSize(size);
  lat->SetLineWidth(line_wd);
}

void Format_typeC(TPolyLine *pl, Float_t line_wd, Color_t line_fill_col, Float_t alpha) {
  pl->SetLineColorAlpha(line_fill_col, alpha);
  pl->SetLineWidth(line_wd);
  pl->SetFillColorAlpha(line_fill_col, alpha);
}

void Format_Marker(TMarker *mk, Float_t mark_size=0, Color_t mark_col=0, Float_t alpha=0) {
  mk->SetMarkerSize(mark_size);
  mk->SetMarkerColorAlpha(mark_col, alpha);
}


void Format_Canvas(TCanvas *c2, int divide_x, int divide_y, float space) 
{
  c2->Divide(divide_x, divide_y,0,0);

  char name1[200];
  TPad *c2_1;
  gStyle->SetOptStat(0000);

  for (int centr = 0; centr < divide_x * divide_y; centr++) 
  {
      c2->cd(centr + 1);
      gStyle->SetOptStat(1);
      c2->Range(-0.8012048, -9.907216, 8.532129, -0.9670103);

      c2->SetBorderMode(1);
      c2->SetBorderSize(1);

      sprintf(name1,"c2_%d",centr + 1);
      c2_1= (TPad*) c2->GetListOfPrimitives()->FindObject(name1);
      c2_1->SetTickx();
      c2_1->SetTicky();
      c2_1->SetLogy();

      if ((centr % divide_x) == 0 )
          c2_1->SetLeftMargin(0.2);

      if (((centr+1) % divide_x) == 0)
          c2_1->SetRightMargin(0.1);

      if (centr < divide_x ){
          c2_1->SetTopMargin(0);
      }
      if (centr >= divide_x*(divide_y - 1) ){
          c2_1->SetBottomMargin(0.2);
      }

      if (space != 0){
          c2_1->SetLeftMargin(0.1354148);
          c2_1->SetRightMargin(0.02);
          c2_1->SetTopMargin(0.02);
          c2_1->SetBottomMargin(0.1594148);
      }
  }
}

void Format_tex(float x = 0.3, float y = 1.4, float tex_size = 0.05402776, TString title = "part , |y|<0.35, #sqrt{s_{NN}} = 200 GeV"){
    TLatex *tex;
    tex = new TLatex(x,y, title);
    tex->SetTextFont(2);
    tex->SetTextSize(tex_size);
    tex->SetLineWidth(2);
    tex->Draw();
}
void Format_Pad(double_t left, double_t right, double_t min, double_t max, const char *title_x, const char *title_y, double_t offset_x, double_t offset_y, double_t Tsize, double_t Lsize, const char *title, int NdivisionsY = 4, int NdivisionsX = 9 ) 
{
  TH1F *second = new TH1F("", "", 100, left, right);

  second->SetMinimum(min);
  second->SetMaximum(max);
  second->SetStats(0);
  second->SetTitle(title);

  second->GetXaxis()->SetTitle(title_x);
  second->GetXaxis()->SetLabelFont(42);
  second->GetXaxis()->SetTitleFont(42);
  second->GetXaxis()->SetLabelSize(Lsize);
  second->GetXaxis()->SetTickSize(0.03);
  second->GetXaxis()->SetTitleSize(Tsize);
  second->GetXaxis()->SetTitleOffset(offset_x);
  second->GetXaxis()->SetNdivisions(NdivisionsX);

  second->GetYaxis()->SetTitle(title_y);
  second->GetYaxis()->SetLabelFont(42);
  second->GetYaxis()->SetTitleFont(42);
  second->GetYaxis()->SetLabelSize(Lsize);
  second->GetYaxis()->SetTickSize(0.02);
  second->GetYaxis()->SetTitleSize(Tsize);
  second->GetYaxis()->SetTitleOffset(offset_y);
  second->GetYaxis()->SetNdivisions(NdivisionsY);

  second->Draw();
}

double PAD_MIN = 0.00011;
double PAD_MAX = 299989;

void FormatSpectraPad( double texScale = 1, TString pad_title_x = "m_{T} - m_{0} [GeV/c]" )
{
    cout<< "Format pad" << endl;
    double ll = 0.01, rl = 2.49, pad_min = PAD_MIN, pad_max = PAD_MAX, 
            pad_offset_x = 1., pad_offset_y = 1., 
            pad_tsize = 0.09 * texScale, pad_lsize=0.08 * texScale;
    TString pad_title_y = "d^{2}N/(p_{T}dydp_{T})";
  
    Format_Pad(ll, rl, pad_min, pad_max, pad_title_x, pad_title_y, pad_offset_x, pad_offset_y, pad_tsize, pad_lsize, "");        
}


#endif //WORKUU_FORMATOFEVERYTHING_H
