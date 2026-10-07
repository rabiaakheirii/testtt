 
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "WriteReadFiles.h"
#include "def.h"

int ipar0[3] = {2, 0, 1};
int ipar2[3] = {3, 0, 1};
int ipar4[3] = {4, 0, 1};
 
bool isParamsFileExist = false;
struct GlobalChi2 {
   GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,
                ROOT::Math::IMultiGenFunction & f2,
                ROOT::Math::IMultiGenFunction & f3) :
      fChi2_1(&f1), fChi2_2(&f2), fChi2_3(&f3) {}

   const  ROOT::Math::IMultiGenFunction * fChi2_1;
   const  ROOT::Math::IMultiGenFunction * fChi2_2;
   const  ROOT::Math::IMultiGenFunction * fChi2_3;

   // par[0] = constant;
   // par[1] = Tf;
   // par[2] = beta;
   // par[3] = mass	pi
   // par[4] = mass	K
   // par[5] = mass	p
   double operator() (const double *par) const 
   {
      const int Nparams = 4, Nfunc = 3;
      double p[Nfunc][Nparams];

      for (int i = 0; i < Nfunc; i++)
       {
         p[i][0] = par[2 + i];
         p[i][1] = par[0]; // T
         p[i][2] = par[1]; // beta
         p[i][3] = masses[2 * i];
      }

      double chi2 = (*fChi2_1)(p[0]) + (*fChi2_2)(p[1]) + (*fChi2_3)(p[2]);

      double beta = par[1];
      double T = par[0];
      double T_pred = 0.129557 + beta * 0.232774 - beta * beta * 0.339709;
      double penalty = (T - T_pred) * (T - T_pred) / (0.01 * 0.01);

      return chi2 + penalty;
   }
};
 
void GlobalFitCentr( int centr, int charge = 0 ) 
{
   cout << " GlobalFitCentr " << endl;
   double xmin = 0.3, xmax = 1.;

   ROOT::Math::WrappedMultiTF1 wf0(*ifuncxGlobal[0 + charge][centr], 1);
   ROOT::Math::WrappedMultiTF1 wf2(*ifuncxGlobal[2 + charge][centr], 1);
   ROOT::Math::WrappedMultiTF1 wf4(*ifuncxGlobal[4 + charge][centr], 1);

   ROOT::Fit::DataOptions opt;
   ROOT::Fit::DataRange range0, range2, range4;

   xmin = 0.05, xmax = 0.8;
   range0.SetRange(xmin, xmax);
   ROOT::Fit::BinData data0(opt, range0);
   ROOT::Fit::FillData(data0, grSpectra[0 + charge][centr]);

   xmin = 0.12, xmax = 0.8;
   range2.SetRange(xmin, xmax);
   ROOT::Fit::BinData data2(opt, range2);
   ROOT::Fit::FillData(data2, grSpectra[2 + charge][centr]);

   xmin = 0.2, xmax = 1.;
   range4.SetRange(xmin, xmax);
   ROOT::Fit::BinData data4(opt, range4);
   ROOT::Fit::FillData(data4, grSpectra[4 + charge][centr]);

   ROOT::Fit::Chi2Function chi2_0(data0, wf0);
   ROOT::Fit::Chi2Function chi2_2(data2, wf2);
   ROOT::Fit::Chi2Function chi2_4(data4, wf4);

   GlobalChi2 globalChi2(chi2_0, chi2_2, chi2_4);

   ROOT::Fit::Fitter fitter;

   double parResults[6][4];
   for (int part : {0 + charge, 2 + charge, 4 + charge})
   {
      ReadParams(2, centr, parResults[part]);
   }

   const int Npar = 5;
   // double par0[5] = 
   // {
   //    GetT(handBeta[centr]), handBeta[centr],
   //    parResults[0 + charge][0] * 8,
   //    parResults[2 + charge][0] * 8,
   //    parResults[4 + charge][0] * 10
   // };
   double par0[5];
   par0[0] = GetT(handBeta[centr]);
   par0[1] = handBeta[centr]; 
   par0[2] = handConst[0 + charge][centr];
   par0[3] = handConst[2 + charge][centr];
   par0[4] = handConst[4 + charge][centr];

   cout << "!!!!!!!!!!!!  " << par0[0] << "  " << par0[1] << "  " << par0[2] << "  " << par0[3] << "  "<< par0[4] << endl;
   fitter.Config().SetParamsSettings(Npar, par0);
   fitter.Config().ParSettings(0).SetLimits(par0[0] * 0.9, par0[0] * 1.1);
   fitter.Config().ParSettings(1).SetLimits(par0[1] * 0.9, par0[1] * 1.1);
   fitter.Config().ParSettings(2).SetLimits(par0[2] * 0.001, par0[2] * 100.9);
   fitter.Config().ParSettings(3).SetLimits(par0[3] * 0.001, par0[3] * 100.9);
   fitter.Config().ParSettings(4).SetLimits(par0[4] * 0.001, par0[4] * 100.9);
   
   fitter.Config().MinimizerOptions().SetPrintLevel(0);
   fitter.Config().SetMinimizer("Minuit2", "Migrad");

   fitter.FitFCN(5, globalChi2, 0, data0.Size() + data2.Size() + data4.Size(), true);

   ROOT::Fit::FitResult result = fitter.Result();
   result.Print(std::cout);

   if(result.Status() != 0 || result.Chi2() < 1e-6) {
       cout << "BAD GLOBAL FIT: centr " << centr 
            << " charge " << charge
            << " status: " << result.Status() 
            << " chi2: " << result.Chi2() << endl;
       return;
   }

   double chi2 = result.Chi2();
   int ndf = result.Ndf();
   double chi2_ndf = (ndf > 0) ? chi2 / ndf : -1;
   cout << "Global fit results for centr " << centr << " charge " << charge << ": "
        << "Chi2/NDF = " << chi2_ndf 
        << " (Chi2 = " << chi2 
        << ", NDF = " << ndf << ")" << endl;

   const double *fitResults = result.GetParams();
   for (int i = 0; i < 5; i++ ) paramsGlobal[charge][centr][i] = fitResults[i];

   string chargeFlag = (charge == 0) ? "pos" : "neg";
   cout<< "Result " << paramsGlobal[charge][centr][0] << "  " << paramsGlobal[charge][centr][1] << "  " << paramsGlobal[charge][centr][2] << "  " << paramsGlobal[charge][centr][3] << "  " << paramsGlobal[charge][centr][4] << endl;
}

void DrawFitSpectra( string chargeFlag = "all" )
{
   TCanvas *c2 = new TCanvas("c2", "c2", 29, 30, 2000, 1200);
   Format_Canvas(c2, 3, 2, 0);

   int padN = 1;   
   double shiftXarr[6] = {0, 0.1, 0.18, 0, 0.1, 0.18};
   for (int part: {0, 2, 4, 1, 3, 5})
   {
      double shiftX = shiftXarr[padN - 1];
      double texScale = (padN <= 3) ? 1 : 0.9;
      
      c2->cd(padN++);
      FormatSpectraPad(1);

      if (chargeFlag == "pos" && part % 2 == 1) continue;
      if (chargeFlag == "neg" && part % 2 == 0) continue;



      TLegend *legend = new TLegend(0.45 - shiftX, 0.7, 0.98 - shiftX, 0.9); //1 column
      legend->SetNColumns(2);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      legend->SetTextSize(0.07 * texScale);

      TLatex *titleTex = new TLatex(0.4, 500, partTitles[part].c_str());
      titleTex->SetTextFont(42);
      titleTex->SetTextSize(0.08);
      titleTex->SetLineWidth(2 * texScale);

      for (int centr: CENTR_MALAEV)
      {
         double parResults[5];
         getGlobalParams(part, centr, parResults);

         ifuncxGlobal[part][centr]->SetParameters(parResults);
         ifuncxGlobal[part][centr]->SetLineColor(centrColors[centr]);
         ifuncxGlobal[part][centr]->Draw("SAME");
         grSpectra[part][centr]->GetListOfFunctions()->Add(ifuncxGlobal[part][centr]);
         grSpectra[part][centr]->SetMarkerStyle(8);
         grSpectra[part][centr]->SetMarkerSize(0.5);
         grSpectra[part][centr]->Draw("P SAME");

         legend->AddEntry(ifuncxGlobal[part][centr], centrTitlesMalaev[centr].c_str(), "l");        
      }

      legend->Draw();
      titleTex->Draw(); 
   }
 
   c2->SaveAs("output/BlastWaveGlobalFit_Malaev.pdf");
   gROOT->ProcessLine(".q");
}

void BlastWaveGlobal(string chargeFlag = "all") 
{
   // ++++++ Read data ++++++++++++++++++++++++++++++++++++

    // string inputFileName = "postprocess_mpdpid10";
    // SetSpectra(inputFileName, "mt");
      SetSpectraMalaev("mt");
   // +++++++++ Fit +++++++++++++++++++++++++++++++++++++++

   double parResults[6][4];
   for (int centr: CENTR_MALAEV)
   {
      TVirtualFitter::SetDefaultFitter("Minuit");  
      // TMinuit* minuit = new TMinuit(5); 
      TF1 *funcx;
      MyIntegFunc *integ;
      double xmin = 0.2, xmax = 0.8;

      funcx = new TF1("funcx", bwfitfunc, 0.01, 10, 5);
      funcx->SetParameters(2,1);
      funcx->SetParNames("constant", "T", "beta", "mass", "pt");
      integ = new MyIntegFunc(funcx);

      for (int part: PARTS_ALL)
      {
         string ifuncxName = "BW_" + to_string(part);
         ifuncxGlobal[part][centr] = new TF1("ifuncx", integ, xmin, xmax, 4, ifuncxName.c_str());
      }

      if (chargeFlag != "neg") GlobalFitCentr(centr, 0); // positive charged
      if (chargeFlag != "pos") GlobalFitCentr(centr, 1); // negative charged
   }

   if (chargeFlag != "neg") WriteGlobalParams(&isParamsFileExist, 0);
   if (chargeFlag != "pos") WriteGlobalParams(&isParamsFileExist, 1);

   DrawFitSpectra();
}