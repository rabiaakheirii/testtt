 
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "def.h"

int ipar0[3] = {2, 0, 1};
int ipar2[3] = {3, 0, 1};
int ipar4[3] = {4, 0, 1};
 
bool isParamsFileExist = false;

void WriteParams( const char filename[30] = "output/txtParams/GlobalBWparams.txt" )
{
   cout << " WriteParams " << endl;
   ofstream txtFile;
   if (isParamsFileExist) 
   {
      txtFile.open(filename, ios::app);
      txtFile << endl;
   }
   else txtFile.open(filename);

   for (int centr: CENTR)
   {
      txtFile  << centr << "  " 
               << paramsGlobalAllParts[centr][0]  << "  " << paramsGlobalAllParts[centr][1] << "  "
               << paramsGlobalAllParts[centr][2] << "   " << paramsGlobalAllParts[centr][3] << "   " << paramsGlobalAllParts[centr][4] << endl;
   }
    
   txtFile.close();
   isParamsFileExist = true;
}

// Create the GlobalCHi2 structure
struct GlobalChi2 {
   GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,
                ROOT::Math::IMultiGenFunction & f2,
                ROOT::Math::IMultiGenFunction & f3,
                ROOT::Math::IMultiGenFunction & f4,
                ROOT::Math::IMultiGenFunction & f5,
                ROOT::Math::IMultiGenFunction & f6) :
      fChi2_1(&f1), fChi2_2(&f2), fChi2_3(&f3), fChi2_4(&f4), fChi2_5(&f5), fChi2_6(&f6) {}

   const  ROOT::Math::IMultiGenFunction * fChi2_1;
   const  ROOT::Math::IMultiGenFunction * fChi2_2;
   const  ROOT::Math::IMultiGenFunction * fChi2_3;
   const  ROOT::Math::IMultiGenFunction * fChi2_4;
   const  ROOT::Math::IMultiGenFunction * fChi2_5;
   const  ROOT::Math::IMultiGenFunction * fChi2_6;

   // par[0] = constant;
   // par[1] = Tf;
   // par[2] = beta;
   // par[3] = mass	pi
   // par[4] = mass	K
   // par[5] = mass	p
   double operator() (const double *par) const 
   {
      const int Nparams = 4, Nfunc = 6;
      double p[Nfunc][Nparams];
   
      for (int i = 0; i < Nfunc; i++)
      {
         p[i][0] = par[2 + i];
         p[i][1] = par[0]; 
         p[i][2] = par[1]; 
         p[i][3] = masses[i];
      }
 
      return (*fChi2_1)(p[0]) + (*fChi2_2)(p[1]) + (*fChi2_3)(p[2]) + (*fChi2_4)(p[3]) + (*fChi2_5)(p[4]) + (*fChi2_6)(p[5]);
   }
};
 
void GlobalFitCentr( int centr) 
{
   cout << " GlobalFitCentr " << endl;
   double xmin = 0.3, xmax = 1.2;

   // perform now global fit
   ROOT::Math::WrappedMultiTF1 wf0(*ifuncxGlobal[0][centr], 1);
   ROOT::Math::WrappedMultiTF1 wf1(*ifuncxGlobal[1][centr], 1);
   ROOT::Math::WrappedMultiTF1 wf2(*ifuncxGlobal[2][centr], 1);
   ROOT::Math::WrappedMultiTF1 wf3(*ifuncxGlobal[3][centr], 1);
   ROOT::Math::WrappedMultiTF1 wf4(*ifuncxGlobal[4][centr], 1);
   ROOT::Math::WrappedMultiTF1 wf5(*ifuncxGlobal[5][centr], 1);
 
   ROOT::Fit::DataOptions opt;
   ROOT::Fit::DataRange range;
   range.SetRange(xmin, xmax);

   ROOT::Fit::BinData data0(opt, range);
   ROOT::Fit::FillData(data0, grSpectra[0][centr]);
 
   ROOT::Fit::BinData data1(opt, range);
   ROOT::Fit::FillData(data1, grSpectra[1][centr]);

   ROOT::Fit::BinData data2(opt, range);
   ROOT::Fit::FillData(data2, grSpectra[2][centr]);
 
   ROOT::Fit::BinData data3(opt, range);
   ROOT::Fit::FillData(data3, grSpectra[3][centr]);

   ROOT::Fit::BinData data4(opt, range);
   ROOT::Fit::FillData(data4, grSpectra[4][centr]);

   ROOT::Fit::BinData data5(opt, range);
   ROOT::Fit::FillData(data5, grSpectra[5][centr]);
   
   ROOT::Fit::Chi2Function chi2_0(data0, wf0);
   ROOT::Fit::Chi2Function chi2_2(data2, wf2);
   ROOT::Fit::Chi2Function chi2_4(data4, wf4);
   ROOT::Fit::Chi2Function chi2_1(data1, wf1);
   ROOT::Fit::Chi2Function chi2_3(data3, wf3);
   ROOT::Fit::Chi2Function chi2_5(data4, wf5);
   GlobalChi2 globalChi2(chi2_0, chi2_1, chi2_2, chi2_3, chi2_4, chi2_5);
 
   ROOT::Fit::Fitter fitter;
 
   const int Npar = 8;
   double par0[Npar] = { 0.10, 0.7, con[0], con[1], con[2], con[3], con[4], con[5]};
 
   // create before the parameter settings in order to fix or set range on them
   fitter.Config().SetParamsSettings(Npar, par0);

   fitter.Config().ParSettings(0).SetLimits(0.08, 0.2);
   fitter.Config().ParSettings(1).SetLimits(0.5, 0.95);
   fitter.Config().ParSettings(2).SetLimits(conminGlobal[0], conmaxGlobal[0]);
   fitter.Config().ParSettings(3).SetLimits(conminGlobal[1], conmaxGlobal[1]);
   fitter.Config().ParSettings(4).SetLimits(conminGlobal[2], conmaxGlobal[2]);
   fitter.Config().ParSettings(5).SetLimits(conminGlobal[3], conmaxGlobal[3]);
   fitter.Config().ParSettings(6).SetLimits(conminGlobal[4], conmaxGlobal[4]);
   fitter.Config().ParSettings(7).SetLimits(conminGlobal[5], conmaxGlobal[5]);

   fitter.Config().MinimizerOptions().SetPrintLevel(0);
   fitter.Config().SetMinimizer("Minuit2","Migrad");
 
   // fit FCN function directly
   // (specify optionally data size and flag to indicate that is a chi2 fit)
   fitter.FitFCN(8, globalChi2, 0, data0.Size() + data1.Size() + data2.Size() + data3.Size() + data4.Size() + data5.Size(), true);
   ROOT::Fit::FitResult result = fitter.Result();
   result.Print(std::cout);

   const double *fitResults = result.GetParams();
   for (int i = 0; i < 8; i++ ) paramsGlobalAllParts[centr][i] = fitResults[i];

   cout<< "Result " << paramsGlobalAllParts[centr][0] << "  " << paramsGlobalAllParts[centr][1] << "  " << paramsGlobalAllParts[centr][2] << "  " << paramsGlobalAllParts[centr][3] << "  " << paramsGlobalAllParts[centr][4] << endl;
}

void DrawFitSpectra( string chargeFlag = "all" )
{
   TCanvas *c2 = new TCanvas("c2", "c2", 29, 30, 1100, 1200);
   Format_Canvas(c2, 2, 3, 0);

   int padN = 1;   
   for (int part: PARTS_ALL)
   {
      c2->cd(padN++);
      FormatSpectraPad(1);

      double shiftX = (part % 2 == 0) ? 0 : 0.1;
      double texScale = (part < 3) ? 1 : 0.9;

      TLegend *legend = new TLegend(0.55 - shiftX, 0.7, 0.98 - shiftX, 0.9); //1 column
      legend->SetNColumns(2);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      legend->SetTextSize(0.07 * texScale);

      TLatex *titleTex = new TLatex(0.6, 500, partTitles[part].c_str());
      titleTex->SetTextFont(42);
      titleTex->SetTextSize(0.08);
      titleTex->SetLineWidth(2 * texScale);

      for (int centr: CENTR)
      {
         double parResults[5];
         getGlobalParamsAllParts(part, centr, parResults);

         ifuncxGlobal[part][centr]->SetParameters(parResults);
         ifuncxGlobal[part][centr]->SetLineColor(centrColors[centr]);
         ifuncxGlobal[part][centr]->Draw("SAME");
         grSpectra[part][centr]->GetListOfFunctions()->Add(ifuncxGlobal[part][centr]);
         grSpectra[part][centr]->SetMarkerStyle(8);
         grSpectra[part][centr]->SetMarkerSize(1);
         grSpectra[part][centr]->Draw("P SAME");

         legend->AddEntry(ifuncxGlobal[part][centr], centrTitles[centr].c_str(), "l");        
      }

      legend->Draw();
      titleTex->Draw(); 
   }
 
   c2->SaveAs("output/BlastWaveGlobalFitAllParts.pdf");
   gROOT->ProcessLine(".q");
}

void BlastWaveGlobalAllParts(string chargeFlag = "all") 
{
   // ++++++ Read data ++++++++++++++++++++++++++++++++++++

    string inputFileName = "postprocess_mpdpid10";
    SetSpectra(inputFileName, "mt");

   // +++++++++ Fit +++++++++++++++++++++++++++++++++++++++

   for (int centr: CENTR)
   {
      TVirtualFitter::SetDefaultFitter("Minuit");  
      // TMinuit* minuit = new TMinuit(5); 
      TF1 *funcx;
      MyIntegFunc *integ;
      double xmin = 0.3, xmax = 1.2;

      funcx = new TF1("funcx", bwfitfunc, 0.01, 10, 5);
      funcx->SetParameters(2,1);
      funcx->SetParNames("constant", "T", "beta", "mass", "pt");
      integ = new MyIntegFunc(funcx);

      for (int part: PARTS_ALL)
      {
         string ifuncxName = "BW_" + to_string(part);
         ifuncxGlobal[part][centr] = new TF1("ifuncx", integ, xmin, xmax, 4, ifuncxName.c_str());
         ifuncxGlobal[part][centr]->FixParameter(3, masses[part]); 	    //	mass
         ifuncxGlobal[part][centr]->SetParameter(0, con[part]);	        //	constant
         ifuncxGlobal[part][centr]->SetParameter(1, 0.1);	            //	temp.
         ifuncxGlobal[part][centr]->SetParameter(2, 0.75);	            //	beta

         ifuncxGlobal[part][centr]->SetParLimits(0, conmin[part], conmax[part]);	//	constant.
         ifuncxGlobal[part][centr]->SetParLimits(1, 0.06, 0.2);          	        //	temp.
         ifuncxGlobal[part][centr]->SetParLimits(2, 0.5, 0.95);	                    //	beta
      }

      GlobalFitCentr(centr); // positive charged
   }

   WriteParams("output/txtParams/GlobalBWparamsAllParts.txt");

   DrawFitSpectra();
}
