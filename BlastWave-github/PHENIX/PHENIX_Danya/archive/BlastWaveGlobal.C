/* Фитируем глобально спектры */

// Подключение необходимых заголовочных файлов

#include "input/headers/def.h"
#include "input/headers/WriteReadFiles.h"
#include "input/headers/BlastWaveFit.h"

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"


// Глобальные переменные и параметры
// Индексы параметров для разных частиц
int ipar0[3] = {2, 0, 1};
int ipar2[3] = {3, 0, 1};
int ipar4[3] = {4, 0, 1};

// Флаг существования файла параметров
bool isParamsFileExist = false;

// Структура для расчета глобального хи-квадрат
struct GlobalChi2 
{
   GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,
                ROOT::Math::IMultiGenFunction & f2,
                ROOT::Math::IMultiGenFunction & f3) :
   fChi2_1(&f1), fChi2_2(&f2), fChi2_3(&f3) {}

   // Указатели на хи-квадрат функции для трех частиц
   const  ROOT::Math::IMultiGenFunction * fChi2_1;
   const  ROOT::Math::IMultiGenFunction * fChi2_2;
   const  ROOT::Math::IMultiGenFunction * fChi2_3;
   
   // Оператор расчета общего хи-квадрат
   double operator() (const double *par) const 
   {
      // par[0] = constant;
      // par[1] = Tf;
      // par[2] = beta;
      // par[3] = mass	pi
      // par[4] = mass	K
      // par[5] = mass	p

      const int Nparams = 4, Nfunc = 3;
      double p[Nfunc][Nparams];
      int commonParams = 3;

      for (int i = 0; i < Nfunc; i++)
      {
         p[i][0] = par[2 + i];
         p[i][1] = par[0]; 
         p[i][2] = par[1]; 
         p[i][3] = masses[2 * i];
      }
      
      // Сумма хи-квадрат для трех частиц
      return (*fChi2_1)(p[0]) + (*fChi2_2)(p[1]) + (*fChi2_3)(p[2]);
   }
};


/* ---------------------- Функция фитирования ---------------------- */


// Основная функция фитирования для определенной центральности
void GlobalFitCentr( int centr, int charge = 0 ) 
{
   cout << " ==================== GlobalFitCentr " << centr << " ==================== " << endl;

   double xmin, xmax;
   if (systN == 0) {
      xmin = 0.2;
      xmax = 2.0;
   } else {
      xmin = 0.3;
      xmax = 1.0;
   }

   // 1. Подготовка функций
   ROOT::Math::WrappedMultiTF1 wf0(*ifuncxGlobal[0 + charge][centr], 1);
   ROOT::Math::WrappedMultiTF1 wf2(*ifuncxGlobal[2 + charge][centr], 1);
   ROOT::Math::WrappedMultiTF1 wf4(*ifuncxGlobal[4 + charge][centr], 1);

   // 2. Настройка данных
   ROOT::Fit::DataOptions opt;
   ROOT::Fit::DataRange range0, range2, range4;
   
   // 3. Загрузка данных из графиков (TGraphErrors)
   range0.SetRange(xmin, xmax);
   ROOT::Fit::BinData data0(opt, range0);
   ROOT::Fit::FillData(data0, grSpectra[0 + charge][centr]);
 
   range2.SetRange(xmin, xmax);
   ROOT::Fit::BinData data2(opt, range2);
   ROOT::Fit::FillData(data2, grSpectra[2 + charge][centr]);
 
   range4.SetRange(xmin, xmax);
   ROOT::Fit::BinData data4(opt, range4);
   ROOT::Fit::FillData(data4, grSpectra[4 + charge][centr]);
   
   // 4. Создание хи-квадрат функций
   ROOT::Fit::Chi2Function chi2_0(data0, wf0);
   ROOT::Fit::Chi2Function chi2_2(data2, wf2);
   ROOT::Fit::Chi2Function chi2_4(data4, wf4);

   // 5. Инициализация глобального хи-квадрат
   GlobalChi2 globalChi2(chi2_0, chi2_2, chi2_4);

   // 6. Настройка фиттера
   ROOT::Fit::Fitter fitter;
   const int Npar = 5;
   double par0[5] = {handT[centr], handBeta[centr], 
                     handConst[0 + charge][centr], 
                     handConst[2 + charge][centr], 
                     handConst[4 + charge][centr]};
 
   // create before the parameter settings in order to fix or set range on them
   fitter.Config().SetParamsSettings(Npar, par0); 

   // 7. Установка ограничений на параметры
   if (centr < 10) {
      fitter.Config().ParSettings(0).SetLimits(0.08, 0.18);
      fitter.Config().ParSettings(1).SetLimits(0.30, 0.80);
      fitter.Config().ParSettings(2).SetLimits(handConst[0 + charge][centr] * 0, handConst[0 + charge][centr] * 3);
      fitter.Config().ParSettings(3).SetLimits(handConst[2 + charge][centr] * 0, handConst[2 + charge][centr] * 3);
      fitter.Config().ParSettings(4).SetLimits(handConst[4 + charge][centr] * 0, handConst[4 + charge][centr] * 3);
   }
   else if (centr < 11) {
      fitter.Config().ParSettings(0).SetLimits(0.165, 0.20);
      fitter.Config().ParSettings(1).SetLimits(0.30, 0.55);
      fitter.Config().ParSettings(2).SetLimits(handConst[0 + charge][centr] * 0, handConst[0 + charge][centr] * 0.0002);
      fitter.Config().ParSettings(3).SetLimits(handConst[2 + charge][centr] * 0, handConst[2 + charge][centr] * 0.1);
      fitter.Config().ParSettings(4).SetLimits(handConst[4 + charge][centr] * 0, handConst[4 + charge][centr] * 0.0003);
   }
   else if (centr < 12) {
      fitter.Config().ParSettings(0).SetLimits(0.165, 0.20);
      fitter.Config().ParSettings(1).SetLimits(0.30, 0.41);
      fitter.Config().ParSettings(2).SetLimits(handConst[0 + charge][centr] * 0, handConst[0 + charge][centr] * 0.0001);
      fitter.Config().ParSettings(3).SetLimits(handConst[2 + charge][centr] * 0, handConst[2 + charge][centr] * 0.1);
      fitter.Config().ParSettings(4).SetLimits(handConst[4 + charge][centr] * 0.00005, handConst[4 + charge][centr] * 0.00009);
   }
      
   // 8. Выполнение фита
   fitter.Config().MinimizerOptions().SetPrintLevel(0);
   
   // Первый проход
   // fitter.Config().ParSettings(0).Fix();
   // fitter.Config().ParSettings(1).Fix(); // Фиксируем beta
   // // fitter.Config().ParSettings(2).Fix();
   // // fitter.Config().ParSettings(3).Fix();
   // // fitter.Config().ParSettings(4).Fix();
   // // fitter.Config().SetMinimizer("Minuit2", "Samplex");
   // // fitter.Config().SetMinimizer("GSLSimAn");
   // fitter.Config().SetMinimizer("Genetic");  // Глобальный поиск
   // fitter.Config().SetMinimizer("Minuit2", "Samplex");  // Точная локальная минимизация
   // fitter.FitFCN(5, globalChi2, 0, data0.Size() + data2.Size() + data4.Size(), true);

   // fitter.Config().ParSettings(0).Release();
   // fitter.Config().ParSettings(1).Release();
   // fitter.Config().ParSettings(2).Release();
   // fitter.Config().ParSettings(3).Release();
   // fitter.Config().ParSettings(4).Release();
   // fitter.Config().SetMinimizer("Minuit2", "Samplex");
   // fitter.FitFCN(5, globalChi2, 0, data0.Size() + data2.Size() + data4.Size(), true);

   // // Второй проход (разблокируем параметры)
   fitter.Config().ParSettings(0).Fix();
   fitter.Config().ParSettings(1).Fix();
   fitter.Config().SetMinimizer("Minuit2", "Migrad");
   fitter.FitFCN(5, globalChi2, 0, data0.Size() + data2.Size() + data4.Size(), true);

   fitter.Config().ParSettings(0).Release();
   fitter.Config().ParSettings(1).Release();
   fitter.Config().SetMinimizer("GSLSimAn"); // Глобальный поиск
   fitter.Config().SetMinimizer("Genetic");  // Глобальный поиск
   fitter.Config().SetMinimizer("Minuit2", "Migrad");  // Точная локальная минимизация
   fitter.FitFCN(5, globalChi2, 0, data0.Size() + data2.Size() + data4.Size(), true);

   ROOT::Fit::FitResult result = fitter.Result();
   result.Print(std::cout);

   double chi2 = result.MinFcnValue();
   int total_points = data0.Size() + data2.Size() + data4.Size();
   int n_free_params = result.NFreeParameters();
   int ndf = total_points - n_free_params;
   double chi2_ndf = chi2 / ndf;
   
   cout << "Chi2/NDF = " << chi2_ndf 
        << " (Chi2 = " << chi2 
        << ", NDF = " << ndf << ")" << endl;

   // 9. Сохранение результатов
   const double *fitResults = result.GetParams();
   for (int i = 0; i < 5; i++ ) paramsGlobal[charge][centr][i] = fitResults[i];

   string chargeFlag = (charge == 0) ? "pos" : "neg";
   cout << "Result " << paramsGlobal[charge][centr][0] << "  " 
                     << paramsGlobal[charge][centr][1] << "  " 
                     << paramsGlobal[charge][centr][2] << "  " 
                     << paramsGlobal[charge][centr][3] << "  " 
                     << paramsGlobal[charge][centr][4] << endl;
}


/* ---------------------- Построение спектров частиц ---------------------- */


// Функция визуализации результатов
void DrawFitSpectra( int systN, string chargeFlag = "all" )
{
   TCanvas *c2 = new TCanvas("c2", "c2", 30, 30, 1440, 2160);
   Format_Canvas(c2, 2, 3, 0);

   // Цикл по частицам и центральностям
   int padN = 1;   
   for (int part: PARTS_ALL)
   {
      c2->cd(padN++);
      FormatSpectraPad(1);

      if (chargeFlag == "pos" && part % 2 == 1) continue;
      if (chargeFlag == "neg" && part % 2 == 0) continue;

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

      for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
         int centr = CENTR_SYST[systN][j];
         double parResults[5];
         getGlobalParams(part, centr, parResults);

         ifuncxGlobal[part][centr]->SetParameters(parResults);
         ifuncxGlobal[part][centr]->SetLineColor(centrColors[centr]);
         ifuncxGlobal[part][centr]->Draw("SAME");

         grSpectra[part][centr]->GetListOfFunctions()->Add(ifuncxGlobal[part][centr]);
         grSpectra[part][centr]->SetMarkerStyle(8);
         grSpectra[part][centr]->SetMarkerSize(1);
         grSpectra[part][centr]->Draw("P SAME");

         legend->AddEntry(ifuncxGlobal[part][centr], centrTitlesAuAu[centr].c_str(), "l");        
      }

      legend->Draw();
      titleTex->Draw(); 
   }
 
   c2->SaveAs("output/pics/BlastWaveGlobal_" + systNamesT[systN] + ".png");
   gROOT->ProcessLine(".q");
}


/* ---------------------- Главная функция ---------------------- */


void BlastWaveGlobal(string chargeFlag = "all") 
{
   // Чтение данных
   if (systN == 0) ReadFromFileAuAu();                    // Для системы AuAu
   else for (int part: PARTS) ReadFromFile(part, systN);  // Для других систем 

   // +++++++++ Fit +++++++++++++++++++++++++++++++++++++++

   for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
      int centr = CENTR_SYST[systN][j];

      TVirtualFitter::SetDefaultFitter("Minuit");  

      TF1 *funcx;
      MyIntegFunc *integ;

      double xmin, xmax;
      if (systN == 0) {
         xmin = 0.2;
         xmax = 2.0;
      } else {
         xmin = 0.3;
         xmax = 1.2;
      }
      
      funcx = new TF1("funcx", bwfitfunc, 0.01, 10, 5);
      funcx->SetParameters(2, 1);
      funcx->SetParNames("constant", "T", "beta", "mass", "pt");
      integ = new MyIntegFunc(funcx);

      for (int part: PARTS_ALL)
      {
         string ifuncxName = "BW_" + to_string(part);
         ifuncxGlobal[part][centr] = new TF1("ifuncx", integ, xmin, xmax, 4, ifuncxName.c_str());
         double handParams[4] = {handConst[part][centr], handT[centr], handBeta[centr], masses[part]};

         ifuncxGlobal[part][centr]->SetParameters(handParams);
         ifuncxGlobal[part][centr]->SetParLimits(0, handConst[part][centr], handConst[part][centr]); 
         ifuncxGlobal[part][centr]->SetParameter(1, systN == 0 ? TAuAu[centr] : handT[centr]);	    // temp.
         ifuncxGlobal[part][centr]->SetParameter(1, 0.118);	       // temp.
         ifuncxGlobal[part][centr]->SetParameter(2, systN == 0 ? betaAuAu[centr] : beta_[centr]);	 // beta
         ifuncxGlobal[part][centr]->SetParLimits(2, 0.3, 0.88);	 // beta
         ifuncxGlobal[part][centr]->FixParameter(3, masses[part]); // mass

         // // double xmin = GetMt(part, 0.5), xmax = GetMt(part, 1.1);
         // string ifuncxName = "BW_" + to_string(part);
         // ifuncxGlobal[part][centr] = new TF1("ifuncx", integ, xmin, xmax, 4, ifuncxName.c_str());
         // ifuncxGlobal[part][centr]->FixParameter(3, masses[part]); 	    //	mass
         // ifuncxGlobal[part][centr]->SetParameter(0, con[part]);	        //	constant
         // ifuncxGlobal[part][centr]->SetParameter(1, 0.118);	            //	temp.
         // ifuncxGlobal[part][centr]->SetParameter(2, systN == 0 ? betaAuAu[centr] : beta_[centr]);	 // beta

         // ifuncxGlobal[part][centr]->SetParLimits(0, conmin[part], conmax[part]);	//	constants
         // ifuncxGlobal[part][centr]->SetParLimits(1, 0.06, 0.2);          	        //	temp.
         // ifuncxGlobal[part][centr]->SetParLimits(2, 0.1, 0.99);	                    //	beta
      }

      if (chargeFlag != "neg") GlobalFitCentr(centr, 0); // positive charged
      if (chargeFlag != "pos") GlobalFitCentr(centr, 1); // negative charged
   }

   if (chargeFlag != "neg") WriteGlobalParams(&isParamsFileExist, 0, systN, "output/parameters/GlobalBWparams_" + systNamesT[systN] + ".txt");
   if (chargeFlag != "pos") WriteGlobalParams(&isParamsFileExist, 1, systN, "output/parameters/GlobalBWparams_" + systNamesT[systN] + ".txt");

   DrawFitSpectra(systN, chargeFlag);
} 
