#include "def.h"
#include "WriteReadFiles.h"


using namespace std;


class BlastWaveFit {
public:
    int currentSystematicType = 0; 

    bool isContour = true;
    bool isDraw = true;
    
    double outParams[N_PARTS][N_CENTR][4];
    double outParamsErr[N_PARTS][N_CENTR][4];
    double paramsSystematics[N_PARTS][N_CENTR][4];
    double lLimitMult = 0.5, rLimitMult = 1.5; // for parLimits in case 4 (Systematic)
    double lLimitMultPi = 0.5, rLimitMultPi = 1.; // for parLimits in case 4 (Systematic Pi meson)
    

    void Fit( int initParamsType = 0 )
    {    
        // ++++++ Read data +++++++++++++++++++++++++++++++++++++

        // Чтение данных в зависимости от системы
        if (systN == 0) 
            ReadFromFileAuAu(); // Для системы AuAu
        else                    // Для других систем
            for (int part: PARTS) ReadFromFile(part, systN);

        // +++++++++ Fit +++++++++++++++++++++++++++++++++++++++

        gMinuit = new TMinuit(5);  // Инициализация глобального Minuit
        auto funcx = new TF1("funcx", bwfitfunc, 0.01, 10, 5);

        funcx->SetParameters(2,1);
        funcx->SetParNames("constant", "T", "beta", "mass", "pt");
        MyIntegFunc *integ = new MyIntegFunc(funcx);

        for (int part: PARTS)
        {  
            for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
                int centr = CENTR_SYST[systN][j];
   
                // cout << "PART: " << part << "   CENTR: " << centr << endl;
                string ifuncxName = "MyIntegFunc_" + to_string(part) + "_" + to_string(centr);
                ifuncx[part][centr] = new TF1("ifuncx", integ, xmin[part], xmax[part], 4, ifuncxName.c_str());

                switch(initParamsType)
                {
                    case 0: { /* DEFAULT */
                        // ================== version1 Params from Global fit ============================
                        double parResults[5];

                        std::string filename = "output/parameters/ALL_GlobalBWparams_" + std::string(systNamesT[systN]) + ".txt";
                        // std::string filename = "output/parameters/ALL_FinalBWparams_" + std::string(systNamesT[systN]) + ".txt";
                        ReadGlobalParams(systN, paramsGlobal, filename.c_str());
                        getGlobalParams(part, centr, parResults);
                        if (parResults[0] == 0) continue;
                            
                        // Установка начальных параметров
                        ifuncx[part][centr]->SetParameters(parResults);
                        
                        if (systN == 0) {
                            // ifuncx[part][centr]->FixParameter(1, parResults[1]); // Фиксируем T
                            // ifuncx[part][centr]->FixParameter(2, parResults[2]); // Фиксируем beta
                            ifuncx[part][centr]->SetParLimits(1, parResults[1] * 0.99, parResults[1] * 1.285);
                            ifuncx[part][centr]->SetParLimits(2, parResults[2] * 0.99, parResults[2] * 1.285);
                            ifuncx[part][centr]->SetParLimits(0, parResults[0] * 0, parResults[0] * 100000); // Широкие границы для константы
                        } else if (systN == 1) {
                            ifuncx[part][centr]->SetParLimits(1, parResults[1] * 1.05, parResults[1] * 1.1); // T ±1%
                            ifuncx[part][centr]->SetParLimits(2, parResults[2] * 1.05, parResults[2] * 1.1); // beta ±1%
                            ifuncx[part][centr]->SetParLimits(0, parResults[0] * 0, parResults[0] * 100000); // Широкие границы для константы
                        } else if (systN == 2) {
                            ifuncx[part][centr]->SetParLimits(1, parResults[1] * 0.99, parResults[1] * 1.01); // T ±1%
                            ifuncx[part][centr]->SetParLimits(2, parResults[2] * 1.1, parResults[2] * 1.2); // beta ±1%
                            ifuncx[part][centr]->SetParLimits(0, parResults[0] * 0, parResults[0] * 100000); // Широкие границы для константы
                        } else if (systN == 3) {
                            // ifuncx[part][centr]->FixParameter(1, parResults[1]); // Фиксируем T
                            ifuncx[part][centr]->SetParLimits(1, parResults[1] * 0.95, parResults[1] * 1.05); // T ±1%
                            ifuncx[part][centr]->SetParLimits(2, parResults[2] * 1.1, parResults[2] * 1.2); // beta ±1%
                            ifuncx[part][centr]->SetParLimits(0, parResults[0] * 0, parResults[0] * 100000); // Широкие границы для константы
                        } else if (systN == 4) {
                            // ifuncx[part][centr]->FixParameter(1, parResults[1]); // Фиксируем T
                            ifuncx[part][centr]->SetParLimits(1, parResults[1] * 0.85, parResults[1] * 0.9); // T ±1%
                            ifuncx[part][centr]->SetParLimits(2, parResults[2] * 1.1, parResults[2] * 1.2); // beta ±1%
                            ifuncx[part][centr]->SetParLimits(0, parResults[0] * 0, parResults[0] * 100000); // Широкие границы для константы
                        }
                        
                        ifuncx[part][centr]->FixParameter(3, masses[part]); // masses

                        TFitResultPtr fitResult = grSpectra[part][centr]->Fit(ifuncx[part][centr], "QR+SEX0", "", xmin[part], xmax[part]);

                        if(fitResult->Status() != 0 || fitResult->Chi2() < 1e-6) {
                            cout << "BAD FIT: part " << part << " centr " << centr 
                                << " status: " << fitResult->Status() << endl;
                            continue;
                        }

                        TMatrixDSym corr = fitResult->GetCorrelationMatrix();
                        if(corr.GetNrows() >= 3) { // Добавляем проверку размерности
                            double T_beta_corr = corr(1,2);
                            if(fabs(T_beta_corr) > 0.8) {
                                outParamsErr[part][centr][1] *= 1.2;
                                outParamsErr[part][centr][2] *= 1.2;
                            }
                        }

                        // Проверяем валидность результата
                        // if (fitResult->IsValid()) {

                        //     double chi2 = fitResult->Chi2();
                        //     int ndf = fitResult->Ndf();
                        //     double chi2_ndf = (ndf > 0) ? chi2 / ndf : -1;

                        //     std::cout << part << " "
                        //               << centr << " "
                        //               << "Chi2/NDF = " << chi2_ndf 
                        //               << " (Chi2 = " << chi2 
                        //               << ", NDF = " << ndf << ")\n" 
                        //               << std::endl;
                        // }
                        break;

                    } case 1: {
                        // ================== version2 Params from individual fit results =====================
                        double parResults[4];
                        ReadParams(part, centr, parResults);
                        if (parResults[0] == 0)
                            continue;
                            
                        ifuncx[part][centr]->SetParameters(parResults);
                        for (int par = 0; par < 3; par++)
                        {
                            ifuncx[part][centr]->SetParLimits(par, parResults[par] * 0.6, parResults[par] * 1.5);
                        }

                        ifuncx[part][centr]->FixParameter(3, masses[part]);
                        grSpectra[part][centr]->Fit(ifuncx[part][centr],  "QR+SEX0", "", xmin[part], xmax[part]);
                        break;

                    } case 2: {
                        // ================= version 2 Params with limits =================================
                        // ================ Version 2 - Improved individual fit initialization ================
                        double parResults[4];
                        ReadParams(part, centr, parResults);
                        
                        // Добавляем проверку на разумные значения параметров
                        if(parResults[1] < 0.1 || parResults[1] > 0.2) parResults[1] = 0.12; // T
                        if(parResults[2] < 0.3 || parResults[2] > 0.9) parResults[2] = 0.6;   // beta
                        
                        // Ужесточаем границы для параметров
                        const double T_range = 0.15;  // ±15%
                        const double beta_range = 0.2; // ±20%
                        
                        ifuncx[part][centr]->SetParameters(parResults);
                        ifuncx[part][centr]->SetParLimits(0, parResults[0]*0.5, parResults[0]*2.0);  // Константа
                        ifuncx[part][centr]->SetParLimits(1, parResults[1]*(1-T_range), parResults[1]*(1+T_range));
                        ifuncx[part][centr]->SetParLimits(2, parResults[2]*(1-beta_range), parResults[2]*(1+beta_range));
                        
                        // Добавляем регуляризацию
                        ifuncx[part][centr]->SetParameter(1, parResults[1]);
                        ifuncx[part][centr]->SetParameter(2, parResults[2]);
                        
                        // Улучшаем процесс фитирования
                        TFitResultPtr fitResult = grSpectra[part][centr]->Fit(
                            ifuncx[part][centr], "QR+SEX0", "", xmin[part], xmax[part]
                        );
                        
                        // Повторный фит с обновленными параметрами
                        if(fitResult->Chi2()/fitResult->Ndf() > 2.0) {
                            fitResult = grSpectra[part][centr]->Fit(
                                ifuncx[part][centr], "QR+SEX0M", "", xmin[part], xmax[part]
                            );
                        }
                        break;
                    } case 3: {
                        // ================= version 3 hand Params without Fit =============================+==
                        double handParams[4] = {handConst[part][centr], TCuAu[centr], betaCuAu[centr], masses[part]};
                        ifuncx[part][centr]->SetParameters(handParams);
                        break;

                    } case 4: {

                        // ========== Систематические вариации параметров ==========
                        ReadGlobalParams(systN, paramsGlobal, "output/parameters/ALL_GlobalBWparams_" + systNamesT[systN] + ".txt");
                        getGlobalParams(part, centr, paramsSystematics[part][centr]);

                        // Применяем модификации параметров в зависимости от currentSystematicType
                        switch(currentSystematicType) {
                            case 1: paramsSystematics[part][centr][1] *= 0.9; break; // T -10%
                            case 2: paramsSystematics[part][centr][1] *= 1.1; break; // T +10%
                            case 3: paramsSystematics[part][centr][2] *= 0.9; break; // beta -10%
                            case 4: paramsSystematics[part][centr][2] *= 1.1; break; // beta +10%
                            case 5: paramsSystematics[part][centr][0] *= 0.9; break; // Нормировка -10%
                            case 6: paramsSystematics[part][centr][0] *= 1.1; break; // Нормировка +10%
                            case 7: lLimitMult = rLimitMult = 0.9; break;       // Границы множителей
                            case 8: lLimitMult = rLimitMult = 1.1; break;
                            case 9: lLimitMultPi = 1.1; rLimitMultPi = 0.9; break;
                            case 10: lLimitMultPi = 0.9; rLimitMultPi = 1.1; break;
                        }
                        
                        // double minConst, maxConst;
                        // if (currentSystematicType >=7 && currentSystematicType <=10) {
                        //     if (part == 0) { // Для пионов
                        //         minConst = paramsSystematics[part][centr][0] * lLimitMultPi;
                        //         maxConst = paramsSystematics[part][centr][0] * rLimitMultPi;
                        //     } else {
                        //         minConst = paramsSystematics[part][centr][0] * lLimitMult;
                        //         maxConst = paramsSystematics[part][centr][0] * rLimitMult;
                        //     }
                        // } else {
                        //     minConst = paramsSystematics[part][centr][0] * 0.5;
                        //     maxConst = paramsSystematics[part][centr][0] * 1.5;
                        // }
                        // ifuncx[part][centr]->SetParLimits(0, minConst, maxConst);
                        
                        // ================ Version 4 - Improved systematic boundaries ================
                        TVirtualFitter::SetMaxIterations(1000);
                        paramsSystematics[part][centr][2] = std::min(paramsSystematics[part][centr][2], 0.95);

                        // Раздельные настройки для положительных и отрицательных частиц
                        double base_T_range, base_beta_range, const_range;
                        if(part == 0) { // Настройки для ПОЛОЖИТЕЛЬНЫХ частиц
                            base_T_range = 0.30;    // ±15% вместо 10%
                            base_beta_range = 0.30; // ±20% вместо 15%
                            const_range = 0.7;      // ±70% вместо 50%
                        } else { // Настройки для ОТРИЦАТЕЛЬНЫХ частиц
                            base_T_range = 0.30;
                            base_beta_range = 0.30;
                            const_range = 0.7;
                        }

                        // Адаптивные границы 
                        double ranges[3] = {
                            const_range,
                            base_T_range * (1.0 + 2.0*(paramsSystematics[part][centr][1] - 0.1)),
                            base_beta_range * (1.0 + (paramsSystematics[part][centr][2] - 0.5))
                        };

                        // Модифицированные ограничения для положительных частиц
                        for(int par = 0; par < 3; par++) {
                            double central = paramsSystematics[part][centr][par];
                            double min_val = central*(1.0 - ranges[par]);
                            double max_val = central*(1.0 + ranges[par]);

                            // Ослабление ограничений для положительных частиц
                            if(part == 0) {
                                if(par == 1) { // T
                                    min_val = std::max(min_val, 0.03);  // Нижний предел 0.03 вместо 0.05
                                    max_val = std::min(max_val, 0.30);   // Верхний предел 0.30 вместо 0.25
                                }
                                if(par == 2) { // beta
                                    min_val = std::max(min_val, 0.2);    // Нижний предел 0.2 вместо 0.3
                                    max_val = std::min(max_val, 1.0);    // Верхний предел 1.0 вместо 0.9
                                }
                            } else {
                                // Старые ограничения для отрицательных
                                if(par == 1) {
                                    min_val = std::max(min_val, 0.03);
                                    max_val = std::min(max_val, 0.30);
                                }
                                if(par == 2) {
                                    min_val = std::max(min_val, 0.2);
                                    max_val = std::min(max_val, 1.0);
                                }
                            }
                            
                            ifuncx[part][centr]->SetParLimits(par, min_val, max_val);
                        }

                        // Усиленный фит для положительных частиц
                        TFitResultPtr fitResult;
                        for(int refit = 0; refit < (part == 0 ? 3 : 2); ++refit) { // 3 итерации для положительных
                            fitResult = grSpectra[part][centr]->Fit(
                                ifuncx[part][centr], 
                                refit ? "QR+SEX0M" : "QR+SEX0", 
                                "", xmin[part], xmax[part]
                            );
                            if(fitResult->Chi2()/fitResult->Ndf() < (part == 0 ? 3.0 : 2.0)) break;
                        }

                        // Корректировка ошибок с учетом заряда частицы
                        TMatrixDSym corr = fitResult->GetCorrelationMatrix();
                        if(corr.GetNrows() >= 3) {
                            double T_beta_corr = corr(1,2);
                            if(part == 0) {
                                // Более агрессивная коррекция для положительных
                                if(fabs(T_beta_corr) > 0.6) {
                                    double scale = 1.0 + 0.8*fabs(T_beta_corr);
                                    outParamsErr[part][centr][1] *= scale;
                                    outParamsErr[part][centr][2] *= scale;
                                }
                            } else {
                                if(fabs(T_beta_corr) > 0.7) {
                                    outParamsErr[part][centr][1] *= 1.0 + 0.5*fabs(T_beta_corr);
                                    outParamsErr[part][centr][2] *= 1.0 + 0.5*fabs(T_beta_corr);
                                }
                            }
                        }

                        break;
                    }
                }
                
                ifuncx[part][centr]->SetLineColor(centrColors[centr]);
                

        // +++++++++ Metrics ++++++++++++++++++++++++++++++++++++

                double *params = ifuncx[part][centr]->GetParameters();
                const double *paramsErr = ifuncx[part][centr]->GetParErrors();
                std::copy(params, params + 4, outParams[part][centr]);
                std::copy(paramsErr, paramsErr + 4, outParamsErr[part][centr]);

                // NormalizeErrors on chi2/NDF
                double chi2 = ifuncx[part][centr]->GetChisquare();
                double ndf = ifuncx[part][centr]->GetNDF();
                double chi2Ndf = chi2 / ndf;
                for (int i = 0; i < 3; i++ )
                {
                    // outParamsErr[part][centr][i] *= sqrt(chi2Ndf);
                    double scale = std::min(std::max(sqrt(chi2Ndf), 0.8), 1.2);
                    outParamsErr[part][centr][i] *= scale; 
                }

                int N = 0, fitN = 0;
                double *x, *y, d = 0;
                x = grSpectra[part][centr]->GetX();
                y = grSpectra[part][centr]->GetY();
                N = grSpectra[part][centr]->GetN();

                for (int i = 0; i < N; i++)
                {
                    if (x[i] >= xmin[part] && x[i] <= xmax[part])
                    {
                        d += pow((y[i] - ifuncx[part][centr]->Eval(x[i])) / y[i], 2);
                        fitN++;
                    }
                }
                d = sqrt(d) / fitN;

                cout << part << " " << centr << " d: " << d << " chi2Ndf: " << chi2Ndf << endl;
            }
        }
    }
};