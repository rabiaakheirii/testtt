#include "input/headers/def.h"
#include "input/headers/WriteReadFiles.h"


// Определяем двумерный массив цветов: первая ось – заряд, вторая – система
// Новые цвета
Color_t systColors[6] = {kGray+2, kAzure+4, kSpring-5, kOrange+7, kViolet-6};
int chargeSystColors[2][5] = {
    // charge = 0 (положительные)          charge = 1 (отрицательные)
    {kAzure+4,  kTeal+3,   kOrange+7,  kViolet-6,  kPink+6},    // основные цвета
    {kAzure-2,  kTeal-7,   kOrange-3,  kViolet+2,  kPink-3}     // производные оттенки
};

// Стили маркеров
int markerStyles[2][5] = {
    {20, 21, 22, 23, 24},  // charge = 0
    {25, 26, 27, 28, 29}   // charge = 1
};

const int NSYST = 5; // Количество систем, например pAl, HeAu, CuAu, UU, и т.д.
double T_syst_err[NSYST][2][MAX_CENTR] = {0};
double beta_syst_err[NSYST][2][MAX_CENTR] = {0};

void DrawParam(TString paramName = "T")
{
    TCanvas *c2 = new TCanvas("", "", 30, 30, 1200, 1000);
    c2->cd();
    c2->SetGrid();
    c2->SetLogx();
    c2->SetLeftMargin(0.18);
    c2->SetRightMargin(0.05);
    c2->SetTopMargin(0.05);
    c2->SetBottomMargin(0.18);

    double ll = 1, rl = 10000., 
           pad_min = (paramName == "T") ? 0.05 : 0., 
           pad_max = (paramName == "T") ? 0.25 : 1., 
           pad_offset_x = 1, pad_offset_y = 1.1, 
           pad_tsize = 0.055, pad_lsize = 0.055;

    TString pad_title_x = "N_{part}";
    TString pad_title_y;
    if (paramName == "T") { 
        pad_title_y = "T [GeV]"; 
    } else if (paramName == "beta") {
        pad_title_y = "#beta [GeV]";
    } 

    Format_Pad(ll, rl, pad_min, pad_max, pad_title_x, pad_title_y, pad_offset_x, pad_offset_y, pad_tsize, pad_lsize, "", 8);        

    // Положение легенды в зависимости от параметра
    double leg_x1, leg_y1, leg_x2, leg_y2, 
           leg_x11, leg_y11, leg_x22, leg_y22;
    if (paramName == "T") {
        leg_x1 = 0.62; leg_y1 = 0.70; 
        leg_x2 = 0.92; leg_y2 = 0.90; 
        leg_x11 = 0.62; leg_y11 = 0.60; 
        leg_x22 = 0.92; leg_y22 = 0.69; 
    } else if (paramName == "beta") {
        leg_x1 = 0.62; leg_y1 = 0.40; 
        leg_x2 = 0.92; leg_y2 = 0.59; 
        leg_x11 = 0.62; leg_y11 = 0.30; 
        leg_x22 = 0.92; leg_y22 = 0.39; 
    }

    TLegend *legend = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetHeader("PHENIX, #sqrt{s_{NN}} = 200 GeV","C");
    legend->SetNColumns(2);
    legend->SetTextSize(0.04);

    TLegend *legend1 = new TLegend(0.195, 0.19, 0.32, 0.32); 
    legend1->SetBorderSize(0);
    legend1->SetFillStyle(0);
    legend1->SetNColumns(1); 
    legend1->SetTextSize(0.04);

    TLegend *legend2 = new TLegend(leg_x11, leg_y11, leg_x22, leg_y22); 
    legend2->SetBorderSize(0);
    legend2->SetFillStyle(0);
    legend2->SetHeader("PHENIX, #sqrt{s_{NN}} = 192 GeV","C");
    legend2->SetTextSize(0.04);
    legend2->SetNColumns(2);

    for (int systN: SYSTS)
    {
        for (int charge: {0, 1}) // БЫЛО: {0, 1}
        {
            // Рисуем систематические ошибки в виде TBox
            for (int i = 0; i < gr[charge][systN]->GetN(); i++) {
                Double_t x, y;
                gr[charge][systN]->GetPoint(i, x, y);

                // Получаем индекс центральности из массива CENTR_SYST
                int centr_index = CENTR_SYST[systN][i];
                
                // Проверка на выход за границы массива
                if(centr_index >= MAX_CENTR || centr_index < 0) {
                    cerr << "Invalid centr index: " << centr_index << endl;
                    continue;
                }
                
                double ey_syst = (paramName == "T") ? T_syst_err[systN][charge][centr_index] 
                                                    : beta_syst_err[systN][charge][centr_index]; 

                // Добавляем ограничение на максимальный размер ошибки
                // if (paramName == "T") {
                //     ey_syst = std::min(ey_syst, 0.05); // Максимум 0.05 GeV для T
                // }
                
                // Отладочный вывод 
                if (paramName == "T") {
                    cout << "SYST T: " << systNamesTT[systN] 
                        << " charge = " << charge 
                        << " centr = " << centr_index 
                        << " T = " << y << " T_syst " << ey_syst 
                        << endl;
                } else { 
                    cout << "SYST BETA: " << systNamesTT[systN] 
                        << " charge = " << charge 
                        << " centr = " << centr_index 
                        << " BETA = " << y << " beta_syst " << ey_syst 
                        << endl;
                }
                
                double ex_width = x * 0.08;
                double box_alpha = 0.3; // Прозрачность заливки

                TBox *box = new TBox(x - ex_width, y - ey_syst, x + ex_width, y + ey_syst);
                box->SetFillColorAlpha(chargeSystColors[charge][systN], box_alpha);
                box->SetLineColor(chargeSystColors[charge][systN]);
                box->SetLineWidth(2);
                box->SetFillStyle(1001); // Сплошная заливка с прозрачностью
                box->Draw("SAME");
            }

            gr[charge][systN]->Draw("P E SAME");
            if (systN == 4) {
                legend2->AddEntry(gr[charge][systN], systNamesTT[systN], "P");
            } else {
                legend->AddEntry(gr[charge][systN], systNamesTT[systN], "P");
            }
        }
    }

    if (paramName == "T" || paramName == "beta") 
    {
        const int nPoints = 13; 
    
        double yValues[nPoints];
    
        // Добавляем точки ARTICLE для AuAu
        if (paramName == "T") {
            for (int i = 0; i < nPoints; i++) yValues[i] = T_AuAu_ART[i];
        } else if (paramName == "beta") {
            for (int i = 0; i < nPoints; i++) yValues[i] = beta_AuAu_ART[i];
        }
    
        TGraph *theoryGraph = new TGraph(nPoints, Npart[0], yValues);
        theoryGraph->SetMarkerStyle(43); 
        theoryGraph->SetMarkerColor(kMagenta);
        theoryGraph->SetMarkerSize(2.5);
        theoryGraph->SetLineColor(kMagenta);
        theoryGraph->SetLineStyle(2); 
        theoryGraph->Draw("P SAME");
    
        legend1->AddEntry(theoryGraph, "Au+Au, #sqrt{s_{NN}} = 200 GeV (C72:014903, 2005)", "P");

        // // Добавляем точки STAR для AuAu
        // TGraph *starAuAu;
        // if (paramName == "T") {
        //     starAuAu = new TGraph(9, Npart_AuAu_STAR, T_AuAu_STAR);
        // } else { // для "beta"
        //     starAuAu = new TGraph(9, Npart_AuAu_STAR, beta_AuAu_STAR);
        // }
        // starAuAu->SetMarkerStyle(20);  // можно выбрать другой стиль
        // starAuAu->SetMarkerColor(kBlue);  // или другой цвет
        // starAuAu->SetMarkerSize(1.5);
        // starAuAu->Draw("P SAME");
        // legend->AddEntry(starAuAu, "AuAu_{STAR}", "P");

        // // Добавляем точки STAR для UU
        // TGraph *starUU;
        // if (paramName == "T") {
        //     starUU = new TGraph(9, Npart_UU_STAR, T_UU_STAR);
        // } else {
        //     starUU = new TGraph(9, Npart_UU_STAR, beta_UU_STAR);
        // }
        // starUU->SetMarkerStyle(21); 
        // starUU->SetMarkerColor(kGreen);
        // starUU->SetMarkerSize(1.5);
        // starUU->Draw("P SAME");
        // legend->AddEntry(starUU, "UU_{STAR}", "P");
    }
    legend->Draw();
    legend1->Draw();
    legend2->Draw();

    // double avgValue = (paramName == "T") ? gAvgT : gAvgUt;
    
    // TLine *avgLine = new TLine(ll, avgValue, rl, avgValue);
    // avgLine->SetLineColor(kBlack);
    // avgLine->SetLineStyle(9); 
    // avgLine->SetLineWidth(2);
    // avgLine->Draw("SAME");

    // double avgLeg_x1 = 0.2, avgLeg_y1 = 0.21; 
    // double avgLeg_x2 = 0.4, avgLeg_y2 = 0.31;

    // TLegend *avgLegend = new TLegend(avgLeg_x1, avgLeg_y1, avgLeg_x2, avgLeg_y2);
    // avgLegend->SetBorderSize(0);
    // avgLegend->SetFillStyle(0);
    // avgLegend->SetTextSize(0.05);

    // TString avgLabel = (paramName == "T") 
    //     ? Form("T_{avr} = %.3f GeV", avgValue) 
    //     : Form("#beta_{avr} = %.3f GeV", avgValue);

    // avgLegend->AddEntry(avgLine, avgLabel, "L");
    // avgLegend->Draw();

    c2->SaveAs("output/pics/SYST_BWparamFinal_" + paramName + ".png");
    c2->SaveAs("output/pics/SYST_BWparamFinal_" + paramName + ".pdf");

    // delete c2;
    // delete legend;
    // delete legend1;
    // delete avgLine;
    // delete avgLegend;
}


void SetGraphs( int systN, TString paramName, TString fitType = "FINAL" ) {
    double xerr[MAX_CENTR];
    for (int i: CENTR_SYST[systN]) xerr[i] = 0.;

    TString filename;
    if( fitType == "GLOBAL" )
    {
        filename = "output/parameters/ALL_GlobalBWparams_" + systNamesT[systN] + ".txt";
    }
    else if( fitType == "FINAL" )
    {
        filename = "output/parameters/ALL_FinalBWparams_" + systNamesT[systN] + ".txt";
    }

    cout << "Reading file: " << filename << endl;
    
    ReadGlobalParams(systN, paramsGlobal, filename);

    TString systFilename = "output/parameters/SYSTErr_" + systNamesT[systN] + ".txt";
    ifstream systFile(systFilename.Data());

    cout << "Reading file: " << systFilename << endl;

    if (systFile.is_open()) {
        int part, centr;
        double constanta, T, T_stat, T_syst, beta, beta_stat, beta_syst;
        while (systFile >> part >> centr >> constanta >> T >> T_stat 
                        >> T_syst >> beta >> beta_stat >> beta_syst) {
            
            // cout << part << " " << centr << " " << constanta << " " 
            //      << T << " " << T_stat << " " << T_syst << " " 
            //      << beta << " " << beta_stat << " " << beta_syst << endl;

            // cout << "T " << T << " T_syst " << T_syst << " beta " 
            //      << beta << " beta_syst " << beta_syst << endl;
            
            T_syst_err[systN][part][centr] = T_syst;
            beta_syst_err[systN][part][centr] = beta_syst;

            cout << part << " " << centr << " T_syst " << T_syst_err[systN][part][centr] 
                 << " beta_syst " << beta_syst_err[systN][part][centr] << endl;
        }
        systFile.close();
    } else {
        cerr << "Error opening systematic file: " << systFilename << endl;
    }

    cout << "СЧИТАЛИ ФАЙЛ: " << filename << " ГРАФИКОВ: " << N_CENTR_SYST[systN] << endl;
    cout << "СЧИТАЛИ ФАЙЛ: " << "\n" << endl;

    // Заполнение массивов Tpar и utPar в зависимости от типа
    for (int charge: {0, 1}) // БЫЛО: {0, 1}
    {
        // Добавляем массивы для ошибок
        double T_err[MAX_CENTR], beta_err[MAX_CENTR];

        for (int centr = 0; centr < N_CENTR_SYST[systN]; centr++)
        {
            if( fitType == "GLOBAL" )
            {
                Tpar[charge][centr] = paramsGlobal[charge][centr][1];  // T value
                T_err[centr] = paramsGlobal[charge][centr][2];         // T error (индекс 2)
                utPar[charge][centr] = paramsGlobal[charge][centr][3];  // beta value
                beta_err[centr] = paramsGlobal[charge][centr][4];       // beta error (индекс 4)
            }
            else if( fitType == "FINAL" )
            {
                Tpar[charge][centr] = paramsGlobal[charge][centr][1];   // T value
                T_err[centr] = paramsGlobal[charge][centr][2];          // T error
                utPar[charge][centr] = paramsGlobal[charge][centr][3]; // beta value
                beta_err[centr] = paramsGlobal[charge][centr][4];      // beta error
            }
        }

        // Создаем графики ВНУТРИ цикла по charge
        if(paramName == "T") {
            gr[charge][systN] = new TGraphErrors(
                N_CENTR_SYST[systN], Npart[systN], 
                Tpar[charge], xerr, T_err
            );
        } else if(paramName == "beta") {
            gr[charge][systN] = new TGraphErrors(
                N_CENTR_SYST[systN], Npart[systN], 
                utPar[charge], xerr, beta_err
            );
        }

        gr[charge][systN]->SetMarkerStyle(markerStyles[charge][systN]);
        gr[charge][systN]->SetMarkerColor(chargeSystColors[charge][systN]);
        gr[charge][systN]->SetLineColor(chargeSystColors[charge][systN]);
        gr[charge][systN]->SetLineWidth(2);
        gr[charge][systN]->SetMarkerSize(2.3);
    }

    // cout << "DEFINE GRAPHS " << N_CENTR_SYST[systN] << endl;
}


void DrawTbeta(TString fitType = "FINAL") 
{
    TCanvas *c3 = new TCanvas("c3", "T vs u_T", 30, 30, 1200, 1000);
    c3->cd();
    c3->SetGrid();
    c3->SetLeftMargin(0.18);
    c3->SetRightMargin(0.05);
    c3->SetTopMargin(0.05);
    c3->SetBottomMargin(0.18);
    
    TString pad_title_x = "#beta [GeV]";
    TString pad_title_y = "T [GeV]";
    double pad_min_x = 0.1, pad_max_x = 1.20;
    double pad_min_y = 0.05, pad_max_y = 0.25;
    
    Format_Pad(pad_min_x, pad_max_x, pad_min_y, pad_max_y, pad_title_x, pad_title_y, 1, 1.1, 0.055, 0.055, "", 8);
    
    TLegend *legend = new TLegend(0.62, 0.70, 0.92, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetHeader("PHENIX, #sqrt{s_{NN}} = 200 GeV","C");
    legend->SetNColumns(2);
    legend->SetTextSize(0.04);

    TLegend *legend1 = new TLegend(0.195, 0.20, 0.32, 0.32); 
    legend1->SetBorderSize(0);
    legend1->SetFillStyle(0);
    legend1->SetNColumns(1); 
    legend1->SetTextSize(0.04); 

    TLegend *legend2 = new TLegend(0.62, 0.60, 0.92, 0.69); 
    legend2->SetBorderSize(0);
    legend2->SetFillStyle(0);
    legend2->SetHeader("PHENIX, #sqrt{s_{NN}} = 192 GeV","C");
    legend2->SetTextSize(0.04);
    legend2->SetNColumns(2);
    
    // Считываем данные из файлов для каждого systN и charge
    for (int systN : SYSTS)
    {
        for (int charge: {0, 1}) // БЫЛО: {0, 1}
        {
            double T_err[MAX_CENTR], beta_err[MAX_CENTR];

            TString filename;
            if (fitType == "GLOBAL")
            {
                filename = "output/parameters/ALL_GlobalBWparams_" + systNamesT[systN] + ".txt";
            }
            else if (fitType == "FINAL")
            {
                filename = "output/parameters/ALL_FinalBWparams_" + systNamesT[systN] + ".txt";
            }
            
            ReadGlobalParams(systN, paramsGlobal, filename);
            
            // Заполняем массивы Tpar и utPar
            for (int centr = 0; centr < N_CENTR_SYST[systN]; centr++)
            {
                if (fitType == "GLOBAL")
                {
                    Tpar[charge][centr] = paramsGlobal[charge][centr][0];
                    utPar[charge][centr] = paramsGlobal[charge][centr][1];
                }
                else if (fitType == "FINAL") ///////////// ТУТ ВКЛЮЧАЮТСЯ ВЫСЧИТАННЫЕ АВТОМАТИЧЕСКИ ОШИБКИ ////////////
                {
                    Tpar[charge][centr] = paramsGlobal[charge][centr][1];
                    T_err[centr] = paramsGlobal[charge][centr][2];    // Ошибка T (индекс 2)
                    utPar[charge][centr] = paramsGlobal[charge][centr][3];
                    beta_err[centr] = paramsGlobal[charge][centr][4]; // Ошибка beta (индекс 4)
                }
            }
            
            // Создаем график T vs beta
            TGraphErrors *gr_TvsUt = new TGraphErrors(N_CENTR_SYST[systN], utPar[charge], Tpar[charge], beta_err, T_err);
            gr_TvsUt->SetMarkerStyle(markerStyles[charge][systN]);
            gr_TvsUt->SetMarkerSize(2.5);
            gr_TvsUt->SetMarkerColor(chargeSystColors[charge][systN]);
            gr_TvsUt->SetLineColor(chargeSystColors[charge][systN]); // Цвет усов
            gr_TvsUt->SetMarkerColor(chargeSystColors[charge][systN]); // Цвет точек
            gr_TvsUt->SetLineWidth(2); // Толщина линий ошибок

            // // Создаем отдельный график для 10% усов
            // const int nPoints = gr_TvsUt->GetN();
            // std::vector<double> x_10(nPoints), y_10(nPoints);
            // std::vector<double> ex_10(nPoints), ey_10(nPoints);

            // for(int i = 0; i < nPoints; i++) 
            // {
            //     Double_t x, y;
            //     gr_TvsUt->GetPoint(i, x, y);
                
            //     // 10% ошибки от текущих значений
            //     ex_10[i] = x * 0.10;
            //     ey_10[i] = y * 0.10;
            //     x_10[i] = x;
            //     y_10[i] = y;
            // }

            // TGraphErrors *gr_10perc = new TGraphErrors(nPoints, x_10.data(), y_10.data(), 
            //                                          ex_10.data(), ey_10.data());
            
            // // Настройка стиля усов
            // gr_10perc->SetLineColor(chargeSystColors[charge][systN]);
            // gr_10perc->SetLineStyle(1);    // Пунктирная линия
            // gr_10perc->SetLineWidth(2);    // Толщина линии
            // gr_10perc->SetMarkerSize(0);   // Скрываем маркеры
            // gr_10perc->Draw("E same");     // Рисуем только усы

            for (int i = 0; i < gr_TvsUt->GetN(); ++i) {
                Double_t x, y;
                gr_TvsUt->GetPoint(i, x, y);
                int centr_index = CENTR_SYST[systN][i];
                if (centr_index < 0 || centr_index >= MAX_CENTR) continue;
            
                double ex_syst = beta_syst_err[systN][charge][centr_index];
                double ey_syst = T_syst_err[systN][charge][centr_index];
            
                TBox *box = new TBox(x - ex_syst, y - ey_syst,
                                    x + ex_syst, y + ey_syst);
                box->SetFillColorAlpha(chargeSystColors[charge][systN], 0.3);
                box->SetLineColor(chargeSystColors[charge][systN]);
                box->SetLineWidth(1);
                box->SetFillStyle(1001);
                box->Draw("SAME");
            }

            gr_TvsUt->Draw("P E SAME"); // Основные точки поверх усов
            gr_TvsUt->Draw("P SAME"); // Основные точки поверх усов


        if (systN == 4) {
            legend2->AddEntry(gr_TvsUt, systNamesTT[systN], "P");

        } else {
            legend->AddEntry(gr_TvsUt, systNamesTT[systN], "P");
        }
    }
}
    
    // Добавляем точки ARTICLE для AuAu
    const int nPoints = 13;
    TGraph *articleAuAu = new TGraph(nPoints, beta_AuAu_ART, T_AuAu_ART);
    articleAuAu->SetMarkerStyle(43);
    articleAuAu->SetMarkerColor(kMagenta);
    articleAuAu->SetMarkerSize(2.5);
    articleAuAu->SetLineColor(kMagenta);
    articleAuAu->SetLineStyle(2);
    articleAuAu->Draw("P SAME");
    legend1->AddEntry(articleAuAu, "Au+Au, #sqrt{s_{NN}} = 200 GeV (C72:014903, 2005)", "P");
    
    // +++ ДОБАВЛЯЕМ НОВУЮ ТОЧКУ MPD: Bi+Bi, 9.2 GeV +++
    const int nMPD = 1;
    double T_mpd[nMPD] = {0.109};
    double beta_mpd[nMPD] = {0.6};
    double T_err_mpd[nMPD] = {0.002};
    double beta_err_mpd[nMPD] = {0.11};
    
    TGraphErrors *grMPD = new TGraphErrors(nMPD, beta_mpd, T_mpd, beta_err_mpd, T_err_mpd);
    grMPD->SetMarkerStyle(29); // Звездочка
    grMPD->SetMarkerColor(kRed);
    grMPD->SetMarkerSize(5.0);
    grMPD->SetLineColor(kRed);
    grMPD->SetLineWidth(2);
    grMPD->Draw("P E SAME");
    legend1->AddEntry(grMPD, "MPD: Bi+Bi, #sqrt{s_{NN}} = 9.2 GeV", "P");
    
    // // Добавляем точки STAR для AuAu
    // TGraph *starAuAu = new TGraph(9, beta_AuAu_STAR, T_AuAu_STAR);
    // starAuAu->SetMarkerStyle(20);
    // starAuAu->SetMarkerColor(kBlue);
    // starAuAu->SetMarkerSize(1.5);
    // starAuAu->Draw("P SAME");
    // legend->AddEntry(starAuAu, "AuAu_{STAR}", "P");
    
    // // Добавляем точки STAR для UU
    // TGraph *starUU = new TGraph(9, beta_UU_STAR, T_UU_STAR);
    // starUU->SetMarkerStyle(21);
    // starUU->SetMarkerColor(kGreen);
    // starUU->SetMarkerSize(1.5);
    // starUU->Draw("P SAME");
    // legend->AddEntry(starUU, "UU_{STAR}", "P");
    
    legend->Draw();
    legend1->Draw();
    legend2->Draw();
    
    // Добавляем линии средних значений
    // Горизонтальная линия для среднего T
    // TLine *avgLineT = new TLine(pad_min_x, gAvgT, pad_max_x, gAvgT);
    // avgLineT->SetLineColor(kBlack);
    // avgLineT->SetLineStyle(9); 
    // avgLineT->SetLineWidth(2);
    // avgLineT->Draw("SAME");

    // // Вертикальная линия для среднего beta
    // TLine *avgLineBeta = new TLine(gAvgUt, pad_min_y, gAvgUt, pad_max_y);
    // avgLineBeta->SetLineColor(kBlack);
    // avgLineBeta->SetLineStyle(9);
    // avgLineBeta->SetLineWidth(2);
    // avgLineBeta->Draw("SAME");

    // // Обновляем легенду
    // TLegend *avgLegend = new TLegend(0.2, 0.22, 0.41, 0.32);
    // avgLegend->SetBorderSize(0);
    // avgLegend->SetFillStyle(0);
    // avgLegend->SetTextSize(0.045);
    // avgLegend->AddEntry(avgLineT, Form("T_{avr} = %.3f GeV", gAvgT), "L");
    // avgLegend->AddEntry(avgLineBeta, Form("#beta_{avr} = %.3f GeV", gAvgUt), "L");
    // avgLegend->Draw();

    // Сохраняем график
    c3->SaveAs("output/pics/SYST_BWparamFinal_T_beta.png");
    c3->SaveAs("output/pics/SYST_BWparamFinal_T_beta.pdf");

    // delete c3;
    // delete legend;
    // delete legend1;
}


void NpartDrawParams_syst (void) 
{
    // cout << "=== Loaded systematics for " << systNamesTT[systN] << " ===" << endl;
    // for (int charge: {0, 1}) // БЫЛО: {0, 1}
    // {
    //     for (int centr : CENTR_SYST[systN]) {
    //         cout << "System " << systNamesTT[systN] 
    //             << " charge " << charge 
    //             << " centr " << centr
    //             << ": T_syst = " << T_syst_err[charge][centr]
    //             << ", beta_syst = " << beta_syst_err[charge][centr] << endl;
    //     }
    // }

    // Для параметра T
    for (int systN: SYSTS) 
    {
        cout << systN << " " << systNamesTT[systN] << endl;
        SetGraphs(systN, "T");
    }
    // CalculateAverage("T"); 
    DrawParam("T");
    
    // Для параметра beta
    for (int systN: SYSTS) 
    {
        cout << systN << " " << systNamesTT[systN] << endl;
        SetGraphs(systN, "beta");
    }
    // CalculateAverage("beta"); 
    DrawParam("beta");
    DrawTbeta();

    WriteAveragesToFile();
    // DrawTbeta();
}