#include "input/headers/def.h"
#include "input/headers/WriteReadFiles.h"
#include "input/headers/BlastWaveFit.h"

/* ОШИБКИ РАБОТАЮТ ОТВРАТИТЕЛЬНО, ПОКА ЧТО ОНИ УБРАНЫ В 249 СТРОЧКЕ */

using namespace std;

void GetMtMinusM0Range(int particle, double& xmin, double& xmax) {
    if (particle == 0 || particle == 1) { // Пионы
        xmin = 0.35;
        xmax = 1.2;  // Оптимально для пионов
    }
    else if (particle == 2 || particle == 3) { // Каоны
        xmin = 0.05;
        xmax = 1.2;
    }
    else { // Протоны и антипротоны
        xmin = 0.05;
        xmax = 1.2;
    }
}

// Преобразование pT -> mT - m0
double ConvertPtToMtMinusM0(double pt, double mass) {
    double mt = sqrt(pt*pt + mass*mass);
    return mt - mass; 
}

// Создание преобразованного графика
TGraphErrors* CreateMtMinusM0Graph(TGraphErrors* gr, double mass) {
    TGraphErrors* grMT = new TGraphErrors();
    int npoints = gr->GetN();
    for (int i = 0; i < npoints; i++) {
        double pt, yield;
        gr->GetPoint(i, pt, yield);
        double errX = gr->GetErrorX(i);
        double errY = gr->GetErrorY(i);
        
        double mt_minus_m0 = ConvertPtToMtMinusM0(pt, mass);
        double mt = mt_minus_m0 + mass;
        double errX_transformed = (pt / mt) * errX;  // Преобразование ошибки
        
        grMT->SetPoint(i, mt_minus_m0, yield);
        grMT->SetPointError(i, errX_transformed, errY);
    }
    return grMT;
}

// Создаёт TGraph (без ошибок) из TGraphErrors
TGraph* CreateTGraphWithoutErrors(TGraphErrors* gr) {
    TGraph* grSimple = new TGraph();
    int npoints = gr->GetN();
    
    for (int i = 0; i < npoints; i++) {
        double x, y;
        gr->GetPoint(i, x, y);
        grSimple->SetPoint(grSimple->GetN(), x, y);
    }
    
    // Копируем стили
    grSimple->SetMarkerStyle(gr->GetMarkerStyle());
    grSimple->SetMarkerSize(gr->GetMarkerSize());
    grSimple->SetMarkerColor(gr->GetMarkerColor());
    grSimple->SetLineColor(gr->GetLineColor());
    grSimple->SetLineWidth(gr->GetLineWidth());
    
    return grSimple;
}

// Создание преобразованной функции
TF1* CreateMtMinusM0Function(TF1* originalFunc, double mass, double xmin, double xmax) {
    return new TF1(
        Form("%s_MT", originalFunc->GetName()), 
        [=](double* x, double* par) {
            double mt_minus_m0 = x[0];
            double mT = mt_minus_m0 + mass;
            double pt = (mT > mass) ? sqrt(mT*mT - mass*mass) : 0.0;
            return originalFunc->Eval(pt);
        },
        xmin, xmax, 0
    );
}

// Функция для получения строки с энергией столкновения
TString GetCollisionEnergyLabel(int systN) {
    if (systN == 4) { // U+U
        return "#sqrt{s_{NN}} = 193 GeV";
    } else { // Все остальные системы
        return "#sqrt{s_{NN}} = 200 GeV";
    }
}

// Эта функция рисует спектры двух заданных частиц и сохраняет график
void DrawSpectraPart( TString partName, int part1, int part2 )
{
    TCanvas *c4 = new TCanvas("c2", "c2", 30, 30, 1200, 1200);
    int NpadX = 1, NpadY = 2;
    Format_Canvas(c4, NpadX, NpadY, 0);

    int padN = 1;

    std::vector<TGraphErrors*> tempGraphs;
    std::vector<TF1*> tempFunctions;
    
    // Проходим по двум переданным номерам частиц
    for (int i: {part1, part2})
    {
        c4->SetLogy();
        c4->cd(padN);  

        double shiftX = (i % NpadX == 0) ? 0 : 0.1;
        double texScale = (padN == 1 ) ? 1 : 0.9;
        TLegend *legend = new TLegend(0.5 - shiftX, 0.65, 0.9 - shiftX, 0.95); //1 column
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetNColumns(1);
        legend->SetTextSize(0.073 * texScale);

        TLatex *titleTex = new TLatex(0.4, 2000, partTitles[i].c_str());
        titleTex->SetTextFont(42);
        titleTex->SetTextSize(0.09);
        titleTex->SetLineWidth(2 * texScale);

        FormatSpectraPad(texScale);

        // Проходим по разным классам центральности, 
        // если данные для центральности существуют
        for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
            int centr = CENTR_SYST[systN][j];
            if (!ifuncx[i][centr]) continue;

            // Определяем границы для текущей частицы
            double xmin_mt, xmax_mt;
            GetMtMinusM0Range(i, xmin_mt, xmax_mt);

            // Преобразование данных
            TGraphErrors* grMT = CreateMtMinusM0Graph(grSpectra[i][centr], masses[i]);
            
            // Преобразование функции с индивидуальными границами
            TF1* funcMT = CreateMtMinusM0Function(ifuncx[i][centr], masses[i], xmin_mt, xmax_mt);
            
            grMT->Draw("P SAME");
            funcMT->Draw("SAME");

            legend->AddEntry(grMT, centrTitles[centr].c_str(), "p");

            // grSpectra[i][centr]->SetMarkerColor(centrColors[centr]);
            // grSpectra[i][centr]->SetMarkerSize(2);
            // grSpectra[i][centr]->SetMarkerStyle(8);
            // grSpectra[i][centr]->Draw("P SAME");      
            // ifuncx[i][centr]->Draw("SAME");
            // legend->AddEntry(grSpectra[i][centr], centrTitles[centr].c_str(), "p");  

            tempGraphs.push_back(grMT);
            tempFunctions.push_back(funcMT);      
        }
        legend->Draw();
        titleTex->Draw();    
        padN++;

        for (auto* gr : tempGraphs) delete gr;
        for (auto* func : tempFunctions) delete func;

        gROOT->ProcessLine(".q");
    }

    // c4->SaveAs("output/pics/ALL_BlastWaveFinal_" + systNamesT[systN] + "_" + partName + ".png");
}


// Основная функция анализа
void BlastWaveFinal_all( void )
{
    bool isContour = false;
    bool isDraw = true;

    // Чтение данных в зависимости от системы
    if (systN == 0) 
        ReadFromFileAuAu(); // Для системы AuAu
    else                    // Для других систем
        for (int part: PARTS) ReadFromFile(part, systN);

    // Фитируем определённым кейсом от 0 до 4
    BlastWaveFit *bwFit = new BlastWaveFit();
    bwFit->isContour = true; 
    bwFit->Fit(0);

    WriteParams(systN, bwFit->outParams, bwFit->outParamsErr, true, "output/parameters/ALL_FinalBWparams_" + systNamesT[systN] + ".txt");
    // WriteParams(systN, bwFit->outParams, bwFit->outParamsErr, false, "output/parameters/ALL_FinalBWparams_" + systNamesT[systN] + ".txt");
   
    if (!isDraw)
        return;

    // ++++++ Draw spectra All +++++++++++++++++++++++++++++++++++++

    TCanvas *c2 = new TCanvas("c2", "c2", 29, 30, 1840, 2160);
    Format_Canvas(c2, 2, 3, 0);

    std::vector<TGraphErrors*> tempGraphs;
    std::vector<TF1*> tempFunctions;

    // Цикл по всем частицам 
    for (int i: {0, 1, 2, 3, 4, 5})
    {
        c2->SetLogy();
        c2->cd(i + 1);  
        
        double shiftX = (i % 2 == 0) ? 0 : 0.1;
        double texScale = (i < 3) ? 1 : 0.95;

        TLegend *legend = new TLegend(0.5 - shiftX, 0.55, 0.95 - shiftX, 0.88); 
        if (systN == 0) {
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            legend->SetNColumns(2); // Устанавливаем два столбца
            legend->SetTextSize(0.075 * texScale);
            legend->SetColumnSeparation(0.05); // Разделение между столбцами
        } else {
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            legend->SetNColumns(1);
            legend->SetTextSize(0.075 * texScale);
        }

        TLatex *titleTex = new TLatex(0.3, 5000000, partTitles[i].c_str());
        titleTex->SetTextFont(42);
        titleTex->SetTextSize(0.1);
        titleTex->SetLineWidth(2 * texScale);

        FormatSpectraPad(texScale);
        for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
            int centr = CENTR_SYST[systN][j];
            if (!ifuncx[i][centr]) continue;

            // Определяем границы для текущей частицы
            double xmin_mt, xmax_mt;
            GetMtMinusM0Range(i, xmin_mt, xmax_mt);

            // Преобразование данных
            TGraphErrors* grMT = CreateMtMinusM0Graph(grSpectra[i][centr], masses[i]);
            
            // Преобразование функции с индивидуальными границами
            TF1* funcMT = CreateMtMinusM0Function(ifuncx[i][centr], masses[i], xmin_mt, xmax_mt);
            
            // Настройка стилей
            grMT->SetMarkerColor(centrColors[centr]);
            grMT->SetMarkerSize(2);
            grMT->SetMarkerStyle(8);
            grMT->SetLineColor(centrColors[centr]);
            grMT->SetLineWidth(0);
            
            funcMT->SetLineColor(centrColors[centr]);
            
            // Отрисовка
            grMT->Draw("P SAME");
            funcMT->Draw("L SAME");
            
            // Сохраняем для последующего удаления
            // tempGraphs.push_back(grMT);
            // tempFunctions.push_back(funcMT);

            // grSpectra[i][centr]->SetMarkerColor(centrColors[centr]);
            // grSpectra[i][centr]->SetMarkerSize(2);
            // grSpectra[i][centr]->SetMarkerStyle(8);
            // grSpectra[i][centr]->SetLineColor(centrColors[centr]);
            // grSpectra[i][centr]->SetLineWidth(1);
            // grSpectra[i][centr]->Draw("P SAME");      
            
            // ifuncx[i][centr]->SetLineColor(centrColors[centr]);
            // ifuncx[i][centr]->Draw("L SAME");
            
            // Выбираем подпись центральности из нужного подмассива
            TString centrLabel;
            if (systN == 0)
                centrLabel = centrStr[4][centr];  // AuAu – последний подмассив
            else if (systN == 1)
                centrLabel = centrStr[0][centr];  // pAl – первый подмассив
            else if (systN == 2)
                centrLabel = centrStr[1][centr];  // HeAu – второй подмассив
            else if (systN == 3)
                centrLabel = centrStr[2][centr];  // CuAu – третий подмассив
            else if (systN == 4)
                centrLabel = centrStr[3][centr];  // UU – четвёртый подмассив
            
            // legend->AddEntry(grSpectra[i][centr], centrLabel.Data(), "pf");     
            // if (systN == 0) {
            //     legend->AddEntry(grSpectra[i][centr], Form("%s", centrLabel.Data()), "p");
            // } else {
            //     legend->AddEntry(grSpectra[i][centr], Form("%s, x10^{%d}", centrLabel.Data(), j), "p");
            // }

            if (systN == 0) {
                legend->AddEntry(grMT, Form("%s", centrLabel.Data()), "p");
            } else {
                legend->AddEntry(grMT, Form("%s, x10^{%d}", centrLabel.Data(), j), "p");
            }
        }
        legend->Draw();
        titleTex->Draw();    

        // Добавляем информацию об энергии только для π⁺ (i=0)
        if (i == 0) { // Только для π⁺
            TString energyLabel;
            if (systN == 4) { // U+U система
                energyLabel = "#sqrt{s_{NN}} = 193 GeV";
            } else {
                energyLabel = "#sqrt{s_{NN}} = 200 GeV";
            }
            
            TLatex* energyTex = new TLatex(0.3, 0.15, energyLabel);
            energyTex->SetNDC();
            energyTex->SetTextFont(42);
            energyTex->SetTextSize(0.1 * texScale);
            energyTex->SetTextAlign(12); // Выравнивание по левому краю и центру по вертикали
            energyTex->Draw();
        }

        if (i == 0) {
            TString collisionLabelAll = Form("%s", systNamesTT[systN].Data());
            TLatex *systemLabelAll = new TLatex(0.3, 0.8, collisionLabelAll);
            systemLabelAll->SetNDC();
            systemLabelAll->SetTextFont(42);
            systemLabelAll->SetTextSize(0.1);
            systemLabelAll->Draw();
        }

    }

    c2->SaveAs("output/pics/spectra/ALL_BlastWaveFinal_" + systNamesT[systN] + ".png");
    c2->SaveAs("output/pics/spectra/ALL_BlastWaveFinal_" + systNamesT[systN] + ".pdf");

    delete c2;
    for (auto* gr : tempGraphs) delete gr;
    for (auto* func : tempFunctions) delete func;

    DrawSpectraPart("pi", 0, 1);
    DrawSpectraPart("K", 2, 3);
    DrawSpectraPart("p", 4, 5);










    

    // //++++++++ Draw Contour plots ++++++++++++++++++++++++++++++

    // if (!isContour) {
    //     gROOT->ProcessLine(".q");
    //     return;
    // }

    // TCanvas *c3 = new TCanvas("c3", "c3", 29, 30, 1840, 2160);
    // c3->cd();
    // c3->SetGrid();

    // // 1. Добавим проверку инициализации контуров
    // bool contoursExist = false;

    // for (int part : PARTS) {
    //     for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
    //         int centr = CENTR_SYST[systN][j];
    //         for (int nsigma = 1; nsigma < N_SIGMA; nsigma++) {
    //             if (contour[part][centr][nsigma]) {
    //                 contoursExist = true;
    //                 break;
    //             }
    //         }
    //         if (contoursExist) break;
    //     }
    //     if (contoursExist) break;
    // }

    // if (!contoursExist) {
    //     cerr << "Error: No contour plots found!" << endl;
    //     return;
    // }

    // // 2. Форматирование pad вынесем после проверок
    // double ll = 0.00001, rl = 2, pad_min = 0.00001, pad_max = 0.2, 
    //     pad_offset_x = 1., pad_offset_y = 1., 
    //     pad_tsize = 0.05, pad_lsize = 0.05;
    // TString pad_title_y = "T";
    // TString pad_title_x = "#beta";
    // Format_Pad(ll, rl, pad_min, pad_max, pad_title_x, pad_title_y, 
    //         pad_offset_x, pad_offset_y, pad_tsize, pad_lsize, "", 8);

    // TLegend *legendContour = new TLegend(0.6, 0.35, 0.85, 0.85);
    // legendContour->SetBorderSize(0);
    // legendContour->SetFillStyle(0);
    // legendContour->SetTextSize(0.04);

    // // 3. Модифицированный цикл отрисовки
    // for (int part : PARTS) {
    //     for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
    //         int centr = CENTR_SYST[systN][j];
            
    //         // Проверка существования хотя бы одного контура
    //         bool hasAnyContour = false;
    //         for (int nsigma = 1; nsigma < N_SIGMA; nsigma++) {
    //             if (contour[part][centr][nsigma]) {
    //                 hasAnyContour = true;
    //                 break;
    //             }
    //         }
    //         if (!hasAnyContour) continue;

    //         // Добавление в легенду
    //         string legendText = partTitles[part] + ", " + centrTitles[centr];
            
    //         // Используем первый существующий контур для легенды
    //         TGraph* firstValidContour = nullptr;
    //         for (int nsigma = 1; nsigma < N_SIGMA; nsigma++) {
    //             if (contour[part][centr][nsigma]) {
    //                 firstValidContour = contour[part][centr][nsigma];
    //                 break;
    //             }
    //         }
            
    //         if (firstValidContour) {
    //             legendContour->AddEntry(firstValidContour, legendText.c_str(), "l");
    //         }

    //         // Отрисовка контуров с разными стилями
    //         for (int nsigma = 1; nsigma < N_SIGMA; nsigma++) {
    //             if (contour[part][centr][nsigma]) {
    //                 // Настройка стилей для разных сигм
    //                 contour[part][centr][nsigma]->SetLineColor(centrColors[centr]);
    //                 contour[part][centr][nsigma]->SetLineStyle(nsigma);
    //                 contour[part][centr][nsigma]->SetLineWidth(2);
    //                 contour[part][centr][nsigma]->Draw("lf same");
    //             }
    //         }
    //     }
    // }

    // // 4. Добавим заголовок и дополнительную разметку
    // TLatex* header = new TLatex(0.4, 0.95, Form("System: %s", systNamesT[systN].Data()));
    // header->SetNDC();
    // header->SetTextSize(0.045);
    // header->Draw();

    // legendContour->Draw();
    // c3->Update();

    // c3->SaveAs("output/pics/spectra/ALL_BlastWave_contour_" + systNamesT[systN] + ".png");
    // c3->SaveAs("output/pics/spectra/ALL_BlastWave_contour_" + systNamesT[systN] + ".pdf");
    gROOT->ProcessLine(".q");
}

