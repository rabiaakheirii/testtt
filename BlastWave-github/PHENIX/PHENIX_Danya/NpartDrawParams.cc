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

void DrawParam(TString paramName = "T")
{
    TCanvas *c2 = new TCanvas("c2", "c2", 30, 30, 1200, 1000);
    c2->cd();
    c2->SetGrid();
    c2->SetLogx();
    c2->SetLeftMargin(0.18);
    c2->SetRightMargin(0.05);
    c2->SetTopMargin(0.05);
    c2->SetBottomMargin(0.18);

    double ll = 1, rl = 10000., pad_min = (paramName == "T") ? 0.05 : 0., pad_max = (paramName == "T") ? 0.25 : 1., 
            pad_offset_x = 1, pad_offset_y = 1.1, 
            pad_tsize = 0.055, pad_lsize = 0.055;
    // TString pad_title_y = "T [GeV]";
    TString pad_title_x = "N_{part}";
    TString pad_title_y;
    if (paramName == "T") 
    {
        pad_title_y = "T [GeV]";
    } else if (paramName == "beta") 
    {
        pad_title_y = "#beta [GeV]";
    } 

    Format_Pad(ll, rl, pad_min, pad_max, pad_title_x, pad_title_y, pad_offset_x, pad_offset_y, pad_tsize, pad_lsize, "", 8);        

    // Положение легенды в зависимости от параметра
    double leg_x1, leg_y1, leg_x2, leg_y2;
    if (paramName == "T") 
    {
        leg_x1 = 0.60; leg_y1 = 0.60; 
        leg_x2 = 0.90; leg_y2 = 0.90; 
    } else if (paramName == "beta") 
    {
        leg_x1 = 0.60; leg_y1 = 0.24; 
        leg_x2 = 0.90; leg_y2 = 0.54; 
    }

    TLegend *legend = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetNColumns(2);
    legend->SetTextSize(0.05);

    for (int systN: SYSTS)
    {
        for (int charge: {0, 1})
        {
            gr[charge][systN]->Draw("P SAME");
            legend->AddEntry(gr[charge][systN], systNamesTT[systN], "P");
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
        theoryGraph->SetMarkerSize(1.5);
        theoryGraph->SetLineColor(kMagenta);
        theoryGraph->SetLineStyle(2); 
        theoryGraph->Draw("P SAME");
    
        legend->AddEntry(theoryGraph, "AuAu_{Th}", "P");

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

    double avgValue = (paramName == "T") ? gAvgT : gAvgUt;
    
    TLine *avgLine = new TLine(ll, avgValue, rl, avgValue);
    avgLine->SetLineColor(kBlack);
    avgLine->SetLineStyle(9); 
    avgLine->SetLineWidth(2);
    avgLine->Draw("SAME");

    double avgLeg_x1 = 0.2, avgLeg_y1 = 0.21; 
    double avgLeg_x2 = 0.4, avgLeg_y2 = 0.31;

    TLegend *avgLegend = new TLegend(avgLeg_x1, avgLeg_y1, avgLeg_x2, avgLeg_y2);
    avgLegend->SetBorderSize(0);
    avgLegend->SetFillStyle(0);
    avgLegend->SetTextSize(0.05);

    TString avgLabel = (paramName == "T") 
    ? Form("T_{avr} = %.1f MeV", avgValue * 1000) 
    : Form("#beta_{avr} = %.1f MeV", avgValue * 1000);

    avgLegend->AddEntry(avgLine, avgLabel, "L");
    avgLegend->Draw();

    // TString name = "output/pics/BWparamGlobal_" + paramName + ".png";
    c2->SaveAs("output/pics/ALL_BWparamFinal_" + paramName + ".png");
    c2->SaveAs("output/pics/ALL_BWparamFinal_" + paramName + ".pdf"); 

    delete c2;
    delete legend;
    delete avgLine;
    delete avgLegend;
}


void SetGraphs( int systN, TString paramName, TString fitType = "FINAL" )
{
    cout << systNamesTT[systN] << endl;

    // Объявляем xerr ДО использования
    double xerr[MAX_CENTR];
    for (int i: CENTR_SYST[systN]) xerr[i] = 0.;

    // Выбор имени файла в зависимости от типа фитирования
    TString filename;
    if( fitType == "GLOBAL" )
    {
        filename = "output/parameters/ALL_GlobalBWparams_" + systNamesT[systN] + ".txt";
    }
    else if( fitType == "FINAL" )
    {
        filename = "output/parameters/ALL_FinalBWparams_" + systNamesT[systN] + ".txt";
        // filename = "output/parameters/ALL_FinalBWparams_AuAu.txt";
        // filename = "output/parameters/ALL_FinalBWparams_CuAu.txt";
    }
    cout << "Reading file: " << filename << endl;
    
    ReadGlobalParams(systN, paramsGlobal, filename);

    cout << filename << "  " <<  N_CENTR_SYST[systN] << endl;
    
    // Заполнение массивов Tpar и utPar в зависимости от типа
    for (int charge: {0, 1})
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
                N_CENTR_SYST[systN], 
                Npart[systN], 
                Tpar[charge], 
                xerr, 
                T_err
            );
        }
        else if(paramName == "beta") {
            gr[charge][systN] = new TGraphErrors(
                N_CENTR_SYST[systN], 
                Npart[systN], 
                utPar[charge], 
                xerr, 
                beta_err
            );
        }

        // Настройки стиля
        gr[charge][systN]->SetMarkerStyle(markerStyles[charge][systN]);
        gr[charge][systN]->SetMarkerColor(chargeSystColors[charge][systN]);
        gr[charge][systN]->SetLineWidth(2);
        gr[charge][systN]->SetMarkerSize(2.3);
    }
    cout << "DEFINE GRAPHS " << N_CENTR_SYST[systN] << endl;

    for (int charge: {0, 1})
    {
        if (paramName == "T")
            gr[charge][systN] = new TGraphErrors(N_CENTR_SYST[systN], Npart[systN], Tpar[charge], xerr, xerr); // Tpar_err[part]);
        else if (paramName == "beta")
            gr[charge][systN] = new TGraphErrors(N_CENTR_SYST[systN], Npart[systN], utPar[charge], xerr, xerr); // utPar_err[part]);
        
            gr[charge][systN]->SetMarkerStyle(markerStyles[charge][systN]);
            gr[charge][systN]->SetMarkerColor(chargeSystColors[charge][systN]);
            gr[charge][systN]->SetLineWidth(5);
            gr[charge][systN]->SetMarkerSize(3);
    }
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
    
    TLegend *legend = new TLegend(0.55, 0.60, 0.95, 0.90);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetNColumns(2);
    legend->SetTextSize(0.05);
    
    
    // Считываем данные из файлов для каждого systN и charge
    for (int systN : SYSTS)
    {
        for (int charge : {0, 1})
        {
            // Считываем параметры с ошибками
            double T_err[MAX_CENTR], beta_err[MAX_CENTR];

            // Выбор имени файла в зависимости от типа фитирования
            TString filename;
            if (fitType == "GLOBAL")
            {
                filename = "output/parameters/ALL_GlobalBWparams_" + systNamesT[systN] + ".txt";
            }
            else if (fitType == "FINAL")
            {
                filename = "output/parameters/ALL_FinalBWparams_" + systNamesT[systN] + ".txt";
                // filename = "output/parameters/ALL_FinalBWparams_CuAu.txt";
            }
            
            // Считываем параметры из файла
            cout << "Reading file: " << filename << endl;
            ReadGlobalParams(systN, paramsGlobal, filename);
            cout << "File read successfully." << endl;

            
            // Заполняем массивы Tpar и utPar
            for (int centr = 0; centr < N_CENTR_SYST[systN]; centr++)
            {
                if (fitType == "GLOBAL")
                {
                    Tpar[charge][centr] = paramsGlobal[charge][centr][0];
                    utPar[charge][centr] = paramsGlobal[charge][centr][1];
                }
                else if (fitType == "FINAL")
                {
                    Tpar[charge][centr] = paramsGlobal[charge][centr][1];
                    T_err[centr] = paramsGlobal[charge][centr][2];    // Ошибка T (индекс 2)
                    utPar[charge][centr] = paramsGlobal[charge][centr][3];
                    beta_err[centr] = paramsGlobal[charge][centr][4]; // Ошибка beta (индекс 4)
                }
            }
            
            // Создаем график T vs beta
            TGraphErrors *gr_TvsUt = new TGraphErrors(N_CENTR_SYST[systN], utPar[charge], Tpar[charge], 0, 0);
            gr_TvsUt->SetMarkerStyle(markerStyles[charge][systN]);
            gr_TvsUt->SetMarkerSize(2.5);
            gr_TvsUt->SetMarkerColor(chargeSystColors[charge][systN]);
            gr_TvsUt->Draw("P SAME");
            legend->AddEntry(gr_TvsUt, systNamesTT[systN], "P");
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
    legend->AddEntry(articleAuAu, "AuAu_{Th}", "P");
    
    // Добавляем точки STAR для AuAu
    TGraph *starAuAu = new TGraph(9, beta_AuAu_STAR, T_AuAu_STAR);
    starAuAu->SetMarkerStyle(20);
    starAuAu->SetMarkerColor(kCyan);
    starAuAu->SetMarkerSize(2.5);
    starAuAu->Draw("P SAME");
    legend->AddEntry(starAuAu, "AuAu_{STAR}", "P");
    
    // Добавляем точки STAR для UU
    TGraph *starUU = new TGraph(9, beta_UU_STAR, T_UU_STAR);
    starUU->SetMarkerStyle(21);
    starUU->SetMarkerColor(kTeal);
    starUU->SetMarkerSize(2.5);
    starUU->Draw("P SAME");
    legend->AddEntry(starUU, "UU_{STAR}", "P");
    
    legend->Draw();
    
        // Добавляем линии средних значений
    // Горизонтальная линия для среднего T
    TLine *avgLineT = new TLine(pad_min_x, gAvgT, pad_max_x, gAvgT);
    avgLineT->SetLineColor(kBlack);
    avgLineT->SetLineStyle(9); 
    avgLineT->SetLineWidth(2);
    avgLineT->Draw("SAME");

    // Вертикальная линия для среднего beta
    TLine *avgLineBeta = new TLine(gAvgUt, pad_min_y, gAvgUt, pad_max_y);
    avgLineBeta->SetLineColor(kBlack);
    avgLineBeta->SetLineStyle(9);
    avgLineBeta->SetLineWidth(2);
    avgLineBeta->Draw("SAME");

    // Обновляем легенду
    TLegend *avgLegend = new TLegend(0.2, 0.22, 0.41, 0.32);
    avgLegend->SetBorderSize(0);
    avgLegend->SetFillStyle(0);
    avgLegend->SetTextSize(0.045);
    avgLegend->AddEntry(avgLineT, Form("T_{avr} = %.1f MeV", gAvgT * 1000), "L");
    avgLegend->AddEntry(avgLineBeta, Form("#beta_{avr} = %.1f MeV", gAvgUt * 1000), "L");
    avgLegend->Draw();

    // Сохраняем график
    c3->SaveAs("output/pics/ALL_BWparamFinal_T_beta.png");
    c3->SaveAs("output/pics/ALL_BWparamFinal_T_beta.pdf");

    delete c3;
    delete legend;
}


void NpartDrawParams (void) 
{
    // Для параметра T
    for (int systN: SYSTS) 
    {
        cout << systN << " " << systNamesTT[systN] << endl;
        SetGraphs(systN, "T");
    }
    CalculateAverage("T"); 
    DrawParam("T");
    
    // Для параметра beta
    for (int systN: SYSTS) 
    {
        cout << systN << " " << systNamesTT[systN] << endl;
        SetGraphs(systN, "beta");
    }
    CalculateAverage("beta"); 
    DrawParam("beta");

    WriteAveragesToFile();
    DrawTbeta();
} 

int main() {
    NpartDrawParams();
}
