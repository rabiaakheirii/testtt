#include "input/headers/def.h"
#include "input/headers/WriteReadFiles.h"
#include <cstdlib>


// Определяем двумерный массив цветов: первая ось – заряд, вторая – система
Color_t systColors[6] = {kRed, kBlue, kGreen+2, kMagenta, kOrange+1, kViolet-6};
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

const double CENTR_MEAN[5][MAX_CENTR] = {
    // AuAu (12 элементов)
    {-1000, 2.5, 7.5, 12.5, 17.5, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 86.0},
    // pAl (4 элемента)
    {-1000, 10.0, 30.0, 56.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    // HeAu (5 элементов)
    {-1000, 10.0, 30.0, 50.0, 74.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    // CuAu (5 элементов)
    {-1000, 10.0, 30.0, 50.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    // UU (4 элемента)
    {-1000, 30.0, 50.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
};


int systMarkerStyles[5] = {20, 21, 22, 23, 29};

const int NSYST = 5; // Количество систем: AuAu, pAl, HeAu, CuAu, UU, и т.д.
const int NPART = 6;
double T_syst_err[NSYST][NPART][MAX_CENTR] = {0};
double beta_syst_err[NSYST][NPART][MAX_CENTR] = {0};

void DrawParam(TString paramName = "T_{kin}")
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
           pad_min = (paramName == "T_{kin}") ? 0.05 : 0., 
           pad_max = (paramName == "T_{kin}") ? 0.25 : 1., 
           pad_offset_x = 1, pad_offset_y = 1.1, 
           pad_tsize = 0.055, pad_lsize = 0.055;

    TString pad_title_x = "N_{part}";
    TString pad_title_y;
    if (paramName == "T_{kin}") { 
        pad_title_y = "T [GeV]"; 
    } else if (paramName == "beta") {
        pad_title_y = "#beta [GeV]";
    } 

    Format_Pad(ll, rl, pad_min, pad_max, pad_title_x, pad_title_y, pad_offset_x, pad_offset_y, pad_tsize, pad_lsize, "", 8);        

    // Положение легенды в зависимости от параметра
    double leg_x1, leg_y1, leg_x2, leg_y2, 
           leg_x11, leg_y11, leg_x22, leg_y22;
    if (paramName == "T_{kin}") {
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
                
                double ey_syst = (paramName == "T_{kin}") ? T_syst_err[systN][charge][centr_index] 
                                                    : beta_syst_err[systN][charge][centr_index]; 
                
                // Отладочный вывод 
                if (paramName == "T_{kin}") {
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

    if (paramName == "T_{kin}" || paramName == "beta") 
    {
        const int nPoints = 13; 
    
        double yValues[nPoints];
    
        // Добавляем точки ARTICLE для AuAu
        if (paramName == "T_{kin}") {
            for (int i = 0; i < nPoints; i++) yValues[i] = T_AuAu_ART[i];
        } else if (paramName == "beta") {
            for (int i = 0; i < nPoints; i++) yValues[i] = beta_AuAu_ART[i];
        }
    
        TGraph *theoryGraph = new TGraph(nPoints, Npart[0], yValues);
        theoryGraph->SetMarkerStyle(43); 
        theoryGraph->SetMarkerColor(kMagenta);
        theoryGraph->SetMarkerSize(3);
        theoryGraph->SetLineColor(kMagenta);
        theoryGraph->SetLineStyle(2); 
        theoryGraph->Draw("P SAME");
    
        legend1->AddEntry(theoryGraph, "Au+Au, #sqrt{s_{NN}} = 200 GeV (C72:014903, 2005)", "P");
    }
    legend->Draw();
    legend1->Draw();
    legend2->Draw();
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
                            
            // // Ограничение систематики для UU (systN=4) и π+ (part=0)
            // if(systN == 4 && part == 0) {  // UU и π+
            //     // Применяем коэффициент 0.5 к систематическим ошибкам
            //     T_syst *= 0.3;
            //     beta_syst *= 0.3;
            // }
            
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
    cout << "СЧИТАЛИ ФАЙЛ: " << systFilename << "\n" << endl;

    // Заполнение массивов Tpar и utPar в зависимости от типа
    for (int part = 0; part < NPART; part++) // БЫЛО: {0, 1}
    {
        // Добавляем массивы для ошибок
        double T_err[MAX_CENTR], beta_err[MAX_CENTR];

        for (int centr = 0; centr < N_CENTR_SYST[systN]; centr++)
        {
            if( fitType == "GLOBAL" )
            {
                Tpar[part][centr] = paramsGlobal[part][centr][1];  // T value
                T_err[centr] = paramsGlobal[part][centr][2];         // T error (индекс 2)
                utPar[part][centr] = paramsGlobal[part][centr][3];  // beta value
                beta_err[centr] = paramsGlobal[part][centr][4];       // beta error (индекс 4)
            }
            else if( fitType == "FINAL" )
            {
                Tpar[part][centr] = paramsGlobal[part][centr][1];   // T value
                T_err[centr] = paramsGlobal[part][centr][2];          // T error
                utPar[part][centr] = paramsGlobal[part][centr][3]; // beta value
                beta_err[centr] = paramsGlobal[part][centr][4];      // beta error
            }
        }

        // Создаем графики ВНУТРИ цикла по charge
        if(paramName == "T_{kin}") {
            gr[part][systN] = new TGraphErrors(
                N_CENTR_SYST[systN], Npart[systN], 
                Tpar[part], xerr, T_err
            );
        } else if(paramName == "beta") {
            gr[part][systN] = new TGraphErrors(
                N_CENTR_SYST[systN], Npart[systN], 
                utPar[part], xerr, beta_err
            );
        }

        gr[part][systN]->SetMarkerStyle(markerStyles[part][systN]);
        gr[part][systN]->SetMarkerColor(chargeSystColors[part][systN]);
        gr[part][systN]->SetLineColor(chargeSystColors[part][systN]);
        gr[part][systN]->SetLineWidth(2);
        gr[part][systN]->SetMarkerSize(2.3);
    }
}

void DrawParticleParams(int part, TString paramName = "T_{kin}") {
    // part: 0 - pi+, 1 - pi-, 2 - K+, 3 - K-, 4 - p, 5 - antip
    TString partNames[6] = {"#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}", "p", "#bar{p}"};
    
    TCanvas *c = new TCanvas(Form("c_%s_%d", paramName.Data(), part), "", 30, 30, 1440, 1440);
    c->SetGrid();
    c->SetLeftMargin(0.18);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.05);
    c->SetBottomMargin(0.18);

    // Настройка осей
    TH1F *frame = c->DrawFrame(0, 0, 100, (paramName == "T_{kin}") ? 0.3 : 1.2);
    frame->SetTitle(Form(";Centrality, %%;%s", (paramName == "T_{kin}") ? "T [GeV]" : "#beta [GeV]"));
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->GetXaxis()->SetNdivisions(505);

    TLegend *leg = new TLegend(0.7, 0.65, 0.95, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetHeader(Form("%s parameters", partNames[part].Data()));

    for (int systN : SYSTS) {
        
        // Чтение параметров
        TString filename = Form("output/parameters/ALL_FinalBWparams_%s.txt", systNamesT[systN].Data());
        ReadGlobalParams(systN, paramsGlobal, filename);

        // Подготовка данных
        const int nPoints = N_CENTR_SYST[systN];
        double x[nPoints], y[nPoints], ex[nPoints], ey_stat[nPoints], ey_syst[nPoints];

        for (int i = 0; i < nPoints; ++i) {
            int centrIndex = CENTR_SYST[systN][i];
            x[i] = CENTR_MEAN[systN][i]; // Среднее значение центральности
            ex[i] = 2.0; // Ширина бина центральности
            
            if (paramName == "T_{kin}") {
                y[i] = paramsGlobal[part][centrIndex][1];
                ey_stat[i] = paramsGlobal[part][centrIndex][2];
                ey_syst[i] = T_syst_err[systN][part][centrIndex];
            } else {
                y[i] = paramsGlobal[part][centrIndex][3];
                ey_stat[i] = paramsGlobal[part][centrIndex][4];
                ey_syst[i] = beta_syst_err[systN][part][centrIndex];
            }
        }

        // График со статистическими ошибками
        TGraphErrors *gr_stat = new TGraphErrors(nPoints, x, y, ex, ey_stat);
        gr_stat->SetMarkerStyle(systMarkerStyles[systN]);
        gr_stat->SetMarkerColor(systColors[systN]);
        gr_stat->SetLineColor(systColors[systN]);
        gr_stat->SetLineWidth(2);
        gr_stat->SetMarkerSize(2.0);
        gr_stat->Draw("P same");

        // Систематические ошибки (TBox)
        for (int i = 0; i < nPoints; ++i) {
            TBox *box = new TBox(x[i]-ex[i], y[i]-ey_syst[i], x[i]+ex[i], y[i]+ey_syst[i]);
            box->SetFillColorAlpha(systColors[systN], 0.3);
            box->SetLineColor(systColors[systN]);
            box->SetLineWidth(1);
            box->Draw("same");
        }

        leg->AddEntry(gr_stat, systNamesTT[systN], "P");
    }

    leg->Draw();
    
    // Сохранение графиков
    TString dir = Form("output/pics/%s/", partNames[part].Data());
    gSystem->Exec(Form("mkdir -p %s", dir.Data()));
    c->SaveAs(Form("%s/%s_vs_cent_%s.png", dir.Data(), paramName.Data(), partNames[part].Data()));
    c->SaveAs(Form("%s/%s_vs_cent_%s.pdf", dir.Data(), paramName.Data(), partNames[part].Data()));
}

void DrawCombinedParticleParamsCentr(TString paramName = "T_{kin}") {
    TString partNames[6] = {"#pi^{+}", "#pi^{#minus}", "K^{+}", "K^{#minus}", "p", "#bar{p}"};
    
    TCanvas *c2 = new TCanvas("c2", "c2", 1440, 2160);
    // c2->Divide(2, 3, 0, 0);
    Format_Canvas(c2, 2, 3, 0);

    // Цвета и стили для систем
    Color_t systColors[5] = {kRed, kBlue, kGreen+2, kMagenta, kOrange+1};
    int systMarkerStyles[5] = {20, 21, 22, 23, 29};

    for (int part = 0; part < 6; ++part) {
        c2->cd(part+1);
        gPad->SetGrid();

        // Создаем фрейм
        TH1F *frame = gPad->DrawFrame(-5, (paramName == "T_{kin}") ? -0.01 : -0.1, 
                                   108, (paramName == "T_{kin}") ? 0.36 : 1.3);
        frame->SetTitle(Form("%s;Centrality, %% ;%s", 
                           " ",
                           (paramName == "T_{kin}") ? "T_{kin} [GeV]" : "#beta_{T} [GeV]"));
        
        if (part == 0 || part == 2) {
            TLatex *tex = new TLatex();
            tex->SetNDC(); // Используем нормализованные координаты
            tex->SetTextFont(42); // Полужирный шрифт
            tex->SetTextSize(0.12); // Размер текста
            tex->SetTextAlign(11); // Выравнивание: левый верхний угол
            tex->DrawLatex(0.35, 0.81, partNames[part]); // (x, y) в NDC координатах
        }

        if (part == 1 || part == 3) {
            TLatex *tex = new TLatex();
            tex->SetNDC(); // Используем нормализованные координаты
            tex->SetTextFont(42); // Полужирный шрифт
            tex->SetTextSize(0.12); // Размер текста
            tex->SetTextAlign(11); // Выравнивание: левый верхний угол
            tex->DrawLatex(0.15, 0.81, partNames[part]); // (x, y) в NDC координатах
        }

        if (part == 4) {
            TLatex *tex = new TLatex();
            tex->SetNDC(); // Используем нормализованные координаты
            tex->SetTextFont(42); // Полужирный шрифт
            tex->SetTextSize(0.12); // Размер текста
            tex->SetTextAlign(11); // Выравнивание: левый верхний угол
            tex->DrawLatex(0.35, 0.85, partNames[part]); // (x, y) в NDC координатах
        }

        if (part == 5) {
            TLatex *tex = new TLatex();
            tex->SetNDC(); // Используем нормализованные координаты
            tex->SetTextFont(42); // Полужирный шрифт
            tex->SetTextSize(0.12); // Размер текста
            tex->SetTextAlign(11); // Выравнивание: левый верхний угол
            tex->DrawLatex(0.15, 0.85, partNames[part]); // (x, y) в NDC координатах
        }

        frame->GetXaxis()->SetTitleSize(0.09);
        frame->GetXaxis()->SetLabelSize(0.07);
        frame->GetXaxis()->SetTitleOffset(1.1);
        frame->GetYaxis()->SetTitleSize(0.09);
        frame->GetYaxis()->SetLabelSize(0.07);
        frame->GetYaxis()->SetTitleOffset(1.2);

        // // Добавляем отладочный вывод
        // cout << "\n=== Drawing " << partNames[part] << " ===" << endl;

        for (int systN : SYSTS) {
            // !!! КРИТИЧЕСКИ ВАЖНО: читаем параметры для каждой системы !!!
            TString filename = Form("output/parameters/ALL_FinalBWparams_%s.txt", 
                                   systNamesT[systN].Data());
            ReadGlobalParams(systN, paramsGlobal, filename); 

            const int nPoints = N_CENTR_SYST[systN];
            double x[nPoints], y[nPoints], ex[nPoints], ey_stat[nPoints], ey_syst[nPoints];

            for (int i = 0; i < nPoints; ++i) {
                int centrIndex = CENTR_SYST[systN][i];
                x[i] = CENTR_MEAN[systN][i];
                ex[i] = 2.0;

                // Отладочный вывод
                // cout << "System: " << systNamesTT[systN] 
                //      << " Centr: " << centrIndex 
                //      << " X: " << x[i];
                
                if (paramName == "T_{kin}") {
                    y[i] = paramsGlobal[part][centrIndex][1];
                    ey_stat[i] = paramsGlobal[part][centrIndex][2];
                    ey_syst[i] = T_syst_err[systN][part][centrIndex];
                    // cout << " T: " << y[i] << " ± " << ey_stat[i] 
                    //      << " (syst: " << ey_syst[i] << ")" << endl;
                } else {
                    y[i] = paramsGlobal[part][centrIndex][3];
                    ey_stat[i] = paramsGlobal[part][centrIndex][4];
                    ey_syst[i] = beta_syst_err[systN][part][centrIndex];
                    // cout << " beta: " << y[i] << " ± " << ey_stat[i] 
                    //      << " (syst: " << ey_syst[i] << ")" << endl;
                }
            }

            // График со статистическими ошибками
            TGraphErrors *gr = new TGraphErrors(nPoints, x, y, ex, ey_stat);
            gr->SetMarkerStyle(systMarkerStyles[systN]);
            gr->SetMarkerSize(3);
            gr->SetMarkerColor(systColors[systN]);
            gr->SetLineColor(systColors[systN]);
            gr->Draw("P same");

            // Систематические ошибки
            for (int i = 0; i < nPoints; ++i) {
                TBox *box = new TBox(x[i]-ex[i], y[i]-ey_syst[i], x[i]+ex[i], y[i]+ey_syst[i]);
                box->SetFillColorAlpha(systColors[systN], 0.3);
                box->SetLineColor(systColors[systN]);
                box->Draw("same");
            }
        }

        // Добавляем легенду для первого pad'а
        if (part == 1) {
            TLegend *leg = new TLegend(0.54, 0.55, 0.93, 0.90);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.09);
            for (int systN : SYSTS) {
                // Создаем временный график для легенды
                TGraph *dummy = new TGraph();
                dummy->SetMarkerStyle(systMarkerStyles[systN]);
                dummy->SetMarkerColor(systColors[systN]);
                dummy->SetMarkerSize(3);  // Устанавливаем размер маркера
                
                // Добавляем запись и настраиваем её
                TLegendEntry *entry = leg->AddEntry(dummy, systNamesTT[systN], "P");
                entry->SetMarkerSize(3);  // Явно задаем размер маркера в легенде
            }
            leg->Draw();
        }

        c2->cd(0); // Переключаемся на основной холст (pad 0)
        TLatex* texEnergy = new TLatex();
        texEnergy->SetNDC(kTRUE); // Используем нормализованные координаты
        texEnergy->SetTextSize(0.03); // Размер текста
        texEnergy->SetTextFont(42);  // Обычный шрифт (не жирный)
        texEnergy->DrawLatex(0.15, 0.125, "p+Al, He+Au, Cu+Au: #sqrt{s_{NN}} = 200 GeV");
        texEnergy->DrawLatex(0.15, 0.10, "U+U: #sqrt{s_{NN}} = 193 GeV");
    }

    TString dir = "output/pics/centr/";
    gSystem->Exec("mkdir -p " + dir);
    c2->SaveAs(dir + Form("centr_%s.png", paramName.Data()));
    c2->SaveAs(dir + Form("centr_%s.pdf", paramName.Data()));
}

void DrawCombinedParticleParamsCentrNpart(TString paramName = "T_{kin}") {
    TString partNames[6] = {"#pi^{+}", "#pi^{#minus}", "K^{+}", "K^{#minus}", "p", "#bar{p}"};
    
    TCanvas *c2 = new TCanvas("c2", "c2", 1440, 2160);
    // c->Divide(2, 3, 0, 0);
    Format_Canvas(c2, 2, 3, 0);

    // Цвета и стили для систем
    Color_t systColors[5] = {kRed, kBlue, kGreen+2, kMagenta, kOrange+1};
    int systMarkerStyles[5] = {20, 21, 22, 23, 29};
    

    for (int part = 0; part < 6; ++part) {
        c2->cd(part+1);
        gPad->SetGrid();
        gPad->SetLogx(); // Логарифмическая шкала для Npart

        // Определяем диапазоны осей
        double xmin = 1, xmax = 2000; // Npart диапазон
        double ymin = (paramName == "T_{kin}") ? -0.02 : 0.005;
        double ymax = (paramName == "T_{kin}") ? 0.27 : 1.25;

        // Создаем фрейм
        TH1F *frame = gPad->DrawFrame(xmin, ymin, xmax, ymax);
        frame->SetTitle(Form("%s;N_{part};%s", 
                           " ",
                           (paramName == "T_{kin}") ? "T_{kin} [GeV]" : "#beta_{T} [GeV]"));
        
        if (part == 0 || part == 2) {
            TLatex *tex = new TLatex();
            tex->SetNDC(); // Используем нормализованные координаты
            tex->SetTextFont(42); // Полужирный шрифт
            tex->SetTextSize(0.12); // Размер текста
            tex->SetTextAlign(11); // Выравнивание: левый верхний угол
            tex->DrawLatex(0.35, 0.81, partNames[part]); // (x, y) в NDC координатах
        }

        if (part == 1 || part == 3) {
            TLatex *tex = new TLatex();
            tex->SetNDC(); // Используем нормализованные координаты
            tex->SetTextFont(42); // Полужирный шрифт
            tex->SetTextSize(0.12); // Размер текста
            tex->SetTextAlign(11); // Выравнивание: левый верхний угол
            tex->DrawLatex(0.15, 0.81, partNames[part]); // (x, y) в NDC координатах
        }

        if (part == 4) {
            TLatex *tex = new TLatex();
            tex->SetNDC(); // Используем нормализованные координаты
            tex->SetTextFont(42); // Полужирный шрифт
            tex->SetTextSize(0.12); // Размер текста
            tex->SetTextAlign(11); // Выравнивание: левый верхний угол
            tex->DrawLatex(0.35, 0.85, partNames[part]); // (x, y) в NDC координатах
        }

        if (part == 5) {
            TLatex *tex = new TLatex();
            tex->SetNDC(); // Используем нормализованные координаты
            tex->SetTextFont(42); // Полужирный шрифт
            tex->SetTextSize(0.12); // Размер текста
            tex->SetTextAlign(11); // Выравнивание: левый верхний угол
            tex->DrawLatex(0.15, 0.85, partNames[part]); // (x, y) в NDC координатах
        }
        
        // Настройка осей
        frame->GetXaxis()->SetTitleSize(0.09);
        frame->GetXaxis()->SetLabelSize(0.07);
        frame->GetXaxis()->SetTitleOffset(1.1);
        frame->GetYaxis()->SetTitleSize(0.09);
        frame->GetYaxis()->SetLabelSize(0.07);
        frame->GetYaxis()->SetTitleOffset(1.2);

        // Отладочный вывод
        cout << "\n=== Drawing " << partNames[part] << " vs Npart ===" << endl;

        for (int systN : SYSTS) {
            TString filename = Form("output/parameters/ALL_FinalBWparams_%s.txt", 
                                   systNamesT[systN].Data());
            ReadGlobalParams(systN, paramsGlobal, filename); 

            const int nPoints = N_CENTR_SYST[systN];
            double x[nPoints], y[nPoints], ex[nPoints], ey_stat[nPoints], ey_syst[nPoints];

            for (int i = 0; i < nPoints; ++i) {
                int centrIndex = CENTR_SYST[systN][i];
                
                // Берем Npart из массива
                x[i] = Npart[systN][centrIndex];
                ex[i] = 0; // Ошибка по X не используется

                // Рассчитываем ошибку X как 10% от Npart
                if (x[i] > 0) ex[i] = x[i] * 0.10;

                if (paramName == "T_{kin}") {
                    y[i] = paramsGlobal[part][centrIndex][1];
                    ey_stat[i] = paramsGlobal[part][centrIndex][2];
                    ey_syst[i] = T_syst_err[systN][part][centrIndex];
                } else {
                    y[i] = paramsGlobal[part][centrIndex][3];
                    ey_stat[i] = paramsGlobal[part][centrIndex][4];
                    ey_syst[i] = beta_syst_err[systN][part][centrIndex];
                }

                cout << "System: " << systNamesTT[systN] 
                     << " Npart: " << x[i] 
                     << " Param: " << y[i] << endl;
            }

            // График со статистическими ошибками
            TGraphErrors *gr = new TGraphErrors(nPoints, x, y, ex, ey_stat);
            gr->SetMarkerStyle(systMarkerStyles[systN]);
            gr->SetMarkerSize(3);
            gr->SetMarkerColor(systColors[systN]);
            gr->SetLineColor(systColors[systN]);
            gr->Draw("P same");

            // Систематические ошибки (TBox)
            for (int i = 0; i < nPoints; ++i) {
                TBox *box = new TBox(x[i]-ex[i], y[i]-ey_syst[i], 
                                    x[i]+ex[i], y[i]+ey_syst[i]);
                box->SetFillColorAlpha(systColors[systN], 0.3);
                box->SetLineColor(systColors[systN]);
                box->Draw("same");
            }
        }

        // Легенда для первого pad'а
        if (part == 1) {
            TLegend *leg = new TLegend(0.54, 0.55, 0.93, 0.90);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.09);
            for (int systN : SYSTS) {
                TGraph *dummy = new TGraph();
                dummy->SetMarkerStyle(systMarkerStyles[systN]);
                dummy->SetMarkerColor(systColors[systN]);
                dummy->SetMarkerSize(3);
                leg->AddEntry(dummy, systNamesTT[systN], "P");
            }
            leg->Draw();
        }

        c2->cd(0); // Переключаемся на основной холст (pad 0)
        TLatex* texEnergy = new TLatex();
        texEnergy->SetNDC(kTRUE); // Используем нормализованные координаты
        texEnergy->SetTextSize(0.03); // Размер текста
        texEnergy->SetTextFont(42);  // Обычный шрифт (не жирный)
        texEnergy->DrawLatex(0.15, 0.125, "p+Al, He+Au, Cu+Au: #sqrt{s_{NN}} = 200 GeV");
        texEnergy->DrawLatex(0.15, 0.10, "U+U: #sqrt{s_{NN}} = 193 GeV");
    }

    // Сохранение
    TString dir = "output/pics/Npart/";
    gSystem->Exec("mkdir -p " + dir);
    c2->SaveAs(dir + Form("%s_vs_Npart.png", paramName.Data()));
    c2->SaveAs(dir + Form("%s_vs_Npart.pdf", paramName.Data()));
}

void DrawAllParticleParams() {
    DrawCombinedParticleParamsCentr("T_{kin}");
    DrawCombinedParticleParamsCentr("beta_{T}");

    DrawCombinedParticleParamsCentrNpart("T_{kin}");
    DrawCombinedParticleParamsCentrNpart("beta_{T}");
}

void DrawAveragedParamsCentr(TString paramName = "T_{kin}") {
    TCanvas *c = new TCanvas("c2","с2", 1200, 1000);
    c->SetGrid();
    c->SetLeftMargin(0.17);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.05);
    c->SetBottomMargin(0.15);

    // Настройка осей
    TH1F *frame = c->DrawFrame(0, (paramName == "T_{kin}") ? 0.0 : 0.0, 
                             100, (paramName == "T_{kin}") ? 0.3 : 1.2);
    frame->SetTitle(Form(";Centrality, %% ;%s [GeV]", paramName.Data()));
    frame->GetXaxis()->SetTitleSize(0.06);
    frame->GetYaxis()->SetTitleSize(0.06);
    frame->GetXaxis()->SetLabelSize(0.05);
    frame->GetYaxis()->SetLabelSize(0.05);

    TLegend *legend = new TLegend(0.71, 0.70, 0.91, 0.89);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.05);

    // Добавляем новую легенду с энергиями
    TLegend *legEnergy = new TLegend(0.03, 0.2, 0.8, 0.32);
    legEnergy->SetBorderSize(0);
    legEnergy->SetFillStyle(0);
    legEnergy->SetTextSize(0.05);
    legEnergy->AddEntry((TObject*)0, "p+Al, He+Au, Cu+Au: #sqrt{s_{NN}} = 200 GeV", "");
    legEnergy->AddEntry((TObject*)0, "U+U: #sqrt{s_{NN}} = 193 GeV", "");

    // Цвета и стили для систем
    Color_t systColors[5] = {kRed, kBlue, kGreen+2, kMagenta, kOrange+1};
    int systMarkerStyles[5] = {20, 21, 22, 23, 29};

    for (int systN : SYSTS) {
        // Пропускаем UU (systN=4), если нужно, или оставляем все
        TString filename = Form("output/parameters/ALL_FinalBWparams_%s.txt", 
                              systNamesT[systN].Data());
        ReadGlobalParams(systN, paramsGlobal, filename);

        vector<double> x, y, ex, ey_stat, ey_syst;

        for (int i = 0; i < N_CENTR_SYST[systN]; ++i) {
            double centrVal = CENTR_MEAN[systN][i];
            // Пропускаем невалидные значения центральности
            if (centrVal < 0) continue;

            double sumParam = 0, sumStat2 = 0, sumSyst2 = 0;
            int count = 0;

            for (int part = 0; part < 6; ++part) {
                double val, stat_err, syst_err;
                if (paramName == "T_{kin}") {
                    val = paramsGlobal[part][i][1];
                    stat_err = paramsGlobal[part][i][2];
                    syst_err = T_syst_err[systN][part][i];
                } else {
                    val = paramsGlobal[part][i][3];
                    stat_err = paramsGlobal[part][i][4];
                    syst_err = beta_syst_err[systN][part][i];
                }

                sumParam += val;
                sumStat2 += stat_err * stat_err;
                sumSyst2 += syst_err * syst_err;
                count++;
            }

            if (count == 0) continue;

            double mean = sumParam / count;
            double meanStatErr = sqrt(sumStat2) / count;
            double meanSystErr = sqrt(sumSyst2) / count;

            x.push_back(centrVal);
            y.push_back(mean);
            ex.push_back(2.0); // Ширина бина центральности
            ey_stat.push_back(meanStatErr);
            ey_syst.push_back(meanSystErr);
        }

        // Создание графика
        TGraphErrors *gr = new TGraphErrors(x.size(), x.data(), y.data(), 
                                           ex.data(), ey_stat.data());
        gr->SetMarkerStyle(systMarkerStyles[systN]);
        gr->SetMarkerColor(systColors[systN]);
        gr->SetMarkerSize(3);
        gr->SetLineColor(systColors[systN]);
        gr->SetLineWidth(2);
        gr->Draw("P SAME");

        // Систематические ошибки (TBox)
        for (size_t j = 0; j < x.size(); ++j) {
            TBox *box = new TBox(x[j] - ex[j], y[j] - ey_syst[j],
                                x[j] + ex[j], y[j] + ey_syst[j]);
            box->SetFillColorAlpha(systColors[systN], 0.3);
            box->SetLineColor(systColors[systN]);
            box->SetLineWidth(1);
            box->Draw("SAME");
        }

        legend->AddEntry(gr, systNamesTT[systN], "P");
    }

    legend->Draw();
    legEnergy->Draw();
    c->SaveAs(Form("output/pics/centr/Averaged_%s_vs_Centrality.png", paramName.Data()));
    c->SaveAs(Form("output/pics/centr/Averaged_%s_vs_Centrality.pdf", paramName.Data()));
}

// Функция для усреднения параметров по всем частицам и рисования vs Npart
void DrawAveragedParamsNpart(TString paramName = "T_{kin}") {
    TCanvas *c = new TCanvas("cAveragedNpart", "Averaged Parameters vs Npart", 1200, 1000);
    c->SetGrid();
    c->SetLogx();
    c->SetLeftMargin(0.17);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.05);
    c->SetBottomMargin(0.15);

    // Цвета и стили для систем
    Color_t systColors[5] = {kRed, kBlue, kGreen+2, kMagenta, kOrange+1};
    int systMarkerStyles[5] = {20, 21, 22, 23, 29};

    // Настройка осей
    TH1F *frame = c->DrawFrame(1, (paramName == "T_{kin}") ? 0.02: -0.15, 2000, 
                              (paramName == "T_{kin}") ? 0.22 : 1.2);
    frame->SetTitle(Form(";N_{part}; %s [GeV]", paramName.Data()));
    frame->GetXaxis()->SetTitleSize(0.06);
    frame->GetYaxis()->SetTitleSize(0.06);
    frame->GetXaxis()->SetLabelSize(0.05);
    frame->GetYaxis()->SetLabelSize(0.05);

    TLegend *leg = new TLegend(0.75, 0.75, 0.95, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.042);

    TLegend *leg1 = new TLegend(0.18, 0.17, 0.4, 0.28);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.042);

    // Объявляем данные из статьи (должны быть определены в def.h)
    const int nPointsART = 13;
    double Npart_AuAu_ART[nPointsART] = {
        109.1, 351.4, 299.0, 253.9, 215.3, 166.6, 
        114.2, 74.4, 45.5, 25.7, 13.4, 6.3 
    };

    // Добавляем точки из статьи для AuAu
    TGraph* grART = nullptr;
    if (paramName == "T_{kin}") {
        grART = new TGraph(nPointsART, Npart_AuAu_ART, T_AuAu_ART);
        grART->SetMarkerStyle(43); 
        grART->SetMarkerColor(kBlack);
        grART->SetMarkerSize(3.5);
        grART->SetLineColor(kBlack);
        grART->SetLineStyle(2);
        grART->Draw("P SAME");
    }
    else if (paramName == "#beta_{T}") {
        grART = new TGraph(nPointsART, Npart_AuAu_ART, beta_AuAu_ART);
        grART->SetMarkerStyle(43);
        grART->SetMarkerColor(kBlack);
        grART->SetMarkerSize(3.5);
        grART->SetLineColor(kBlack);
        grART->SetLineStyle(2);
        grART->Draw("P SAME");
    }

    const int nPoints_STAR = 9;
    double Npart_AuAu_STAR_ch[nPoints_STAR] = {
        324.84, 237.99, 169.90, 118.09, 64.31, 24.34, 15.0, 8.0, 3.0
    };

    TGraph* grSTAR = nullptr;
    if (paramName == "T_{kin}") {
        grSTAR = new TGraph(nPoints_STAR, Npart_AuAu_STAR_ch, T_AuAu_STAR_ch);
    }
    else if (paramName == "#beta_{T}") {
        grSTAR = new TGraph(nPoints_STAR, Npart_AuAu_STAR_ch, beta_AuAu_STAR_ch);
    }
    
    if (grSTAR) {
        grSTAR->SetMarkerStyle(42);  // Крестики
        grSTAR->SetMarkerColor(kBlack);
        grSTAR->SetMarkerSize(3.5);
        grSTAR->SetLineColor(kBlack);
        grSTAR->SetLineWidth(2);     // Жирные линии
        grSTAR->Draw("P SAME");
    }

    for (int systN : SYSTS) {
        TString filename = Form("output/parameters/ALL_FinalBWparams_%s.txt", 
                               systNamesT[systN].Data());
        ReadGlobalParams(systN, paramsGlobal, filename);

        const int nPoints = N_CENTR_SYST[systN];
        vector<double> x, y, ex, ey_stat, ey_syst;

        for (int i = 0; i < nPoints; ++i) {
            int centrIndex = CENTR_SYST[systN][i];
            double sum = 0, sumErr2_stat = 0, sumErr2_syst = 0;
            int count = 0;

            // Усреднение по всем частицам
            for (int part = 0; part < 6; ++part) {
                double val, err_stat, err_syst;
                if (paramName == "T_{kin}") {
                    val = paramsGlobal[part][centrIndex][1];
                    err_stat = paramsGlobal[part][centrIndex][2];
                    err_syst = T_syst_err[systN][part][centrIndex];
                } else {
                    val = paramsGlobal[part][centrIndex][3];
                    err_stat = paramsGlobal[part][centrIndex][4];
                    err_syst = beta_syst_err[systN][part][centrIndex];
                }
                
                sum += val;
                sumErr2_stat += pow(err_stat, 2);
                sumErr2_syst += pow(err_syst, 2);
                count++;
            }

            if (count == 0) continue;

            double mean = sum / count;
            double meanErr_stat = sqrt(sumErr2_stat) / count;
            double meanErr_syst = sqrt(sumErr2_syst) / count;

            x.push_back(Npart[systN][centrIndex]);
            y.push_back(mean);
            ex.push_back(Npart[systN][centrIndex] * 0.10); // 10% ошибка Npart
            ey_stat.push_back(meanErr_stat);
            ey_syst.push_back(meanErr_syst);
        }

        // Создание графиков
        TGraphErrors *gr = new TGraphErrors(x.size(), x.data(), y.data(), 
                                           ex.data(), ey_stat.data());
        gr->SetMarkerStyle(systMarkerStyles[systN]);
        gr->SetMarkerSize(3);
        gr->SetMarkerColor(systColors[systN]);
        gr->SetLineColor(systColors[systN]);
        gr->Draw("P same");

        // Систематические ошибки
        for (size_t i = 0; i < x.size(); ++i) {
            TBox *box = new TBox(x[i]-ex[i], y[i]-ey_syst[i], 
                                x[i]+ex[i], y[i]+ey_syst[i]);
            box->SetFillColorAlpha(systColors[systN], 0.3);
            box->SetLineColor(systColors[systN]);
            box->Draw("same");
        }

        leg->AddEntry(gr, systNamesTT[systN], "P");

        // Добавляем новую легенду с энергиями
        TLegend *legEnergy = new TLegend(0.05, 0.28, 0.8, 0.38);
        legEnergy->SetBorderSize(0);
        legEnergy->SetFillStyle(0);
        legEnergy->SetTextSize(0.04);
        legEnergy->SetNColumns(1);
        legEnergy->AddEntry((TObject*)0, "p+Al, He+Au, Cu+Au: #sqrt{s_{NN}} = 200 GeV", "");
        legEnergy->AddEntry((TObject*)0, "U+U: #sqrt{s_{NN}} = 193 GeV", "");
        legEnergy->Draw();
    }

    if (grART) {
        leg1->AddEntry(grART, "Au+Au, #sqrt{s_{NN}} = 200 GeV (PHENIX, C72:014903)", "P");
    }
    if (grSTAR) {
        leg1->AddEntry(grSTAR, "Au+Au, #sqrt{s_{NN}} = 39 GeV (STAR, Hot Quarks 2025)", "P");
    }

    leg1->Draw();
    leg->Draw();
    
    c->SaveAs(Form("output/pics/Npart/Averaged_%s_vs_Npart.png", paramName.Data()));
    c->SaveAs(Form("output/pics/Npart/Averaged_%s_vs_Npart.pdf", paramName.Data()));
}

void DrawAveragedTvsBeta() {
    TCanvas *c = new TCanvas("cAvgTbeta", "Averaged T vs #beta", 1200, 1000);
    c->SetGrid();
    c->SetLeftMargin(0.18);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.05);
    c->SetBottomMargin(0.15);

    // Настройка осей
    TH1F *frame = c->DrawFrame(0.1, 0.02, 1.2, 0.22);
    frame->SetTitle(";#beta_{T} [GeV]; T_{kin} [GeV]");
    frame->GetXaxis()->SetTitleSize(0.06);
    frame->GetYaxis()->SetTitleSize(0.06);
    frame->GetXaxis()->SetLabelSize(0.05);
    frame->GetYaxis()->SetLabelSize(0.05);

    TLegend *legend = new TLegend(0.71, 0.70, 0.91, 0.89);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.05);

    // Добавляем новую легенду с энергиями
    TLegend *legEnergy = new TLegend(0.05, 0.28, 0.8, 0.38);
    legEnergy->SetBorderSize(0);
    legEnergy->SetFillStyle(0);
    legEnergy->SetTextSize(0.042);
    legEnergy->AddEntry((TObject*)0, "p+Al, He+Au, Cu+Au: #sqrt{s_{NN}} = 200 GeV", "");
    legEnergy->AddEntry((TObject*)0, "U+U: #sqrt{s_{NN}} = 193 GeV", "");

    TLegend *leg = new TLegend(0.18, 0.17, 0.4, 0.28);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.042);

    for (int systN : SYSTS) {
        // Читаем параметры для системы
        TString filename = Form("output/parameters/ALL_FinalBWparams_%s.txt", 
                              systNamesT[systN].Data());
        ReadGlobalParams(systN, paramsGlobal, filename);

        const int nPoints = N_CENTR_SYST[systN];
        vector<double> avgBeta, avgT, betaStatErr, TStatErr, betaSystErr, TSystErr;

        for (int i = 0; i < nPoints; ++i) {
            int centrIndex = CENTR_SYST[systN][i];
            
            double sumT = 0, sumBeta = 0;
            double sumTerr2 = 0, sumBetaerr2 = 0;
            double sumTsyst2 = 0, sumBetasyst2 = 0;
            int count = 0;

            // Усреднение по всем частицам
            for (int part = 0; part < 6; ++part) {
                sumT += paramsGlobal[part][centrIndex][1];
                sumBeta += paramsGlobal[part][centrIndex][3];
                
                sumTerr2 += pow(paramsGlobal[part][centrIndex][2], 2);
                sumBetaerr2 += pow(paramsGlobal[part][centrIndex][4], 2);
                
                sumTsyst2 += pow(T_syst_err[systN][part][centrIndex], 2);
                sumBetasyst2 += pow(beta_syst_err[systN][part][centrIndex], 2);
                
                count++;
            }

            if (count == 0) continue;

            avgT.push_back(sumT / count);
            avgBeta.push_back(sumBeta / count);
            
            TStatErr.push_back(sqrt(sumTerr2) / count);
            betaStatErr.push_back(sqrt(sumBetaerr2) / count);
            
            TSystErr.push_back(sqrt(sumTsyst2) / count);
            betaSystErr.push_back(sqrt(sumBetasyst2) / count);
        }

        // Создаем график со статистическими ошибками
        TGraphErrors *gr = new TGraphErrors(avgBeta.size(), 
                                           avgBeta.data(), avgT.data(),
                                           betaStatErr.data(), TStatErr.data());
        gr->SetMarkerStyle(systMarkerStyles[systN]);
        gr->SetMarkerColor(systColors[systN]);
        gr->SetMarkerSize(4);
        gr->SetLineColor(systColors[systN]);
        gr->SetLineWidth(2);
        gr->Draw("P SAME");

        // Рисуем систематические ошибки
        for (size_t i = 0; i < avgBeta.size(); ++i) {
            TBox *box = new TBox(avgBeta[i] - betaSystErr[i],
                                avgT[i] - TSystErr[i],
                                avgBeta[i] + betaSystErr[i],
                                avgT[i] + TSystErr[i]);
            box->SetFillColorAlpha(systColors[systN], 0.3);
            box->SetLineColor(systColors[systN]);
            box->SetLineWidth(1);
            box->Draw("SAME");
        }

        legend->AddEntry(gr, systNamesTT[systN], "P");
        legEnergy->Draw();
    }

    // Добавляем теоретические данные (пример)
    TGraph *theory = new TGraph(13, beta_AuAu_ART, T_AuAu_ART);
    theory->SetMarkerStyle(43);
    theory->SetMarkerColor(kBlack);
    theory->SetMarkerSize(3.5);
    theory->Draw("P SAME"); 

    const int nPoints_STAR = 9; 
    TGraph *gr_STAR_ch = new TGraph(nPoints_STAR, beta_AuAu_STAR_ch, T_AuAu_STAR_ch);
    gr_STAR_ch->SetMarkerStyle(42);  
    gr_STAR_ch->SetMarkerColor(kBlack);
    gr_STAR_ch->SetMarkerSize(3.5);
    gr_STAR_ch->Draw("P SAME");

    leg->AddEntry(theory, "Au+Au, #sqrt{s_{NN}} = 200 GeV (PHENIX, C72:014903)", "P");
    leg->AddEntry(gr_STAR_ch, "Au+Au, #sqrt{s_{NN}} = 39 GeV (STAR, Hot Quarks 2025)", "P"); 

    leg->Draw();
    legend->Draw();

    c->SaveAs("output/pics/Tbeta/Averaged_T_vs_Beta.png");
    c->SaveAs("output/pics/Tbeta/Averaged_T_vs_Beta.pdf");
}

// Обновлённая функция для вызова всех графиков
void DrawAllAveragedParams() {
    DrawAveragedParamsNpart("T_{kin}");
    DrawAveragedParamsNpart("#beta_{T}");

    DrawAveragedParamsCentr("T_{kin}");
    DrawAveragedParamsCentr("#beta_{T}");

    DrawAveragedTvsBeta();
}

void NpartDrawParams_all(void) 
{
    // Для параметра T
    for (int systN: SYSTS) 
    {
        cout << systN << " " << systNamesTT[systN] << endl;
        SetGraphs(systN, "T_{kin}");
    }
    DrawParam("T_{kin}");
    
    // Для параметра beta
    for (int systN: SYSTS) 
    {
        cout << systN << " " << systNamesTT[systN] << endl;
        SetGraphs(systN, "beta_{T}");
    }
    DrawParam("beta_{T}");

    DrawAllParticleParams();
    DrawAllAveragedParams();
    WriteAveragesToFile();

    // Явная очистка холстов
    if (gROOT->GetListOfCanvases()) {
        gROOT->GetListOfCanvases()->Delete();
    }

    // Принудительный выход, минуя деструкторы ROOT
    gROOT->ProcessLine(".q");
}