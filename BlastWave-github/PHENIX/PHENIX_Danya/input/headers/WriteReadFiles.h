#ifndef __WRITEREADFILES_H_
#define __WRITEREADFILES_H_

#include "def.h"


/* ================================ BlastWaveGlobal.C ================================ */

void ReadFromFile(int part, int systN)
{
    int N; 
    // mT – поперечная масса, pT – поперечный импульс, s – спектральные данные, 
    // s_e – ошибки, s_s – систематические ошибки, x_e – горизонтальные ошибки
    double mT[30], pT[30], s[30], s_e[30], s_s[30], x_e[30];

    for (int i = 0; i < 30; i++) x_e[i] = 0.05;

    ifstream f;
    string fileName = "input/PHENIX/" + systNames[systN] + "/Spectra_particle_" + to_string(part) + "_" + to_string(part) + ".txt";
    cout << fileName << endl;
    f.open(fileName.c_str());

    f >> N;
    for (int centr = 0; centr < Ncentr[systN]; centr++)
    {   
        for (int i = 0; i < N; i++)
        {
            // Read the data: pT, s, s_e, and s_s. These will be read in a way that C++ correctly interprets any format.
            f >> pT[i] >> s[i] >> s_e[i] >> s_s[i];
          
            // Calculate mT with respect to the given pT and particle mass
            mT[i] = sqrt(pT[i] * pT[i] + masses[part] * masses[part]) - masses[part];
        }

        // If systN == 1, multiply the entire spectrum for the current centrality by 10^centr
        if (systN == 0 || systN == 1 || systN == 2 || systN == 3 || systN == 4)
        {
            double scaleFactor = pow(10, centr);  // 10^centr для текущей центральности
            for (int i = 0; i < N; i++)
            {
                // Scale all the spectral values (s[i]) including those in scientific notation
                s[i] *= scaleFactor;  
            }
        }        

        // Create a graph with errors (TGraphErrors) for the current part and centrality, and store it in the global array grSpectra
        grSpectra[part][centr] = new TGraphErrors(N, mT, s, s_e, x_e);
    }
}


// Функция для чтения спектральных данных для системы AuAu
void ReadFromFileAuAu( void )
{
    int N;
    // mT и pT – поперечные масса и импульс, s и s_e – двумерные массивы (по центральностям и точкам), x_e – ошибки по оси X
    double mT[30], pT[30], s[12][30], s_e[12][30], x_e[30];

    for (int i = 0; i < 30; i++) x_e[i] = 0.05;

    ifstream f;
    f.open("input/PHENIX/AuAu/spectra.txt");

    for (int part = 0; part < 6; part++)
    {
        f >> N;
        cout << N << endl; 

        for (int i = 0; i < N; i++)
        {   
            f >> pT[i];
            mT[i]  = sqrt(pT[i] * pT[i] + masses[part] * masses[part]) - masses[part];

            for (int centr = 0; centr < Ncentr[0]; centr++)
            {
                f >> s[centr][i] >> s_e[centr][i];
            }
        }

        // Создаём графики с ошибками для каждой центральности и сохраняем в глобальный массив grSpectra
        for (int centr = 0; centr < Ncentr[0]; centr++)
        {
            grSpectra[part][centr] = new TGraphErrors(N, mT, s[centr], x_e, s_e[centr]);
        }
    }
}


// Извлечение глобальных параметры модели 
void getGlobalParams( int part, int centr, double parResults[4] )
{
    int charge = part % 2; 

    parResults[0] = paramsGlobal[charge][centr][2 + part / 2];
    parResults[1] = paramsGlobal[charge][centr][0]; 
    parResults[2] = paramsGlobal[charge][centr][1]; 
    parResults[3] = masses[part];     
}


// Запись параметров в нужный файл (если не задано другое название - используется стандартное из функции)
void WriteGlobalParams( bool *isParamsFileExist, int charge, int systN, const char filename[30] = "output/txtParams/GlobalBWparams.txt" )
{
   cout << " WriteGlobalParams " << endl;
   ofstream txtFile;
   if (*isParamsFileExist) 
   {
      txtFile.open(filename, ios::app);
      txtFile << endl;
   }
   else txtFile.open(filename);

   for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
    int centr = CENTR_SYST[systN][j];
      txtFile  << charge << "  "  << centr << "  " 
               << paramsGlobal[charge][centr][0] << "  " 
               << paramsGlobal[charge][centr][1] << "  "
               << paramsGlobal[charge][centr][2] << "  " 
               << paramsGlobal[charge][centr][3] << "  " 
               << paramsGlobal[charge][centr][4] << endl;
   }
    
   txtFile.close();
   *isParamsFileExist = true;
}


/* ================================ BlastWaveFinal.C ================================ */


// Запись параметров от финального перефита с использованием кейса из input/BlastWaveFit.h
// Запись параметров производится в файл с определённым названием, если не задано другого 
void WriteParams( int systN, double par[N_PARTS][N_CENTR][4], double parErr[N_PARTS][N_CENTR][4],
                  int printAll = true, const char filename[30] = "output/parameters/BWparams_pAl.txt")
{
    ofstream txtFile;

    if (printAll)
    {
        txtFile.open(filename);
        for (int part: PARTS)
        {
            for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
                int centr = CENTR_SYST[systN][j];
                txtFile << part << "  " 
                        << centr << "  " 
                        << par[part][centr][0] << "  "
                        << par[part][centr][1] << "  " 
                        << parErr[part][centr][1] << "  "
                        << par[part][centr][2] << "  " 
                        << parErr[part][centr][2] << endl;
            }
        }
    }
    else 
    {
        txtFile.open("output/parameters/BWparams_noErrs.txt");
        string centrTitle[6] = {"0-10\\%", "10-20\\%", "20-30\\%", "30-40\\%", "40-60\\%", "60-80\\%"};
        for (int part: PARTS)
        {
            for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
                int centr = CENTR_SYST[systN][j];
                txtFile << centrTitle[centr] << " & " << 
                int(par[part][centr][1] * 1000) << " & " << 
                floorf(par[part][centr][2] * 100) / 100. << " \\\\ " << endl;
            }
            txtFile << endl;
        }
    }
    
    txtFile.close();
}


// Запись параметров с систематикой в файл (добавлен параметр systN)
void WriteParamsSyst( int systN, double par[N_PARTS][N_CENTR][4], double parErr[N_PARTS][N_CENTR][4], double parSyst[N_PARTS][N_CENTR][4],
    const char filename[30] = "output/parameters/BWparamsSyst.txt" )
{
    ofstream txtFile;
    txtFile.open(filename);

    for (int part : PARTS)
    {
        for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
            int centr = CENTR_SYST[systN][j];
            txtFile << part << "  " << centr << "  " 
                << par[part][centr][0] << "  "
                << par[part][centr][1] << "   " << parErr[part][centr][1] << "   "  << parSyst[part][centr][1] << "   "
                << par[part][centr][2] << "   " << parErr[part][centr][2] << "   "  << parSyst[part][centr][2] << endl;
        }
    }

    txtFile.close();
}


/* ================================ BWDrawParams.C ================================ */


// Чтение финальных параметров от BlastWaveFinal.C в зависимости от bool isSyst = true or false
void ReadParam( int parN, double par[N_PARTS][N_CENTR], double parErr[N_PARTS][N_CENTR],
                 const char filename[30] = "output/parameters/FinalBWparams_AuAu.txt" )
{
    ifstream f;
    f.open(filename);

    int p, c;
    double tmpPar[4] = {0, 0, 0, 0}, tmpParErr[4] = {0, 0, 0, 0};

    for (int part: PARTS)
    {
        for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
            int centr = CENTR_SYST[systN][j];
            f >> p >> c 
              >> tmpPar[0] 
              >> tmpPar[1] 
              >> tmpParErr[1] 
              >> tmpPar[2] 
              >> tmpParErr[2];

            par[part][centr] = tmpPar[parN];
            parErr[part][centr] = tmpParErr[parN];
        }
    }
    f.close();
}

void ReadParam( int parN, double par[N_PARTS][N_CENTR], double parErr[N_PARTS][N_CENTR], double parSyst[N_PARTS][N_CENTR],
                 const char filename[30] = "output/parameters/ALL_FinalBWparams_HeAu.txt" )
{
    ifstream f;
    f.open(filename);

    int p, c;
    double tmpPar[4] = {0, 0, 0, 0}, tmpParErr[4] = {0, 0, 0, 0}, tmpParSyst[4] = {0, 0, 0, 0};

    for (int part: PARTS)
    {
        for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
            int centr = CENTR_SYST[systN][j];
            f >> p >> c 
              >> tmpPar[0] 
              >> tmpPar[1] 
              >> tmpParErr[1] 
              >> tmpParSyst[1] 
              >> tmpPar[2] 
              >> tmpParErr[2]
              >> tmpParSyst[2];

            //     >> tmpPar[0] 
            //     >> tmpPar[1] 
            //     >> tmpParErr[1] 
            // //   >> tmpParSyst[1] 
            //     >> tmpPar[2] 
            //     >> tmpParErr[2]; 
            // //   >> tmpParSyst[2];
              
            par[part][centr] = tmpPar[parN];
            parErr[part][centr] = tmpParErr[parN];
            // parSyst[part][centr] = tmpParSyst[parN] * par[part][centr];
        }
    }
    f.close();
}



/* ================================ GlobalDrawParams.C ================================ */


// Чтение параметров из глобального фита (можно брать другие данные, но с тем же форматом)
void ReadGlobalParams(int systN, double paramsGlobal[2][N_CENTR][5], const char filename[30] = "output/parameters/ALL_GlobalBWparams_AuAu.txt")
{
    std::ifstream txtFile;
    txtFile.open(filename);

    int charge;
    while (true)
    {
        for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
            int centr = CENTR_SYST[systN][j];
            
            // Чтение данных из файла
            txtFile 
            >> charge 
            >> centr 
            >> paramsGlobal[charge][centr][0] // Константа
            >> paramsGlobal[charge][centr][1]  // T value
            >> paramsGlobal[charge][centr][2]  // T error
            >> paramsGlobal[charge][centr][3]  // beta value
            >> paramsGlobal[charge][centr][4]; // beta error

            // // Проверка на NaN и пропуск, если найдено
            // if (std::isnan(paramsGlobal[charge][centr][0]) ||
            //     std::isnan(paramsGlobal[charge][centr][1]) ||
            //     std::isnan(paramsGlobal[charge][centr][2]) ||
            //     std::isnan(paramsGlobal[charge][centr][3]) ||
            //     std::isnan(paramsGlobal[charge][centr][4])) {
            //     std::cerr << "Warning: NaN found in line " << charge << " " << centr << ", skipping this entry." << std::endl;
            //     continue; // Пропустить эту строку, если хотя бы одно значение NaN
            // }

            // cout << paramsGlobal[charge][centr][0] << endl;
        }
        
        if (txtFile.eof()) break;
    }
    
    txtFile.close();
}


// Расчёт средней температуры и средней скорости радиального потока наших данных
void CalculateAverage(TString paramName) {
    double sum = 0.0;
    int count = 0;

    for (int systN : SYSTS) {
        // Пропускаем системные номера, соответствующие ARTICLE и STAR
        if(systNamesT[systN].Contains("ARTICLE") || systNamesT[systN].Contains("STAR")) 
        continue;

        for (int charge : {0, 1}) 
        {
            for (int centr = 0; centr < N_CENTR_SYST[systN]; centr++) 
            {
                double value = (paramName == "T") ? Tpar[charge][centr] : utPar[charge][centr];
                sum += value;
                count++;
            }
        }
    }

    if (count > 0) 
    {
        if(paramName == "T") 
        {
            gAvgT = sum / count;
            cout << "Average " << paramName << " = " << gAvgT 
            << " (from " << count << " points)" << endl;
        } 
        else 
        {
            gAvgUt = sum / count;
            cout << "Average " << paramName << " = " << gAvgUt 
             << " (from " << count << " points)" << endl;
        }
        
    }
}


// Поиск и запись среднего значения всех точек на графиках для параметров T и u_t 
void WriteAveragesToFile() 
{
    ofstream outFile("output/parameters/AVR_T_beta.txt", ios::out);

    outFile << "T_avr = " << gAvgT * 1000 << " MeV\n";
    outFile << "beta_avr = " << gAvgUt * 1000 << " MeV\n";
    
    outFile.close();
}













void WriteParams( double par[N_PARTS][N_CENTR][5], double parErr[N_PARTS][N_CENTR][5],
                  const char filename[30] = "output/txtParams/HagedornParams.txt" )
{
    ofstream txtFile;
    txtFile.open(filename);

    for (int part: PARTS)
    {
        for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
            int centr = CENTR_SYST[systN][j];
            txtFile << part << "  " << centr << "  " << par[part][centr][0] << "  "
                    << par[part][centr][1] << "   " << parErr[part][centr][1] << "   "
                    << par[part][centr][2] << "   " << parErr[part][centr][2] << "   "
                    << par[part][centr][3] << "   " << parErr[part][centr][3] << endl;
        }
    }
    
    txtFile.close();
}



void ReadParams( double par[N_PARTS][N_CENTR][4], double parErr[N_PARTS][N_CENTR][4],
                 const char filename[30] = "output/txtParams/BWparams.txt" )
{
    ifstream f;
    f.open(filename);

    int p, c;

    for (int part: PARTS)
    {
        for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
            int centr = CENTR_SYST[systN][j];
            f >> p >> c >> par[p][c][0]
              >> par[p][c][1] >> parErr[p][c][1] >> par[p][c][2] >> parErr[p][c][2];
        }
    }
    f.close();
}

void ReadParams( int part, int centr, double par[4],
                 const char filename[30] = "output/parameters/BWparams.txt" )
{
    ifstream f;
    f.open(filename);

    int p, c;
    double dummyErr;

    for (int part_: PARTS)
    {
        for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
            int centr = CENTR_SYST[systN][j];

            f >> p >> c >> par[0] >> par[1] >> dummyErr >> par[2] >> dummyErr;
            if (j == centr) break;
        }

        par[3] = masses[part];
        if (part_ == part) break;
    }
    f.close();
}



#endif /* __WRITEREADFILES_H_ */