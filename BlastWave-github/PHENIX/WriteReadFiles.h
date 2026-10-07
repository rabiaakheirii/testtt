#ifndef __WRITEREADFILES_H_
#define __WRITEREADFILES_H_

#include "def.h"

void ReadFromFile( int part, int systN )
{
    int N;
    double mT[30], pT[30], s[30], s_e[30], s_s[30], x_e[30];
    for (int i = 0; i < 30; i++) x_e[i] = 0.05;

    ifstream f;
    string fileName = "input/PHENIX_spectra/" + systNames[systN] + "/Spectra_particle_" + to_string(part) + "_"+ to_string(part) +".txt";
    cout << fileName << endl;
    f.open(fileName.c_str());

    f >> N;
    for (int centr = 0; centr < Ncentr[systN]; centr++)
    {   
        for (int i = 0; i < N; i++)
        {
            f >> pT[i] >> s[i] >> s_e[i] >> s_s[i];
            mT[i]  = sqrt(pT[i] * pT[i] + masses[part] * masses[part]) - masses[part];
        }
        grSpectra[part][centr] = new TGraphErrors(N, mT, s, s_e, x_e);
    }
}

void ReadFromFileAuAu( void )
{
    int N;
    double mT[30], pT[30], s[12][30], s_e[12][30], x_e[30];
    for (int i = 0; i < 30; i++) x_e[i] = 0.05;

    ifstream f;
    f.open("input/PHENIX_spectra/AuAu/spectra.txt");

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
                if (part == 2 || part == 3) s[centr][i] *= 10;
            }
        }
        for (int centr = 0; centr < Ncentr[0]; centr++)
        {
            grSpectra[part][centr] = new TGraphErrors(N, mT, s[centr], x_e, s_e[centr]);
        }
    }
}


void WriteGlobalParams( bool *isParamsFileExist, int charge, const char filename[30] = "output/txtParams/GlobalBWparams.txt" )
{
   cout << " WriteParams " << endl;
   ofstream txtFile;
   if (*isParamsFileExist) 
   {
      txtFile.open(filename, ios::app);
      txtFile << endl;
   }
   else txtFile.open(filename);

   for (int centr: CENTR)
   {
      txtFile  << charge << "  "  << centr << "  " 
               << paramsGlobal[charge][centr][0]  << "  " << paramsGlobal[charge][centr][1] << "  "
               << paramsGlobal[charge][centr][2] << "   " << paramsGlobal[charge][centr][3] << "   " << paramsGlobal[charge][centr][4] << endl;
   }
    
   txtFile.close();
   *isParamsFileExist = true;
}

void ReadGlobalParams( double paramsGlobal[2][N_CENTR][5], const char filename[30] = "output/txtParams/GlobalBWparams.txt" )
{
    cout << " ReadParams: " << filename << endl;
    ifstream txtFile;
    txtFile.open(filename);

    int charge;
    while (true)
    {
        for (int centr: CENTR)
        {
            txtFile >> charge >> centr 
                    >> paramsGlobal[charge][centr][0] >> paramsGlobal[charge][centr][1]
                    >> paramsGlobal[charge][centr][2] >> paramsGlobal[charge][centr][3] >> paramsGlobal[charge][centr][4];

            // cout << paramsGlobal[charge][centr][0] << endl;
        }
        
        if( txtFile.eof() ) break;
    }
    
    
    txtFile.close();
}


void WriteParams( double par[N_PARTS][N_CENTR][4], double parErr[N_PARTS][N_CENTR][4],
                  int printAll = true, const char filename[30] = "output/txtParams/BWparams.txt")
{
    ofstream txtFile;

    if (printAll)
    {
        txtFile.open(filename);
        for (int part: PARTS)
        {
            for (int centr: CENTR)
            {
                txtFile << part << "  " << centr << "  " << par[part][centr][0] << "  "
                        << par[part][centr][1] << "   " << parErr[part][centr][1] << "   "
                        << par[part][centr][2] << "   " << parErr[part][centr][2] << endl;
            }
        }
    }
    else 
    {
        txtFile.open("output/txtParams/BWparams_noErrs.txt");
        string centrTitle[6] = {"0-10\\%", "10-20\\%", "20-30\\%", "30-40\\%", "40-60\\%", "60-80\\%"};
        for (int part: PARTS)
        {
            for (int centr: CENTR)
            {
                txtFile << centrTitle[centr] << " & " << int(par[part][centr][1] * 1000) << " & "  << floorf(par[part][centr][2] * 100) / 100. << " \\\\ " << endl;
            }
            txtFile << endl;
        }
    }
    
    
    
    txtFile.close();
}

void WriteParams( double par[N_PARTS][N_CENTR][5], double parErr[N_PARTS][N_CENTR][5],
                  const char filename[30] = "output/txtParams/HagedornParams.txt" )
{
    ofstream txtFile;
    txtFile.open(filename);

    for (int part: PARTS)
    {
        for (int centr: CENTR)
        {
            txtFile << part << "  " << centr << "  " << par[part][centr][0] << "  "
                    << par[part][centr][1] << "   " << parErr[part][centr][1] << "   "
                    << par[part][centr][2] << "   " << parErr[part][centr][2] << "   "
                    << par[part][centr][3] << "   " << parErr[part][centr][3] << endl;
        }
    }
    
    txtFile.close();
}

void WriteParamsSyst( double par[N_PARTS][N_CENTR][4], double parErr[N_PARTS][N_CENTR][4], double parSyst[N_PARTS][N_CENTR][4],
                  const char filename[30] = "output/txtParams/BWparamsSyst.txt" )
{
    ofstream txtFile;
    txtFile.open(filename);

    for (int part: PARTS)
    {
        for (int centr: CENTR)
        {
            txtFile << part << "  " << centr << "  " << par[part][centr][0] << "  "
                    << par[part][centr][1] << "   " << parErr[part][centr][1] << "   "  << parSyst[part][centr][1] << "   "
                    << par[part][centr][2] << "   " << parErr[part][centr][2] << "   "  << parSyst[part][centr][2] << endl;
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
        for (int centr: CENTR)
        {
            f >> p >> c >> par[p][c][0]
              >> par[p][c][1] >> parErr[p][c][1] >> par[p][c][2] >> parErr[p][c][2];
        }
    }
    f.close();
}

void ReadParams( int part, int centr, double par[4],
                 const char filename[30] = "output/txtParams/BWparams.txt" )
{
    ifstream f;
    f.open(filename);

    int p, c;
    double dummyErr;

    for (int part_: PARTS)
    {
        for (int centr_: CENTR)
        {
            f >> p >> c >> par[0] >> par[1] >> dummyErr >> par[2] >> dummyErr;
            if (centr_ == centr) break;
        }
        par[3] = masses[part];
        if (part_ == part) break;
    }
    f.close();
}


void ReadParam( int parN, double par[N_PARTS][N_CENTR], double parErr[N_PARTS][N_CENTR],
                 const char filename[30] = "output/txtParams/BWparams.txt" )
{
    ifstream f;
    f.open(filename);

    int p, c;
    double tmpPar[4] = {0, 0, 0, 0}, tmpParErr[4] = {0, 0, 0, 0};

    for (int part: PARTS)
    {
        for (int centr: CENTR)
        {
            f >> p >> c >> tmpPar[0] >> tmpPar[1] >> tmpParErr[1] >> tmpPar[2] >> tmpParErr[2];
            par[part][centr] = tmpPar[parN];
            parErr[part][centr] = tmpParErr[parN];
        }
    }
    f.close();
}

void ReadParam( int parN, double par[N_PARTS][N_CENTR], double parErr[N_PARTS][N_CENTR], double parSyst[N_PARTS][N_CENTR],
                 const char filename[30] = "output/txtParams/BWparamsSyst.txt" )
{
    ifstream f;
    f.open(filename);

    int p, c;
    double tmpPar[4] = {0, 0, 0, 0}, tmpParErr[4] = {0, 0, 0, 0}, tmpParSyst[4] = {0, 0, 0, 0};

    for (int part: PARTS)
    {
        for (int centr: CENTR)
        {
            f >> p >> c 
              >> tmpPar[0] >> tmpPar[1] >> tmpParErr[1] >> tmpParSyst[1] 
              >> tmpPar[2] >> tmpParErr[2] >> tmpParSyst[2];
            par[part][centr] = tmpPar[parN];
            parErr[part][centr] = tmpParErr[parN];
            parSyst[part][centr] = tmpParSyst[parN] * par[part][centr];
        }
    }
    f.close();
}

#endif /* __WRITEREADFILES_H_ */