#include "TF1.h"
#include "TLegend.h"
#include "input/FormatOfEverything.h"
#include "TMath.h"
using namespace std;

const double PI = 3.1415;
int N = 23;

void ReadSpectra( int systN, int part, float pt[30], float y[5][30], float y_e[5][30], float y_s[5][30] )
{
    float tmp;
    int Ntemp;

    ifstream myfile_p;
    char name_file[100];
    if (systN == 0) sprintf(name_file, "input/PHENIX_spectra/spectrapAl/Spectra_particle_%d_%d.txt", part, part);
    if (systN == 1) sprintf(name_file, "input/PHENIX_spectra/spectraHeAu/Spectra_particle_%d_%d.txt", part, part);
    if (systN == 2) sprintf(name_file, "input/PHENIX_spectra/spectraCuAu/Spectra_particle_%d_%d.txt", part, part);
    if (systN == 3) sprintf(name_file, "input/PHENIX_spectra/spectraUU/Spectra_particle_%d_%d.txt", part, part);

    myfile_p.open(name_file);
    myfile_p >> Ntemp;

    for (int centr = 0; centr < 5; centr++)
    {
        for (int i = 0; i < Ntemp; i++)
        {
            myfile_p >> pt[i] >> y[centr][i] >>y_e[centr][i]>>y_s[centr][i];
        }
        for (int i = Ntemp; i < N; i++)
        {
            y[centr][i] = 0;
        }
    }
        
    myfile_p.close();
}

Double_t func(Double_t *x, Double_t *par)
{
    double pT = x[0], mT, T, m, mu, g = 2;

    T = par[0];
    m = par[1];
    mu = par[2];
    mT = sqrt(m * m + pT * pT);
   
    double R = 14; // Au
    double V = 4. / 3. * PI * pow(R, 3);
    double Const = g * V / (2 * PI * PI);

    Double_t f = Const * m * m  * T * TMath::BesselK(2, m / T) * exp(mu / T);
    return f;
}

double Ngc (double T, double mu, double m)
{
    double g = 2.;
    double R = 14; // Au
    double V = 4. / 3. * PI * pow(R, 3);
    double Const = g * V / (2 * PI * PI);

    double beta = 0.55;
    double gamma = 1.0 / sqrt(1 - beta * beta);
    double Teff = T * gamma * (1 + beta);

    T= Teff;
    cout << T << endl;
    Double_t f = Const * m * m  * T * TMath::BesselK(2, m / T) * exp(mu / T);
    return f;
}

void ChemicalPotential_p_pi_PHENIX ( void )
{
    double mPi = 0.139, mP = 0.938, Tch = 0.177;

    double mu = 0.029;
    TF1 *funcP = new TF1("funcP", func, 0, 1.05, 3);
    funcP->SetParameters(Tch, mP, mu);

    int isyst = 3; // Cu+Au

    float pT[30], y[5][30], y_e[5][30], y_s[5][30];
    ReadSpectra(isyst, 4, pT, y, y_e, y_s);

    double Rat[5];

    for (int centr: {1})
    {
        double sumP = 0, maxPt = 0;

        for (int i = 0; i < 30; i++)
        {
            if (y[centr][i] > 1e4) continue;
            if (y[centr][i] < 1e-10) continue;
            // if (pT[i] > 2.1)  continue;
            double dpT = pT[i] > 2. ? 0.2 : 0.1;
        
            sumP += y[centr][i] * pT[i]  * pT[i]* dpT * 2 * PI; 
            cout << pT[i] << " " << dpT << " " << y[centr][i] * pT[i] * pT[i] * dpT * 2 * PI * 0.35 * 2 << endl;
            maxPt = pT[i];
        }
        cout << maxPt << endl;
        
        //double Ip = funcP->Integral(0.1, maxPt);
        double Ip =Ngc(Tch, mu, mP);
        //funcP->Draw();
        cout << "Ip: " << Ip << "  sumP: " << sumP << endl;     
    }
        
}
