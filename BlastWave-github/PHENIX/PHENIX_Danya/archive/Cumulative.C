#include "Cumulative.h"


void Cumulative ( void ) 
{
    double sNN = 9.2 * 9.2;
    double W = calculateW(sNN);

    auto Efunc = new TF1("funcx", maxE, 0.1, 179, 1);
    Efunc->SetParameters(2,1);
    Efunc->FixParameter(0, W);
    Efunc->Draw();

    double theta1 = getTheta(0.5); // быстрота = 3 -> theta = 5.7
    double m1 = mN; // рождается кумулятивный протон
    double params[2] = {W, m1}; 
    double Emax = maxE(&theta1, params);
    double pmax = sqrt(Emax * Emax - m1 * m1);

    cout << pmax << endl;
}