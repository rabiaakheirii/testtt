#ifndef __CUMULATIVE_H_
#define __CUMULATIVE_H_

#include "math.h"
#include "TF1.h"
#include "def.h"

const double mN = 0.931494;
const double mPI = 0.13957061;
const double PI = acos(-1);

double mI = mN, mII = mN;  // - сталкиваются нуклон ядра и нуклон мишени
double m2 = 2 * mN; // - рождается 1 частица с массой m1 и остаются два нуклона с массой m2
double Minit = mI + mII; // масса покоя начальной системы сталкивающихся нуклонов
    
double maxE(double *x, double *params)
{
    double Etot, W, m1, E1, E2, A1, pI, theta1;

    W = params[0]; // - кинетическая энергия частицы I
    m1 = params[1]; // - масса рождающейся кумулятивной частицы
    //theta1 = params[1]; // угол вылета пи-мезона
    theta1 = x[0] * PI / 180.; 

    Etot = mI + mII + W; // - полная энергия
    pI = sqrt(W * W + 2 * W * mI);

    A1 = Minit * Minit + 2 * mII * W + m1 * m1 - m2 * m2;
    
    double sq = A1 * A1 - 4 * m1 * m1 * (Etot * Etot - pow( pI * cos(theta1), 2));
    double E1_nom = A1 * Etot + pI * cos(theta1) * sqrt(sq);
    double E1_denom = 2 * (Etot * Etot - pow( pI * cos(theta1), 2));

    E1 = E1_nom / E1_denom;

    return E1;
}

double calculateW( double sNN ) // считает кин энергию налетающей частицы в лск
{
    return (sNN - mI * mI - mII * mII) / 2 * mII - mI;
}

double getTheta( double y ) // считает угол по быстроте
{
    return 2 * atan(exp(y)) * 180 / PI;
}


void DrawCumulativeBorder( int part, double pad_min, double pad_max )
{
    double sNN = 9.2 * 9.2;
    double W = calculateW(sNN);
    double params[2] = {W, masses[part]};

    double theta1 = getTheta(0.5); // быстрота = 3 -> theta = 5.7
    double Emax = maxE(&theta1, params);
    double pmax = sqrt(Emax * Emax - masses[part] * masses[part]);

    TLine *line= new TLine(pmax, pad_min, pmax, pad_max); 
    line->SetLineWidth(2);
    line->SetLineStyle(kDashed);
    line->Draw("same");

}

#endif /* __CUMULATIVE_H_ */