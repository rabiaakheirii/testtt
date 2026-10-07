#include <iostream>
#include <fstream>
#include "TF1.h"
#include "TMath.h"

using namespace std;

//	blastwave function from E.J.Kim "flow_c12.C"
double bwfitfunc(double *x, double *par)
{
    //  *x is r in this case
    //  pt is par 3
	//	par[0] = constant;
	//	par[1] = Tf;
	//	par[2] = beta;
	//  par[3] = mass	
	//  par[4] = pt
	double out; 
	
	double con = par[0]; 
	double mass = par[3];
    double mt = par[4] + mass;
	double pt = sqrt( mt * mt - mass * mass );

	double rho = TMath::ATanH(par[2]) * pow(x[0] / 13.0, 1);
	out = con * x[0] * mt 
			* TMath::BesselI0(pt * TMath::SinH(rho) / par[1])
			* TMath::BesselK1(mt * TMath::CosH(rho) / par[1]);

	return out;
}


//	structure representing the integral of a function between 0 and radius
struct MyIntegFunc
{
	//	constructor using the TF1 pointer
	MyIntegFunc(TF1 *f):
		fFunc(f) {}

	TF1 *fFunc; // pointer to the integral function
   	double param[5]; 

	//	evaluate the integral of fFunc (pt, r) in r
	double operator() (double *x, double *p) 
	{
		// *x is pt in this case
		// p[] is Tf, alpha and beta
		double radius = 13.0;	//	radius = 13.0 fm (Rmax)
		std::copy(p, p + 4, param); 
		param[4] = *x;    // set value of pt for integrand function (fFunc)
		// cout << param[0] << " " << param[1] << " " << param[2] << " " << param[3] << " " << endl;
		fFunc->SetParameters(param);
		return fFunc->Integral(0.0001, radius, 1.e-10);
	}
};