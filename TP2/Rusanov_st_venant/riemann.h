#ifndef _RIEM_H
#define _RIEM_H


#define _SMALL 1e-6
#define _GPES 9.81

void riem_stvenant(double *wL, double *wR, double xi,	double *w);

double Z(double hs, double h);

double dZ(double hs, double h);

void plot_riem(double *wL, double *wR);

double Heaviside(double x);

double Dirac(double x);

#endif
