#ifndef _RIEM_H
#define _RIEM_H

void riem_stvenant(double *wL, double *wR, double xi,	double *w);

void flux_riem_2d(double *wL, double *wR, double *vnorm, double *flux);

double Z(double hs, double h);

double dZ(double hs, double h);

void plot_riem(double *wL, double *wR);

double Heaviside(double x);

double Dirac(double x);

#endif
