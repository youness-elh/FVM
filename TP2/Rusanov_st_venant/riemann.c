#include "riemann.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//
// void main(void) {
//   double hL = 1.;
//   double uL = 1.;
//
//   double hR = 1.;
//   double uR = -1.;
//
//   double wL[2] = {hL, hL * uL};
//   double wR[2] = {hR, hR * uR};
//
//   double xi;
//   double w[2];
//
//   riem_stvenant(wL, wR, xi, w);
//   plot_riem(wL,wR);
//
// }

void riem_stvenant(double *wL, double *wR, double xi, double *w) {

  double g = _GPES;

  double hL = wL[0];
  double uL = wL[1]/wL[0];

  double hR = wR[0];
  double uR = wR[1]/wR[0];

  double hs = _SMALL;
  int itermax = 10;

  for (int i = 0; i < itermax; ++i){
    double f = uL - (hs - hL) * Z(hs, hL) - uR - (hs - hR) * Z(hs, hR);
    double df = -(hs - hL) * dZ(hs, hL) - Z(hs, hL) - (hs - hR) * dZ(hs,hR) - Z(hs,hR);
    double dhs = -f / df;
    hs += dhs;

    // printf("i=%d f=%e df=%e hs=%e dhs=%e \n", i, f, df, hs, dhs);
  }

  double us = uL - (hs - hL) * Z(hs, hL);

  double v1m, v1p, v2m, v2p;

  // 1 - onde
  if (hs < hL){ // détente

    v1m = uL - sqrt(g*hL);
    v1p = us - sqrt(g*hs);

  } else { // choc

    double a = sqrt(hs) / (sqrt(hs) + sqrt(hL));
    double u = a * us + (1 - a) * uL;
    double h = (hs + hL) / 2;

    v1m = u - sqrt(g*h);
    v1p = v1m;

  }

  // 2 - onde
  if (hs < hR){ // détente

    v2m = us + sqrt(g*hs);
    v2p = uR + sqrt(g*hR);

  } else { // choc

    double a = sqrt(hs) / (sqrt(hs) + sqrt(hR));
    double u = a * us + (1 - a) * uR;
    double h = (hs + hR) / 2;

    v2m = u + sqrt(g*h);
    v2p = v2m;

  }

  // printf("v = %f %f %f %f\nhs=%f us=%f\n", v1m, v1p, v2m, v2p, hs, us);

  if (xi < v1m){

    w[0] = wL[0];
    w[1] = wL[1];

  } else if (xi < v1p) {

    double u = (uL + 2*xi + 2*sqrt(g * hL)) / 3;
    double h = (u - xi) * (u - xi) / g;

    w[0] = h;
    w[1] = h * u;

  } else if (xi < v2m) {

    w[0] = hs;
    w[1] = hs * us;

  } else if (xi < v2p) {

    double u = (uR + 2*xi - 2*sqrt(g*hR)) / 3;
    double h = (u - xi) * (u - xi) / g;

    w[0] = h;
    w[1] = h * u;

  } else {

    w[0] = wR[0];
    w[1] = wR[1];

  }

}


void plot_riem(double *wL, double *wR){

  FILE *plotfile;

  double xmin = -10;
  double xmax = 10;

  int n = 300;
  double dx = (xmax - xmin)/n;

  plotfile = fopen("riem.dat","w");

  for (int i = 0; i < n; ++i){
    double xi = xmin + i*dx;
    double w[2];
    riem_stvenant(wL,wR,xi,w);
    double h = w[0];
    double u = w[1]/w[0];
    printf("h=%f\n",h);
    printf("u = %f",u);
    fprintf(plotfile, "%f %f %f\n", xi, h, u);
  }
  fclose(plotfile);

  system("gnuplot riemcom");
}

double Z(double hs, double h){

  double g = _GPES;

  if (hs < h)
    return 2*sqrt(g)/(sqrt(hs)+sqrt(h));

  return sqrt(g*(hs+h)/(2*hs*h));

}

double dZ(double hs, double h){

  double g = _GPES;

  if (hs < h) /* choc */
    return - sqrt(g) / sqrt(hs) / (sqrt(hs)+sqrt(h)) / (sqrt(hs)+sqrt(h));
  /* détente */
  return - sqrt(g * hs * h / 2 / (hs + h)) / 2 / hs / hs;

}

double Heaviside(double x){
  return (x>0)? 1 : 0;
}

double Dirac(double x){
  return 0;
}
