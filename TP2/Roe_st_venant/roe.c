#include "roe.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void roe_stvenant(double *wL, double *wR, double xi, double *w) {

  double g = _GPES;

  double hL = wL[0];
  double uL = wL[1]/wL[0];

  double hR = wR[0];
  double uR = wR[1]/wR[0];


  double uc = (uR + uL)/2;
  double hc = (hL + hR)/2;

  double cc = sqrt(g*hc);

  double lbd1 = uc - cc;
  double lbd2 = uc + cc;

  double h, u;
  if (lbd1 < 0 && lbd2 < 0){
      h = hR;
      u = uR;
  }
  else if (lbd1 > 0 && lbd2 > 0) {
      h = hL;
      u = uL;
  } else {
      h = hc*(1 - (uR - uL)/cc);
      u = uc - g*(hR - hL)/2/cc;
  }
  w[0] = h;
  w[1] = h*u;

}
