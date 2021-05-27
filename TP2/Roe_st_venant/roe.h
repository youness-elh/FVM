#ifndef _ROE_H
#define _ROE_H


#define _SMALL 1e-6
#define _GPES 9.81

void roe_stvenant(double *wL, double *wR, double xi, double *w);

#endif
