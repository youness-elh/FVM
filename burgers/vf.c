#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
//#include <gc.h> // conservative garbage collector

#include "vf.h"

//#define GC_MALLOC_ATOMIC malloc // without gc


int main(void){

  double xmin = -1;
  double xmax = 2;

  int m = 1;
  int n = 100;

  double cfl = 0.5;

  vf svf = {0};
  vf_init(&svf, m, n, cfl, xmin, xmax);

  double tmax = 2;
  svf.tmax = tmax;
  vf_solve(&svf, tmax);

  vf_plot(&svf);

  vf_error(&svf);
  godunov_clean(&svf);


  return 0;
}

/******************************************************************************/

void vf_init(vf *svf, int m, int n, double cfl, double xmin, double xmax){

  svf->m = m;
  svf->n = n;
  svf->cfl = cfl;
  svf->xmin = xmin;
  svf->xmax = xmax;

  svf->xi = malloc(m*(n + 2) * sizeof(double));//#GC_MALLOC_ATOMIC((n+2) * sizeof(double));

  svf->dx = (xmax - xmin) / n;

  for (int i = 0; i < n+2; ++i){
    svf->xi[i] = xmin + i * svf->dx - svf->dx/2;
  }

  double t = 0;
  svf->tnow = t;

  svf->wn = malloc(m*(n + 2) * sizeof(double));//GC_MALLOC_ATOMIC( m * (n+2) * sizeof(double));
  svf->wnp1 = malloc(m*(n + 2) * sizeof(double));//GC_MALLOC_ATOMIC( m * (n+2) * sizeof(double));

  for(int i = 0; i < n+2; ++i){
    double *w = svf->wn + i * m;
    double x = svf->xi[i];
    solexact(x,t,w);
    w = svf->wnp1 + i * m;
    solexact(x,t,w);
  }
}


void riemann(double *wL, double *wR, double xi, double *w){

	if (*wL > *wR){ // choc

		if (xi < (*wL + *wR)/2.) {
		  *w = *wL;
		} else {
		  *w = *wR;
		}

	} else { // détente

        if (xi < *wL) {
		  *w = *wL;
      } else if (xi > *wR) {
		  *w = *wR;
      } else {
          *w = xi;
      }
	}
  *w = *wL;

}

void solexact(double x, double t, double *w) // eq.Carac. [Burger]
{
  if (t < 1){ // avant le choc
    if (x < t){
      *w = 1;
    }
    else if (x > 1){
      *w = 0;
    }
    else {
      *w = (1-x)/(1-t);
    }
  }
  else{ // choc
    if (x < 0.5*(1+t)){
      *w = 1;
    }
    else{
      *w = 0;
    }
  }
}

void fphy(double *w, double *flux){

  *flux = (*w) * (*w) /2.;

}

void fnum(double *wL, double *wR, double *flux){

  double xi = 0;
  double w[20]; // constante magique (marge)

  riemann(wL, wR, xi, w);
  fphy(w, flux); // godounov

}


void vf_solve(vf *svf, double tmax){


  int m = svf->m;
  int n = svf->n;

  while (svf->tnow < tmax){

    // calcul du pas de temps
    double vmax = 0;
    for (int i=0; i < n + 2; ++i){
      vmax = vmax > *(svf->wn + i*m) ? vmax : *(svf->wn + i*m);
    }
    svf->dt = svf->dx / vmax * svf->cfl;

    double beta = svf->dt / svf->dx;

    // cells loop
    for(int i = 1; i < n + 1; ++i){

      double flux[m];
      double *wm = svf->wnp1 + i * m;

      // flux en i+1/2
      double *wL = svf->wn + i / m;
      double *wR = svf->wn + (i+1) * m;
      fnum(wL, wR, flux);
      for(int k = 0; k < m; ++k){
				wm[k] += -beta * flux[k];
      }
      // flux en i-1/2
      wL = svf->wn + (i-1) / m;
      wR = svf->wn + i * m;
      fnum(wL, wR, flux);
      for(int k = 0; k < m; ++k){
				wm[k] += beta * flux[k];
      }
    }

    svf->tnow += svf->dt;

    // valeurs aux bords
    int i = 0; // bord gauche
    double* wm = svf->wnp1 + i * m;
    solexact(svf->xi[i], svf->tnow, wm);
    i = svf->n + 1; // bord droite
    wm = svf->wnp1 + i * m;
    solexact(svf->xi[i], svf->tnow, wm);

    for (int j = 0; j < svf->m * (svf->n + 2); ++j){
      svf->wn[j] = svf->wnp1[j];
    }
  }
}

void vf_plot(vf *svf){

  int numvar = 0;

  FILE *plotfile;

  plotfile = fopen("plot.dat","w");

  for(int i = 0; i < svf->n + 2; ++i){
    double x = svf->xi[i];
    double t = svf->tnow;
    double wex[svf->m];

    solexact(x, t, wex);

    double *wnum = svf->wn + i * svf->m;

    fprintf(plotfile, "%f %f %f\n", x, wnum[numvar], wex[numvar]);

  }

  fclose(plotfile);

  system("gnuplot plotcom");

}

void vf_error(vf *svf){

	double errL1;
  int numvar = 0;

  double wex[svf->m];

	FILE *errorfile;
  errorfile = fopen("error.dat","w");

  for(double dx = 0.01; dx > 1e-5; dx = dx/2.){
  	svf->dx = dx;
  	printf("dx = %f\n",svf->dx);
  	svf->tnow = 0;
  	errL1 = 0;
		vf_solve(svf, svf->tmax);
		for(int i = 0; i < svf->n + 2; ++i){
			double x = svf->xi[i];
		  double *wnum = svf->wn + i * svf->m;
		  solexact(x, svf->tmax, wex);
		  errL1 += fabs( wnum[numvar] - wex[numvar] );
		}
		fprintf(errorfile, "%f %f\n", svf->dx, errL1*dx);
  }

	fclose(errorfile);

  system("gnuplot errcom");
}

void godunov_clean(vf *gd)
{
        free(gd->xi);
        free(gd->wn);
        free(gd->wnp1);
}