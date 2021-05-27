#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
//#include <gc.h> // conservative garbage collector

#include "riemann.h"
#include "riemann.c"
#include "VFRoe.h"
#include "roe.h"
#include "roe.c"

#define GC_MALLOC_ATOMIC malloc // without gc

double g = _GPES;

int main(void){

    double xmin = -10;
    double xmax = 10.;

    int m = 2;
    int n = 1000;

    double cfl = 0.5;

    vf svf = {0};
    vf_init(&svf, m, n, cfl, xmin, xmax);

    double tmax = 0.5;
    svf.tmax = tmax;
    vf_solve(&svf, tmax);

    vf_plot(&svf);

    //vf_error(&svf);

  return 0;
}

/******************************************************************************/

void vf_init(vf *svf, int m, int n, double cfl, double xmin, double xmax){

    svf->m = m;
    svf->n = n;
    svf->cfl = cfl;
    svf->xmin = xmin;
    svf->xmax = xmax;

    svf->xi = GC_MALLOC_ATOMIC((n+2) * sizeof(double));

    svf->dx = (xmax - xmin) / n;

    for (int i = 0; i < n+2; ++i){
    svf->xi[i] = xmin + i * svf->dx - svf->dx/2;
    }

    double t = 0;
    svf->tnow = t;

    svf->wn = GC_MALLOC_ATOMIC( m * (n+2) * sizeof(double));
    svf->wnp1 = GC_MALLOC_ATOMIC( m * (n+2) * sizeof(double));

    for(int i = 0; i < n+2; ++i){
        double *w = svf->wn + i * m;
        double x = svf->xi[i];
        solexact(x,t,w);
        w = svf->wnp1 + i * m;
        solexact(x,t,w);
    }
}

void solexact(double x, double t, double *w) {

    double hL = 2;
    double uL = 0;

    double hR = 1;
    double uR = 0;


    // double hL = 1.;
    // double uL = -1.;

    // double hR = 1./4;
    // double uR = uL + 2*sqrt(_GPES) * (sqrt(hL)-sqrt(hR));

    double wL[2] = {hL, hL * uL};
    double wR[2] = {hR, hR * uR};

    double xi = x / (t + _SMALL);

    riem_stvenant(wL, wR, xi, w);
}

double compute_vmax(double *w){

    double g = _GPES;

    double h = w[0];
    double u = w[1]/w[0];

    return fabs(u) + sqrt(g*h);
}

void fphy(double *w, double *flux){
    double g = _GPES;

    double h = w[0];
    double u = w[1] / w[0];

    flux[0] = h*u;
    flux[1] = h*u*u + g*h*h/2;
}

void fnum(double *wL, double *wR, double *flux){

    double xi = 0;
    double w[20]; // 20 = constante magique (marge)

    roe_stvenant(wL, wR, xi, w);
    fphy(w, flux); // VFRoe Scheme

}


void vf_solve(vf *svf, double tmax){

    int m = svf->m;
    int n = svf->n;

    while (svf->tnow < tmax){

        // calcul du pas de temps
        double vmax = 0;
        double v_i = 0;
        for (int i = 0; i < n + 2; ++i){
            double *w = svf->wn + i*m;
            v_i = compute_vmax(w);
            vmax = vmax > v_i ? vmax : v_i;
        }
        svf->dt = svf->dx / vmax * svf->cfl;

        // valeurs aux bords
        int i;
        double wb[2];
        double *wm;

        i = 0; // bord gauche
        wm = svf->wn + (i+1)*m;
        // état miroir
        wb[0] = wm[0];
        wb[1] = - wm[1];
        wm = svf->wn + i*m;
        for(int iv = 0; iv < m; ++iv) wm[iv] = wb[iv];

        i = n + 1; // bord droite
        wm = svf->wn + (i-1) * m;
        wb[0] = wm[0];
        wb[1] = - wm[1];
        wm = svf->wn + i * m;
        for(int iv=0; iv<m; ++iv) wm[iv] = wb[iv];


        double beta = svf->dt / svf->dx;
        // cells loop
        for(int i = 1; i < n + 1; ++i){

            double flux[m];
            double *wm = svf->wnp1 + i*m;

            // flux en i+1/2
            double *wL = svf->wn + i*m;
            double *wR = svf->wn + (i+1) * m;
            fnum(wL, wR, flux);
            for(int k = 0; k < m; ++k){
        	       wm[k] += -beta * flux[k];
            }
            // flux en i-1/2
            wL = svf->wn + (i-1)*m;
            wR = svf->wn + i*m;
            fnum(wL, wR, flux);
            for(int k = 0; k < m; ++k){
    	           wm[k] += beta * flux[k];
            }
        }

        printf("tmax= %f tnow= %f vmax= %f dt= %f \n", tmax, svf->tnow, vmax, svf->dt);

        svf->tnow += svf->dt;

        for (int j = 0; j < (n + 2)*m ; ++j){
            svf->wn[j] = svf->wnp1[j];
        }
    }
}

void vf_plot(vf *svf){

    int numvar = 1;

    FILE *plotfile;

    plotfile = fopen("plot.dat","w");

    for(int i = 0; i < svf->n + 2; ++i){
        double x = svf->xi[i];
        double t = svf->tnow;
        double wex[svf->m];

        solexact(x, t, wex);

        double *wnum = svf->wn + i * svf->m;

    //     fprintf(plotfile, "%f %f %f\n", x, wnum[numvar], wex[numvar]);
    // }

    // fclose(plotfile);

    // system("gnuplot plotcom");
    fprintf(plotfile, "%f %f %f %f %f \n", x, wnum[0], wex[0],wnum[numvar]/ wnum[0], wex[numvar]/wex[0]);
    }

    fclose(plotfile);

    system("gnuplot plotcom_multi");

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
