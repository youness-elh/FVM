#include "godunov.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
// #include <omp.h>

#define _C 1 // vitesse de transport

#define _G 9.81

#define MAX(a, b) (((a) > (b)) ? (a) : (b))

/******************************************************************************
 *                              OPTIONS
 *****************************************************************************/

// #define _POOL

// #define _MUSCL

// #define _VFROE

// #define _CORR_VFROE

// #define _RUSANOV



/******************************************************************************
 *                              MAIN
 *****************************************************************************/

int main(int argc, char *argv[])
{
        double xmin = -10;
        double xmax = 10;
        double cfl = 0.5;
        int m = 2;
        int N = 1000;
        double tmax = 2.0;

        extern char * optarg;
        extern int optind;
	int c;
	while( ( c = getopt(argc, argv, "+t:N:") ) != -1 ){
		switch( c ){
			case 't' :
				tmax = atof(optarg);
				break;
			case 'N' :
				N = atoi(optarg);
				break;
			default :
                                fprintf (stderr, "stop \n") ;
				exit(1);
		}
	}

        // printf("start simulation for N = %d and tmax = %.2f \n", N, tmax);

        godunov gd = {0};
        godunov_init(&gd, xmin, xmax, cfl, m, N);
        godunov_solve(&gd, tmax);
        godunov_plot(&gd);

        // double errors[2] = {0.0, 0.0};
        // godunov_error(&gd, errors);
        // printf("erreur sur h : %.6f \nerreur sur q : %.6f \n", errors[0], errors[1]);
        // FILE *file = fopen("vfroe_error.dat", "a");
        // fprintf(file, "%f %f %f\n", gd.dx, errors[0], errors[1]);
        // fclose(file);


        godunov_clean(&gd);
        // printf(">>> done simulation for N = %d and tmax = %.2f \n", N, tmax);

        return 0;
}

void solexacte(double x, double t, double *w) 
{
        // double hL = 2;
        // double uL = 0;

        // double hR = 1;
        // double uR = 0;

        double hL = 1.0; // q7
        double uL = -1.0;

        double hR = 0.25;
        double uR = uL + 2.0 * sqrt(_G) * (sqrt(hL) - sqrt(hR));

        double wL[2] = {hL, hL * uL};
        double wR[2] = {hR, hR * uR};

        double xi = x / (t + 1.e-12);
        
        riem_stvenant(wL, wR, xi, w);

}

double minmod(double a, double b, double c)
{
        double res = 0.0;
        if ((a > 0) && (b > 0) && (c > 0))
                res = a < b ? (a < c ? a : c) : (b < c ? b : c);
        else if ((a < 0) && (b < 0) && (c < 0))
                res = a > b ? (a > c ? a : c) : (b > c ? b : c);
        return res;
}

double Heaviside(double x)
{
        if (x > 0)
                return 1;
        else
                return 0;
}

double Dirac(double x)
{
        return 0;
}

double Z(double hs, double h)
{
        double t0 = 2.0 * sqrt(_G) / (sqrt(hs) + sqrt(h)) * Heaviside(h - hs) + sqrt(2.0) * sqrt(_G) * sqrt(h + hs) / sqrt(h * hs) / 2.0 - sqrt(2.0) * sqrt(_G) * sqrt(h + hs) / sqrt(h * hs) * Heaviside(h - hs) / 2.0;
        return t0;
}

double dZ(double hs, double h)
{
        double t0 = -sqrt(_G) / pow(sqrt(hs) + sqrt(h), 2.0) * Heaviside(h - hs) / sqrt(hs) - 2.0 * sqrt(_G) / (sqrt(hs) + sqrt(h)) * Dirac(-h + hs) + sqrt(2.0) * sqrt(_G) / sqrt(h + hs) / sqrt(h * hs) / 4.0 - sqrt(2.0) * sqrt(_G) * sqrt(h + hs) / sqrt(h * h * h * hs * hs * hs) * h / 4.0 - sqrt(2.0) * sqrt(_G) / sqrt(h + hs) / sqrt(h * hs) * Heaviside(h - hs) / 4.0 + sqrt(2.0) * sqrt(_G) * sqrt(h + hs) / sqrt(h * h * h * hs * hs * hs) * Heaviside(h - hs) * h / 4.0 + sqrt(2.0) * sqrt(_G) * sqrt(h + hs) / sqrt(h * hs) * Dirac(-h + hs) / 2.0;
        return t0;
}

double lambda_max(double *w) 
{
        // w = (h,hu)
        double h = w[0];
        double u = w[1] / w[0];
        return fabs(u) + sqrt(_G * h);
}

void fluxphy(double *w, double *flux) // flux de St Venant
{
        flux[0] = w[1];
        flux[1] = pow(w[1], 2.0) / w[0] + _G * pow(w[0], 2.0) / 2.0;
}

void riem_stvenant(double *wL, double *wR, double xi, double *w)
{

        double hL = wL[0];
        double uL = wL[1] / wL[0];
        double hR = wR[0];
        double uR = wR[1] / wR[0];

        double hs = 1e-6;
        int itermax = 10;

        for (int it = 0; it < itermax; it++)
        {
                double f = uL - (hs - hL) * Z(hs, hL) -
                           uR - (hs - hR) * Z(hs, hR);
                double df = -(hs - hL) * dZ(hs, hL) -
                            Z(hs, hL) -
                            (hs - hR) * dZ(hs, hR) -
                            Z(hs, hR);
                double dhs = -f / df;

                hs += dhs;

                // printf("it=%d f=%e df=%e hs=%e dhs=%e\n", it, f, df, hs, dhs);
        }

        double us = uL - (hs - hL) * Z(hs, hL);

        double v1m, v1p, v2m, v2p;

        // 1-onde
        if (hs < hL)
        {
                v1m = uL - sqrt(_G * hL);
                v1p = us - sqrt(_G * hs);
        }
        else
        {
                double a = sqrt(hs) / (sqrt(hs) + sqrt(hL));
                double u = a * us + (1 - a) * uL;
                double h = (hs + hL) / 2;
                v1m = u - sqrt(_G * h);
                v1p = v1m;
        }

        // 2 onde
        if (hs < hR)
        {
                v2m = us + sqrt(_G * hs);
                v2p = uR + sqrt(_G * hR);
        }
        else
        {
                double a = sqrt(hs) / (sqrt(hs) + sqrt(hR));
                double u = a * us + (1 - a) * uR;
                double h = (hs + hR) / 2;
                v2m = u + sqrt(_G * h);
                v2p = v2m;
        }

        //printf("v=%f %f %f %f\n hs=%f us=%f\n", v1m,v1p,v2m,v2p, hs,us);

        if (xi < v1m)
        {
                w[0] = wL[0];
                w[1] = wL[1];
        }
        else if (xi < v1p)
        {
                double u = (uL + 2 * xi + 2 * sqrt(_G * hL)) / 3;
                double h = (u - xi) * (u - xi) / _G;
                w[0] = h;
                w[1] = h * u;
        }
        else if (xi < v2m)
        {
                w[0] = hs;
                w[1] = hs * us;
        }
        else if (xi < v2p)
        {
                double u = (uR + 2 * xi - 2 * sqrt(_G * hR)) / 3;
                double h = (u - xi) * (u - xi) / _G;
                w[0] = h;
                w[1] = h * u;
        }
        else
        {
                w[0] = wR[0];
                w[1] = wR[1];
        }
}

void riem_vfroe(double *wL, double * wR, double xi, double *w)
{
        double hL = wL[0];
        double uL = wL[1] / wL[0];
        double hR = wR[0];
        double uR = wR[1] / wR[0];

        double hb = (hL + hR) / 2.0;
        double ub = (uL + uR) / 2.0;
        double cb = sqrt(_G * hb);

        double lambda_1 = ub - cb; 
        double lambda_2 = ub + cb; 

        // double lambda_1L = uL - sqrt(_G * hL);
        // double lambda_2L = uL + sqrt(_G * hL);
        // double lambda_1R = uR - sqrt(_G * hR);
        // double lambda_2R = uR + sqrt(_G * hR);

        
        double h, u;
        if ((lambda_1 > 0.0) && (lambda_2 > 0.0)) // y(0,t) = y_L
        {
                h = hL;
                u = uL;
        }
        else if((lambda_1 < 0.0) && (lambda_2 < 0.0)) // y(0,t) = y_R
        {
                h = hR;
                u = uR;
        }
        else // y(0,t) = y^*
        {
                h = (hL + hR) / 2.0 - hb * (uR - uL) / 2.0 / cb;
                u = (uL + uR) / 2.0 - _G * (hR - hL) / 2.0 / cb;
        }
        // double epsL = 0.0;
        // double epsR = 0.0;
        // if (lambda_1L < 0.0 && lambda_1R > 0.0)
        // {
                // epsL = lambda_1R - lambda_1L;
                // printf("%f %f \n", lambda_1L, lambda_1R);
        // }
        // if (lambda_2L < 0.0 && lambda_2R > 0.0)
        // {
                // epsR = lambda_2R - lambda_2L;
        // }
        // double eps = 0.0; // 1.0; //MAX(epsL, epsR);
        // h = h + eps * (hR - hL) / 2.0;
        // u = u + eps * (uR - uL) / 2.0;
        w[0] = h;
        w[1] = h * u;
}

#if !defined(_VFROE) && !defined(_RUSANOV) && !defined(_CORR_VFROE)

void fluxnum(double *wL, double *wR, double *flux)
{
        double w[2]; // (h, hu)
        riem_stvenant(wL, wR, 0., w);
        fluxphy(w, flux);
}

#elif defined(_VFROE)

void fluxnum(double *wL, double *wR, double *flux)
{
        double w[2]; // (h, hu)
        riem_vfroe(wL, wR, 0., w);
        fluxphy(w, flux);
}

#elif defined(_RUSANOV)

void fluxnum(double *wL, double *wR, double *flux)
{ 
        double lamb = MAX(lambda_max(wL), lambda_max(wR));
        flux[0] = 0.5 * (wL[1] + wR[1]) - 0.5 * lamb * (wR[0] - wL[0]);
        flux[1] = 0.5 * (pow(wL[1], 2.0) / wL[0] + _G * pow(wL[0], 2.0) / 2.0 + pow(wR[1], 2.0) / wR[0] + _G * pow(wR[0], 2.0) / 2.0) - 0.5 * lamb * (wR[1] - wL[1]);
}

#elif defined(_CORR_VFROE)

void fluxnum(double *wL, double *wR, double *flux)
{
        double w[2]; 
        double lamb;

        double hL = wL[0];
        double uL = wL[1] / wL[0];
        double hR = wR[0];
        double uR = wR[1] / wR[0];

        double lambda_1L = uL - sqrt(_G * hL);
        double lambda_2L = uL + sqrt(_G * hL);
        double lambda_1R = uR - sqrt(_G * hR);
        double lambda_2R = uR + sqrt(_G * hR);

        if (((lambda_1L < 0.0) && (lambda_1R > 0.0)) || ((lambda_2L < 0.0) && (lambda_2R > 0.0)))
        {
                lamb = MAX(lambda_max(wL), lambda_max(wR));
                flux[0] = 0.5 * (wL[1] + wR[1]) - 0.5 * lamb * (wR[0] - wL[0]);
                flux[1] = 0.5 * (pow(wL[1], 2.0) / wL[0] + _G * pow(wL[0], 2.0) / 2.0 + pow(wR[1], 2.0) / wR[0] + _G * pow(wR[0], 2.0) / 2.0) - 0.5 * lamb * (wR[1] - wL[1]);
        }
        else
        {
                riem_vfroe(wL, wR, 0., w);
                fluxphy(w, flux);
        }
}

#else
#error Unknown choice
#endif

void godunov_init(godunov *gd, double xmin, double xmax, double cfl, int m, int N)
{
        gd->xmin = xmin;
        gd->xmax = xmax;
        gd->m = m;
        gd->N = N;
        gd->cfl = cfl;
        gd->dx = (xmax - xmin) / N;
        gd->dt = 0; // provisoire

        gd->xi = malloc((N + 2) * sizeof(double));
        gd->un = malloc((N + 2) * sizeof(double) * m);
        gd->unp1 = malloc((N + 2) * sizeof(double) * m);

        double t = 0;
        for (int i = 0; i < N + 2; i++)
        {
                gd->xi[i] = xmin + gd->dx / 2 + (i - 1) * gd->dx; // ctrl K ctrl F met au propre
                // un = [..., wi,0 , wi,1 , ... , wi,m-1,
                //            wi+1,0,   ...   , wi+1,m-1 ,
                //                              ...]
                solexacte(gd->xi[i], t, gd->un + i * m);
        }
}

void godunov_solve(godunov *gd, double tmax)
{
        double tnow = 0;
        int m = gd->m;
        assert(m == 2);
#ifdef _MUSCL
        double *si = calloc(m * (gd->N + 2), sizeof(double));
        double *ri = calloc(m * (gd->N + 2), sizeof(double));
#endif
        while (tnow < tmax)
        {
                double vmax = 0;
                // calcul vitesse max
                for (int i = 0; i < gd->N + 2; i++)
                {
                        double vloc = lambda_max(gd->un + m * i);
                        vmax = vmax > vloc ? vmax : vloc;
                }
                gd->dt = gd->cfl * gd->dx / vmax;
#ifdef _MUSCL
                // calcul des pentes
                for(int i = 1; i < gd->N + 1; ++i)
                {
                        for(int iv = 0; iv < m; ++iv)
                        {
                                double a = (gd->un[i * m + iv] - gd->un[(i - 1) * m + iv]) / gd->dx;
                                double b = (gd->un[(i + 1) * m + iv] - gd->un[i * m + iv]) / gd->dx;
                                double c = (gd->un[(i + 1) * m + iv] - gd->un[(i - 1) * m + iv]) / 2.0 / gd->dx;
                                si[i * m + iv] = minmod(a, b, c);
                        }
                        ri[i * m] = -si[i * m + 1];
                        ri[i * m + 1] = (pow(gd->un[i * m + 1], 2.0) / pow(gd->un[i * m], 2.0) - _G * gd->un[i * m]) * si[i * m] 
                                        - 2.0 * gd->un[i * m + 1] / gd->un[i * m] * si[i * m + 1];
                }
#endif
                for (int i = 1; i < gd->N + 1; i++)
                {
                        
                        double flux[gd->m];
                        double wL[m];
                        double wR[m];
                        for (int iv = 0; iv < m; ++iv)
                        {
#ifdef _MUSCL
                                wL[iv] = gd->un[i * m + iv] + si[i * m + iv] * gd->dx / 2.0 + ri[i * m + iv] * gd->dt / 2.0;
                                wR[iv] = gd->un[(i + 1)*m + iv] - si[(i + 1) * m + iv] * gd->dx / 2.0 + ri[(i + 1) * m + iv] * gd->dt / 2.0;
#else
                                wL[iv] = gd->un[i * m + iv];
                                wR[iv] = gd->un[(i + 1)*m + iv];
#endif
                        }
                        
                        fluxnum(wL, wR, flux); // flux droite
                        for (int iv = 0; iv < m; iv++)
                        {
                                gd->unp1[i * m + iv] = gd->un[i * m + iv] - gd->dt * flux[iv] / gd->dx;
                        }

                        for (int iv = 0; iv < m; ++iv)
                        {
#ifdef _MUSCL
                                wL[iv] = gd->un[(i - 1) * m + iv] + si[(i - 1)*m + iv] * gd->dx / 2.0 + ri[(i - 1) * m + iv] * gd->dt / 2.0;
                                wR[iv] = gd->un[i * m + iv] - si[i * m + iv] * gd->dx / 2.0 + ri[i*m+iv] * gd->dt / 2.0;
#else
                                wL[iv] = gd->un[(i - 1) * m + iv];
                                wR[iv] = gd->un[i * m + iv];
#endif
                        }
                        
                        fluxnum(wL, wR, flux); // flux gauche
                        for (int iv = 0; iv < m; iv++)
                        {
                                gd->unp1[i * m + iv] += gd->dt * flux[iv] / gd->dx;
                        }
                        
                }

                // mise à jour
                tnow += gd->dt;
#ifdef _POOL
                // conditions aux limites
                int i = 0; // à gauche
                gd->unp1[i * m] = gd->unp1[(i + 1) * m]; 
                gd->unp1[i * m + 1] = - gd->unp1[(i + 1) * m + 1]; 

                i = gd->N + 1; // à droite
                gd->unp1[i * m] = gd->unp1[(i - 1) * m];
                gd->unp1[i * m + 1] = - gd->unp1[(i - 1) * m  + 1]; 
#else
                // conditions aux limites
                int i = 0; // à gauche
                solexacte(gd->xi[i], tnow, gd->unp1 + i * m);

                i = gd->N + 1; // à droite
                solexacte(gd->xi[i], tnow, gd->unp1 + i * m);
#endif
                memcpy(gd->un, gd->unp1, (gd->N + 2) * m * sizeof(double));
        }
        gd->tfin = tnow;
#ifdef _MUSCL
        free(si);
        free(ri);
#endif
}

void godunov_error(godunov *gd)
{
        double errors[2] = {0.0, 0.0};
        errors[0] = errors[1] = 0.0;
        for (int i = 1; i < gd->N + 1; i++)
        {
                double uex[gd->m];
                solexacte(gd->xi[i], gd->tfin, uex);
                errors[0] += gd->dx * fabs(uex[0] - gd->un[i * gd->m + 0]);
                errors[1] += gd->dx * fabs(uex[1] - gd->un[i * gd->m + 1]);
        }
}

void godunov_plot(godunov *gd)
{
        FILE *fic = fopen("godu.dat", "w");
        for (int i = 0; i < gd->N + 2; i++)
        {
                double h = gd->un[i * gd->m];
                double u = gd->un[i * gd->m + 1] / gd->un[i * gd->m];

// #ifndef _POOL
                double uex[gd->m]; // (h, hu)^T
                solexacte(gd->xi[i], gd->tfin, uex);
                double h_ex = uex[0];
                double u_ex = uex[1] / uex[0];
                fprintf(fic, "%f %f %f %f %f\n", gd->xi[i], h_ex, h, u_ex, u);
                // fprintf(fic, "%f %f %f %f %f\n", t, h_ex, h, u_ex, u);
// #else
                // fprintf(fic, "%f %f %f\n", gd->xi[i], h, u);
// #endif
        }

        fclose(fic);

// #ifdef _POOL
//         int status = system("gnuplot plotcom_pool -persist");
// #else
        int status = system("gnuplot plotcom -persist");
// #endif
        assert(status == EXIT_SUCCESS);
}

void godunov_clean(godunov *gd)
{
        free(gd->xi);
        free(gd->un);
        free(gd->unp1);
}
