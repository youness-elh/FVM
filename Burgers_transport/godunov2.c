#include "godunov.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
// #include <omp.h>

#define _C 1 // vitesse de transport

/******************************************************************************
 *                              OPTIONS
 *****************************************************************************/

// #define _TRANSPORT // exercice 1 TP 1

#define _BURGERS // exercice 2 TP 1

//#define _ERROR_FILE

// #define _ORDER_EST

//#define _MUSCL

/******************************************************************************
 *                              MAIN
 *****************************************************************************/

int main(int argc, char *argv[])
{
#ifdef _TRANSPORT
        double xmin = 0;
        double xmax = 1;

        double cfl = 0.99;
        int m = 1;
        int N = 1000;

        if (argc > 1)
                N = atoi(argv[1]);

        double tmax = 0.5;
#endif

#ifdef _BURGERS
        double xmin = -1;
        double xmax = 2;
        double cfl = 0.5;
        int m = 1;
        int N = 10;

        if (argc > 1)
                N = atoi(argv[1]);

        double tmax = 2;
#endif

        godunov gd = {0};
        godunov_init(&gd, xmin, xmax, cfl, m, N);
        godunov_solve(&gd, tmax);
        double erreurN = godunov_error(&gd);

#ifdef _ERROR_FILE
        FILE *fic2 = fopen("godu_error.dat", "a");
        fprintf(fic2, "%f %f\n", gd.dx, erreurN);
        fclose(fic2);
#else
        godunov_plot(&gd);
        printf("Erreur L1 = %f\n", erreurN);
#endif

        godunov_clean(&gd);

#ifdef _ORDER_EST
        N *= 2;
        godunov_init(&gd, xmin, xmax, cfl, m, N);
        godunov_solve(&gd, tmax);
        double erreur2N = godunov_error(&gd);
        godunov_clean(&gd);

        double order = log(erreurN / erreur2N) / log(2);
        printf("ordre estimé: %f\n", order);
#endif
      
        return 0;
}

// gcc godunov.c -lm //* -lm pour avoir la librairie math si besoin

/******************************************************************************
 *                              TRANSPORT
 *****************************************************************************/

#ifdef _TRANSPORT
void fluxphy(double *w, double *flux)
{
        flux[0] = _C * w[0];
}

void riemann(double *a, double *b, double z, double *w)
{
        if (z < _C)
        {
                w[0] = a[0];
        }
        else
        {
                w[0] = b[0];
        }
}

void solexacte(double x, double t, double *w)
{
        double uL = exp(-t + x / _C);
        double uR = 0;

        if (x < _C * t)
        {
                w[0] = uL;
        }
        else
        {
                w[0] = uR;
        }
}
#endif

/******************************************************************************
 *                              BURGERS
 *****************************************************************************/

#ifdef _BURGERS
void fluxphy(double *w, double *flux)
{
        flux[0] = pow(w[0], 2.0) / 2.0;
}

void riemann(double *a, double *b, double z, double *w) // z = x/t, a = u_L, b = u_R
{
        if (a[0] < b[0]) // u_L < u_R
        {
                if (z < a[0]) // x/t < u_L
                {
                        w[0] = a[0];
                }
                else if (z > b[0]) // x/t > u_R
                {
                        w[0] = b[0];
                }
                else
                {
                        w[0] = z;
                }
        }
        else
        {
                if (z < (a[0] + b[0]) / 2.0) // x/t < (u_L+u_R)/2
                {
                        w[0] = a[0];
                }
                else
                {
                        w[0] = b[0];
                }
        }
}

void solexacte(double x, double t, double *w)
{
        if (t < 1.0)
        {
                if (x < t)
                {
                        w[0] = 1.0;
                }
                else if (x > 1.0)
                {
                        w[0] = 0.0;
                }
                else
                {
                        w[0] = (1.0 - x) / (1.0 - t);
                }
        }
        else
        {
                if (x < 0.5 * t + 0.5)
                {
                        w[0] = 1.0;
                }
                else
                {
                        w[0] = 0.0;
                }
        }
}
#endif

/******************************************************************************
 *                              COMMUN
 *****************************************************************************/

double minmod(double a, double b, double c)
{
        double res = 0.0;
        if ((a > 0) && (b > 0) && (c > 0))
        {
                res = a < b ? (a < c ? a : c) : (b < c ? b : c);
        }
        else if ((a < 0) && (b < 0) && (c < 0))
        {
                res = a > b ? (a > c ? a : c) : (b > c ? b : c);
        }
        return res;
        // if ((a > 0) && (b > 0) && (c > 0))
        // {
        //         double ab = a < b ? a : b;
        //         return ab < c ? ab : c;
        // }
        // else if ((a < 0) && (b < 0) && (c < 0))
        // {
        //         double ab = a > b ? a : b;
        //         return ab > c ? ab : c;
        // }
        // else
        // {
        //         return 0;
        // }
}

double lambda_max(double *u)
{
        return _C;
        // return u[0];
}

void fluxnum(double *a, double *b, double *flux)
{
        double w[10];
        riemann(a, b, 0., w);
        fluxphy(w, flux);
}

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
        assert(m == 1);
        double *si = calloc(gd->N + 2, sizeof(double));
        double *ri = calloc(gd->N + 2, sizeof(double));

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

                // calcul de pentes
                for (int i = 1; i < gd->N + 1; i++)
                {
                        // pour MUSCL
                        double a = (gd->un[i] - gd->un[i - 1]) / gd->dx;
                        double b = (gd->un[i + 1] - gd->un[i]) / gd->dx;
                        double c = (gd->un[i + 1] - gd->un[i - 1]) / 2.0 / gd->dx;

#ifdef _MUSCL
                        si[i] = minmod(a,b,c);
#else
                        si[i] = 0.0; 
#endif
                        ri[i] = -gd->un[i] * si[i]; // Burgers
                }
                

                // fluxs
                for (int i = 1; i < gd->N + 1; i++)
                {
                        double flux[gd->m];
                        double wL[m];
                        double wR[m];

                        wL[0] = gd->un[i] + si[i] * gd->dx / 2.0 + ri[i] * gd->dt / 2.0;
                        wR[0] = gd->un[i + 1] - si[i + 1] * gd->dx / 2.0 + ri[i + 1] * gd->dt / 2.0;

                        fluxnum(wL, wR, flux); // flux droite
                        for (int iv = 0; iv < m; iv++)
                        {
                                gd->unp1[i * m + iv] = gd->un[i * m + iv] - gd->dt * flux[iv] / gd->dx;
                        }

                        wL[0] = gd->un[i - 1] + si[i - 1] * gd->dx / 2.0 + ri[i - 1] * gd->dt / 2.0;
                        wR[0] = gd->un[i] - si[i] * gd->dx / 2.0 + ri[i] * gd->dt / 2.0;


                        fluxnum(wL, wR, flux); // flux gauche
                        for (int iv = 0; iv < m; iv++)
                        {
                                gd->unp1[i * m + iv] += gd->dt * flux[iv] / gd->dx;
                        }
                }

                // mise à jour
                tnow += gd->dt;

                // conditions aux limites
                int i = 0; // à gauche
                solexacte(gd->xi[i], tnow, gd->unp1 + i * m);

                i = gd->N + 1; // à droite
                solexacte(gd->xi[i], tnow, gd->unp1 + i * m);

                memcpy(gd->un, gd->unp1, (gd->N + 2) * m * sizeof(double));
        }
        gd->tfin = tnow;

        free(si);
        free(ri);
}

double godunov_error(godunov *gd)
{
        double erreur = 0;
        for (int i = 1; i < gd->N + 1; i++)
        {
                double uex[gd->m];
                solexacte(gd->xi[i], gd->tfin, uex);
                int iv = 0;
                erreur += gd->dx * fabs(uex[iv] - gd->un[i * gd->m + iv]);
        }
        return erreur;
}

void godunov_plot(godunov *gd)
{
        FILE *fic = fopen("godu.dat", "w");
        for (int i = 0; i < gd->N + 2; i++)
        {
                double uex[gd->m];
                solexacte(gd->xi[i], gd->tfin, uex);
                int iv = 0; // nombre de variables

                fprintf(fic, "%f %f %f\n", gd->xi[i], uex[iv], gd->un[i * gd->m + iv]);
        }

        fclose(fic);
        int status = system("gnuplot plotcom -persist");
        assert(status == EXIT_SUCCESS);
}

void godunov_clean(godunov *gd)
{
        free(gd->xi);
        free(gd->un);
        free(gd->unp1);
}