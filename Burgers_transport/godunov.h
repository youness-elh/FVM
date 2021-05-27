#ifndef _GODUNOV_H
#define _GODUNOV_H

// on voudra pouvoir résoudre en général un système de m équations avec un shéma de type Godunov

typedef struct godunov
{
        int m; // nombre de varaibles conservatives

        int N; // nombre de cellules

        double dt, dx; // pas de temps et d'espace
        double tfin; // temps final 
        
        double cfl; // condition de CFL = vmax . dt / dx

        double xmin, xmax; // bornes de l'intervalle

        double *xi; // tableau des centres des cellules

        double *un;   // solution à l'instant n
        double *unp1; // solution à l'instant n+1
} godunov;

// flux numérique
void fluxnum(double *a, double *b, double *flux);

// calcul vmax
double lambda_max(double *u);

// f
void fluxphy(double *w, double *flux);

// solveur de Riemann
void riemann(double *a, double *b, double z, double *w);

void godunov_init(godunov *gd, double xmin, double xmax, double cfl, int m, int N);

void godunov_solve(godunov *gd, double tmax);

void godunov_plot(godunov *gd);

void godunov_error(godunov *gd);

void godunov_clean(godunov *gd);

// solution exacte
void solexacte(double x, double t, double *w);

// St Venant
double Z(double h1, double h2);

double dZ(double h1, double h2);

void riem_stvenant(double *wL, double *wR, double xi, double *w);

void flux_riem(double *wL, double *wR, double *flux);

#endif
