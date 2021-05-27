#ifndef _VF_H
#define _VF_H


typedef struct vf{

	// nombre de variables conservatives
	int m;

	// nombre de cellules internes (2 ghost cells)
	int n;

	double cfl;

	double dt, dx;

	double xmin, xmax;

	// temps courant et final
	double tnow, tmax;

	// centres des cellules
	double *xi;

	// solution au pas de temps courant
	// et au pas de temps précédent
	double *wn;
	double *wnp1;

} vf;

/******************************************************************************/



// solution exacte - condition initiale
void solexact(double x, double t, double *w);

// solveur de riemann
void riemann(double *wL, double *wR, double xi, double *w);

double vmax(double *w);

// flux physique
void fphy(double *w, double *flux);

// flux numérique
void fnum(double *wL, double *wR, double* flux);

void vf_init(vf *svf, int m, int n, double cfl, double xmin, double xmax);

void vf_solve(vf *svf, double tmax);

void vf_plot(vf *svf);

void vf_error(vf *svf);

void godunov_free(vf *svf);

#endif
