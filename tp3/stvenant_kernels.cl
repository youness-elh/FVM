#define _NX _nx_
#define _NY _ny_
#define _DX _dx_
#define _DY _dy_
#define _DT _dt_
#define _M _m_

#define double float

#define _LAMBDA _lambda_
#define SMALL 1e-6

#define _G (9.81f)

#define _VOL (_DX * _DY)

__constant int dir[4][2] = { {1, 0}, {-1, 0},
			     {0, 1}, {0, -1}}; // E, W, N, S

__constant float ds[4] = { _DY, _DY, _DX, _DX };

///////////////////////////////////////////////////////

double Z(double hs, double h){

  double g = _G;

  if (hs < h)
    return 2*sqrt(g)/(sqrt(hs)+sqrt(h));

  return sqrt(g*(hs+h)/(2*hs*h));

}

double dZ(double hs, double h){

  double g = _G;

  if (hs < h) /* choc */
    return - sqrt(g) / sqrt(hs) / (sqrt(hs)+sqrt(h)) / (sqrt(hs)+sqrt(h));
  /* détente */
  return - sqrt(g * hs * h / 2 / (hs + h)) / 2 / hs / hs;

}

void fluxphy(float *w, float* n, float *flux){

  float h = w[0];
  float u = w[1] / h;
  float v = w[2] / h;

  float un = u * n[0] + v * n[1];

  flux[0] = h * un;
  flux[1] = h * u * un + _G * h * h / 2 * n[0];
  flux[2] = h * v * un + _G * h * h / 2 * n[1];

}

//Rusanov
void fluxnum(float *wL, float *wR, float* vnorm, float* flux){
  float fL[_M];
  float fR[_M];

  fluxphy(wL, vnorm, fL);
  fluxphy(wR, vnorm, fR);
  // Rusanov, TODO : Godunov (comparer les résulats)
  for(int iv = 0; iv < _M; iv++){
    flux[iv] = 0.5f * (fL[iv] + fR[iv]) - 0.5f * _LAMBDA * (wR[iv] - wL[iv]);
  }

}

// condition initiale
void exact_sol(float* xy, float t, float* w){
  float h = 1;
  float x = xy[0];// - 0.5f;
  float y = xy[1];// - 0.5f;

  if (x < 0.75f && x > 0.25f && y < 0.75f && y > 0.25f) h = 2;
  //if (x * x + y * y < 0.2f * 0.2f) h = 2;
  w[0] = h;
  w[1] = 0;
  w[2] = 0;
}

// ---

void riem_stvenant(double *wL, double *wR, double xi, double *w)
{
  double hL = wL[0];
  double uL = wL[1] / wL[0];

  double hR = wR[0];
  double uR = wR[1] / wR[0];

  double hs = SMALL;
  int itermax = 10;

  for (int i = 0; i < itermax; ++i)
  {
    double f = uL - (hs - hL) * Z(hs, hL) - uR - (hs - hR) * Z(hs, hR);
    double df = -(hs - hL) * dZ(hs, hL) - Z(hs, hL) - (hs - hR) * dZ(hs,hR) - Z(hs,hR);
    double dhs = -f / df;
    hs += dhs;
    // printf("i=%d f=%e df=%e hs=%e dhs=%e \n", i, f, df, hs, dhs);
  }

  double us = uL - (hs - hL) * Z(hs, hL);

  double v1m, v1p, v2m, v2p;

  // 1 - onde
  if (hs < hL) // détente
  {
    v1m = uL - sqrt(_G*hL);
    v1p = us - sqrt(_G*hs);

  } else { // choc

    double a = sqrt(hs) / (sqrt(hs) + sqrt(hL));
    double u = a * us + (1 - a) * uL;
    double h = (hs + hL) / 2;

    v1m = u - sqrt(_G*h);
    v1p = v1m;

  }

  // 2 - onde
  if (hs < hR) // détente
  {
    v2m = us + sqrt(_G*hs);
    v2p = uR + sqrt(_G*hR);

  } else { // choc

    double a = sqrt(hs) / (sqrt(hs) + sqrt(hR));
    double u = a * us + (1 - a) * uR;
    double h = (hs + hR) / 2;

    v2m = u + sqrt(_G*h);
    v2p = v2m;

  }

  // printf("v = %f %f %f %f\nhs=%f us=%f\n", v1m, v1p, v2m, v2p, hs, us);

  if (xi < v1m){

    w[0] = wL[0];
    w[1] = wL[1];

  } else if (xi < v1p) {

    double u = (uL + 2*xi + 2*sqrt(_G * hL)) / 3;
    double h = (u - xi) * (u - xi) / _G;

    w[0] = h;
    w[1] = h * u;

  } else if (xi < v2m) {

    w[0] = hs;
    w[1] = hs * us;

  } else if (xi < v2p) {

    double u = (uR + 2*xi - 2*sqrt(_G*hR)) / 3;
    double h = (u - xi) * (u - xi) / _G;

    w[0] = h;
    w[1] = h * u;

  } else {

    w[0] = wR[0];
    w[1] = wR[1];

  }

}

// Godunov
// void fluxnum(float *wL, float *wR, float* vnorm, float* flux){
//   float fL[_M];
//   float fR[_M];
// 	float* w;
// 	riem_stvenant(wL, wR, 0, w);
//   fluxphy(w, vnorm,flux);
// }

void flux_riem_2d(double *wL, double *wR, double *vnorm, double *flux)
{
    double qnL = wL [1] * vnorm [0] + wL [2] * vnorm [1];
    double qnR = wR [1] * vnorm [0] + wR [2] * vnorm [1];

    double qtL = -wL [1] * vnorm [1] + wL [2] * vnorm [0];
    double qtR = -wR [1] * vnorm [1] + wR [2] * vnorm [0];

    double vL [2] = {wL [0], qnL};
    double vR [2] = {wR [0], qnR};

    double v [2];
    double xi = 0;

    riem_stvenant (vL, vR, xi, v);

    double un = v [1] / v [0];
    double ut;

    ut = (un > 0) ? qtL / wL [0] : qtR / wR [0];

    double qn = v [1];
    double qt = ut * v [0];

    double w [3];
    w [0] = v [0];
    w [1] = qn * vnorm [0] - qt * vnorm [1];
    w [2] = qn * vnorm [1] + qt * vnorm [0];

    // puis calcul du flux physique (appliqué à w)
	fluxphy (w, vnorm, flux);
}

// ---

// ----- OPENCL -----
// WARNING : access to global memory (wn) takes time !

__kernel void init_sol(__global  float *wn){

  int id = get_global_id(0);

  int i = id % _NX;
  int j = id / _NX;

  int ngrid = _NX * _NY;

  float wnow[_M];

  float t = 0;
  float xy[2] = {i * _DX + _DX / 2, j * _DY + _DY / 2};

  exact_sol(xy, t, wnow);
  //printf("x=%f, y=%f \n",xy[0],xy[1]);
  // load middle value
  for(int iv = 0; iv < _M; iv++){ // iv-th macrocopic variable
    int imem = i + j * _NX + iv * ngrid; // flattened
    wn[imem] =  wnow[iv];
    // boundary values
    //if (i == 0 || i == _NX - 1 || j == 0 || j == _NY - 1)
    // wn[imem] = _WBORD;
  }

}

__kernel void time_step(__global  float *wn, __global float *wnp1){

  int id = get_global_id(0);

  int i = id % _NX;
  int j = id / _NX;


  int ngrid = _NX * _NY;

  float wnow[_M];
  float wnext[_M];

  // load middle value
  for(int iv = 0; iv < _M; iv++){
    int imem = i + j * _NX + iv * ngrid;
    wnow[iv] = wn[imem];
    wnext[iv] = wnow[iv];
   }

  // only compute internal cells
  if (i > 0 && i < _NX - 1 && j > 0 && j < _NY - 1){
    float flux[_M];

    // loop on all directions:
    // idir = 0 (east)
    // idir = 1 (west)
    // idir = 2 (north)
    // idir = 3 (south)
    for(int idir = 0; idir < 4; idir++){
      float vn[2];
      vn[0] = dir[idir][0];
      vn[1] = dir[idir][1];
      int iR = i + dir[idir][0];
      int jR = j + dir[idir][1];
      // load neighbour values
      float wR[_M];
      for(int iv = 0; iv < _M; iv++){
	int imem = iv * ngrid + iR + jR * _NX;
	wR[iv] = wn[imem];
      }
      if (iR == 0 || iR == _NX - 1 || jR == 0 || jR == _NY - 1){
	wR[0] = wnow[0];
	float qn = wnow[1] * vn[0] + wnow[2] * vn[1]; // qte de mvt
	wR[1] = wnow[1] - 2 * qn * vn[0]; // bounce
	wR[2] = wnow[2] - 2 * qn * vn[1]; // bounce
      }
      fluxnum(wnow, wR, vn, flux);

      // time evolution
      for(int iv = 0; iv < _M; iv++){
	wnext[iv] -= _DT * ds[idir] / _VOL * flux[iv]; // VF
      }
    }
  }   // end test for internal cells

  // copy wnext into global memory
  // (including boundary cells)

  for(int iv = 0; iv < _M; iv++){
    int imem = iv * ngrid + i + j * _NX;
    wnp1[imem] = wnext[iv];
  }


}
