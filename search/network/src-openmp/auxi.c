#ifndef __SETTINGS_H__
#define __SETTINGS_H__

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <complex.h>
#include <float.h>

#include "auxi.h"

#include <omp.h>


// lin2ast described in Phys. Rev. D 82, 022005 (2010) (arXiv:1003.0844)
void lin2ast (double be1, double be2, int pm, double sepsm, double cepsm, 
	      double *sinal, double *cosal, double *sindel, double *cosdel) {

  *sindel = be1*sepsm-(2*pm-3)*sqrt(1.-sqr(be1)-sqr(be2))*cepsm;
  *cosdel = sqrt(1.-sqr(*sindel));
  *sinal = (be1-sepsm*(*sindel))/(cepsm*(*cosdel));
  *cosal = be2/(*cosdel);

} /* lin2ast() */


int ast2lin (FLOAT_TYPE alfa, FLOAT_TYPE delta, double epsm, double *be) {

  /* alfa - right ascension [rad]
     delta - declination [rad]
     Mean obliquity of the equator with respect to the ecliptic at J2000.0:
     epsm =  84381.448*pi/(3600*180)
  */

    be[0] = cos(epsm)*sin(alfa)*cos(delta)+sin(epsm)*sin(delta);
    be[1] = cos(alfa)*cos(delta);

    //#mb this is not needed at the moment 
 
    double d1 = asin(be[0]*sin(epsm) 
            + sqrt(1. - be[0]*be[0] - be[1]*be[1])*cos(epsm)) - delta;

//  double d2 = asin(be[0]*sin(epsm) 
//          - sqrt(1. - be[0]*be[0] - be[1]*be[1])*cos(epsm)) - delta;

    int pm; 

    if(fabs(d1)  < 10.*DBL_EPSILON)
        pm = 1;
    else
        pm = 2;

    return pm; 
} /* ast2lin */



inline void spline(complex double *y, int n, complex double *y2) {

  int i, k;
  complex double invp, qn, un;
  static complex double *u = NULL;

  if (!u) u = (complex double *)malloc((n-1)*sizeof(complex double));
  y2[0] = u[0] = 0.;

  for (i=1; i<n-1; ++i) {
    //p = .5*y2[i-1]+2.;
    //y2[i] = -.5/p;
    //u[i] = y[i+1]-2.*y[i]+y[i-1];
    //u[i] = (3.*u[i]-.5*u[i-1])/p;
    invp = 2./(y2[i-1]+4.);
    y2[i] = -.5*invp;
    u[i] = y[i-1]-2.*y[i]+y[i+1];
    u[i] = (-.5*u[i-1]+3.*u[i])*invp;
  }
  qn = un = 0.;
  y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.);
  for (k=n-2; k>=0; --k)
    y2[k] = y2[k]*y2[k+1]+u[k];
} /* spline() */


inline complex double splint (complex double *ya, complex double *y2a, int n, double x) {

  int klo, khi;
  double b, a;

  if (x<0 || x>n-1)
    return 0.;
  klo = floor (x);
  khi = klo+1;
  a = khi - x;
  b = x - klo;
  return a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])/6.0;
} /* splint() */


void splintpad (complex double *ya, double *shftf, int N, int interpftpad,
	   complex double *out) {
  /* Cubic spline with "natural" boundary conditions.
     Input:
     ya[i] - value of the function being interpolated in x_i = i,
     for i = 0 .. (interpftpad*N-1)	(changed on exit);
     Interpolating spline will be calculated at the points
     interpftpad*(i-shftf[i]), for i = 0 .. (N-1);
     N - number of output data points.
     Output:
     out[i] - value of the interpolating function
     at interpftpad*(i-shftf[i]).
  */
  complex double *y2;
  double x;
  int i;

  y2 = (complex double *) malloc (interpftpad*N*sizeof (complex double)); //vector twice-size of N
  spline (ya, interpftpad*N, y2);
#if defined(_OPENMP)
#pragma omp parallel default(shared) private(x)
#endif
  {
#if defined(_OPENMP)         
#pragma omp for schedule(static)
#endif       
    for (i=0; i<N; ++i) {
      x = interpftpad*(i-shftf[i]);
      out[i] = splint (ya, y2, interpftpad*N, x);
    } /* for i */
  }
  free (y2);
} /* splintpad */


// pci test
void linterp (complex double *ya, double *shftf, int N, int interpftpad,
	      complex double *out) {
     /* linear interpolation
	Input:
	ya[i] - value of the function being interpolated in x_i = i,
	for i = 0 .. (interpftpad*N-1)
	Interpolating spline will be calculated at the points interpftpad*(i-shftf[i]), for i = 0 .. (N-1);
	N - number of output data points.
	Output:
	out[i] - value of the interpolating function at interpftpad*(i-shftf[i]).
     */
     
     double x;
     int i;
     
     for (i=0; i<N; ++i) {
	  x = interpftpad*(i-shftf[i]);
	  int i1 = (int) x;
	  int i2 = i1+1;
	  double dx = x-i1;
	  double dya_abs = cabs(ya[i2]) - cabs(ya[i1]);
	  double dya_arg = carg(ya[i2]) - carg(ya[i1]);
	  //out[i] = (ya[i2]-ya[i1])*dx;
	  out[i] = dya_abs*dx*cexp(dya_arg*dx*I);
     }
     
}



// test version
void triginterp (complex double *ya, complex double *yb, double *shftf, int N, int nfft, complex double *outa, complex double *outb) {
     /* trigonometric interpolation - direct sum of fourier modes
	Input:
	ya[i], yb[i] - value of the function being interpolated in x_i = i,
	for i = 0 .. (N-1)
	Interpolating spline will be calculated at the points (i-shftf[i]), for i = 0 .. (N-1);
	N - number of output data points.
	Output:
	outa[i], outb[i] - value of the interpolating function at (i-shftf[i]).
     */
     
     double x;
     int i, k;
     complex double ef, ff, oua, oub;
     double invnfft = 1./nfft;

     ef = 2*M_PI*invnfft*I;
     
     for (i=0; i<N; ++i) {
	  x = (double)i-shftf[i];
	  ff = ef*x;
	  oua = 0.;
	  oub = 0.;
#if defined(_OPENMP)
#pragma omp parallel for default(shared) schedule(static) reduction(+:oua,oub)
#endif
	  for (k=0; k<nfft; ++k) {
	       oua += ya[k] * cexp(ff*k);
	       oub += yb[k] * cexp(ff*k);
	  }

	  outa[i] = oua*invnfft;
	  outb[i] = oub*invnfft;
	  //printf("i=%d   ya = %f %f   out = %f %f\n", i, creal(ya[i]), cimag(ya[i]), creal(out[i]), cimag(out[i]) );
     }
     //i=2100;
     //printf("i=%d   ya = %f %f   out = %f %f\n", i, creal(ya[i]), cimag(ya[i]), creal(outa[i]), cimag(outa[i]) );
}



double var (double *x, int n) {
  /* var(x, n) returns the variance (square of the standard deviation)
     of a given vector x of length n.
  */

  int i;
  double mean=0., variance=0.;

  for (i=0; i<n; i++)
    mean += x[i];
  mean /= n;
  for (i=0; i<n; i++)
    variance += sqr (x[i]-mean);
  variance /= (n-1);
  return variance;
} /* var() */



void gridr (double *M, int *spndr, int *nr, int *mr, double oms, double Smax) {

  double cof, Mp[16], smx[64], d, Ob;
  int i, j, indx[4];

  /* Grid range */

  // input:
  // *M - pointer to the array that generates the grid
  // *spndr - pointer to the range of spindowns in grid units
  // i.e., integer numbers
  // *nr and *mr - pointers to the range of sky positions
  // in grid units i.e., integer numbers

  // from settings() :
  // maximal value of the spindown:
  // Smax = 2.*M_PI*(fpo+B)*dt*dt/(2.*tau_min)
  //
  // oms equals 2.*M_PI*fpo*dt

  Ob = M_PI;
  cof = oms + Ob;

  //Mp - macierz transponowana
  for (i=0; i<4; i++)
    for (j=0; j<4; j++)
      Mp[4*i+j] = M[4*j+i];
  ludcmp (Mp, 4, indx, &d);

  for (i=0; i<8; i++) {
    smx[8*i+2] = cof;
    smx[8*i+6] = -cof;
  }
  for (i=0; i<4; i++) {
    smx[16*i+3] = smx[16*i+7] = cof;
    smx[16*i+11] = smx[16*i+15] = -cof;
  }
  for (i=0; i<8; i++) {
    smx[4*i] = Ob;
    smx[4*i+32] = -Ob;
  }
  for (i=0; i<2; i++)
    for (j=0; j<4; j++) {
      smx[32*i+4*j+1] = -Smax;
      smx[32*i+4*j+17] = 0.;
    }
  for (i=0; i<16; i++)
    lubksb (Mp, 4, indx, smx+4*i);

  spndr[0] = nr[0] = mr[0] = 16384;
  spndr[1] = nr[1] = mr[1] = -16384;

  for (i=0; i<16; i++) {
    if (floor(smx[4*i+1]) < spndr[0])
      spndr[0] = floor(smx[4*i+1]);
    if (ceil(smx[4*i+1]) > spndr[1])
      spndr[1] = ceil(smx[4*i+1]);

    if (floor(smx[4*i+2]) < nr[0])
      nr[0] = floor(smx[4*i+2]);
    if (ceil(smx[4*i+2]) > nr[1])
      nr[1] = ceil(smx[4*i+2]);

    if (floor(smx[4*i+3]) < mr[0])
      mr[0] = floor(smx[4*i+3]);
    if (ceil(smx[4*i+3]) > mr[1])
      mr[1] = ceil(smx[4*i+3]);
  }
} /* gridr() */

double FStat (double *F, int nfft, int nav, int indx) {
  /* FStat Smoothed F-statistic */

  // input:
  // *F - pointer to the value of F statistic
  // nfft - the length of the FFT data
  // nav	- length of the block (nfft/nav is the number of blocks)
  // indx - block index

  int i, j;
  double mu, *fr, pxout=100.;

  indx /= nav;
  fr = F;
  for (j=0; j<nfft/nav; j++) {
    mu = 0.;
    for (i=0; i<nav; i++)
      mu += *fr++;
    mu /= 2.*nav;
    if (mu<0.6) pxout=mu; // pci test
    //if (j == indx) pxout = mu;
    fr -= nav;
    for (i=0; i<nav; i++)
      *fr++ /= mu;
  } /* for j */
  return pxout;
} /* FStat() */

int ludcmp (double *a, int n, int *indx, double *d) {
/*	LU decomposition of a given real matrix a[0..n-1][0..n-1]
	Input:
	a		- an array containing elements of matrix a
	(changed on exit)
	n		- number of rows and columns of a
	Output:
	indx - row permutation effected by the partial pivoting
	d		- +-1 depending on whether the number of rows
	interchanged was even or odd, respectively
*/

  int i, imax = -1, j, k;
  double big, dum, sum, temp;
  double *vv;

  vv = (double *) calloc (n, sizeof (double));
  *d = 1.0;
  for (i=0; i<n; i++) {
    big = 0.0;
    for (j=0; j<n; j++)
      if ((temp=fabs (a[n*i+j])) > big)
	big = temp;
    if (big == 0.0)
      return 1;
    vv[i] = 1.0/big;
  }
  for (j=0; j<n; j++) {
    for (i=0; i<j; i++) {
      sum = a[n*i+j];
      for (k=0; k<i; k++)
	sum -= a[n*i+k]*a[n*k+j];
      a[n*i+j] = sum;
    }
    big = 0.0;
    for (i=j; i<n; i++) {
      sum = a[n*i+j];
      for (k=0; k<j; k++)
	sum -= a[n*i+k]*a[n*k+j];
      a[n*i+j] = sum;
      if ((dum = vv[i]*fabs (sum)) >= big) {
	big = dum;
	imax = i;
      }
    }
    if (j != imax) {
      for (k=0; k<n; k++) {
	dum = a[n*imax+k];
	a[n*imax+k] = a[n*j+k];
	a[n*j+k] = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (a[n*j+j] == 0.0)
      a[n*j+j] = TINY;
    if (j != n) {
      dum = 1.0/(a[n*j+j]);
      for (i=j+1; i<n; i++)
	a[n*i+j] *= dum;
    }
  }
  free (vv);
  return 0;
} /* ludcmp() */

int lubksb (double *a, int n, int *indx, double *b) {
/* Solves the set of n linear equations A X=B.
   Input:
   a[0..n-1][0..n-1] - LU decomposition af a matrix A,
   determined by ludcmp()
   n				- number of rows and columns of a
   indx[0..n-1]		- permutation vector returned by ludcmp
   b[0..n-1]			 - right-hand side vector B
   (changed on exit)
   Output:
   b[0..n-1]			- solution vector X
*/

  int i, ii=-1, ip, j;
  double sum;

  for (i=0; i<n; i++) {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if (ii>=0)
      for (j=ii; j<=i-1; j++)
	sum -= a[n*i+j]*b[j];
    else if (sum)
      ii = i;
    b[i] = sum;
  }
  for (i=n-1; i>=0; i--) {
    sum = b[i];
    for (j=i+1; j<n; j++)
      sum -= a[n*i+j]*b[j];
    b[i] = sum/a[n*i+i];
  }
  return 0;
} /* lubksb() */

int invm (const double *a, int N, double *y) {
     /* Inverse of a real matrix a[0..N-1][0..N-1].
	Input:
		a[0..N-1][0..N-1] - given matrix (saved on exit)
		N	      - number of rows and columns of a
        Output:
		y[0..N-1][0..N-1] - inverse of a
     */

  double d, *col, *al;
  int i, j, *indx;

  al = (double *) calloc (sqr(N), sizeof (double));
  indx = (int *) calloc (N, sizeof (int));
  col = (double *) calloc (N, sizeof (double));
  for (i=0; i<sqr(N); i++)
    al[i] = a[i];
  if (ludcmp (al, N, indx, &d))
    return 1;
  for (j=0; j<N; j++) {
    for (i=0; i<N; i++)
      col[i] = 0.0;
    col[j] = 1.0;
    lubksb (al, N, indx, col);
    for (i=0; i<N; i++)
      y[N*i+j] = col[i];
  }
  free (col);
  free (indx);
  free (al);
  return 0;
} /* invm() */

double det (const double *a, int N) {
  /* determinant of a real matrix a[0..N-1][0..N-1] */

  double d, *al;
  int j, *indx;

  al = (double *) calloc (sqr(N), sizeof (double));
  indx = (int *) calloc (N, sizeof (int));
  for (j=0; j<sqr(N); j++)
    al[j] = a[j];
  ludcmp (al, N, indx, &d);
  for (j=0; j<N; j++)
    d *= al[N*j+j];
  free (indx);
  free (al);
  return d;
} /* det() */

int compared2c(const void *a, const void *b) {

  double* da = (double*)a;
  double* db = (double*)b;
  
  int diff1 = (da[0] > db[0]) - (da[0] < db[0]);
  if (diff1 != 0) return diff1;
  return (da[1] > db[1]) - (da[1] < db[1]);

}


#endif
