#include <math.h>
//#include <fftw3.h>
#include <complex.h>
#include <cuda.h>
#include <cufft.h>
#include "auxi.h"
#include "settings.h"

extern double *t2, oms;

double FStat (double *F, int nfft, int nav, int indx) {
  /* FStat Smoothed F-statistic */

// input: 
// *F - pointer to the value of F statistic 
// nfft - the length of the FFT data
// nav  - length of the block (nfft/nav is the number of blocks)
// indx - block index 

  int i, j;
  double mu, *fr, pxout=0.;

  indx /= nav;
  fr = F;
  for (j=0; j<nfft/nav; j++) {
    mu = 0.;
    for (i=0; i<nav; i++)
      mu += *fr++;
    mu /= 2.*nav;
    if (j == indx)
      pxout = mu;
    fr -= nav;
    for (i=0; i<nav; i++)
      *fr++ /= mu;
  } /* for j */
  return pxout;
} /* FStat() */

