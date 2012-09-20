/* lvcvirgo.h */

#ifndef LVCVIRGO_H
#define LVCVIRGO_H

#include <fftw3.h>

#ifdef __cplusplus
extern "C"
{
#endif

#define NPAR 5 		/* no. of trigger parameters */

#define INT 1		/* simplest interpolation */
#define FFT 2		/* refined (fft) interpolation */

  extern int nod, N, Nv, nfft, s, nd, fftpad, interpftpad;
  extern double dt, B, oms, deg, c, AU, epsma, a, f, b, Omegar,		\
    omr, SIDday, TAIday, ephi, elam, eheight, egam, epsi, r, yr,	\
    tau_min, alfa, Smax, c1, c2, c3, c4, c5, c6, c7, c8, c9;

  void lvcvirgo (double);
  int rogcvir (void);
  void lin2ast (double, double, int, double, double,
		double *, double *, double *, double *);
  void modvir (double, double, double, double, double, double,
		double *, double *, int);
  double FStat (double *, int, int, int);

  double *JobCore (int, int, int, int, int, double *,			\
		   double *, double *, int, int,			\
		   int, fftw_complex *, fftw_complex *,			\
		   fftw_complex *, fftw_complex *,			\
		   fftw_plan, fftw_plan, fftw_plan,			\
		   int, int, double, double, double, double,	 	\
		   int *, int, int, int *, double, double,		\
		   double, int);

  void gridr (double *, int *, int *, int *);

#ifdef __cplusplus
}
#endif

#endif                          /* LVCVIRGO_H */
