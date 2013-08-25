#ifndef SETTINGS_H
#define SETTINGS_H

#define NPAR 5      // no. of trigger parameters

#define INT 1       // simplest interpolation
#define FFT 2       // refined (fft) interpolation

#ifdef __cplusplus
extern "C"
{
#endif

void set_search_parameters (double, char *);

typedef struct search_parameters_ {

    int nod, N, Nv, nd, s, nfft, fftpad, interpftpad;
	int spndr[2], nr[2], mr[2], pmr[2]; 
    double oms, epsma, omr, Smax, sepsm, cepsm, sphir, cphir;
	double c1, c2, c3, c4, c5, c6, c7, c8, c9;

} search_parameters;

extern search_parameters pars;

void gridr (double *, int *, int *, int *);

#ifdef __cplusplus
}
#endif

#endif // settings_h
