#ifndef __AUXI_H__
#define __AUXI_H__

#include <complex.h>

#define sqr(x) ((x)*(x))
#define TOSTRA(x) #x
#define TOSTR(x) TOSTRA(x)

#define TINY 1.0e-20
#define NINTERP 3			 /* degree of the interpolation polynomial */
							 /* Do not change this value!!! */
#define NAVFSTAT 4096
#define round(x) floor((x)+0.5)

// Define COMP_FLOAT this to change the double/single precision of triggers 
//
// #define COMP_FLOAT
#ifdef COMP_FLOAT // if single-precision
    #define FLOAT_TYPE float
#else             // if double-precision
    #define FLOAT_TYPE double
#endif

void lin2ast(double be1, double be2, int pm, double sepsm, double cepsm,	\
         double *sinal, double *cosal, double *sindel, double *cosdel);

int ast2lin(FLOAT_TYPE alfa, FLOAT_TYPE delta, double epsm, double *be); 
 
void spline(complex double *, int, complex double *);
complex double splint (complex double *, complex double *, int, double);
void splintpad (complex double *, double *, int, int, complex double*);
void linterp (complex double *, double *, int, int, complex double*);
void triginterp (complex double *ya, complex double *yb, double *shftf, int N, int nfft, complex double *outa, complex double *outb);
double var (double *, int);

void gridr (double *, int *, int *, int *, double, double);
double FStat (double *, int, int, int);

int ludcmp (double *, int, int *, double *);
int lubksb (double *, int, int *, double *);

// gridopt 
int invm (const double *, int, double *);
double det (const double *, int);

// for qsorting the lines
int compared2c (const void *, const void *);

// signal handler
static void sig_handler(int signo);

#endif
