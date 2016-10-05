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
//#define round(x) floor((x)+0.5)

#include "floats.h"

void lin2ast(double be1, double be2, int pm, double sepsm, double cepsm, \
	     double *sinal, double *cosal, double *sindel, double *cosdel);

double var (double *, int);

void gridr (double *, int *, int *, int *, double, double);

int ludcmp (double *, int, int *, double *);
int lubksb (double *, int, int *, double *);

// gridopt 
int invm (const double *, int, double *);
double det (const double *, int);

#endif
