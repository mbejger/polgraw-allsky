#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxi.h"
#include "settings.h"

#ifndef DTAPREFIX
#define DTAPREFIX ../data
#endif

int
lineph (double epsm, double *r1, double *r2, char *dtaprefix, \
  char *ident, int N) {
  FILE *f1, *f2;
  double rSSB[3], rDet[3], ce;
  int i;
  char filename1[512], filename2[512];

  ce = cos (epsm);

  sprintf (filename1, "%s/rSSB%s.bin", dtaprefix, ident);
  sprintf (filename2, "%s/rDet%s.bin", dtaprefix, ident);

  if ((f1 = fopen (filename1, "r")) == NULL) {
    perror (filename1);
    abort ();
  }
  if ((f2 = fopen (filename2, "r")) == NULL) {
    perror (filename2);
    abort ();
  }

  for (i=0; i<N; i++) {
    fread ((void *)(rSSB), sizeof (double), 3, f1);
    fread ((void *)(rDet), sizeof (double), 3, f2);
    r1[i] = rSSB[1]/ce+rDet[1]*ce;
    r2[i] = rSSB[0]+rDet[0];
  }
  fclose (f1);
  fclose (f2);

  return 0;
} /* lineph() */
