/* test LAL pulsar library */

#include <stdio.h>
#include <stdlib.h>

//#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/LALBarycenter.h>

#define efile "earth00-19-DE405.dat"
#define sfile "sun00-19-DE405.dat"
#define N 258492 //47127
#define EPSILON 0.40909280422232891

#include "EphemerisDetector.h"

int get_barycenter (double, Detectors, EphemerisData *,	\
		    double *, double *, double, int);


int main (int argc, char *argv[]) {
  double gps1 = 1.1260846080e+09;  /* 010 O1 */
  double dt = 2.0;
  double bandwidth;
  char name[] = "L1";
  double position[4], mjd1, phir, elam;
  Detectors detector;
  EphemerisData *edat;
  double *DetSSB, *rDet, *rSSB;
  FILE *f;
  int j;

  DetSSB = (double *) calloc (3*N+2, sizeof(double));
  rDet = (double *) calloc (3*N, sizeof(double));
  rSSB = (double *) calloc (3*N, sizeof(double));
  
  detector = get_detector(name);
  fprintf (stderr, "Generating ephemeris for %s detector\n",
	   names[detector]);
  get_position (detector, position);
  fprintf (stderr, "ephi = %f\nelam = %f\neheight = %f\negam = %f\n",
	   position[0], position[1], position[2], position[3]);
  elam = position[1];
  mjd1 = gps2mjd (gps1);
  phir = sid (mjd1, elam);
  fprintf (stderr, "mjd = %f\nphir = %f\n", mjd1, phir);
  edat = XLALInitBarycenter (efile, sfile);
  bandwidth = 1/(2*dt);
  get_barycenter (gps1, L1, edat, DetSSB, rDet, dt, N); //here also change detector name
  DetSSB[3*N] = phir;
  DetSSB[3*N+1] = EPSILON;
  for (j=0; j<3*N; j++)
    rSSB[j] = DetSSB[j] - rDet[j];

  if ((f=fopen ("DetSSB.bin", "w")) != NULL) {
    fwrite ((void *)DetSSB, sizeof(double), 3*N+2, f);
    fclose (f);
  }

  if ((f=fopen ("rDet.bin", "w")) != NULL) {
    fwrite ((void *)rDet, sizeof(double), 3*N, f);
    fclose (f);
  }

  if ((f=fopen ("rSSB.bin", "w")) != NULL) {
    fwrite ((void *)rSSB, sizeof(double), 3*N, f);
    fclose (f);
  }

  printf ("success\n");
  return 0;
}
