/* test LAL pulsar library */

#include <stdio.h>
#include <stdlib.h>

#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/LALBarycenter.h>

#define efile "earth00-19-DE405.dat"
#define sfile "sun00-19-DE405.dat"
#define N 344656

#include "EphemerisDetector.h"

int get_barycenter (double, Detectors, EphemerisData *,	\
		    double *, double *, double, int);


int main (int argc, char *argv[]) {
  double gps1 = 9.437876e+08;  /* 021 */
  char name[] = "H1";
  double position[4], mjd1, phir, elam;
  Detectors detector;
  EphemerisData *edat;
  double *DetSSB, *rDet;
  FILE *f;

  DetSSB = (double *) calloc (3*N, sizeof(double));
  rDet = (double *) calloc (3*N, sizeof(double));
  
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
  get_barycenter (gps1, H1, edat, DetSSB, rDet, 0.5, N);

  if ((f=fopen ("DetSSB.bin", "w")) != NULL) {
    fwrite ((void *)DetSSB, sizeof(double), 3*N, f);
    fclose (f);
  }

  if ((f=fopen ("rDet.bin", "w")) != NULL) {
    fwrite ((void *)rDet, sizeof(double), 3*N, f);
    fclose (f);
  }

  printf ("success\n");
  return 0;
}
