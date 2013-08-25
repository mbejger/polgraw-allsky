#include <math.h>

#include "auxi.h"
#include "settings.h"

void
gridr (double *M, int *spndr, int *nr, int *mr) {
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
  cof = pars.oms + Ob;

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
      smx[32*i+4*j+1] = -pars.Smax;
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


