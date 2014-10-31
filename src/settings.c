#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
					
#include "settings.h"
#include "auxi.h"
#include "struct.h"

/* Initializes settings structs: 
 * first creates local variables, initializes them 
 * and copies them to struct of which pointer is passed as the third argument
 */

void settings(
	Search_settings* sett, 
	Command_line_opts *opts, 
	Aux_arrays *aux) {

  double dt, B, oms, omr, alfa, Smin, Smax;
  int i=0, nod, N, nfft, s, nd, fftpad, interpftpad;

  /* Network of detectors' discovery: 
   * subdirectories in the main input directory
   * by convention should be named like V1, L1, H1 
	 * and contain input data and ephemerids. 
	 */ 

  char dirname[512];
  // Main input directory name 
  sprintf (dirname, "%s/%03d", opts->dtaprefix, opts->ident); 

  DIR *dp;
  struct dirent *ep;
  char **detnames = malloc(MAX_DETECTORS*sizeof(char*));   

  dp = opendir (dirname);
  if (dp != NULL) {
		while ((ep = readdir (dp))) { 

			// Subdirectory names: 2 char long
			if((ep->d_type == DT_DIR) && 
		  	(strlen(ep->d_name)==2) && 
		  	strncmp(&ep->d_name[0],".",1)) { 

			  	detnames[i] = malloc(DETNAME_LENGTH); 
			  	strncpy(detnames[i], ep->d_name, DETNAME_LENGTH);
  			  i++; 
			}
    } 
      
		(void) closedir (dp);

  } else perror ("Couldn't open the input directory");

	sett->nifo=i;      // number of detectors  
  printf("Settings - number of detectors: %d\n", sett->nifo); 

	/* Search settings: 
	 * lenghts, bandwidth and Earth details
   */

  dt = 0.5;												// data sampling time
  B = 0.5/dt;											// Bandwidth
  oms = 2.*M_PI*(sett->fpo)*dt;		// Dimensionless angular frequency

  omr = C_OMEGA_R*dt;

  nod = 2;												// Observation time in days
  N = round (nod*C_SIDDAY/dt);		// No. of data points

  nfft = 1 << (int)ceil(log(N)/log(2.));	// length of FFT
  s = 1;																	// No. of spindowns

  Smin = 1000.*C_YEARSEC;					// Minimum spindown time 
																	// [sec.]

  // Maximum spindown (1000 years) [angular, dimensionless]
  Smax = 2.*M_PI*(sett->fpo + B)*dt*dt/(2.*Smin);   

  alfa = .01;			// False alarm probability
  nd = 2;					// Degree of freedom, 
									// (2*nd = deg. no ofrees of freedom for chi^2)

  fftpad = 2;			// Zero padding (original grid: 2, new grids: 1)
  interpftpad = 2;

  sett->dt=dt;        	// sampling time
  sett->B=B;          	// bandwidth
  sett->oms=oms;      	// dimensionless angular frequency
  sett->omr=omr;      	// C_OMEGA_R * dt
  sett->nod=nod;      	// number of days of observation
  sett->N=N;          	// number of data points
  sett->nfft=nfft;    	// length of fft
  sett->s=s;          	// number of spindowns
  sett->Smin=Smin;    	// minimum spindown
  sett->Smax=Smax;    	// maximum spindown
  sett->alfa=alfa;    	// false alarm probability
  sett->nd=nd;        	// degrees of freedom
  sett->fftpad=fftpad;  // zero padding
  sett->interpftpad=interpftpad;

  // Because of frequency-domain filters, we search
  // F-statistic in range (nmin+1, nmax) of data points
  sett->nmin = sett->fftpad*NAV;
  sett->nmax = (sett->nfft/2 - NAV)*sett->fftpad;

  for(i=0; i<sett->nifo; i++) { 

    // Virgo detector
    if(!strcmp("V1", detnames[i])) {

      strncpy(ifo[i].name, detnames[i], DETNAME_LENGTH);
      // Geographical latitude phi in radians
      ifo[i].ephi = (43.+37./60.+53.0880/3600.)/RAD_TO_DEG;
      // Geographical longitude in radians
      ifo[i].elam = (10.+30./60.+16.1885/3600.)/RAD_TO_DEG;
      // Height h above the Earth ellipsoid in meters
      ifo[i].eheight = 53.238;
      // Orientation of the detector gamma
      ifo[i].egam = (135. - (19.0+25./60.0+57.96/3600.))/RAD_TO_DEG;

      printf("Using %s IFO (as detector #%d)...\n", ifo[i].name, i);

    // Hanford H1 detector
    } else if(!strcmp("H1", detnames[i])) {

      strncpy(ifo[i].name, detnames[i], DETNAME_LENGTH);
      // Geographical latitude phi in radians
      ifo[i].ephi = (46+(27+18.528/60.)/60.)/RAD_TO_DEG;
      // Geographical longitude in radians
      ifo[i].elam = -(119+(24+27.5657/60.)/60.)/RAD_TO_DEG;
      // Height h above the Earth ellipsoid in meters
      ifo[i].eheight = 142.554;
      // Orientation of the detector gamma
      ifo[i].egam	= 170.9994/RAD_TO_DEG;

      printf("Using %s IFO (as detector #%d)...\n", ifo[i].name, i);
  
    // Livingston L1 detector
    } else if(!strcmp("L1", detnames[i])) {

      strncpy(ifo[i].name, detnames[i], DETNAME_LENGTH);
      // Geographical latitude phi in radians
      ifo[i].ephi = (30+(33+46.4196/60.)/60.)/RAD_TO_DEG;
      // Geographical longitude in radians
      ifo[i].elam = -(90+(46+27.2654/60.)/60.)/RAD_TO_DEG;
      // Height h above the Earth ellipsoid in meters
      ifo[i].eheight = -6.574;
      // Orientation of the detector gamma
      ifo[i].egam = 242.7165/RAD_TO_DEG;

      printf("Using %s IFO (as detector #%d)...\n", ifo[i].name, i);

    } else {
      printf("Meh, unknown detector %s (see settings.c) Exiting...\n", detnames[i]);
      exit(EXIT_FAILURE);
    }

  } 

  // memory free for detnames 
  for(i=0; i<MAX_DETECTORS; i++) free(detnames[i]); 
  free(detnames); 

  // Auxiliary arrays, Earth's rotation
  aux->t2 = (double *) calloc (N, sizeof (double));
  aux->cosmodf = (double *) calloc (N, sizeof (double));
  aux->sinmodf = (double *) calloc (N, sizeof (double));
  double omrt;

  for (i=0; i<N; i++) {
    omrt = omr*i;     //Earth angular velocity * dt * i
    aux->t2[i] = sqr((double)i);
    aux->cosmodf[i] = cos (omrt);
    aux->sinmodf[i] = sin (omrt);
  }

} // settings


  /* Coefficients of the amplitude modulation functions
   * of the Virgo detector
   */ 

void rogcvir(
	Ampl_mod_coeff* amod, 
	Search_settings* sett) {

  // In the notation of Phys. Rev. D 58, 063001 (1998):
  // ephi = lambda (geographical latitude phi in radians)
  // egam = gamma (orientation of the detector)
  //
  // (see modvir function in jobcore.c
  // for full Eqs. 12 and 13)

  amod->c1 = .25*sin(2.*ifo[0].egam)*(1+sqr(sin(ifo[0].ephi)));
  amod->c2 = -.5*cos(2.*ifo[0].egam)*sin(ifo[0].ephi);
  amod->c3 = .5*sin(2.*ifo[0].egam)*sin(2.*ifo[0].ephi);
  amod->c4 = -cos(2.*ifo[0].egam)*cos(ifo[0].ephi);
  amod->c5 = .75*sin(2.*ifo[0].egam)*sqr(cos(ifo[0].ephi));
  amod->c6 = cos(2.*ifo[0].egam)*sin(ifo[0].ephi);
  amod->c7 = .5*sin(2.*ifo[0].egam)*(1.+sqr(sin(ifo[0].ephi)));
  amod->c8 = cos(2.*ifo[0].egam)*cos(ifo[0].ephi);
  amod->c9 = .5*sin(2.*ifo[0].egam)*sin(2.*ifo[0].ephi);

} // rogcvir

