#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "settings.h"
#include "auxi.h"
#include "struct.h"

#include "cuda_error.h"


/* initializes settings struct
 * first creates local variables, initializes them 
 * and copies them to struct of which pointer is passed as the third argument
 */
void settings (Detector_settings* sett, Command_line_opts *opts, Arrays *arr)
{
  double fpo = sett->fpo;
  char* ifo_choice = opts->ifo_choice;

  int nod, N, nfft, s, nd, fftpad, interpftpad;
  double dt, B, oms;
  double  omr;
  double ephi, elam, eheight, egam;//, epsi;
  double alfa, Smin, Smax;

  dt = 0.5;				 // data sampling time
  B = 0.5/dt;				 // Bandwidth
  oms = 2.*M_PI*fpo*dt;			 // Dimensionless angular frequency

  omr = C_OMEGA_R*dt;
  nod = 2;				 // Observation time in days
  N = round (nod*C_SIDDAY/dt);           // No. of data points

  nfft = 1 << (int)ceil(log(N)/log(2.)); // length of FFT
  s = 1;				 // No. of spindowns
  Smin = 1000.*C_YEARSEC;		 // Minimum spindown time in sec.
  Smax = 2.*M_PI*(fpo+B)*dt*dt/(2.*Smin); // Maximum spindown (1000 years)
                                              // [angular, dimensionless]


  alfa = .01;		// False alarm probability
  nd = 2;		// Degree of freedom, 2*nd = deg. no ofrees of freedom for Chi^2
  fftpad = 2;	        // Zero padding (original grid: 2, new grids: 1)
  interpftpad = 2;

  if(!strcmp("V1", ifo_choice)) {
    // Geographical location of the Virgo detector
    //
    // Geographical latitude phi in radians
    ephi = (43.+37./60.+53.0880/3600.)/RAD_TO_DEG;
    // Geographical longitude in radians
    elam = (10.+30./60.+16.1885/3600.)/RAD_TO_DEG;
    // Height h above the Earth ellipsoid in meters
    eheight = 53.238;
    // Orientation of the detector gamma
    egam = (135. - (19.0+25./60.0+57.96/3600.))/RAD_TO_DEG;

    printf("Using the Virgo detector...\n");

  } else if(!strcmp("H1", ifo_choice)) {

    // Geographical location of the Hanford H1 detector
    //
    // Geographical latitude phi in radians
    ephi = (46+(27+18.528/60.)/60.)/RAD_TO_DEG;
    // Geographical longitude in radians
    elam = -(119+(24+27.5657/60.)/60.)/RAD_TO_DEG;
    // Height h above the Earth ellipsoid in meters
    eheight = 142.554;
    // Orientation of the detector gamma
    egam	= 170.9994/RAD_TO_DEG;

    printf("Using the LIGO Hanford detector...\n");

  } else if(!strcmp("L1", ifo_choice)) {

    // Geographical location of the Livingston L1 detector
    //
    // Geographical latitude phi in radians
    ephi = (30+(33+46.4196/60.)/60.)/RAD_TO_DEG;
    // Geographical longitude in radians
    elam = - (90+(46+27.2654/60.)/60.)/RAD_TO_DEG;
    // Height h above the Earth ellipsoid in meters
    eheight = - 6.574;
    // Orientation of the detector gamma
    egam = 242.7165/RAD_TO_DEG;

    printf("Using the LIGO Livingston detector...\n");

  } else {

    printf("Meh, unknown detector. Exiting...\n");
    exit(EXIT_FAILURE);

  }


  /*
    ############ Earth's rotation ################
  */
	
  //	int i;

  /*
  //>
  double omrt;
  for (i=0; i<N; i++) {
  omrt = omr*i; //Earth angular velocity * dt * i
  aux->t2[i] = sqr ((double)i);
  aux->cosmodf[i] = cos (omrt);
  aux->sinmodf[i] = sin (omrt);
  }
  */
  CudaSafeCall ( cudaMalloc((void**)&arr->cu_sinmodf, sizeof(double)*N));
  CudaSafeCall ( cudaMalloc((void**)&arr->cu_cosmodf, sizeof(double)*N));
	
  arr->cosmodf = (double *) calloc (N, sizeof (double));
  arr->sinmodf = (double *) calloc (N, sizeof (double));
  double omrt;
  int i;
  for (i=0; i<N; i++) {
    omrt = omr*i; //Earth angular velocity * dt * i
    arr->cosmodf[i] = cos (omrt);
    arr->sinmodf[i] = sin (omrt);
  }
  CudaSafeCall ( cudaMemcpy(arr->cu_sinmodf, arr->sinmodf, sizeof(double)*N, cudaMemcpyHostToDevice));
  CudaSafeCall ( cudaMemcpy(arr->cu_cosmodf, arr->cosmodf, sizeof(double)*N, cudaMemcpyHostToDevice));

	
  //	compute_sincosmodf_wrapper(arr->c_sinmodf, arr->c_cosmodf, omr, N);
  //	compute_sincosmodf<<<N/256+1,256>>>(arr->c_sinmodf, arr->c_cosmodf, omr, N);
  // change to kernel


  sett->dt=dt;           // sampling time
  sett->B=B;             // bandwidth
  sett->oms=oms;         // dimensionless angular frequency
  sett->omr=omr;         // C_OMEGA_R * dt
  sett->nod=nod;         // number of days of observation
  sett->N=N;             // number of data points
  sett->nfft=nfft;       // length of fft
  sett->s=s;             // number of spindowns
  sett->Smin=Smin;       // minimum spindown
  sett->Smax=Smax;       // maximum spindown
  sett->alfa=alfa;       // false alarm probability
  sett->nd=nd;           // degrees of freedom
  sett->fftpad=fftpad;   // zero padding
  sett->interpftpad=interpftpad;
  sett->ephi=ephi;
  sett->elam=elam;
  sett->eheight=eheight;
  sett->egam=egam;

  sett->Ninterp = sett->interpftpad*sett->nfft;


  // Because of frequency-domain filters, we search
  // F-statistic in range (nmin+1, nmax) of data points
  sett->nmin = sett->fftpad*NAV;
  sett->nmax = (sett->nfft/2 - NAV)*sett->fftpad;


} // settings()



void rogcvir (Ampl_mod_coeff* amod, Detector_settings* sett) {

  /* Coefficients of the amplitude modulation functions
     for VIRGO detector
  */

  // In the notation of Phys. Rev. D 58, 063001 (1998):
  // ephi = lambda (geographical latitude phi in radians)
  // egam = gamma (orientation of the detector)
  //
  // (see modvir function in JobNAllSky-common.c
  // for full Eqs. 12 and 13)

  amod->c1 = .25*sin(2.*sett->egam)*(1+sqr(sin(sett->ephi)));
  amod->c2 = -.5*cos(2.*sett->egam)*sin(sett->ephi);
  amod->c3 = .5*sin(2.*sett->egam)*sin(2.*sett->ephi);
  amod->c4 = -cos(2.*sett->egam)*cos(sett->ephi);
  amod->c5 = .75*sin(2.*sett->egam)*sqr(cos(sett->ephi));
  amod->c6 = cos(2.*sett->egam)*sin(sett->ephi);
  amod->c7 = .5*sin(2.*sett->egam)*(1.+sqr(sin(sett->ephi)));
  amod->c8 = cos(2.*sett->egam)*cos(sett->ephi);
  amod->c9 = .5*sin(2.*sett->egam)*sin(2.*sett->ephi);

} /* rogcvir() */
