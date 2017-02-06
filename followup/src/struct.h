#ifndef __STRUCT_H__
#define __STRUCT_H__

#include <fftw3.h>
#include <complex.h>

#define MAX_DETECTORS 8        // Maximum number of detectors in network 
#define DETNAME_LENGTH 2       // Detector name length (H1, L1, V1...)
#define XDATNAME_LENGTH 512    // Maximum length of input file name xdat*bin 
#define INICANDSIZE 1024       // 1048576? Initial size for array candidates storage; 
                               // realloc called if needed (in coincidences)  

#define MAXL 2048              // Max number of known lines for a detector  

// Command line option struct for search 
typedef struct _comm_line_opts {
  
  int veto_flag,                // veto lines flag 
      simplex_flag,		// Simplex direct maximum search flag
      mads_flag,		// MADS direct maximum search flag
      gauss_flag,		// Generate Gaussian noise instead of reading data
      neigh_flag,		// Area of calculation will be defined as %% from initial value
      help_flag;		
  
  int ident, band, hemi, refr;
  double trl;
  double fpo_val;
  
  char prefix[512], dtaprefix[512], label[512], qname[512], 
    usedet[32], addsig[512], candidates[512], glue[512], 
    gauss[512], *wd;
  
} Command_line_opts;


// input signal arrays
typedef struct _signals {
	
  double *xDat;
  double *DetSSB;       // Ephemeris of the detector
  double *aa, *bb;      // Amplitude modulation functions
  double *shftf, *shft; // Resampling and time-shifting
  
  double epsm, 
         phir, 
         sepsm,	  // sin(epsm)
         cepsm,	  // cos(epsm)
         sphir,	  // sin(phi_r)
         cphir,	  // cos(phi_r)
         crf0,    // number of 0s as: N/(N-Nzeros)
         sig2; 	  // variance of signal

  int Nzeros;
  complex double *xDatma, *xDatmb;

} Signals;

  /* Auxiluary arrays
   */ 

typedef struct _aux_arrays {

  double *sinmodf, *cosmodf; // Earth position
  double *t2;                // time^2

} Aux_arrays;

  /* Search settings 
   */ 

typedef struct _search_settings {

  double fpo,    // Band frequency
         dt,     // Sampling time
         B,      // Bandwidth
         oms,    // Dimensionless angular frequency (fpo)
         omr,    // C_OMEGA_R * dt 
                 // (dimensionless Earth's angular frequency)
         Smin,   // Minimum spindown
         Smax,   // Maximum spindown
         sepsm,	 // sin(epsm)
         cepsm;	 // cos(epsm)
  
  int nod,        // number of days of observation
      N,          // number of data points
      nd,         // degrees of freedom
      fftpad,     // zero padding
      nifo;       // number of detectors

  double *M;      // Grid-generating matrix (or Fisher matrix, 
                  // in case of coincidences) 

  double lines[MAXL][2]; // Array for lines in given band 
  int numlines_band;     // number of lines in band   

} Search_settings;


  /* Amplitude modulation function coefficients
   */ 

typedef struct _ampl_mod_coeff {
	double c1, c2, c3, c4, c5, c6, c7, c8, c9;
} Ampl_mod_coeff;


  /* Detector and its data related settings 
   */ 

typedef struct _detector { 

  char xdatname[XDATNAME_LENGTH]; 
  char name[DETNAME_LENGTH]; 
 
  double ephi, 		// Geographical latitude phi in radians
         elam, 		// Geographical longitude in radians 
         eheight,   // Height h above the Earth ellipsoid in meters
         egam; 		// Orientation of the detector gamma  

  Ampl_mod_coeff amod; 
  Signals sig;  

  double lines[MAXL][2]; // Array for lines: column values 
                         // are beginning and end of line to veto 
  int numlines;                        
 
} Detector_settings; 


  /* Array of detectors (network) 
   */ 

struct _detector ifo[MAX_DETECTORS]; 

#endif 
