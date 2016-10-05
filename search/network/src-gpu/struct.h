#ifndef __STRUCT_H__
#define __STRUCT_H__

// Standard C includes
#include <complex.h>
#include <floats.h>

// clFFT includes
#include <clFFT.h>

#define MAX_DETECTORS 8        // Maximum number of detectors in network 
#define DETNAME_LENGTH 2       // Detector name length (H1, L1, V1...)
#define XDATNAME_LENGTH 512    // Maximum length of input file name xdat*bin 
#define INICANDSIZE 1024       // 1048576? Initial size for array candidates storage; 
                               // realloc called if needed (in coincidences)  

#define MAXL 2048              // Max number of known lines for a detector  

// Command line option struct for search 
typedef struct _comm_line_opts {
  
  int white_flag, 		// white noise flag
      s0_flag,			// no spin-down flag
      checkp_flag,		// checkpointing flag
      veto_flag,                // veto lines flag 
      help_flag;
  
  int ident, band, hemi;
  double trl;
  double fpo_val;
  
  char prefix[512], dtaprefix[512], label[512], 
    range[512], getrange[512], qname[512], usedet[32], addsig[512], *wd;
  
} Command_line_opts;


/* input signal arrays */
typedef struct _signals {
	
  double *xDat, *xDat_d;
  double *DetSSB, *DetSSB_d;   // Ephemeris of the detector
  double *aa_d, *bb_d;         // Amplitude modulation functions
  double *shftf_d, *shft_d;    // Resampling and time-shifting
  
  double epsm, 
         phir, 
         sepsm,	  // sin(epsm)
         cepsm,	  // cos(epsm)
         sphir,	  // sin(phi_r)
         cphir,	  // cos(phi_r)
         crf0,    // number of 0s as: N/(N-Nzeros)
         sig2; 	  // variance of signal

  int Nzeros;
  COMPLEX_TYPE *xDatma_d, *xDatmb_d;

} Signals;


/* fftw arrays */
typedef struct _fft_arrays {

  COMPLEX_TYPE *xa_d, *xb_d;
  int arr_len;

} FFT_arrays;


/* Search range  */ 
typedef struct _search_range {
  int pmr[2], mr[2], nr[2], spndr[2];
  int pst, mst, nst, sst;
} Search_range;


/* FFTW plans  */ 
typedef struct _fft_plans {

  clfftPlanHandle plan,    // main plan
                  pl_int,  // interpolation forward
                  pl_inv;  // interpolation backward
  clfftPlanHandle plan2,   // main plan
                  pl_int2, // interpolation forward
                  pl_inv2; // interpolation backward

} FFT_plans;


/* Auxiluary arrays  */ 
typedef struct _aux_arrays {

  double *sinmodf_d, *cosmodf_d; // Earth position
  double *t2_d  ;                // time^2
  double *tshift_d;
  COMPLEX_TYPE *diag_d, *ldiag_d, *udiag_d, *B_d; //used in spline interpolation
  FLOAT_TYPE *mu_d, *mu_t_d; //arrays for smoothing F-stat

} Aux_arrays;


/* Search settings */ 
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
  
  int nfft,       // length of fft
      nod,        // number of days of observation
      N,          // number of data points
      nfftf,      // nfft * fftpad
      nmax,       // first and last point
      nmin, 	  // of Fstat
      s,          // number of spindowns
      nd,         // degrees of freedom
      interpftpad,
      fftpad,     // zero padding
      Ninterp, 	  // for resampling (set in plan_fftw() init.c)
      nifo;       // number of detectors

  double *M;      // Grid-generating matrix (or Fisher matrix, 
                  // in case of coincidences) 

  double vedva[4][4];   // transformation matrix: its columns are 
                        // eigenvectors, each component multiplied 
                        // by sqrt(eigval), see init.c manage_grid_matrix(): 
                        // sett->vedva[i][j]  = eigvec[i][j]*sqrt(eigval[j])  

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

extern Detector_settings ifo[MAX_DETECTORS]; 

// Command line option struct for coincidences 
typedef struct _comm_line_opts_coinc {
  
  int help_flag; 
  
  int shift, // Cell shifts  (4 digit number corresponding to fsda, e.g. 0101)  
      scale, // Cell scaling (4 digit number corresponding to fsda, e.g. 4824) 
      refr;  // Reference frame 

  // Minimal number of coincidences recorded in the output  
  int mincoin; 

  double fpo, refgps, narrowdown, snrcutoff; 
  
  char prefix[512], dtaprefix[512], trigname[512], refloc[512], *wd;
  
} Command_line_opts_coinc;

typedef struct _triggers { 

  int frameinfo[256][3];    // Info about candidates in frames: 
                            // - [0] frame number, [1] initial number 
                            // of candidates, [2] number of candidates
                            // after sorting    

  int frcount, goodcands; 

} Candidate_triggers; 

#endif
