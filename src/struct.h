#ifndef __STRUCT_H__
#define __STRUCT_H__

#include <fftw3.h>
#include <complex.h>

typedef struct __comm_line_opts {
	int white_flag; 							// white noise flag
	int s0_flag;								// no spin-down flag
	int checkp_flag;							// checkpointing flag
	int help_flag;
	
	int fftinterp;
	int ident, band, hemi;
	double trl;
	double fpo_val;
	
	char prefix[64], dtaprefix[64], label[64], range[64], *wd,
			ifo_choice[3];
	char qname[64];
} Command_line_opts;




//signal arrays
typedef struct _signals {
	
	double *xDat;
	complex double *xDatma, *xDatmb;
} Signals;


//fftw arrays
typedef struct _fftw_arrays {
	fftw_complex *xa, *xb;
					//  *xDa, *xDb,
					 // *rDa, *rDb;
	int arr_len;
	
} FFTW_arrays;

//search range
typedef struct _search_range {
	int pmr[2], mr[2], nr[2], spndr[2];
	int pst, mst, nst, sst;
} Search_range;


//fftw plans
typedef struct _fftw_plans {
	fftw_plan plan, //main plan
				  pl_int, //interpolation forward
				  pl_inv; //interpolation backward
	fftw_plan plan2, //main plan
				  pl_int2, //interpolation forward
				  pl_inv2; //interpolation backward
} FFTW_plans;

//auxiluary arrays
typedef struct _aux_arrays {

	double *sinmodf, *cosmodf; // Earth position
	double *t2; // time^2
	double *aa, *bb; //amplitude modulation functions
	double *shftf, *shft; //used to resample and shift time
	double *DetSSB; //ephemeris of the detector

} Aux_arrays;


typedef struct _detector_settings {
	double fpo, //frequency
			dt, //sampling time
			B, // bandwidth
			oms, //dimensionless angular frequency (fpo)
			omr, //C_OMEGA_R * dt - dimensionless Earth's angular frequency

			//#mb 
			crf0, //number of 0 as: N/(N-Nzeros)
			Smin, //minimum spindown
			Smax, //maximum spindown
			alfa, //false alarm probability
			//#mb
			ephi, 		//position
			elam, 		//of
			eheight, 	//the
			egam, 		//detector

			sig2,		//variance of signal
			sepsm,		// sin(epsm)
			cepsm,		// cos(epsm)
			sphir,		// sin(phi_r)
			cphir;		// cos(phi_r)

	int nfft, // length of fft
		 nod, //number of days of observation
		 N, //number of data points
		 nfftf, //nfft * fftpad
		 nmax,	 	//first and last point
		 nmin, 		//of Fstat
		 s, //number of spindowns
		 nd, //degrees of freedom
		 interpftpad,
		 fftpad, //zero padding
		 Ninterp; 	// for resampling
			 
	double *M; // Grid-generating matrix
} Detector_settings;


//Amplitude modulation function coefficients
typedef struct _ampl_mod_coeff {
	double c1, c2, c3, c4, c5, c6, c7, c8, c9;
} Ampl_mod_coeff;



#endif
