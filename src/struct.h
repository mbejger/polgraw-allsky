#ifndef __STRUCT_H__
#define __STRUCT_H__


#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cufft.h>

#include "floats.h"


#define CHAR_BUFFER_SIZE 1024

typedef struct __comm_line_opts {
	int white_flag; 							// white noise flag
	int s0_flag;								// no spin-down flag
	int checkp_flag;							// checkpointing flag
	int help_flag;
	
	int fftinterp;								// interpolation type
	int ident, band, hemi;
	double trl;									// critical F-stat value
	double fpo_val;							// frequency
	
	char prefix[CHAR_BUFFER_SIZE], dtaprefix[CHAR_BUFFER_SIZE], label[CHAR_BUFFER_SIZE], range[CHAR_BUFFER_SIZE], *wd,
			ifo_choice[3];
	char qname[CHAR_BUFFER_SIZE];
} Command_line_opts;



typedef struct _arrays {
	double *xDat, *cu_xDat; //signal
	cufftDoubleComplex *cu_xa, *cu_xb; //signal times a(t) and b(t) 
	cufftDoubleComplex *cu_xar, *cu_xbr; //resampled and interpolated

	COMPLEX_TYPE *cu_xa_f, *cu_xb_f; //signal times a(t) and b(t) 
	COMPLEX_TYPE *cu_xar_f, *cu_xbr_f; //resampled and interpolated

	COMPLEX_TYPE *cu_xa2_f, *cu_xb2_f; //auxiliary arrays for interbinning	
	
	double *cu_sinmodf, *cu_cosmodf; // Earth position
	double *sinmodf, *cosmodf; // Earth position
	double *aa, *bb;
	double *cu_aa, *cu_bb; //amplitude modulation functions
	double *cu_shftf;//*cu_shft; //used to resample and shift time
	FLOAT_TYPE *cu_shft;
	double *cu_tshift; //shifted time (for spline interpolation)
	double *DetSSB; //ephemeris of the detector
	double *cu_DetSSB; //ephemeris of the detector

//!!!
	FLOAT_TYPE *cu_cand_params; //found candidates (max [>trl])
	FLOAT_TYPE *cu_cand_buffer; //buffer
	FLOAT_TYPE *cand_buffer;
	int *cu_cand_count;
	int cand_params_size, cand_buffer_size;

	cufftDoubleComplex *cu_d, *cu_dl, *cu_du, *cu_B; //used in spline interpolation	

	FLOAT_TYPE *cu_mu, *cu_mu_t; //arrays for smoothing F-stat

	int arr_len;

	//auxiliary arrays used in modvir_gpu
	double *cu_o_aa, *cu_o_bb, *cu_o_aa2, *cu_o_bb2;

} Arrays;


//search range
typedef struct _search_range {
	int pmr[2], mr[2], nr[2], spndr[2];
	int pst, mst, nst, sst;
} Search_range;


//fftw plans
typedef struct _fft_plans {
	cufftHandle plan, //main plan
				   pl_int, //interpolation forward
				   pl_inv; //interpolation backward
} FFT_plans;



typedef struct _detector_settings {
	double fpo, //frequency
			 dt, //sampling time
			 B, // bandwidth
			 oms, //dimensionless angular frequency (fpo)
			 omr, //C_OMEGA_R * dt - dimensionless Earth's angular frequency
			 crf0, //number of 0 as: N/(N-Nzeros)
			 Smin, //minimum spindown
			 Smax, //maximum spindown
			 alfa, //false alarm probability
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
		 nfftf, //nfft * fftpad
		 nod, //number of days of observation
		 N, //number of data points
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
	double c1,c2,c3,c4,c5,c6,c7,c8,c9;

} Ampl_mod_coeff;



#endif
