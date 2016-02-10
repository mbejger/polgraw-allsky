#ifndef __JOBCORE_H__
#define __JOBCORE_H__

#include "auxi.h"
#include "struct.h"
#include "cublas_v2.h"

#define BLOCK_DIM(n, b) ((n)/b + ((n)%b==0 ? 0 : 1))

void search(
	    Search_settings *sett,
	    Command_line_opts *opts,
	    Search_range *s_range,
	    FFT_plans *plans,
	    FFT_arrays *fft_arr,
	    Aux_arrays *aux,
	    int *Fnum,
	    double *F
	    );

/* Main job function
 * The output is stored in single or double precision 
 * (FLOAT_TYPE defined in struct.h)  
 */ 

FLOAT_TYPE* job_core(
		     int pm,                   // hemisphere
		     int mm,                   // grid 'sky position'
		     int nn,                   // other grid 'sky position'
		     Search_settings *sett,    // search settings
		     Command_line_opts *opts,  // cmd opts
		     Search_range *s_range,    // range for searching
		     FFT_plans *plans,        // plans for fftw
		     FFT_arrays *fft_arr,    // arrays for fftw
		     Aux_arrays *aux,          // auxiliary arrays
		     double *F,                // F-statistics array
		     int *sgnlc,               // reference to array with the parameters
		     // of the candidate signal
		     // (used below to write to the file)
		     int *FNum,               // Candidate signal number
		     cublasHandle_t scale      //handle for scaling
		     );
      

#endif
