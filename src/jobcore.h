#ifndef __JOBCORE_H__
#define __JOBCORE_H__

#include "struct.h"

#include "floats.h"

#define BLOCK_SIZE 256
#define BLOCK_SIZE_RED 128
#define BLOCK_DIM(n, b) ((n)/b + ((n)%b==0 ? 0 : 1))
#define NAV_THREADS 16

void search(
				Detector_settings *sett,
				Command_line_opts *opts,
				Search_range *s_range,
				Arrays *arr,
				FFT_plans *plans,
				Ampl_mod_coeff *amod,
				int *Fnum,
				FLOAT_TYPE *cu_F
				);


double* job_core(
			int pm,			// hemisphere
			int mm,			// grid 'sky position'
			int nn,			// other grid 'sky position'
			Detector_settings *sett, // detector settings
			Command_line_opts *opts, // cmd opts
			Search_range *s_range,	 // range for searching
			Arrays *arr,				//arrays
			FFT_plans *plans,       // plans for fftw
			Ampl_mod_coeff *amod,	//amplitude modulation functions coefficients
			FLOAT_TYPE *cu_F,			// F-stat on GPU
			int *FNum,					// Candidate signal number
			int *cand_buffer_count,
			const char* outname
       );


void modvir (double sinal, double cosal, double sindel, double cosdel,	\
        double sphir, double cphir, double *a, double *b, int Np, Ampl_mod_coeff *amod,
        Arrays *arr);


void compute_sincosmodf_wrapper(double* cu_sinmodf, double* cu_cosmodf, double omr, int N);


void save_candidates(FLOAT_TYPE* cu_cand_buffer, FLOAT_TYPE* cand_buffer, int *cand_buffer_count, const char* outname);

void norm_Fstat_nwn(FLOAT_TYPE *cu_F, FLOAT_TYPE *mu, int N, int nav);

void modvir_gpu (double sinal, double cosal, double sindel, double cosdel,
		double sphir, double cphir, double *cu_a, double *cu_b, int N, Arrays *arr);


void FStat_gpu(FLOAT_TYPE *cu_F, int N, int nav, FLOAT_TYPE *cu_mu, FLOAT_TYPE *cu_mu_t);

template<typename T>
struct Square
{
 __host__ __device__ __forceinline__
  T operator()(const T& a) const {
    return a*a;
  }
};

#endif
