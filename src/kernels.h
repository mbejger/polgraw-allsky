#ifndef __KERNELS_H__
#define __KERNELS_H__

#include "floats.h"

void copy_amod_coeff(Ampl_mod_coeff *amod);


__global__ void compute_sincosmodf(double *s, double *c, double omr, int N);

__global__ void shift_time_mod(double shft1, double het0, double ns0, double ns1,
			       double ns2,  double *xDat, cufftDoubleComplex *xa,
			       cufftDoubleComplex *xb, FLOAT_TYPE* shft,
			       double* shftf, double* tshift, double* aa, double* bb,
			       double* DetSSB, double oms, int N, int nfft, int interpftpad);


__global__ void scale_fft(cufftDoubleComplex *xa, cufftDoubleComplex *xb,
			  double ft, int Ninterp);


__global__ void resample_postfft(cufftDoubleComplex *xa, cufftDoubleComplex *xb,
				 int nfft, int Ninterp, int nyqst);


__global__ void double_to_float(cufftDoubleComplex *xa, cufftDoubleComplex *xb,
				COMPLEX_TYPE *xa_f, COMPLEX_TYPE *xb_f, int N);


__global__ void phase_mod_2(COMPLEX_TYPE *xa, COMPLEX_TYPE *xb,
			    COMPLEX_TYPE *xar, COMPLEX_TYPE *xbr,
			    FLOAT_TYPE het1, FLOAT_TYPE sgnlt1, FLOAT_TYPE *shft,
			    int N);


__global__ void compute_Fstat(COMPLEX_TYPE *xa, COMPLEX_TYPE *xb,
			      FLOAT_TYPE *F, FLOAT_TYPE crf0, int N);


__global__ void kernel_norm_Fstat_wn(FLOAT_TYPE *F, FLOAT_TYPE sig2, int N);


__global__ void check_array(double *x, int start, int end);


__global__ void check_array(cufftDoubleComplex *x, int start, int end);


__global__ void find_candidates(FLOAT_TYPE *F, FLOAT_TYPE *params, int *found, FLOAT_TYPE val,
				int nmin, int nmax, double fftpad, double nfft, FLOAT_TYPE sgnl0,
				int ndf, FLOAT_TYPE sgnl1, FLOAT_TYPE sgnl2, FLOAT_TYPE sgnl3);

__global__ void copy_candidates(FLOAT_TYPE *params, FLOAT_TYPE *buffer, int N);

__global__ void pad_zeros(COMPLEX_TYPE *xa, COMPLEX_TYPE *xb, int N);


__global__ void compute_modvir(double *aa, double *bb, double cosalfr, double sinalfr,
			       double c2d, double c2sd, double *sinmodf,
			       double *cosmodf, double sindel, double cosdel,
			       int N);


__global__ void modvir_normalize(double *aa, double *bb, double s_a, double s_b, int N);




//first sum: sum squares of aa and bb
template <unsigned int blockSize>
__global__ void reduce_modvir_sumsq(double *aa, double *bb, double *o_aa, double *o_bb,
				    unsigned int n) {
  __shared__ double sdata_a[blockSize];
  __shared__ double sdata_b[blockSize];

  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*(blockSize*2) + tid;
  unsigned int gridSize = blockSize*2*gridDim.x;
  sdata_a[tid] = 0;
  sdata_b[tid] = 0;

  while (i < n) {
    sdata_a[tid] += aa[i]*aa[i] + aa[i+blockSize]*aa[i+blockSize];
    sdata_b[tid] += bb[i]*bb[i] + bb[i+blockSize]*bb[i+blockSize];
    i += gridSize;
  }
  __syncthreads();

  if (blockSize >= 512) { if (tid < 256) { sdata_a[tid] += sdata_a[tid + 256];
      sdata_b[tid] += sdata_b[tid + 256]; } __syncthreads(); }
  if (blockSize >= 256) { if (tid < 128) { sdata_a[tid] += sdata_a[tid + 128];
      sdata_b[tid] += sdata_b[tid + 128]; } __syncthreads(); }
  if (blockSize >= 128) { if (tid < 64) { sdata_a[tid] += sdata_a[tid + 64];
      sdata_b[tid] += sdata_b[tid + 64]; } __syncthreads(); }

  if (tid < 32) {
    if (blockSize >= 64) { sdata_a[tid] += sdata_a[tid + 32];
      sdata_b[tid] += sdata_b[tid + 32]; }
    if (blockSize >= 32) { sdata_a[tid] += sdata_a[tid + 16];
      sdata_b[tid] += sdata_b[tid + 16]; }
    if (blockSize >= 16) { sdata_a[tid] += sdata_a[tid + 8];
      sdata_b[tid] += sdata_b[tid + 8]; }
    if (blockSize >= 8) { sdata_a[tid] += sdata_a[tid + 4];
      sdata_b[tid] += sdata_b[tid + 4]; }
    if (blockSize >= 4) { sdata_a[tid] += sdata_a[tid + 2];
      sdata_b[tid] += sdata_b[tid + 2]; }
    if (blockSize >= 2) { sdata_a[tid] += sdata_a[tid + 1];
      sdata_b[tid] += sdata_b[tid + 1]; }
  }

  if (tid == 0) {
    o_aa[blockIdx.x] = sdata_a[0];
    o_bb[blockIdx.x] = sdata_b[0];
  }
}

/*
//second sum for modvir
template <unsigned int blockSize>
__global__ void reduce_modvir_sum2(double *aa, double *bb, double *o_aa, double *o_bb, unsigned int n) {
__shared__ double sdata_a[blockSize];
__shared__ double sdata_b[blockSize];

unsigned int tid = threadIdx.x;
unsigned int i = blockIdx.x*(blockSize*2) + tid;
unsigned int gridSize = blockSize*2*gridDim.x;
sdata_a[tid] = 0;
sdata_b[tid] = 0;

while (i < n) {
sdata_a[tid] += aa[i] + aa[i+blockSize];
sdata_b[tid] += bb[i] + bb[i+blockSize];
i += gridSize;
}
__syncthreads();

if (blockSize >= 512) { if (tid < 256) { sdata_a[tid] += sdata_a[tid + 256]; sdata_b[tid] += sdata_b[tid + 256]; } __syncthreads(); }
if (blockSize >= 256) { if (tid < 128) { sdata_a[tid] += sdata_a[tid + 128]; sdata_b[tid] += sdata_b[tid + 128]; } __syncthreads(); }
if (blockSize >= 128) { if (tid < 64) { sdata_a[tid] += sdata_a[tid + 64]; sdata_b[tid] += sdata_b[tid + 64]; } __syncthreads(); }

if (tid < 32) {
if (blockSize >= 64) { sdata_a[tid] += sdata_a[tid + 32]; sdata_b[tid] += sdata_b[tid + 32]; }
if (blockSize >= 32) { sdata_a[tid] += sdata_a[tid + 16]; sdata_b[tid] += sdata_b[tid + 16]; }
if (blockSize >= 16) { sdata_a[tid] += sdata_a[tid + 8]; sdata_b[tid] += sdata_b[tid + 8]; }
if (blockSize >= 8) { sdata_a[tid] += sdata_a[tid + 4]; sdata_b[tid] += sdata_b[tid + 4]; }
if (blockSize >= 4) { sdata_a[tid] += sdata_a[tid + 2]; sdata_b[tid] += sdata_b[tid + 2]; }
if (blockSize >= 2) { sdata_a[tid] += sdata_a[tid + 1]; sdata_b[tid] += sdata_b[tid + 1]; }
}

if (tid == 0) {
o_aa[blockIdx.x] = sdata_a[0];
o_bb[blockIdx.x] = sdata_b[0];
}
}

*/
__global__ void modvir_normalize(double *aa, double *bb, double s_a, double s_b, int N);


__global__ void fstat_sum(float *F, float *out, int N);
__global__ void fstat_norm(float *F, float *mu, int N, int nav);

__global__ void reduction(float * __restrict__ indata, float * __restrict__ outdata, int N);
//__global__ void reduction_sum(float *in, float *out, int N);

__global__ void interbinning_gap(COMPLEX_TYPE *xa, COMPLEX_TYPE *xb,
				 COMPLEX_TYPE *xa_t, COMPLEX_TYPE *xb_t,
				 int nfft);

__global__ void interbinning_interp(COMPLEX_TYPE *xa, COMPLEX_TYPE *xb,
				    int nfft);


#endif
