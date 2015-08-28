#include "struct.h"
#include "settings.h"
#include <stdio.h>

//extern __constant__ Detector_settings *c_sett;
__constant__ double cu_c[9];
//extern double *cu_c;

void copy_amod_coeff(Ampl_mod_coeff *amod) {
  double amod_coeff_tmp[] = { amod->c1,
			      amod->c2,
			      amod->c3,
			      amod->c4,
			      amod->c5,
			      amod->c6,
			      amod->c7,
			      amod->c8,
			      amod->c9 };
  cudaMemcpyToSymbol(cu_c, amod_coeff_tmp, sizeof(double)*9);
}


__global__ void check_array(double *x, int start, int end) {
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx>=start && idx < end) {
    printf("block %d, thread %d, idx %d, value %lf\n",blockIdx.x, threadIdx.x, idx, x[idx]);
  }
}


/* not used */
__global__ void compute_sincosmodf(double *s, double *c, double omr, int N) {
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx<N) {
    s[idx] = sin(omr*idx);
    c[idx] = cos(omr*idx);
  }
}


__global__ void shift_time_mod(double shft1, double het0, double ns0, double ns1,
			       double ns2,  double *xDat, cufftDoubleComplex *xa, 
			       cufftDoubleComplex *xb, FLOAT_TYPE* shft,
			       double* shftf, double* tshift, double* aa, double* bb,
			       double* DetSSB, double oms, int N, int nfft, int interpftpad) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < N) {
    double S = ns0 * DetSSB[i*3+0] + ns1 * DetSSB[i*3+1] + ns2 * DetSSB[i*3+2];
    shft[i] = S;
    shftf[i]= S - shft1;
    
    /* phase mod */
    double phase = -het0*i - oms * S;
    double c = cos(phase), s = sin(phase);
    xa[i].x = xDat[i] * aa[i] * c;
    xa[i].y = xDat[i] * aa[i] * s;
    xb[i].x = xDat[i] * bb[i] * c;
    xb[i].y = xDat[i] * bb[i] * s;
    
    //calculate time positions for spline interpolation
    tshift[i] = interpftpad * ( i - shftf[i] );
  } else if (i < nfft) {
    xa[i].x = xa[i].y = xb[i].x = xb[i].y = 0;
  }
}


__global__ void scale_fft(cufftDoubleComplex *xa, cufftDoubleComplex *xb,
			  double ft, int Ninterp) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (idx < Ninterp) {
    xa[idx].x*=ft;
    xa[idx].y*=ft;
    xb[idx].x*=ft;
    xb[idx].y*=ft;
  }
}


__global__ void resample_postfft(cufftDoubleComplex *xa, cufftDoubleComplex *xb,
				 int nfft, int Ninterp, int nyqst) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  
  //move frequencies from second half of spectrum; loop length: nfft - nyqst =
  // = nfft - nfft/2 - 1 = nfft/2 - 1
  if (idx < nfft/2 - 1) {
    int i = nyqst + Ninterp - nfft + idx;
    int j = nyqst + idx;
    xa[i].x=xa[j].x;
    xa[i].y=xa[j].y;
    xb[i].x=xb[j].x;
    xb[i].y=xb[j].y;
  }
  
  //zero frequencies higher than nyquist, length: Ninterp - nfft
  //loop length: Ninterp - nfft ~ nfft
  if (idx < Ninterp - nfft) {
    xa[nyqst + idx].x = xa[nyqst+idx].y = 0.0;
    xb[nyqst + idx].x = xb[nyqst+idx].y = 0.0;
  }
}


__global__ void double_to_float(cufftDoubleComplex *xa, cufftDoubleComplex *xb,
				COMPLEX_TYPE *xa_f, COMPLEX_TYPE *xb_f, int N)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < N) {
    xa_f[idx].x = xa[idx].x;
    xa_f[idx].y = xa[idx].y;
    xb_f[idx].x = xb[idx].x;
    xb_f[idx].y = xb[idx].y;
  }
}


__global__ void phase_mod_2(COMPLEX_TYPE *xa, COMPLEX_TYPE *xb,
			    COMPLEX_TYPE *xar, COMPLEX_TYPE *xbr,
			    FLOAT_TYPE het1, FLOAT_TYPE sgnlt1, FLOAT_TYPE *shft, 
			    int N)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < N) {
    //		FLOAT_TYPE phase = - ( het1*idx + sgnlt1 * ( (double)idx*idx + 2 * idx * shft[idx] ) );
    FLOAT_TYPE phase = - idx * ( het1 + sgnlt1 * ( idx + 2 * shft[idx] ) );
    FLOAT_TYPE s, c;    
#ifdef COMP_FLOAT
    //c = __cosf(phase); 
    //s = __sinf(phase);
    sincosf(phase, &s, &c);
#else
    //c = cos(phase);
    //s = sin(phase);
    sincos(phase, &s, &c);
#endif
    
    //complex multiplication:
    // ( xar[idx].x + i* xar[idx].y ) * ( c + i * s ) = 
    // = ( xar[idx].x * c - xar[idx].y * s ) + i * ( s * xar[idx].x + c * xar[idx].y )
    xa[idx].x = xar[idx].x*c - xar[idx].y*s;
    xa[idx].y = xar[idx].x*s + xar[idx].y*c;
    
    xb[idx].x = xbr[idx].x*c - xbr[idx].y*s;
    xb[idx].y = xbr[idx].x*s + xbr[idx].y*c;
  }
}


__global__ void pad_zeros(COMPLEX_TYPE *xa, COMPLEX_TYPE *xb, int N) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < N) {
    xa[idx].x=0.f;
    xa[idx].y=0.f;
    xb[idx].x=0.f;
    xb[idx].y=0.f;
  }
}


//rewrite odd frequencies i to 2*i
//i.e. leave a gap between every bin
__global__ void interbinning_gap(COMPLEX_TYPE *xa, COMPLEX_TYPE *xb,
				 COMPLEX_TYPE *xa_t, COMPLEX_TYPE *xb_t,
				 int nfft)
{
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx < nfft/2 ) {
    xa_t[2*idx].x = xa[idx].x;
    xa_t[2*idx].y = xa[idx].y;
    xb_t[2*idx].x = xb[idx].x;
    xb_t[2*idx].y = xb[idx].y;
  }
}


#define _RSQRT_2_ 0.70710678118654746172f
//interpolate with interbinning
__global__ void interbinning_interp(COMPLEX_TYPE *xa, COMPLEX_TYPE *xb,
				    int nfft)
{
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx < (nfft-2)/2 ) {
    int i = 2*idx+1;
    xa[i].x = (xa[i-1].x - xa[i+1].x) * _RSQRT_2_; // ./sqrt(2)
    xa[i].y = (xa[i-1].y - xa[i+1].y) * _RSQRT_2_; // ./sqrt(2)
    xb[i].x = (xb[i-1].x - xb[i+1].x) * _RSQRT_2_; // ./sqrt(2)
    xb[i].y = (xb[i-1].y - xb[i+1].y) * _RSQRT_2_; // ./sqrt(2)
  }
}

__global__ void compute_Fstat(COMPLEX_TYPE *xa, COMPLEX_TYPE *xb,
			      FLOAT_TYPE *F, FLOAT_TYPE crf0, int N) { // N = nmax - nmin
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < N) {
    F[i] = ( xa[i].x*xa[i].x + xa[i].y*xa[i].y + xb[i].x*xb[i].x + xb[i].y*xb[i].y) * (crf0);
  }
}

__global__ void kernel_norm_Fstat_wn(FLOAT_TYPE *F, FLOAT_TYPE sig2, int N) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < N) {
    F[idx] *= sig2;
  }	
}


// This may be confusing, but three "if"s in kernel below
// are faster than one "if" with one, big, combined condition.
// Operations in every "if" are the same, so they are wrapped in macro.

// parameters are:
// [frequency, spindown, position1, position2, snr]
#define ADD_PARAMS_MACRO \
  int p = atomicAdd(found, 1);						\
  params[p*NPAR + 0] = 2.0*M_PI*(idx)*fftpad*nfft+sgnl0;	\
  params[p*NPAR + 1] = sgnl1;	\
  params[p*NPAR + 2] = sgnl2;	\
  params[p*NPAR + 3] = sgnl3;	\
  params[p*NPAR + 4] = sqrt(2*(F[idx]-ndf));


__global__ void find_candidates(FLOAT_TYPE *F, FLOAT_TYPE *params, int *found, FLOAT_TYPE val,
				int nmin, int nmax, double fftpad, double nfft, FLOAT_TYPE sgnl0, int ndf,
				FLOAT_TYPE sgnl1, FLOAT_TYPE sgnl2, FLOAT_TYPE sgnl3) {
  
  int idx = blockIdx.x * blockDim.x + threadIdx.x + nmin;
  
  if (idx > nmin && idx < nmax && F[idx] >= val && F[idx] > F[idx+1] && F[idx] > F[idx-1]) {
    ADD_PARAMS_MACRO
      } else if (idx == nmin && F[idx] >= val && F[idx] > F[idx+1]) {
    ADD_PARAMS_MACRO
      } else if (idx == nmax-1 && F[idx] >= val && F[idx] > F[idx-1]) {
    ADD_PARAMS_MACRO
      } 
}


//temporary
//not used
__global__ void copy_candidates(FLOAT_TYPE *params, FLOAT_TYPE *buffer, int N) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < N) {
    buffer[idx] = params[idx];
  }
}


__global__ void reduction_sumsq(double *in, double *out, int N) {
  extern __shared__ double sdata[];
  
  int tid = threadIdx.x;
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  
  sdata[tid] = (i<N) ? in[i]*in[i] : 0;
  
  __syncthreads();
	
  for (int s = blockDim.x/2; s>0; s>>=1) {
    if (tid < s) {
      sdata[tid] += sdata[tid + s];
    }
    __syncthreads();
  }
  
  if (tid==0) out[blockIdx.x] = sdata[0];
}



__global__ void reduction_sum(float *in, float *out, int N) {
  extern __shared__ float sf_data[];
  
  int tid = threadIdx.x;
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  
  sf_data[tid] = (i<N) ? in[i] : 0;
  
  __syncthreads();
  
  for (int s = blockDim.x/2; s>0; s>>=1) {
    if (tid < s) {
      sf_data[tid] += sf_data[tid + s];
    }
    __syncthreads();
  }
  
  if (tid==0) out[blockIdx.x] = sf_data[0];
}


__global__ void reduction_sum(double *in, double *out, int N) {
  extern __shared__ double sd_data[];
  
  int tid = threadIdx.x;
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  
  sd_data[tid] = (i<N) ? in[i] : 0;
  
  __syncthreads();
  
  for (int s = blockDim.x/2; s>0; s>>=1) {
    if (tid < s) {
      sd_data[tid] += sd_data[tid + s];
    }
    __syncthreads();
  }
  
  if (tid==0) out[blockIdx.x] = sd_data[0];
  
}


__global__ void fstat_norm(float *F, float *mu, int N, int nav) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < N) {
    int block = i/nav; //block index
    F[i] *= 2*nav / mu[block];
  }
}


__global__ void compute_modvir(double *aa, double *bb, double cosalfr, double sinalfr,
			       double c2d, double c2sd, double *sinmodf,
			       double *cosmodf, double sindel, double cosdel,
			       int N) {
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx<N) {
    double c = cosalfr * cosmodf[idx] + sinalfr * sinmodf[idx];
    double s = sinalfr * cosmodf[idx] - cosalfr * sinmodf[idx];
    double c2s = 2.*c*c;
    double cs = c*s;
    aa[idx] = cu_c[0]*(2.-c2d)*c2s + cu_c[1]*(2.-c2d)*2.*cs + cu_c[2]*c2sd*c + cu_c[3]*c2sd*s-
      cu_c[0]*(2.-c2d) + cu_c[4]*c2d;
    bb[idx] = cu_c[5]*sindel*c2s + cu_c[6]*sindel*2.*cs + cu_c[7]*cosdel*c + cu_c[8]*cosdel*s-
      cu_c[5]*sindel;
  }
}


__global__ void modvir_normalize(double *aa, double *bb, double s_a, double s_b, int N) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if ( idx < N ) { 
    aa[idx] *= s_a;
    bb[idx] *= s_b;
  }
}
