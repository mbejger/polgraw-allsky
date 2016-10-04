#include "struct.h"
#include "settings.h"
#include "jobcore.h"
#include <stdio.h>

// maybe it should be dynamic...
__constant__ Ampl_mod_coeff amod_d[MAX_DETECTORS];


void copy_amod_coeff(int nifo) {
  int i;
  Ampl_mod_coeff amod_coeff_tmp[nifo];
  for(i=0; i<nifo; ++i){
    amod_coeff_tmp[i] = ifo[i].amod;
  }
  cudaMemcpyToSymbol(amod_d, amod_coeff_tmp, sizeof(Ampl_mod_coeff)*nifo, 
		     0, cudaMemcpyHostToDevice);
}


__global__ void modvir_kern(double *aa_d, double *bb_d, double cosalfr, double sinalfr,
			    double c2d, double c2sd, 
			    double *sinmodf_d, double *cosmodf_d, 
			    double sindel, double cosdel, int Np, int idet) {

  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx<Np) {
    double c = cosalfr * cosmodf_d[idx] + sinalfr * sinmodf_d[idx];
    double s = sinalfr * cosmodf_d[idx] - cosalfr * sinmodf_d[idx];
    double c2s = 2.*c*c;
    double cs = c*s;

    aa_d[idx] = amod_d[idet].c1*(2.-c2d)*c2s + amod_d[idet].c2*(2.-c2d)*2.*cs + 
      amod_d[idet].c3*c2sd*c + amod_d[idet].c4*c2sd*s - amod_d[idet].c1*(2.-c2d) + amod_d[idet].c5*c2d;
    bb_d[idx] = amod_d[idet].c6*sindel*c2s + amod_d[idet].c7*sindel*2.*cs + 
      amod_d[idet].c8*cosdel*c + amod_d[idet].c9*cosdel*s - amod_d[idet].c6*sindel;
  }
}


__global__ void tshift_pmod_kern(double shft1, double het0, 
				 double ns0, double ns1, double ns2,  
				 double *xDat_d, 
				 cufftDoubleComplex *xa_d, cufftDoubleComplex *xb_d, 
				 FLOAT_TYPE *shft_d, double *shftf_d, 
				 double *tshift_d, 
				 double *aa_d, double *bb_d,
				 double *DetSSB_d, 
				 double oms, int N, int nfft, int interpftpad) {

  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < N) {
    double S = ns0 * DetSSB_d[i*3]
             + ns1 * DetSSB_d[i*3+1] 
             + ns2 * DetSSB_d[i*3+2];
    shft_d[i] = S;
    shftf_d[i]= S - shft1;

    /* phase mod */
    // dlaczego - ?
    double phase = -het0*i - oms * S;
    double c = cos(phase), s = sin(phase);
    xa_d[i].x = xDat_d[i] * aa_d[i] * c;
    xa_d[i].y = xDat_d[i] * aa_d[i] * s;
    xb_d[i].x = xDat_d[i] * bb_d[i] * c;
    xb_d[i].y = xDat_d[i] * bb_d[i] * s;

    //calculate time positions for spline interpolation
    tshift_d[i] = interpftpad * ( i - shftf_d[i] );
    // no need for this on gpu
    //_tmp1[n][i] = aux->t2[i] + (double)(2*i)*ifo[n].sig.shft[i]; 
  } else if (i < nfft) {
    xa_d[i].x = xa_d[i].y = xb_d[i].x = xb_d[i].y = 0.;
  }
}


__global__ void resample_postfft(cufftDoubleComplex *xa_d, cufftDoubleComplex *xb_d,
                                 int nfft, int Ninterp, int nyqst) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  //move frequencies from second half of spectrum; loop length: nfft - nyqst =
  // = nfft - nfft/2 - 1 = nfft/2 - 1
  if (idx < nfft/2 - 1) {
    int i = nyqst + Ninterp - nfft + idx;
    int j = nyqst + idx;
    xa_d[i].x=xa_d[j].x;
    xa_d[i].y=xa_d[j].y;
    xb_d[i].x=xb_d[j].x;
    xb_d[i].y=xb_d[j].y;
  }

  //zero frequencies higher than nyquist, length: Ninterp - nfft
  //loop length: Ninterp - nfft ~ nfft
  if (idx < Ninterp - nfft) {
    xa_d[nyqst+idx].x = xa_d[nyqst+idx].y = 0.;
    xb_d[nyqst+idx].x = xb_d[nyqst+idx].y = 0.;
  }
}



__global__ void compute_sincosmodf(double *s, double *c, double omr, int N) {
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx<N) {
    s[idx] = sin(omr*idx);
    c[idx] = cos(omr*idx);
  }
}



__global__ void phase_mod_1(cufftDoubleComplex *xa, cufftDoubleComplex *xb,
                            cufftDoubleComplex *xar, cufftDoubleComplex *xbr,
                            double het1, double sgnlt1, double *shft,
                            int N){

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < N) {
    //          FLOAT_TYPE phase = - ( het1*idx + sgnlt1 * ( (double)idx*idx + 2 * idx * shft[idx] ) );
    double phase = - idx * ( het1 + sgnlt1 * ( idx + 2 * shft[idx] ) );
    double s, c;

    sincos(phase, &s, &c);

    xa[idx].x = xar[idx].x*c - xar[idx].y*s;
    xa[idx].y = xar[idx].x*s + xar[idx].y*c;

    xb[idx].x = xbr[idx].x*c - xbr[idx].y*s;
    xb[idx].y = xbr[idx].x*s + xbr[idx].y*c;
    //    if (idx==1) printf("xa=(%f, %f)\n", xa[idx].x, xa[idx].y);

  }
}

__global__ void phase_mod_2(cufftDoubleComplex *xa, cufftDoubleComplex *xb,
                            cufftDoubleComplex *xar, cufftDoubleComplex *xbr,
                            double het1, double sgnlt1, double *shft,
                            int N){

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < N) {
    //          FLOAT_TYPE phase = - ( het1*idx + sgnlt1 * ( (double)idx*idx + 2 * idx * shft[idx] ) );
    double phase = - idx * ( het1 + sgnlt1 * ( idx + 2 * shft[idx] ) );
    double s, c;

    sincos(phase, &s, &c);

    xa[idx].x += xar[idx].x*c - xar[idx].y*s;
    xa[idx].y += xar[idx].x*s + xar[idx].y*c;

    xb[idx].x += xbr[idx].x*c - xbr[idx].y*s;
    xb[idx].y += xbr[idx].x*s + xbr[idx].y*c;
    //    if (idx==1) printf("xa=(%f, %f)\n", xa[idx].x, xa[idx].y);

  }
}

extern  __constant__  double maa_d, mbb_d;

__global__ void compute_Fstat(cufftDoubleComplex *xa, cufftDoubleComplex *xb,
                              double *F, int N) { // N = nmax - nmin
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  //  if (i==0) printf("maa_d=%f\n", maa_d);
  if (i < N) {
    F[i] = ( xa[i].x*xa[i].x + xa[i].y*xa[i].y)/maa_d + (xb[i].x*xb[i].x + xb[i].y*xb[i].y)/mbb_d;
  }
}


/*
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
*/

__global__ void fstat_norm_simple(FLOAT_TYPE *F_d, int nav) {
  //  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int i;
  FLOAT_TYPE *fr = F_d + blockIdx.x*nav;
  FLOAT_TYPE mu = 0.;
  for (i=0; i<nav; ++i)
    mu += *fr++;
  mu /= 2.*nav;
  fr = F_d + blockIdx.x*nav;
  for (i=0; i<nav; i++)
    *fr++ /= mu;

}


// parameters are:
// [frequency, spindown, position1, position2, snr]
#define ADD_PARAMS_MACRO						\
  int p = atomicAdd(found, 1);						\
  params[p*NPAR + 0] = 2.0*M_PI*(idx)*fftpad*nfft+sgnl0;		\
  params[p*NPAR + 1] = sgnl1;						\
  params[p*NPAR + 2] = sgnl2;						\
  params[p*NPAR + 3] = sgnl3;						\
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



//---------------------------------------------------------------

//second reduction used in fstat
__global__ void reduction_sum(FLOAT_TYPE *in, FLOAT_TYPE *out, int N) {
  extern __shared__ FLOAT_TYPE sf_data[];

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

  if (tid==0) out[blockIdx.x] = 1.0f/sf_data[0];
}


__global__ void fstat_norm(FLOAT_TYPE *F, FLOAT_TYPE *mu, int N, int nav) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < N) {
    int block = i/nav; //block index
    F[i] *= 2*nav * mu[block];
  }
}
