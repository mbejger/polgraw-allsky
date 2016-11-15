#ifndef __KERNELS_HCL__
#define __KERNELS_HCL__

#include <floats.hcl>

/// <summar>Amplitude modulation function coefficients</summary>
///
typedef struct _ampl_mod_coeff
{
    real_t c1, c2, c3, c4, c5, c6, c7, c8, c9;

} Ampl_mod_coeff;


/// <summary>The purpose of this function was undocumented.</summary>
///
__kernel void modvir_kern(__global real_t* aa_d,
                          __global real_t* bb_d,
                          real_t cosalfr,
                          real_t sinalfr,
                          real_t c2d,
                          real_t c2sd,
                          __global real_t* sinmodf_d,
                          __global real_t* cosmodf_d,
                          real_t sindel,
                          real_t cosdel,
                          int Np,
                          int idet,
                          __constant Ampl_mod_coeff* amod_d);

/// <summary>The purpose of this function was undocumented.</summary>
///
__kernel void tshift_pmod_kern(real_t shft1,
                               real_t het0,
                               real_t ns0,
                               real_t ns1,
                               real_t ns2,
                               __global real_t* xDat_d,
                               __global complex_t* xa_d,
                               __global complex_t* xb_d,
                               __global real_t* shft_d,
                               __global real_t* shftf_d,
                               __global real_t* tshift_d,
                               __global real_t* aa_d,
                               __global real_t* bb_d,
                               __global real_t* DetSSB_d,
                               real_t oms,
                               int N,
                               int nfft,
                               int interpftpad);

/// <summary>Shifts frequencies and remove those over Nyquist.</summary>
///
__kernel void resample_postfft(__global complex_t *xa_d,
                               __global complex_t *xb_d,
                               int nfft,
                               int Ninterp,
                               int nyqst);

/// <summary>Computes sin and cos values and stores them in an array.</summary>
/// <remarks>Most likely a very bad idea. Results are used in modvir and should be computed there in place.</remarks>
///
__kernel void compute_sincosmodf(__global real_t* s,
                                 __global real_t* c,
                                 real_t omr,
                                 int N);

/// <summary>The purpose of this function was undocumented.</summary>
///
__kernel void computeB(__global complex_t* y,
                       __global complex_t* B,
                       int N);

/// <summary>Multiplies the tridiagonal matrix specified by <c>{dl, d, du}</c> with dense vector <c>x</c>.</summary>
///
__kernel void tridiagMul(__global real_t* dl,
                         __global real_t* d,
                         __global real_t* du,
                         __global complex_t* x,
                         __global complex_t* y);

/// <summary>The purpose of this function was undocumented.</summary>
///
__kernel void interpolate(__global real_t* new_x,
                          __global complex_t* new_y,
                          __global complex_t* z,
                          __global complex_t* y,
                          int N,
                          int new_N);


/// <summary>The purpose of this function was undocumented.</summary>
///
__kernel void phase_mod_1(__global complex_t* xa,
                          __global complex_t* xb,
                          __global complex_t* xar,
                          __global complex_t* xbr,
                          real_t het1,
                          real_t sgnlt1,
                          __global real_t* shft,
                          int N);

/// <summary>The purpose of this function was undocumented.</summary>
///
__kernel void phase_mod_2(__global complex_t* xa,
                          __global complex_t* xb,
                          __global complex_t* xar,
                          __global complex_t* xbr,
                          real_t het1,
                          real_t sgnlt1,
                          __global real_t* shft,
                          int N);

/// <summary>Compute F-statistics.</summary>
/// 
__kernel void compute_Fstat(__global complex_t* xa,
                            __global complex_t* xb,
                            __global real_t* F,
                            __constant real_t* maa_d,
                            __constant real_t* mbb_d,
                            int N);

/// <summary>Compute F-statistics.</summary>
/// <precondition>ssi less than or equal to nav</precondition>
/// <precondition>lsi less than or equal to ssi</precondition>
/// <precondition>lsi be an integer power of 2</precondition>
/// 
__kernel void fstat_norm_simple(__global real_t* F,
                                __local real_t* shared,
                                unsigned int ssi,
                                unsigned int nav);

//__global__ void fstat_norm(FLOAT_TYPE *F, FLOAT_TYPE *mu, int N, int nav);


//first reduction used in fstat
//template <unsigned int blockSize>
//__device__ void warpReduce(volatile FLOAT_TYPE *sdata, unsigned int tid) {
//  if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
//  if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
//  if (blockSize >= 16) sdata[tid] += sdata[tid + 8];
//  if (blockSize >= 8) sdata[tid] += sdata[tid + 4];
//  if (blockSize >= 4) sdata[tid] += sdata[tid + 2];
//  if (blockSize >= 2) sdata[tid] += sdata[tid + 1];
//}
//template <unsigned int blockSize>
//__global__ void reduction_sum(FLOAT_TYPE *g_idata, FLOAT_TYPE *g_odata, unsigned int n) {
//  extern __shared__ volatile FLOAT_TYPE sdata[];
//  unsigned int tid = threadIdx.x;
//  unsigned int i = blockIdx.x*(blockSize*2) + tid;
//  unsigned int gridSize = blockSize*2*gridDim.x;
//  sdata[tid] = 0;
//  while (i < n) { sdata[tid] += g_idata[i] + g_idata[i+blockSize]; i += gridSize; }
//  __syncthreads();
//  if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
//  if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
//  if (blockSize >= 128) { if (tid < 64) { sdata[tid] += sdata[tid + 64]; } __syncthreads(); }
//  if (tid < 32) warpReduce<blockSize>(sdata, tid);
//  if (tid == 0) g_odata[blockIdx.x] = sdata[0];
//}
//
//__global__ void reduction_sum(float *in, float *out, int N);
//__global__ void reduction_sum(double *in, double *out, int N);


#endif // __KERNELS_HCL__
