#ifndef __KERNELS_HCL__
#define __KERNELS_HCL__

#include <floats.hcl>


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
__kernel void resample_postfft(complex_t *xa_d,
                               complex_t *xb_d,
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

//__global__ void phase_mod_1(cufftDoubleComplex *xa, cufftDoubleComplex *xb,
//                            cufftDoubleComplex *xar, cufftDoubleComplex *xbr,
//                            double het1, double sgnlt1, double *shft,
//                            int N);
//
//__global__ void phase_mod_2(cufftDoubleComplex *xa, cufftDoubleComplex *xb,
//                            cufftDoubleComplex *xar, cufftDoubleComplex *xbr,
//                            double het1, double sgnlt1, double *shft,
//                            int N);
//
//__global__ void compute_Fstat(cufftDoubleComplex *xa, cufftDoubleComplex *xb,
//                              double *F, 
//			      int N);
//
//__global__ void fstat_norm_simple(FLOAT_TYPE *F_d, int nav);
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
