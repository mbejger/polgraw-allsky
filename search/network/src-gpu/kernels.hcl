#ifndef __KERNELS_H__
#define __KERNELS_H__

//void copy_amod_coeff(int nifo);
//
//__global__ void modvir_kern(double *aa_d, double *bb_d, double cosalfr, double sinalfr,
//			    double c2d, double c2sd, 
//			    double *sinmodf_d, double *cosmodf_d, 
//			    double sindel, double cosdel, int Np, int idet);
//
//__global__ void tshift_pmod_kern(double shft1, double het0, 
//				 double ns0, double ns1, double ns2,  
//				 double *xDat_d, 
//				 cufftDoubleComplex *xa_d, cufftDoubleComplex *xb_d, 
//				 double *shft_d, double *shftf_d, 
//				 double *tshift_d, 
//				 double *aa_d, double *bb_d,
//				 double *DetSSB_d, 
//				 double oms, int N, int nfft, int interpftpad);
//
//__global__ void resample_postfft(cufftDoubleComplex *xa_d, cufftDoubleComplex *xb_d,
//                                 int nfft, int Ninterp, int nyqst);

__kernel void compute_sincosmodf(__global double *s,
                                 __global double *c,
                                 double omr,
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


#endif
