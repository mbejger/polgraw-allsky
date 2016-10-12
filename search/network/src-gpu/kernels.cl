// Polgraw includes
#include <floats.hcl>       // real_t, complex_t
#include <kernels.hcl>      // function declarations


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
                          __constant Ampl_mod_coeff* amod_d)
{
    size_t idx = get_global_id(0);

    real_t c = cosalfr * cosmodf_d[idx] + sinalfr * sinmodf_d[idx];
    real_t s = sinalfr * cosmodf_d[idx] - cosalfr * sinmodf_d[idx];
    real_t c2s = 2.*c*c;
    real_t cs = c*s;

    aa_d[idx] = amod_d[idet].c1*(2. - c2d)*c2s + amod_d[idet].c2*(2. - c2d)*2.*cs +
        amod_d[idet].c3*c2sd*c + amod_d[idet].c4*c2sd*s - amod_d[idet].c1*(2. - c2d) + amod_d[idet].c5*c2d;
    bb_d[idx] = amod_d[idet].c6*sindel*c2s + amod_d[idet].c7*sindel*2.*cs +
        amod_d[idet].c8*cosdel*c + amod_d[idet].c9*cosdel*s - amod_d[idet].c6*sindel;
}

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
                               int interpftpad)
{
    size_t i = get_global_id(0);

    if (i < N)
    {
        real_t S = ns0 * DetSSB_d[i * 3]
            + ns1 * DetSSB_d[i * 3 + 1]
            + ns2 * DetSSB_d[i * 3 + 2];
        shft_d[i] = S;
        shftf_d[i] = S - shft1;

        /* phase mod */
        // dlaczego - ?
        real_t phase = -het0*i - oms * S;
        real_t c = cos(phase), s = sin(phase);
        xa_d[i].x = xDat_d[i] * aa_d[i] * c;
        xa_d[i].y = xDat_d[i] * aa_d[i] * s;
        xb_d[i].x = xDat_d[i] * bb_d[i] * c;
        xb_d[i].y = xDat_d[i] * bb_d[i] * s;

        //calculate time positions for spline interpolation
        tshift_d[i] = interpftpad * (i - shftf_d[i]);
    }
    else if (i < nfft)
    {
        xa_d[i].x = xa_d[i].y = xb_d[i].x = xb_d[i].y = 0.;
    }
}

/// <summary>Shifts frequencies and remove those over Nyquist.</summary>
///
__kernel void resample_postfft(complex_t *xa_d,
                               complex_t *xb_d,
                               int nfft,
                               int Ninterp,
                               int nyqst)
{
    size_t idx = get_global_id(0);

    // move frequencies from second half of spectrum; loop length: nfft - nyqst =
    // = nfft - nfft/2 - 1 = nfft/2 - 1
    if (idx < nfft / 2 - 1)
    {
        int i = nyqst + Ninterp - nfft + idx;
        int j = nyqst + idx;
        xa_d[i].x = xa_d[j].x;
        xa_d[i].y = xa_d[j].y;
        xb_d[i].x = xb_d[j].x;
        xb_d[i].y = xb_d[j].y;
    }

    // zero frequencies higher than nyquist, length: Ninterp - nfft
    // loop length: Ninterp - nfft ~ nfft
    if (idx < Ninterp - nfft)
    {
        xa_d[nyqst + idx].x = xa_d[nyqst + idx].y = 0.;
        xb_d[nyqst + idx].x = xb_d[nyqst + idx].y = 0.;
    }
}

/// <summary>Computes sin and cos values and stores them in an array.</summary>
/// <remarks>Most likely a very bad idea. Results are used in modvir and should be computed there in place.</remarks>
///
__kernel void compute_sincosmodf(__global real_t* s,
                                 __global real_t* c,
                                 real_t omr,
                                 int N)
{
    size_t idx = get_global_id(0);

    if (idx < N)
    {
        s[idx] = sin(omr * idx);
        c[idx] = cos(omr * idx);
    }
}

/// <summary>The purpose of this function was undocumented.</summary>
///
__kernel void computeB(__global complex_t* y,
                       __global complex_t* B,
                       int N)
{
    size_t idx = get_global_id(0);

    if (idx < N - 1)
    {
        B[idx].x = 6 * (y[idx + 2].x - 2 * y[idx + 1].x + y[idx].x);
        B[idx].y = 6 * (y[idx + 2].y - 2 * y[idx + 1].y + y[idx].y);
    }
}

//__global__ void phase_mod_1(cufftDoubleComplex *xa, cufftDoubleComplex *xb,
//                            cufftDoubleComplex *xar, cufftDoubleComplex *xbr,
//                            double het1, double sgnlt1, double *shft,
//                            int N){
//
//  int idx = blockIdx.x * blockDim.x + threadIdx.x;
//  if (idx < N) {
//    //          FLOAT_TYPE phase = - ( het1*idx + sgnlt1 * ( (double)idx*idx + 2 * idx * shft[idx] ) );
//    double phase = - idx * ( het1 + sgnlt1 * ( idx + 2 * shft[idx] ) );
//    double s, c;
//
//    sincos(phase, &s, &c);
//
//    xa[idx].x = xar[idx].x*c - xar[idx].y*s;
//    xa[idx].y = xar[idx].x*s + xar[idx].y*c;
//
//    xb[idx].x = xbr[idx].x*c - xbr[idx].y*s;
//    xb[idx].y = xbr[idx].x*s + xbr[idx].y*c;
//    //    if (idx==1) printf("xa=(%f, %f)\n", xa[idx].x, xa[idx].y);
//
//  }
//}
//
//__global__ void phase_mod_2(cufftDoubleComplex *xa, cufftDoubleComplex *xb,
//                            cufftDoubleComplex *xar, cufftDoubleComplex *xbr,
//                            double het1, double sgnlt1, double *shft,
//                            int N){
//
//  int idx = blockIdx.x * blockDim.x + threadIdx.x;
//  if (idx < N) {
//    //          FLOAT_TYPE phase = - ( het1*idx + sgnlt1 * ( (double)idx*idx + 2 * idx * shft[idx] ) );
//    double phase = - idx * ( het1 + sgnlt1 * ( idx + 2 * shft[idx] ) );
//    double s, c;
//
//    sincos(phase, &s, &c);
//
//    xa[idx].x += xar[idx].x*c - xar[idx].y*s;
//    xa[idx].y += xar[idx].x*s + xar[idx].y*c;
//
//    xb[idx].x += xbr[idx].x*c - xbr[idx].y*s;
//    xb[idx].y += xbr[idx].x*s + xbr[idx].y*c;
//    //    if (idx==1) printf("xa=(%f, %f)\n", xa[idx].x, xa[idx].y);
//
//  }
//}
//
//extern  __constant__  double maa_d, mbb_d;
//
//__global__ void compute_Fstat(cufftDoubleComplex *xa, cufftDoubleComplex *xb,
//                              double *F, int N) { // N = nmax - nmin
//  int i = blockIdx.x * blockDim.x + threadIdx.x;
//  //  if (i==0) printf("maa_d=%f\n", maa_d);
//  if (i < N) {
//    F[i] = ( xa[i].x*xa[i].x + xa[i].y*xa[i].y)/maa_d + (xb[i].x*xb[i].x + xb[i].y*xb[i].y)/mbb_d;
//  }
//}
//
//
///*
//__global__ void reduction_sum(double *in, double *out, int N) {
//  extern __shared__ double sd_data[];
//
//  int tid = threadIdx.x;
//  int i = blockIdx.x * blockDim.x + threadIdx.x;
//
//  sd_data[tid] = (i<N) ? in[i] : 0;
//
//  __syncthreads();
//
//  for (int s = blockDim.x/2; s>0; s>>=1) {
//    if (tid < s) {
//      sd_data[tid] += sd_data[tid + s];
//    }
//    __syncthreads();
//  }
//
//  if (tid==0) out[blockIdx.x] = sd_data[0];
//
//}
//*/
//
//__global__ void fstat_norm_simple(FLOAT_TYPE *F_d, int nav) {
//  //  int idx = blockIdx.x * blockDim.x + threadIdx.x;
//  int i;
//  FLOAT_TYPE *fr = F_d + blockIdx.x*nav;
//  FLOAT_TYPE mu = 0.;
//  for (i=0; i<nav; ++i)
//    mu += *fr++;
//  mu /= 2.*nav;
//  fr = F_d + blockIdx.x*nav;
//  for (i=0; i<nav; i++)
//    *fr++ /= mu;
//
//}
//
//
//// parameters are:
//// [frequency, spindown, position1, position2, snr]
//#define ADD_PARAMS_MACRO						\
//  int p = atomicAdd(found, 1);						\
//  params[p*NPAR + 0] = 2.0*M_PI*(idx)*fftpad*nfft+sgnl0;		\
//  params[p*NPAR + 1] = sgnl1;						\
//  params[p*NPAR + 2] = sgnl2;						\
//  params[p*NPAR + 3] = sgnl3;						\
//  params[p*NPAR + 4] = sqrt(2*(F[idx]-ndf));
//
//
//__global__ void find_candidates(FLOAT_TYPE *F, FLOAT_TYPE *params, int *found, FLOAT_TYPE val,
//                                int nmin, int nmax, double fftpad, double nfft, FLOAT_TYPE sgnl0, int ndf,
//                                FLOAT_TYPE sgnl1, FLOAT_TYPE sgnl2, FLOAT_TYPE sgnl3) {
//
//  int idx = blockIdx.x * blockDim.x + threadIdx.x + nmin;
//
//  if (idx > nmin && idx < nmax && F[idx] >= val && F[idx] > F[idx+1] && F[idx] > F[idx-1]) {
//    ADD_PARAMS_MACRO
//      } else if (idx == nmin && F[idx] >= val && F[idx] > F[idx+1]) {
//    ADD_PARAMS_MACRO
//      } else if (idx == nmax-1 && F[idx] >= val && F[idx] > F[idx-1]) {
//    ADD_PARAMS_MACRO
//      }
//}
//
//
//
////---------------------------------------------------------------
//
////second reduction used in fstat
//__global__ void reduction_sum(FLOAT_TYPE *in, FLOAT_TYPE *out, int N) {
//  extern __shared__ FLOAT_TYPE sf_data[];
//
//  int tid = threadIdx.x;
//  int i = blockIdx.x * blockDim.x + threadIdx.x;
//
//  sf_data[tid] = (i<N) ? in[i] : 0;
//
//  __syncthreads();
//
//  for (int s = blockDim.x/2; s>0; s>>=1) {
//    if (tid < s) {
//      sf_data[tid] += sf_data[tid + s];
//    }
//    __syncthreads();
//  }
//
//  if (tid==0) out[blockIdx.x] = 1.0f/sf_data[0];
//}
//
//
//__global__ void fstat_norm(FLOAT_TYPE *F, FLOAT_TYPE *mu, int N, int nav) {
//  int i = blockIdx.x * blockDim.x + threadIdx.x;
//  if (i < N) {
//    int block = i/nav; //block index
//    F[i] *= 2*nav * mu[block];
//  }
//}
