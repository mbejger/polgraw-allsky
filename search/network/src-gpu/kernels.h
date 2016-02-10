#ifndef __KERNELS_H__
#define __KERNELS_H__

void copy_amod_coeff(int nifo);

__global__ void modvir_kern(double *aa_d, double *bb_d, double cosalfr, double sinalfr,
			    double c2d, double c2sd, 
			    double *sinmodf_d, double *cosmodf_d, 
			    double sindel, double cosdel, int Np, int idet);

__global__ void tshift_pmod_kern(double shft1, double het0, 
				 double ns0, double ns1, double ns2,  
				 double *xDat_d, 
				 cufftDoubleComplex *xa_d, cufftDoubleComplex *xb_d, 
				 double *shft_d, double *shftf_d, 
				 double *tshift_d, 
				 double *aa_d, double *bb_d,
				 double *DetSSB_d, 
				 double oms, int N, int nfft, int interpftpad);

__global__ void resample_postfft(cufftDoubleComplex *xa_d, cufftDoubleComplex *xb_d,
                                 int nfft, int Ninterp, int nyqst);

__global__ void compute_sincosmodf(double *s, double *c, double omr, int N);

__global__ void phase_mod_1(cufftDoubleComplex *xa, cufftDoubleComplex *xb,
                            cufftDoubleComplex *xar, cufftDoubleComplex *xbr,
                            double het1, double sgnlt1, double *shft,
                            int N);

__global__ void phase_mod_2(cufftDoubleComplex *xa, cufftDoubleComplex *xb,
                            cufftDoubleComplex *xar, cufftDoubleComplex *xbr,
                            double het1, double sgnlt1, double *shft,
                            int N);

__global__ void compute_Fstat(cufftDoubleComplex *xa, cufftDoubleComplex *xb,
                              double *F, 
			      int N);


#endif
