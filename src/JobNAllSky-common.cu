#include <math.h>
#include <gsl/gsl_vector.h>
#include <complex.h>
#include <cuda.h>
#include <cufft.h>
#include <assert.h>
#include <time.h>
#include "auxi.h"
#include "settings.h"
#include "cusparse_v2.h"
extern void* memcpy();
extern int nd;
extern double dt, *aa, *bb, *shftf, *shft, *t2, *F;
extern complex double *xDatma, *xDatmb;

void checkCUDAError(const char *msg) {

	cudaError_t err = cudaGetLastError();  
	if( cudaSuccess != err) {
    
		fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString(err));
	    exit(EXIT_FAILURE);

  	}

}


__global__ void crrhs(cufftComplex *funkcja,
					cufftComplex *rhs, 
					int N) {

	int i=threadIdx.x + blockIdx.x*blockDim.x;
  
	if(i<N) {

    	rhs[i].x=6.0*(funkcja[i].x-2.0*funkcja[i+1].x+funkcja[i+2].x);
	    rhs[i].y=6.0*(funkcja[i].y-2.0*funkcja[i+1].y+funkcja[i+2].y);
  	
	}

}


__global__ void crrhs2(cufftComplex *funkcja1,
					cufftComplex *rhs, 
					int N, 
					int Ninterp) {


	int i = threadIdx.x + blockIdx.x*blockDim.x;
	cufftComplex *funkcja = funkcja1 + Ninterp;

	if(i<N) {
      
		rhs[i].x=6.0*(funkcja[i].x-2.0*funkcja[i+1].x+funkcja[i+2].x);
      	rhs[i].y=6.0*(funkcja[i].y-2.0*funkcja[i+1].y+funkcja[i+2].y);
  
	}

}


__global__ void splinterp(cufftDoubleComplex *out, 
						cufftComplex *in, 
						cufftComplex *m, 
						double *shftf, 
						int interpftpad, 
						int Np) {

	int t=threadIdx.x+blockIdx.x*blockDim.x;

  	if(t<Np) {

	    float odleglosc=(interpftpad*((float)t-shftf[t])-floor(interpftpad*((float)t-shftf[t])));
	    int polozenie=(int)(interpftpad*((float)t-shftf[t]));
    	float c=((1-odleglosc)*(1-odleglosc)*(1-odleglosc)-(1-odleglosc))/6.0;
	    float d=(odleglosc*odleglosc*odleglosc-odleglosc)/6.0;
    
		if(polozenie==0) {

    		out[t].x=(1-odleglosc)*in[polozenie].x+odleglosc*in[polozenie+1].x+d*m[polozenie].x;
	      	out[t].y=(1-odleglosc)*in[polozenie].y+odleglosc*in[polozenie+1].y+d*m[polozenie].y;
	    }

	    if(polozenie<0 || polozenie>(interpftpad*Np-1)) {

    		out[t].x=0.0;
	      	out[t].y=0.0;

	    } else {

    		out[t].x=(1-odleglosc)*in[polozenie].x+odleglosc*in[polozenie+1].x+c*m[polozenie-1].x+d*m[polozenie].x;
	      	out[t].y=(1-odleglosc)*in[polozenie].y+odleglosc*in[polozenie+1].y+c*m[polozenie-1].y+d*m[polozenie].y;
    	}

  	}

}

__global__ void splinterp2(cufftDoubleComplex *out, 
						cufftComplex *in2, 
						cufftComplex *m, 
						double *shftf, 
						int interpftpad, 
						int Np, 
						int Ninterp) {

	int t=threadIdx.x+blockIdx.x*blockDim.x;

  	cufftComplex *in=in2+Ninterp;

	if(t<Np) {
	    float odleglosc=(interpftpad*((float)t-shftf[t])-floor(interpftpad*((float)t-shftf[t])));
    	int polozenie=(int)(interpftpad*((float)t-shftf[t]));
	    float c=((1-odleglosc)*(1-odleglosc)*(1-odleglosc)-(1-odleglosc))/6.0;
    	float d=(odleglosc*odleglosc*odleglosc-odleglosc)/6.0;

	    if(polozenie==0) {

    		out[t].x=(1-odleglosc)*in[polozenie].x+odleglosc*in[polozenie+1].x+d*m[polozenie].x;
	        out[t].y=(1-odleglosc)*in[polozenie].y+odleglosc*in[polozenie+1].y+d*m[polozenie].y;

    	}
    
		if( polozenie<0 || polozenie>(interpftpad*Np-1)) {

      		out[t].x=0.0;
		    out[t].y=0.0;

	    } else {

      		out[t].x=(1-odleglosc)*in[polozenie].x+odleglosc*in[polozenie+1].x+c*m[polozenie-1].x+d*m[polozenie].x;
		    out[t].y=(1-odleglosc)*in[polozenie].y+odleglosc*in[polozenie+1].y+c*m[polozenie-1].y+d*m[polozenie].y;
    
		}
  
	}

}


__global__ void smallint(cufftComplex *xa, 
						int nfft, 
						int Npoints) {

	cufftComplex *xb=xa+nfft;
  	int i=threadIdx.x+blockIdx.x*blockDim.x;

    if(i<nfft-Npoints) {
    
		xa[i+Npoints].x =0.0;
		xb[i+Npoints].x =0.0;
	    xa[i+Npoints].y =0.0;
    	xb[i+Npoints].y =0.0;
  
	}

}


__global__ void zeropadd(cufftComplex *xa, 
						int Npoints, 
						int fftpad, 
						int nfft) {

	int i=threadIdx.x+blockIdx.x*blockDim.x;

	if(i<fftpad*nfft-Npoints) {

    	xa[i+Npoints].x=0.0;
	    xa[i+Npoints].y=0.0;
    	xa[i+Npoints+fftpad*nfft].x=0.0;
	    xa[i+Npoints+fftpad*nfft].y=0.0;
  
	}

}


__global__ void 
fillf(double *F, cufftComplex *xao, int nmin, int nmax, int nfft, int fftinterp, int fftpad) {
  int i=threadIdx.x+blockIdx.x*blockDim.x;
  if(i<(nmax-nmin)) {
    if(fftinterp==1) {
      F[i] = (xao[i+nmin].x)*(xao[i+nmin].x) + (xao[i+nmin].y)*(xao[i+nmin].y) +(xao[i+nfft+nmin].x)*(xao[i+nfft+nmin].x) + (xao[i+nfft+nmin].y)*(xao[i+nfft+nmin].y);
    }
    if(fftinterp==2) {
      F[i] = (xao[i+nmin].x)*(xao[i+nmin].x) + (xao[i+nmin].y)*(xao[i+nmin].y) +(xao[i+nfft*fftpad+nmin].x)*(xao[i+nfft*fftpad+nmin].x) + (xao[i+nfft*fftpad+nmin].y)*(xao[i+nfft*fftpad+nmin].y);
    }
  }
}

__global__ void 
fdivsig2(double *F, int nmin, int nmax, double sig2) {
  int i=threadIdx.x+blockIdx.x*blockDim.x;
  if(i<(nmax-nmin)) {
    F[i] = F[i]/sig2;
  }
}

__global__ void 
findsig(double *F, int nmin, int nmax, double sgn1, double sgn2, double sgn3, double trl, int fftpad, int nfft, double sgnl0, int nd, int *trigcount, float *hugemem) {
  int i = threadIdx.x+blockIdx.x*blockDim.x;
  if(i==0) {
    if(F[i]>trl && F[i]>F[i+1]) {
      double frq=2.*M_PI*(i+nmin)/((double)fftpad*nfft)+sgnl0;
      double snr=sqrt (2.*(F[i]-nd));
      int val=5*atomicAdd(trigcount, 1);
      hugemem[val]=frq;
      hugemem[val+1]=sgn1;
      hugemem[val+2]=sgn2;
      hugemem[val+3]=sgn3;
      hugemem[val+4]=snr;
    }
  }
  if(i==(nmax-nmin-1)) {
    if(F[i]>trl && F[i]>F[i-1])	{
      double frq=2.*M_PI*(i+nmin)/((double)fftpad*nfft)+sgnl0;
      double snr=sqrt (2.*(F[i]-nd));
      int val=5*atomicAdd(trigcount, 1);
      hugemem[val]=frq;
      hugemem[val+1]=sgn1;
      hugemem[val+2]=sgn2;
      hugemem[val+3]=sgn3;
      hugemem[val+4]=snr;
    }
  }
  if(i>0 && i<(nmax-nmin-1))  {
    if(F[i]>trl && F[i]>F[i-1] && F[i]>F[i+1]) {
      double frq=2.*M_PI*(i+nmin)/((double)fftpad*nfft)+sgnl0;
      double snr=sqrt (2.*(F[i]-nd));
      int val=5*atomicAdd(trigcount, 1);
      hugemem[val]=frq;
      hugemem[val+1]=sgn1;
      hugemem[val+2]=sgn2;
      hugemem[val+3]=sgn3;
      hugemem[val+4]=snr;
    }
  }
}


__global__ void 
smallinterp1(cufftComplex *xaoout, cufftComplex *xao, int nfft) {
  cufftComplex *xbo=xao+nfft;
  cufftComplex *xboout=xaoout+nfft;
  int i=threadIdx.x+blockIdx.x*blockDim.x;

  /* Interpolation algorithm */
  if (i<nfft/2) {
    xaoout[2*i] = xao[i];
    xboout[2*i] = xbo[i];
  }
}

__global__ void 
smallinterp2(cufftComplex *xaoout, cufftComplex *xao, int nfft) {
  cufftComplex *xbo=xao+nfft;
  cufftComplex *xboout=xaoout+nfft;
  int i=threadIdx.x+blockIdx.x*blockDim.x;
  
  if (i<nfft-1 && i%2==1) {
    xaoout[i].x = (xao[i-1].x-xao[i+1].x)/M_SQRT2;
    xboout[i].x = (xbo[i-1].x-xbo[i+1].x)/M_SQRT2;
    xaoout[i].y = (xao[i-1].y-xao[i+1].y)/M_SQRT2;
    xboout[i].y = (xbo[i-1].y-xbo[i+1].y)/M_SQRT2;
  }
  //  __syncthreads();
}



__global__ void phas2(cufftComplex *xb, 
					int Np, 
					cufftComplex *xa, 
					double sgnlt, 
					double het1, 
					double *shft, 
					cufftDoubleComplex *xDatma, 
					cufftDoubleComplex *xDatmb) {

#if 0
  cufftComplex *xb;

  switch (fftinterp) {
  case 1:
    xb = xa+nfft;
    break;
  case 2:
    xb = xa+fftpad*nfft;
  } 

  xb = xa + xb_shift;
#endif

	int i = threadIdx.x + blockIdx.x*blockDim.x;

	if(i<Np) {
    //double phase2 = het1*(double)i + sgnlt*((double)(i*i)+2.*(double)i*shft[i]);
    //double  cosPH = cos (phase2);
    //double sinPH = sin (phase2);
    //exph = cosPH - I*sinPH;
#if 1
    
		double phase, sp, cp;

	    phase = (het1 + sgnlt*((double)i + 2.*shft[i]))*(double)i;

    	sincos(phase, &sp, &cp);

	    xa[i].x = (float)(xDatma[i].x*cp + xDatma[i].y*sp);
	    xa[i].y = (float)(xDatma[i].y*cp - xDatma[i].x*sp);

	    xb[i].x = (float)(xDatmb[i].x*cp + xDatmb[i].y*sp);
    	xb[i].y = (float)(xDatmb[i].y*cp - xDatmb[i].x*sp);

#else
    xa[i].x =(float)( xDatma[i].x*cos(het1*(double)i + sgnlt*((double)i*(double)i+2.*(double)i*shft[i]))
                    + xDatma[i].y*sin(het1*(double)i + sgnlt*((double)i*(double)i+2.*(double)i*shft[i])));
    xb[i].x = (float)(xDatmb[i].x*cos(het1*(double)i + sgnlt*((double)i*(double)i+2.*(double)i*shft[i]))
		    + xDatmb[i].y*sin(het1*(double)i + sgnlt*((double)i*(double)i+2.*(double)i*shft[i])));
    xa[i].y =(float)( xDatma[i].y*cos(het1*(double)i + sgnlt*((double)i*(double)i+2.*(double)i*shft[i]))
		    - xDatma[i].x*sin(het1*(double)i + sgnlt*((double)i*(double)i+2.*(double)i*shft[i])));
    xb[i].y =(float)(xDatmb[i].y*cos(het1*(double)i + sgnlt*((double)i*(double)i+2.*(double)i*shft[i]))
		    - xDatmb[i].x*sin(het1*(double)i + sgnlt*((double)i*(double)i+2.*(double)i*shft[i])));
#endif
  
	}

}

__global__ void 
resample(cufftComplex *xDa, cufftDoubleComplex *cuxDatma, cufftDoubleComplex *cuxDatmb, int Npoints, int nfft) {
  cufftComplex *xDb=xDa+nfft;
  int i=threadIdx.x+blockIdx.x*blockDim.x;
  if (i<Npoints) {
    xDa[i].x = (float)cuxDatma[i].x;
    xDb[i].x = (float)cuxDatmb[i].x;
    xDa[i].y = (float)cuxDatma[i].y;
    xDb[i].y = (float)cuxDatmb[i].y;
  }
  if (i>=Npoints && i<nfft) {
    xDa[i].x = 0.;
    xDb[i].x = 0.;
    xDa[i].y = 0.;
    xDb[i].y = 0.;
  } else {
  }
}

__global__ void dord1(cufftComplex *xDa, 
					cufftComplex *rDa,
					int Npoints, 
					int nfft, 
					int Ninterp) {

	cufftComplex *xDb=xDa+nfft;
	cufftComplex *rDb=rDa+Ninterp;

	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int nyqst = (nfft+2)/2;

	if(i<nyqst) 
    	rDb[i] = xDb[i];
  
	if (i>=nyqst && i<nyqst + Ninterp - nfft) {

    	rDb[i].x = 0.;
	    rDb[i].y = 0.;
	
	}

	__syncthreads();

}

__global__ void dord2(cufftComplex *xDa, 
					cufftComplex *rDa, 
					int Npoints, 
					int nfft, 
					int Ninterp) {

	int i = threadIdx.x + blockIdx.x*blockDim.x;

	cufftComplex *xDb=xDa+nfft;
	cufftComplex *rDb=rDa+Ninterp;

	int nyqst = (nfft+2)/2;

	if(i>=nyqst+Ninterp-nfft && i<Ninterp) {

	    rDb[i] = xDb[i - Ninterp + nfft];
    	rDa[i] = xDa[i - Ninterp + nfft];
	}

	if(i>=nyqst && i<nyqst + Ninterp - nfft) {

    	rDa[i].x = 0.;
	    rDa[i].y = 0.;

	}

}


__global__ void rdft(cufftComplex *rDa, 
					int interpftpad, 
					int Ninterp) {

	float ft = (float)interpftpad/Ninterp;

  	int i=threadIdx.x+blockIdx.x*blockDim.x;

	cufftComplex *rDb=rDa+Ninterp;

    if(i<Ninterp) {

    	rDa[i].x *=  ft;
	    rDb[i].x *= ft;
    	rDa[i].y *=  ft;
	    rDb[i].y *= ft;
  
	}

}

__global__ void phasfilt(int Npoints, 
						double *cunSource, 
						double het0, 
						double oms, 
						double *cuxDat, 
						double *cuDetSSB, 
						double *aa, 
						double *bb, 
						double *shft, 
						double* shftf, 
						cufftDoubleComplex *cuxDatma, 
						cufftDoubleComplex *cuxDatmb) {

	int i = threadIdx.x + blockIdx.x*blockDim.x;

	if(i<Npoints) {
    
		double shft1 = 0.;
	    int j;

    	shft[i] = 0.;

	    for (j=0; j<3; j++) {

    		shft1 += cunSource[j]*cuDetSSB[j];
	        shft[i] += cunSource[j]*cuDetSSB[i*3+j];
	    }

    	shftf[i] = shft[i] - shft1;

	    cuxDatma[i].x = cuxDat[i]*aa[i]*cos (het0*(double)i+oms*shft[i]);
    	cuxDatmb[i].x = cuxDat[i]*bb[i]*cos (het0*(double)i+oms*shft[i]);
	    cuxDatma[i].y = -cuxDat[i]*aa[i]*sin (het0*(double)i+oms*shft[i]);
    	cuxDatmb[i].y =-cuxDat[i]*bb[i]*sin (het0*(double)i+oms*shft[i]);

	}

  	syncthreads();

}

// lin2ast described in Phys. Rev. D 82, 022005 (2010) (arXiv:1003.0844)
void
lin2ast (double be1, double be2, int pm, double sepsm, double cepsm,	\
	 double *sinal, double *cosal, double *sindel, double *cosdel) {
  *sindel = be1*sepsm-(2*pm-3)*sqrt(1.-sqr(be1)-sqr(be2))*cepsm;
  *cosdel = sqrt(1.-sqr(*sindel));
  *sinal = (be1-sepsm*(*sindel))/(cepsm*(*cosdel));
  *cosal = be2/(*cosdel);
} /* lin2ast() */


__device__ double atomicAddd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                             (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}


template <unsigned int blockSize>
__device__ void warpReduce(volatile double *sdata, unsigned int tid) {

	if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
	if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
	if (blockSize >= 16) sdata[tid] += sdata[tid + 8];
	if (blockSize >= 8) sdata[tid] += sdata[tid + 4];
	if (blockSize >= 4) sdata[tid] += sdata[tid + 2];
	if (blockSize >= 2) sdata[tid] += sdata[tid + 1];

}


__global__ void
modvir2 (double sinal, double cosal, double sindel, double cosdel,	\
	 double sphir, double cphir, double *a, double *b, int Np, double *as, double *bs) {
  /* Amplitude modulation functions */
  //extern double *cosmodf, *sinmodf;
  int t;
  double cosalfr, sinalfr, c2d, c2sd;

  __shared__ double tempa[256];
  __shared__ double tempb[256];


  double omr=0.5*7.2921151467064e-5;
  double deg=180.0/M_PI;
  double ephi = (43.+37./60.+53.0880/3600.)/deg; 
  double egam = (135. - (19.0+25./60.0+57.96/3600.))/deg;
  double  c1 = .25*sin(2.*egam)*(1+sqr(sin(ephi)));
  double  c2 = -.5*cos(2.*egam)*sin(ephi);
  double  c3 = .5*sin(2.*egam)*sin(2.*ephi);
  double  c4 = -cos(2.*egam)*cos(ephi);
  double  c5 = .75*sin(2.*egam)*sqr(cos(ephi));
  double  c6 = cos(2.*egam)*sin(ephi);
  double  c7 = .5*sin(2.*egam)*(1.+sqr(sin(ephi)));
  double  c8 = cos(2.*egam)*cos(ephi);
  double  c9 = .5*sin(2.*egam)*sin(2.*ephi);

  //pci wtf? chyba to można skasować
  /*  if(0==threadIdx.x) {
    tempa[threadIdx.x]=0.0;
    tempb[threadIdx.x]=0.0;
  }
  */

  if(0==threadIdx.x && 0==blockIdx.x) {
    *as=0.0;
    *bs=0.0;
  }
  // pewnie niepotrzebne
  //__syncthreads();
  cosalfr = cosal*cphir+sinal*sphir;
  sinalfr = sinal*cphir-cosal*sphir;
  c2d = sqr(cosdel);
  c2sd = sindel*cosdel;

  t=threadIdx.x+blockIdx.x*blockDim.x;

  if(t<Np) {

    double sot, cot;
    sincos(omr*t, &sot, &cot);

//     a[t]=c1*(2.-c2d)*2.*sqr(cosalfr*cos(omr*t)+sinalfr*sin(omr*t))+c2*(2.-c2d)*2.*(cosalfr*cos(omr*t)+sinalfr*sin(omr*t))*(sinalfr*cos(omr*t)-cosalfr*sin(omr*t))+c3*c2sd*(cosalfr*cos(omr*t)+sinalfr*sin(omr*t))+c4*c2sd*(sinalfr*cos(omr*t)-cosalfr*sin(omr*t))- c1*(2.-c2d)+c5*c2d;
     a[t]=c1*(2.-c2d)*2.*sqr(cosalfr*cot+sinalfr*sot)+c2*(2.-c2d)*2.*(cosalfr*cot+sinalfr*sot)*(sinalfr*cot-cosalfr*sot)+c3*c2sd*(cosalfr*cot+sinalfr*sot)+c4*c2sd*(sinalfr*cot-cosalfr*sot)- c1*(2.-c2d)+c5*c2d;
    
//#mb
//    tempa[threadIdx.x]=(double)a[t]*a[t];

//     b[t] = c6*sindel*2.*sqr(cosalfr*cos(omr*t)+sinalfr*sin(omr*t))+c7*sindel*2.*(cosalfr*cos(omr*t)+sinalfr*sin(omr*t))*(sinalfr*cos(omr*t)-cosalfr*sin(omr*t))+c8*cosdel*(cosalfr*cos(omr*t)+sinalfr*sin(omr*t))+c9*cosdel*(sinalfr*cos(omr*t)-cosalfr*sin(omr*t))- c6*sindel;
      b[t] = c6*sindel*2.*sqr(cosalfr*cot+sinalfr*sot)+c7*sindel*2.*(cosalfr*cot+sinalfr*sot)*(sinalfr*cot-cosalfr*sot)+c8*cosdel*(cosalfr*cot+sinalfr*sot)+c9*cosdel*(sinalfr*cot-cosalfr*sot)- c6*sindel;
    
//#mb
//    tempb[threadIdx.x]=(double)b[t]*b[t];

  }
 
  __syncthreads();

	unsigned int blockSize = 256; 

	unsigned int tid = threadIdx.x;
	unsigned int i = blockIdx.x*(blockSize*2) + tid;
	unsigned int gridSize = blockSize*2*gridDim.x;

	tempa[tid] = 0;
	tempb[tid] = 0;

	while (i < Np) {
 
        tempa[tid] += a[i]*a[i] + a[i+blockSize]*a[i+blockSize];
		tempb[tid] += b[i]*b[i] + b[i+blockSize]*b[i+blockSize];

		i += gridSize; 
		
	}

	__syncthreads();

/*	if (blockSize >= 512) { 
		if (tid < 256) { 
			tempa[tid] += tempa[tid + 256]; 
			tempb[tid] += tempb[tid + 256];

		} 
		__syncthreads(); 

	}
*/

    if (blockSize >= 256) { 
        if (tid < 128) { 
            tempa[tid] += tempa[tid + 128]; 
            tempb[tid] += tempb[tid + 128];

        } 
        __syncthreads(); 

    }


    if (blockSize >= 128) {
        if (tid < 64) {
            tempa[tid] += tempa[tid + 64];
            tempb[tid] += tempb[tid + 64];

        }
        __syncthreads();

    }

	if (tid < 32) { 

    if (blockSize >= 64) tempa[tid] += tempa[tid + 32];
    if (blockSize >= 32) tempa[tid] += tempa[tid + 16];
    if (blockSize >= 16) tempa[tid] += tempa[tid + 8];
    if (blockSize >= 8) tempa[tid] += tempa[tid + 4];
    if (blockSize >= 4) tempa[tid] += tempa[tid + 2];
    if (blockSize >= 2) tempa[tid] += tempa[tid + 1];

    if (blockSize >= 64) tempb[tid] += tempb[tid + 32];
    if (blockSize >= 32) tempb[tid] += tempb[tid + 16];
    if (blockSize >= 16) tempb[tid] += tempb[tid + 8];
    if (blockSize >= 8) tempb[tid] += tempb[tid + 4];
    if (blockSize >= 4) tempb[tid] += tempb[tid + 2];
    if (blockSize >= 2) tempb[tid] += tempb[tid + 1];

//		warpReduce(tempa, tid);
//		warpReduce(tempb, tid);

 	}

	if (tid == 0) { 

			atomicAddd(as,tempa[0]);
		    atomicAddd(bs,tempb[0]);

	}

/*  if(0==threadIdx.x) {
    double sumas=0.0;
    double sumbs=0.0;
    for(int i=0;i<256;i++) {
      if((i+blockIdx.x*blockDim.x)<Np){
	  sumas+=tempa[i];
	  sumbs+=tempb[i];
      }
    }
    atomicAddd(as,sumas);
    atomicAddd(bs,sumbs);
  }
  //__syncthreads();
*/ 

}

__global__ void 
modvirctd( double *a, double *b, int Np, double *as, double *bs) {

  if(0==threadIdx.x && 0==blockIdx.x) {
    *as =*as/(float)Np;
    *bs =*bs/(float)Np;
    *as = sqrt (*as);
    *bs = sqrt (*bs);
  }
  __syncthreads();
}

__global__ void modvirctd2(double *a, 
						double *b, 
						int Np, 
						double *as, 
						double *bs) {

	int t = threadIdx.x + blockIdx.x*blockDim.x;

	if(t<Np) {

	  a[t] /= *as;
	  b[t] /= *bs;

	  //	  if(t==11233)  printf("kernel at[%d] = %f\n", t, a[t]);
      	}
 
} /* modvir() */

double *
JobCore(int pm,			// hemisphere 
	int mm,			// grid 'sky position'
	int nn,			// other grid 'sky position' 
	int smin,		// grid spindown coordinate
	int smax,		// spindown range limit 
	double *M,		// grid generating matrix
	double *cuDetSSB,		// ephemerides array 
	double *cuxDat,		// time-domain input data array 
	int Npoints,		// Number of data points
	int Ninterp,		// interpftpad*nfft (for resampling)
	int nfft,		// size of the FFT transform
	//complex float *xDa,	// Array for resampling
	//complex float *rDa,	// Array for resampling
	complex float *xa,	// input array for plan
	cufftComplex *cuxao,	// output array for plan
	cufftHandle pl_int,	// cufftHandle needed for resampling
	cufftHandle pl_inv,	// cufftHandle needed for resampling
				// (inverse transformation)
	cufftHandle plan,		// main cufftHandle 
	int nmin,		// nmin+1: first point of 
				// the F-statistic to search for a signal
	int nmax,		// nmax: last point of 
				// the F-statistic to search for a signal
	double sepsm,		// sin(epsm)
	double cepsm,		// cos(epsm)
	double sphir,		// sin(phi_r)
	double cphir,		// cos(phi_r)
	int *sgnlc,		// reference to array with the parameters 
				// of the candidate signal 
				// (used below to write to the file) 
	int write_st,		// std output writing flag
	int fftinterp,		// interpolation flag 
				// (INT of FFT, see settings.h) 
	int *FNum,		// Candidate signal number
	double coft,		// = oms
	double trl,		// F-statistic threshold
	double sig2,		// N*(variance of xDat) if white_flag
				// else sig2=-1 
	int s0,			// No-spindown flag
	double *cudaF, 
	cufftComplex *nakarcie1, 
	cufftComplex *curDa, 
	cufftComplex *du, 
	cufftComplex *d, 
	cufftComplex *dl, 
	cufftComplex *splvec, 
	double *cuaa, 
	double *cubb, 
	int *trigs, 
	float *humem) {

  //int j, i;
  double al1, al2, sinalt, cosalt, sindelt, cosdelt, sgnlt[NPAR], \
    nSource[3], het0, sgnl0, *sgnlv;
  double *as, *bs; //sums for modvir

  cudaMalloc((void**)&as, sizeof(double));
  cudaMalloc((void**)&bs, sizeof(double));

  //checkCUDAError("391");
  /* Matrix  M(.,.) (defined on page 22 of PolGrawCWAllSkyReview1.pdf file)
     defines the transformation form integers (bin, ss, nn, mm) determining 
     a grid point to linear coordinates omega, omegadot, alpha_1, alpha_2),
     where bin is the frequency bin number and alpha_1 and alpha_2 are 
     defined on p. 22 of PolGrawCWAllSkyReview1.pdf file. 

     [omega]                           [bin]
     [omegadot]        = M(.,.) \times [ss]
     [alpha_1/omega]                   [nn]
     [alpha_2/omega]                   [mm]

     Array M[.] is related to matrix M(.,.) in the following way;

     [ M[0] M[4] M[8]  M[12] ]
     M(.,.) =  [ M[1] M[5] M[9]  M[13] ]
     [ M[2] M[6] M[10] M[14] ]
     [ M[3] M[7] M[11] M[15] ]

     and 

     M[1] = M[2] = M[3] = M[6] = M[7] = 0

  */ 
 
  al1 = nn*M[10]+mm*M[14];
  al2 = nn*M[11]+mm*M[15];
  sgnlv = NULL;
  *sgnlc = 0;

  // check if the search is in an appropriate region of the grid
  // if not, returns NULL 
  if((sqr(al1)+sqr(al2))/sqr(coft) > 1.) return NULL ; 

  int ss;
  //complex float *xb, *xbo;
    
   
  //xb = xbo = NULL;
  // Switch for different interpolation methods 
  // is has been moved to all the kernels it is needed in  


	lin2ast (al1/coft, al2/coft, pm, sepsm, cepsm, &sinalt, &cosalt, &sindelt, &cosdelt);

    // calculate declination and right ascention
    // written in file as candidate signal sky positions 
	sgnlt[2] = asin (sindelt);
	sgnlt[3] = fmod (atan2 (sinalt, cosalt)+2.*M_PI, 2.*M_PI);

	// sanity check 1  
	printf("%lf %lf\n", sgnlt[2], sgnlt[3]); 

  /* amplitude modulation functions */
 
  modvir2<<<(int)((double)(Npoints)/256.0)+1, 256>>> (sinalt, cosalt, sindelt, cosdelt,sphir, cphir, cuaa, cubb, Npoints, as, bs);
  modvirctd<<<(int)((double)(Npoints)/256.0)+1, 256>>>(cuaa, cubb, Npoints, as, bs);
  modvirctd2<<<(int)((double)(Npoints)/256.0)+1, 256>>>(cuaa, cubb, Npoints, as, bs);

  nSource[0] = cosalt*cosdelt;
  nSource[1] = sinalt*cosdelt;
  nSource[2] = sindelt;

  het0 = fmod (nn*M[8]+mm*M[12], M[0]);

  cufftDoubleComplex *cuxDatma;
  cudaMalloc((void**)&cuxDatma, Npoints*sizeof(cufftDoubleComplex));
  cufftDoubleComplex *cuxDatmb;
  cudaMalloc((void**)&cuxDatmb, Npoints*sizeof(cufftDoubleComplex));

  double *cushft, *cushftf;
  cudaMalloc((void**)&cushft, Npoints*sizeof(double));
  cudaMalloc((void**)&cushftf, Npoints*sizeof(double));

  double *cunSource;
  cudaMalloc((void**)&cunSource, 3*sizeof(double));

  cudaMemcpy(cunSource,nSource, 3*sizeof(double), cudaMemcpyHostToDevice);

  phasfilt<<<(int)((double)(4*Npoints)/256.0)+1, 256>>>(Npoints, cunSource,het0, oms, cuxDat, cuDetSSB, cuaa, cubb, cushft, cushftf, cuxDatma, cuxDatmb);
  resample<<<(int)((double)(4*fftpad*nfft)/256.0)+1, 256>>>(nakarcie1,cuxDatma,cuxDatmb,Npoints,nfft);
  checkCUDAError("545");

  cudaFree(cunSource);
  cudaFree(as);
  cudaFree(bs);

  // Resampling to barycentric time 
  //-------------------------------
	
  cufftExecC2C(pl_int, nakarcie1, nakarcie1, CUFFT_FORWARD);

  cudaMemcpy(curDa,nakarcie1,Ninterp*sizeof(cufftComplex),cudaMemcpyDeviceToDevice);
  dord1<<< (int)((double)(4*fftpad*nfft)/256.0)+1, 256>>>(nakarcie1,curDa,Npoints,nfft,Ninterp);
  dord2<<< (int)((double)(4*fftpad*nfft)/256.0)+1, 256>>>(nakarcie1,curDa,Npoints,nfft,Ninterp);
  cufftExecC2C(pl_inv,curDa,curDa,CUFFT_INVERSE);

  rdft<<<(int)((double)(4*fftpad*nfft)/256.0)+1, 256>>>(curDa,interpftpad,Ninterp);
  checkCUDAError("561");
  //cubic spline interpolation
  crrhs<<<(int)((float)(Ninterp-2)/256.0)+1,256>>>(curDa,splvec,Ninterp-2);
  cusparseHandle_t handle=0;
  cusparseCreate(&handle);
  checkCUDAError("567");
  if(cusparseCgtsv(handle, Ninterp-2,1,dl, d, du, splvec,Ninterp-2)!=CUSPARSE_STATUS_SUCCESS) {
    printf("\n Dupa :( \n");
    abort();
  }
  splinterp<<<(int)((float)(2*Ninterp)/256.0)+1,256>>>(cuxDatma,curDa,splvec,cushftf,interpftpad, Npoints);
  checkCUDAError("614");
  cudaDeviceSynchronize();
  crrhs2<<<(int)((float)(Ninterp-2)/256.0)+1,256>>>(curDa,splvec,Ninterp-2, Ninterp);
  checkCUDAError("616");
  cudaDeviceSynchronize();
  if(cusparseCgtsv(handle, Ninterp-2,1,dl, d, du, splvec,Ninterp-2)!=CUSPARSE_STATUS_SUCCESS) {
    printf("\n Dupa1 :( \n");
    abort();
  }
  checkCUDAError("621");
  splinterp2<<<(int)((float)(2*Ninterp)/256.0)+1,256>>>(cuxDatmb,curDa,splvec,cushftf,interpftpad,Npoints, Ninterp);
  checkCUDAError("623");
  cudaFree(cushftf);
  //End of cubic interpolation
  cudaDeviceSynchronize();


  if (write_st)
    printf ("\n>>%d\t%d\t%d\t[%d..%d]\n", *FNum, mm, nn, smin, smax);

  // if no-spindown
  if (s0) smin = smax;

  // if spindown parameter is taken into account, 
  // smin != smax 

  int xb_shift = nfft;
  if (fftinterp==2) xb_shift = fftpad*nfft;
  cufftComplex *nakarcie_b;
  nakarcie_b = nakarcie1 + xb_shift;

  for (ss=smin; ss<=smax; ss++) {

    // Spindown parameter
    sgnlt[1] = (s0 ? 0. : ss*M[5] + nn*M[9] + mm*M[13]);

    if (sgnlt[1] >= -Smax && sgnlt[1] <= 0.) {
	
      double het1;
      if (write_st) {
	printf (".");
	fflush (stdout);
      }

      het1 = fmod(ss*M[4], M[0]);
      if (het1<0) het1 += M[0];

      sgnl0 = het0+het1;

      phas2<<<(int)((float)Npoints/256.0)+1,256>>>(nakarcie_b, Npoints, nakarcie1, sgnlt[1],het1, cushft, cuxDatma, cuxDatmb);
      checkCUDAError("537");

      /* for i */
      switch (fftinterp) {
      case INT:
	smallint<<<(int)((float)Npoints/256.0)+1,256>>>(nakarcie1,nfft,Npoints);
	cufftExecC2C(plan,nakarcie1,nakarcie1,CUFFT_FORWARD);
	smallinterp1<<<(int)((float)2*nfft/256.0)+1,256>>>(cuxao,nakarcie1,nfft);
	smallinterp2<<<(int)((float)2*nfft/256.0)+1,256>>>(cuxao,cuxao,nfft);
	break;
      case FFT:
	zeropadd<<<(int)((float)(nfft*fftpad-Npoints)/256.0)+1,256>>>(nakarcie1,Npoints,fftpad,nfft);
	cufftExecC2C(plan,nakarcie1,cuxao,CUFFT_FORWARD);
      }

      (*FNum)++;

      fillf<<<(int)((float)(nmax-nmin)/256.0)+1,256>>>(cudaF, cuxao,nmin,nmax,nfft,fftinterp, fftpad);
      //cudaMemcpy(F,cudaF,(nmax-nmin)*sizeof(double),cudaMemcpyDeviceToHost);
      /* Normalize F-statistics */
      // if the noise is not white noise
      //	  if (sig2 < 0.) FStat (F+nmin, nmax-nmin, NAV, 0);
      //	  else for (i=nmin; i<nmax; i++) F[i] /= sig2;
      fdivsig2<<<(int)((float)(nmax-nmin)/256.0)+1,256>>>(cudaF,nmin,nmax,sig2);
      //cudaMemcpy(F,cudaF,(nmax-nmin)*sizeof(double),cudaMemcpyDeviceToHost);


      findsig<<<(int)((float)(nmax-nmin)/256.0)+1,256>>>(cudaF, nmin, nmax, sgnlt[1],sgnlt[2],sgnlt[3], trl, fftpad, nfft, sgnl0, nd, trigs, humem);

    } /* if sgnlt[1] */
  } /* for ss */


  cudaFree(cuxDatma);
  cudaFree(cuxDatmb);
  cudaFree(cushft);

  checkCUDAError("Loop finished");
  return sgnlv;
} /* JobCore() */
