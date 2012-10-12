#include <math.h>
#include <gsl/gsl_vector.h>
#include <complex.h>
#include <fftw3.h>
#include "auxi.h"
#include "lvcvirgo.h"

extern int nd;
extern double dt, *aa, *bb, *shftf, *shft, *t2, *F;
extern complex double *xDatma, *xDatmb;

// lin2ast described in Phys. Rev. D 82, 022005 (2010) (arXiv:1003.0844)
void
lin2ast (double be1, double be2, int pm, double sepsm, double cepsm,	\
	 double *sinal, double *cosal, double *sindel, double *cosdel) {
  *sindel = be1*sepsm-(2*pm-3)*sqrt(1.-sqr(be1)-sqr(be2))*cepsm;
  *cosdel = sqrt(1.-sqr(*sindel));
  *sinal = (be1-sepsm*(*sindel))/(cepsm*(*cosdel));
  *cosal = be2/(*cosdel);
} /* lin2ast() */

void
modvir (double sinal, double cosal, double sindel, double cosdel,	\
	double sphir, double cphir, double *a, double *b, int Np) {
  /* Amplitude modulation functions */
  extern double *cosmodf, *sinmodf;
  int t;
  double cosalfr, sinalfr, c2d, c2sd, c, s, c2s, cs, as, bs;

  cosalfr = cosal*cphir+sinal*sphir;
  sinalfr = sinal*cphir-cosal*sphir;
  c2d = sqr(cosdel);
  c2sd = sindel*cosdel;

  as = bs = 0.;
  for (t=0; t<Np; t++) {
    c = cosalfr*cosmodf[t]+sinalfr*sinmodf[t];
    s = sinalfr*cosmodf[t]-cosalfr*sinmodf[t];
    c2s = 2.*sqr(c);
    cs = c*s;
    /* modulation factors */
    a[t] = c1*(2.-c2d)*c2s+c2*(2.-c2d)*2.*cs+c3*c2sd*c+c4*c2sd*s-	\
      c1*(2.-c2d)+c5*c2d;
    b[t] = c6*sindel*c2s+c7*sindel*2.*cs+c8*cosdel*c+c9*cosdel*s-	\
      c6*sindel;
    as += sqr(a[t]);
    bs += sqr(b[t]);
  } /* for t */
  as /= Np;
  bs /= Np;
  as = sqrt (as);
  bs = sqrt (bs);
  for (t=0; t<Np; t++) {
    a[t] /= as;
    b[t] /= bs;
  } 
} /* modvir() */

double *
JobCore(int pm,			// hemisphere 
	int mm,			// grid 'sky position'
	int nn,			// other grid 'sky position' 
	int smin,		// grid spindown coordinate
	int smax,		// spindown range limit 
	double *M,		// grid generating matrix
	double *DetSSB,		// ephemerides array 
	double *xDat,		// time-domain input data array 
	int Npoints,		// Number of data points
	int Ninterp,		// interpftpad*nfft (for resampling)
	int nfft,		// size of the FFT transform
	fftw_complex *xDa,	// Array for resampling
	fftw_complex *rDa,	// Array for resampling
	fftw_complex *xa,	// input array for plan
	fftw_complex *xao,	// output array for plan
	fftw_plan pl_int,	// fftw_plan needed for resampling
	fftw_plan pl_inv,	// fftw_plan needed for resampling
				// (inverse transformation)
	fftw_plan plan,		// main fftw_plan 
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
				// (INT of FFT, see lvcvirgo.h) 
	int *FNum,		// Candidate signal number
	double coft,		// = oms
	double trl,		// F-statistic threshold
	double sig2,		// N*(variance of xDat) if white_flag
				// else sig2=-1 
	int s0			// No-spindown flag
	) {

  int j, i, nyqst;
  double al1, al2, sinalt, cosalt, sindelt, cosdelt, sgnlt[NPAR], \
    nSource[3], het0, sgnl0, *sgnlv, ft;

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
    double shft1, phase, cp, sp;
    complex double exph;
    fftw_complex *xDb, *rDb, *xb, *xbo;

    xDb = xDa+nfft;
    rDb = rDa+Ninterp;
    xb = xbo = NULL;

    // Switch for different interpolation methods 
    // (set while calling JobCore() from JobNAllSky.c)    
    switch (fftinterp) {
    case INT:
      xb = xa+nfft;
      xbo = xao+nfft;
      break;
    case FFT:
      xb = xa+fftpad*nfft;
    } /* case fftinterp */

    lin2ast (al1/coft, al2/coft, pm, sepsm, cepsm,
	     &sinalt, &cosalt, &sindelt, &cosdelt);

    // calculate declination and right ascention
    // written in file as candidate signal sky positions 
    sgnlt[2] = asin (sindelt);
    sgnlt[3] = fmod (atan2 (sinalt, cosalt)+2.*M_PI, 2.*M_PI);

    /* amplitude modulation functions */
    modvir (sinalt, cosalt, sindelt, cosdelt,
	    sphir, cphir, aa, bb, Npoints);
    nSource[0] = cosalt*cosdelt;
    nSource[1] = sinalt*cosdelt;
    nSource[2] = sindelt;
    shft1 = 0.;
    for (j=0; j<3; j++)
      shft1 += nSource[j]*DetSSB[j];
    het0 = fmod (nn*M[8]+mm*M[12], M[0]);
    for (i=0; i<Npoints; i++) {
      shft[i] = 0.;
      for (j=0; j<3; j++)
	shft[i] += nSource[j]*DetSSB[i*3+j];
      shftf[i] = shft[i]-shft1;
      /* Phase modulation */
      phase = het0*i+oms*shft[i];
#ifdef HAVE_SINCOS
      sincos (phase, &sp, &cp);
#else
      cp = cos (phase);
      sp = sin (phase);
#endif
      exph = cp - I*sp;
      /* Matched filter */
      xDatma[i] = xDat[i]*aa[i]*exph;
      xDatmb[i] = xDat[i]*bb[i]*exph;
    } /* for i */

      /* Resampling */
    for (i=0; i<Npoints; i++) {
      xDa[i] = xDatma[i];
      xDb[i] = xDatmb[i];
    }
    for (i=Npoints; i<nfft; i++) {
      xDa[i] = 0.;
      xDb[i] = 0.;
    }
    fftw_execute (pl_int);
    nyqst = (nfft+2)/2;
    for (i=0; i<nyqst; i++)
      rDb[i] = xDb[i];
    for (i=nyqst; i<nyqst+Ninterp-nfft; i++)
      rDb[i] = 0.;
    for (i=nyqst+Ninterp-nfft,j=nyqst; i<Ninterp; i++,j++) {
      rDb[i] = xDb[j];
      rDa[i] = xDa[j];
    }
    for (i=nyqst; i<nyqst+Ninterp-nfft; i++)
      rDa[i] = 0.;
    fftw_execute (pl_inv);
    ft = (double)interpftpad/Ninterp;
    for (i=0; i<Ninterp; i++) {
      rDa[i] *= ft;
      rDb[i] *= ft;
    }
    splintpad (rDa, shftf, Npoints, interpftpad, xDatma);
    splintpad (rDb, shftf, Npoints, interpftpad, xDatmb);
    /* nearest neighbor interpolation */
    /*	    nnintab (xDatma, shftf, Npoints);
	    nnintab (xDatmb, shftf, Npoints);
    */

    if (write_st)
      printf ("\n>>%d\t%d\t%d\t[%d..%d]\n", *FNum, mm, nn, smin, smax);

    // if no-spindown
    if (s0) smin = smax;

    // if spindown parameter is taken into account, 
    // smin != smax 
    for (ss=smin; ss<=smax; ss++) {

      // Spindown parameter
      sgnlt[1] = (s0 ? 0. : ss*M[5] + nn*M[9] + mm*M[13]);

	if (sgnlt[1] >= -Smax && sgnlt[1] <= 0.) {
	  int ii;
	  double phase2, cosPH, sinPH, Fc, het1;
	  if (write_st) {
	    printf (".");
	    fflush (stdout);
	  }

	  het1 = fmod(ss*M[4], M[0]);
	  if (het1<0) het1 += M[0];

	  sgnl0 = het0+het1;

	  for (i=0; i<Npoints; i++) {
	    phase2 = het1*i + sgnlt[1]*(t2[i]+2.*i*shft[i]);
#ifdef HAVE_SINCOS
	    sincos (phase2, &sinPH, &cosPH);
#else
	    cosPH = cos (phase2);
	    sinPH = sin (phase2);
#endif
	    exph = cosPH - I*sinPH;
	    xa[i] = xDatma[i]*exph;
	    xb[i] = xDatmb[i]*exph;
	  } /* for i */
	  switch (fftinterp) {
	  case INT:
	    for (i=Npoints; i<nfft; i++)
	      xa[i] = xb[i] = 0.;
	    fftw_execute (plan);
	    /* Interpolation algorithm */
	    for (i=nfft/2-1; i>=0; i--) {
	      xao[2*i] = xao[i];
	      xbo[2*i] = xbo[i];
	    }
	    for (i=1; i<nfft-1; i+=2) {
	      xao[i] = (xao[i-1]-xao[i+1])/M_SQRT2;
	      xbo[i] = (xbo[i-1]-xbo[i+1])/M_SQRT2;
	    }
	    break;
	  case FFT:
	    for (i=Npoints; i<fftpad*nfft; i++)
	      xa[i] = xb[i] = 0.;
	    fftw_execute (plan);
	    xao = xa;
	    xbo = xb;
	  }
	  (*FNum)++;

	  for (i=nmin; i<nmax; i++)
	    F[i] = sqr (creal(xao[i])) + sqr (cimag(xao[i])) +
	      sqr (creal(xbo[i])) + sqr (cimag(xbo[i]));

	  /* Normalize F-statistics */
          // if the noise is not white noise
	  if (sig2 < 0.) FStat (F+nmin, nmax-nmin, NAV, 0);
	  else for (i=nmin; i<nmax; i++) F[i] /= sig2;
	  
	  for (i=nmin; i<nmax; i++) {
	    if ((Fc = F[i]) > trl) {
	      /* Find local maximum for neighboring signals */
	      ii = i;
	      while (++i < nmax && F[i] > trl) {
		if (F[i] >= Fc) {
		  ii = i;
		  Fc = F[i];
		} /* if F[i] */
	      } /* while i */
	      // Candidate signal frequency 
	      sgnlt[0] = 2.*M_PI*ii/((double)fftpad*nfft)+sgnl0;
	      // Signal-to-noise ratio
	      sgnlt[4] = sqrt (2.*(Fc-nd));

	      (*sgnlc)++;
	      sgnlv = (double *) realloc (sgnlv,			\
					  NPAR*(*sgnlc)*sizeof (double));
	      for (j=0; j<NPAR; j++)
			  sgnlv[NPAR*(*sgnlc-1)+j] = sgnlt[j];

	      if (write_st)
		printf ("\nSignal %d: %d %d %d %d %d \tsnr=%.2f\n", \
			*sgnlc, pm, mm, nn, ss, ii, sgnlt[4]);
	    } /* if Fc > trl */
	  } /* for i */
	} /* if sgnlt[1] */
    } /* for ss */

  return sgnlv;
} /* JobCore() */
