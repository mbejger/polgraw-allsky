#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <malloc.h>
#include <gsl/gsl_vector.h>
#include <complex.h>
#include <fftw3.h>

/* JobCore file */
#include "jobcore.h"
#include "auxi.h"
#include "settings.h"
#include "timer.h"

#include <assert.h>
#if defined(SLEEF)
//#include "sleef-2.80/purec/sleef.h"
#include <sleefsimd.h>
#elif defined(YEPPP)
#include <yepLibrary.h>
#include <yepCore.h>
#include <yepMath.h>
#endif


void save_array(complex double *arr, int N, const char* file) {
  int i;
  FILE *fc = fopen(file, "w");
  for (i=0; i<N; i++) {
    fprintf(fc, "%d %e + i %e\n", i, creal(arr[i]), cimag(arr[i]));
  }
  fclose(fc);
}

void save_array_double(double *arr, int N, const char* file) {
  int i;
  FILE *fc = fopen(file, "w");
  for (i=0; i<N; i++) {
    fprintf(fc, "%d %e\n", i, arr[i]);
  }
  fclose(fc);
}


// Main searching function (loops inside)
void search(
  Search_settings *sett,
  Command_line_opts *opts,
  Search_range *s_range,
  FFTW_plans *plans,
  FFTW_arrays *fftw_arr,
  Aux_arrays *aux,
  int *FNum,
  double *F) {

  // struct stat buffer;
  struct flock lck;

  int pm, mm, nn;    // hemisphere, sky positions 
  int sgnlc;         // number of candidates
  FLOAT_TYPE *sgnlv;    // array with candidates data

  char outname[512];
  int fd;
  FILE *state;

#ifdef TIMERS
  struct timespec tstart = get_current_time(), tend;
#endif

  state=NULL;
  if(opts->checkp_flag) 
    state = fopen (opts->qname, "w");

  /* Loop over hemispheres
   */ 

  for (pm=s_range->pst; pm<=s_range->pmr[1]; ++pm) {

    sprintf (outname, "%s/triggers_%03d_%03d%s_%d.bin", 
	     opts->prefix, opts->ident, opts->band, opts->label, pm);
    
    /* Two main loops over sky positions 
     */ 

    for (mm=s_range->mst; mm<=s_range->mr[1]; ++mm) {	
      for (nn=s_range->nst; nn<=s_range->nr[1]; ++nn) {	
	
        if(opts->checkp_flag) {
          ftruncate(fileno(state), 0);  
  	      fprintf(state, "%d %d %d %d %d\n", pm, mm, nn, s_range->sst, *FNum);
		      fseek(state, 0, SEEK_SET);
	    }
	
	    /* Loop over spindowns is inside job_core() */
	    sgnlv = job_core(
			 pm,           // hemisphere
			 mm,           // grid 'sky position'
			 nn,           // other grid 'sky position'
			 sett,         // search settings
			 opts,         // cmd opts
			 s_range,      // range for searching
			 plans,        // fftw plans 
			 fftw_arr,     // arrays for fftw
			 aux,          // auxiliary arrays
			 F,            // F-statistics array
			 &sgnlc,       // reference to array with the parameters
                     // of the candidate signal
                     // (used below to write to the file)
			 FNum);        // Candidate signal number
	
	    // Get back to regular spin-down range
	    s_range->sst = s_range->spndr[0];

	    /* Add trigger parameters to a file */

	    // if any signals found (Fstat>Fc)
	    if (sgnlc) {
	      if((fd = open (outname, O_WRONLY|O_CREAT|O_APPEND,
            S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH)) < 0) {
	        perror (outname);
	        return;
	      }

	      lck.l_type = F_WRLCK;
	      lck.l_whence = 0;
	      lck.l_start = 0L;
	      lck.l_len = 0L;

          if (fcntl (fd, F_SETLKW, &lck) < 0) perror ("fcntl()");
          write (fd, (void *)(sgnlv), sgnlc*NPAR*sizeof(FLOAT_TYPE));
          if (close (fd) < 0) perror ("close()");

	    } /* if sgnlc */
	    free (sgnlv);
      } // for nn
      s_range->nst = s_range->nr[0];
    } // for mm
    s_range->mst = s_range->mr[0]; 
  } // for pm

  if(opts->checkp_flag) 
    fclose(state); 

#ifdef TIMERS
  tend = get_current_time();
  // printf("tstart = %d . %d\ntend = %d . %d\n", tstart.tv_sec, tstart.tv_usec, tend.tv_sec, tend.tv_usec);
  double time_elapsed = get_time_difference(tstart, tend);
  printf("Time elapsed: %e s\n", time_elapsed);
#endif

}


  /* Main job 
   */ 

FLOAT_TYPE* job_core(
  int pm,                    // Hemisphere
  int mm,                    // Grid 'sky position'
  int nn,                    // Second grid 'sky position'
  Search_settings *sett,     // Search settings
  Command_line_opts *opts,   // Search options 
  Search_range *s_range,     // Range for searching
  FFTW_plans *plans,         // Plans for fftw
  FFTW_arrays *fftw_arr,     // Arrays for fftw
  Aux_arrays *aux,           // Auxiliary arrays
  double *F,                 // F-statistics array
  int *sgnlc,                // Candidate trigger parameters 
  int *FNum) {               // Candidate signal number

  int i, j, n;
  int smin = s_range->sst, smax = s_range->spndr[1];
  double al1, al2, sinalt, cosalt, sindelt, cosdelt, sgnlt[NPAR], 
    nSource[3], het0, sgnl0, ft;
  double _tmp1[sett->nifo][sett->N];
  FLOAT_TYPE *sgnlv;

#undef NORMTOMAX
#ifdef NORMTOMAX
  double blkavg, threshold = 4.;
  int imax, imax0, iblk, blkstart, ihi;
  int blksize = 1024;
  int nfft = sett->nmax - sett->nmin;
  int *Fmax;
  Fmax = (int *) malloc(2*sett->nfft*sizeof(int));
#endif
  
  /* Matrix	M(.,.) (defined on page 22 of PolGrawCWAllSkyReview1.pdf file)
     defines the transformation form integers (bin, ss, nn, mm) determining
     a grid point to linear coordinates omega, omegadot, alpha_1, alpha_2),
     where bin is the frequency bin number and alpha_1 and alpha_2 are
     defined on p. 22 of PolGrawCWAllSkyReview1.pdf file.

     [omega]                          [bin]
     [omegadot]       = M(.,.) \times [ss]
     [alpha_1/omega]                  [nn]
     [alpha_2/omega]                  [mm]

     Array M[.] is related to matrix M(.,.) in the following way;

                 [ M[0] M[4] M[8]  M[12] ]
      M(.,.) =   [ M[1] M[5] M[9]  M[13] ]
                 [ M[2] M[6] M[10] M[14] ]
                 [ M[3] M[7] M[11] M[15] ]

     and

     M[1] = M[2] = M[3] = M[6] = M[7] = 0
  */

  // Grid positions
  al1 = nn*sett->M[10] + mm*sett->M[14];
  al2 = nn*sett->M[11] + mm*sett->M[15];

  sgnlv = NULL;
  *sgnlc = 0;

  // check if the search is in an appropriate region of the grid
  // if not, returns NULL
  if ((sqr(al1)+sqr(al2))/sqr(sett->oms) > 1.) return NULL ;

  int ss;
  double shft1, phase, cp, sp;
  complex double exph;

  // Change linear (grid) coordinates to real coordinates
  lin2ast(al1/sett->oms, al2/sett->oms, 
	  pm, sett->sepsm, sett->cepsm,
	  &sinalt, &cosalt, &sindelt, &cosdelt);

  // calculate declination and right ascention
  // written in file as candidate signal sky positions
  sgnlt[2] = asin(sindelt);
  sgnlt[3] = fmod(atan2(sinalt, cosalt) + 2.*M_PI, 2.*M_PI);

  het0 = fmod(nn*sett->M[8] + mm*sett->M[12], sett->M[0]);

  // Nyquist frequency 
  int nyqst = (sett->nfft)/2 + 1;

  // Loop for each detector 
  for(n=0; n<sett->nifo; ++n) { 

    /* Amplitude modulation functions aa and bb 
     * for each detector (in signal sub-struct 
     * of _detector, ifo[n].sig.aa, ifo[n].sig.bb) 
     */

    modvir(sinalt, cosalt, sindelt, cosdelt,
	   sett->N, &ifo[n], aux);
	
    // Calculate detector positions with respect to baricenter
    nSource[0] = cosalt*cosdelt;
    nSource[1] = sinalt*cosdelt;
    nSource[2] = sindelt;

    shft1 = nSource[0]*ifo[n].sig.DetSSB[0]
          + nSource[1]*ifo[n].sig.DetSSB[1]
          + nSource[2]*ifo[n].sig.DetSSB[2];

    //    printf("shft1=%f\n",shft1 );
    //printf("DetSSB[0]=%f\n",ifo[n].sig.DetSSB[0] );


    for(i=0; i<sett->N; ++i) {

      ifo[n].sig.shft[i] = nSource[0]*ifo[n].sig.DetSSB[i*3]
	                 + nSource[1]*ifo[n].sig.DetSSB[i*3+1]
	                 + nSource[2]*ifo[n].sig.DetSSB[i*3+2];
    
      ifo[n].sig.shftf[i] = ifo[n].sig.shft[i] - shft1;
      _tmp1[n][i] = aux->t2[i] + (double)(2*i)*ifo[n].sig.shft[i];

      // Phase modulation 
      phase = het0*i + sett->oms*ifo[n].sig.shft[i];
#ifdef NOSINCOS
      cp = cos(phase);
      sp = sin(phase);
#else
      sincos(phase, &sp, &cp);
#endif

      exph = cp - I*sp;

      // Matched filter 
      ifo[n].sig.xDatma[i] = 
        ifo[n].sig.xDat[i]*ifo[n].sig.aa[i]*exph;
      ifo[n].sig.xDatmb[i] = 
        ifo[n].sig.xDat[i]*ifo[n].sig.bb[i]*exph;
  
    }

    /* Resampling using spline interpolation:
     * This will double the sampling rate 
     */ 
  
    for(i=0; i < sett->N; ++i) {
      fftw_arr->xa[i] = ifo[n].sig.xDatma[i];
      fftw_arr->xb[i] = ifo[n].sig.xDatmb[i];
    }
 
    // Zero-padding (filling with 0s up to sett->nfft, 
    // the nearest power of 2)
    for (i=sett->N; i<sett->nfft; ++i) {
      fftw_arr->xa[i] = 0.;
      fftw_arr->xb[i] = 0.;
    }

    fftw_execute(plans->pl_int);  //forward fft (len nfft)
    fftw_execute(plans->pl_int2); //forward fft (len nfft)

    // move frequencies from second half of spectrum; 
    // loop length: nfft - nyqst = nfft - nfft/2 - 1 = nfft/2 - 1
    for(i=nyqst + sett->Ninterp - sett->nfft, j=nyqst; 
        i<sett->Ninterp; ++i, ++j) {
      fftw_arr->xa[i] = fftw_arr->xa[j];
      fftw_arr->xb[i] = fftw_arr->xb[j];
    }
	
    //  zero frequencies higher than nyquist
    for (i=nyqst; i<nyqst + sett->Ninterp - sett->nfft; ++i) {
      fftw_arr->xa[i] = 0.;
      fftw_arr->xb[i] = 0.;
    }

    // Backward fft (len Ninterp = nfft*interpftpad)
    fftw_execute (plans->pl_inv);     
    fftw_execute (plans->pl_inv2); 

    ft = (double)sett->interpftpad / sett->Ninterp; //scale FFT
    for (i=0; i < sett->Ninterp; ++i) {
      fftw_arr->xa[i] *= ft;
      fftw_arr->xb[i] *= ft;
    }

    //  struct timeval tstart = get_current_time(), tend;

    // Spline interpolation to xDatma, xDatmb arrays
    splintpad(fftw_arr->xa, ifo[n].sig.shftf, sett->N, 
      sett->interpftpad, ifo[n].sig.xDatma);   
    splintpad(fftw_arr->xb, ifo[n].sig.shftf, sett->N, 
      sett->interpftpad, ifo[n].sig.xDatmb);

  } // end of detector loop 

  // square sums of modulation factors 
  double aa = 0., bb = 0.; 

  for(n=0; n<sett->nifo; ++n) {

    double aatemp = 0., bbtemp = 0.;
 
    for(i=0; i<sett->N; ++i) {
      aatemp += sqr(ifo[n].sig.aa[i]);
      bbtemp += sqr(ifo[n].sig.bb[i]);
    }

    for(i=0; i<sett->N; ++i) {
      ifo[n].sig.xDatma[i] /= ifo[n].sig.sig2;
      ifo[n].sig.xDatmb[i] /= ifo[n].sig.sig2;
    }

    aa += aatemp/ifo[n].sig.sig2; 
    bb += bbtemp/ifo[n].sig.sig2;   
  }

  //  printf("maa=%f,  mbb=%f\n", aa, bb);
  //  exit(0);

#ifdef YEPPP
#define VLEN 2048
    yepLibrary_Init();
    //printf("npoints=%d, size=%d\n", Npoints, Npoints*sizeof(Yep64f));
    Yep64f _p[VLEN];
    Yep64f _s[VLEN];
    Yep64f _c[VLEN];
    enum YepStatus status;
    int bnd = (sett->N/VLEN)*VLEN;
    //printf("npoints=%d, bnd=%d\n", sett->N, bnd);
#endif
#ifdef SLEEF
    double _p[VECTLENDP], _c[VECTLENDP];
    vdouble2 v;
    vdouble a;
#endif


  /* Spindown loop 
   */

#if TIMERS>2
  struct timespec tstart, tend;
  double spindown_timer = 0;
  int spindown_counter  = 0;
#endif

  printf ("\n>>%d\t%d\t%d\t[%d..%d]\n", *FNum, mm, nn, smin, smax);

  // No-spindown calculations
  if(opts->s0_flag) smin = smax;
  // if spindown parameter is taken into account, smin != smax
  for(ss=smin; ss<=smax; ++ss) {

#if TIMERS>2
    tstart = get_current_time();
#endif 

    // Spindown parameter
    sgnlt[1] = (opts->s0_flag ? 0. : ss*sett->M[5] + nn*sett->M[9] + mm*sett->M[13]);
    
    // Spindown range
    if(sgnlt[1] >= -sett->Smax && sgnlt[1] <= sett->Smax) { 

      int ii;
      double Fc, het1;

#ifdef VERBOSE
      //print a 'dot' every new spindown
      printf ("."); fflush (stdout);
#endif 

      het1 = fmod(ss*sett->M[4], sett->M[0]);
      if(het1<0) het1 += sett->M[0];

      sgnl0 = het0 + het1;

      // phase modulation before fft

#if defined(SLEEF)
      // use simd sincos from the SLEEF library;
      // VECTLENDP is a simd vector length defined in the SLEEF library
      // and it depends on selected instruction set e.g. -DENABLE_AVX
      for (i=0; i<sett->N; i+=VECTLENDP) {
	for(j=0; j<VECTLENDP; j++)
	  _p[j] =  het1*(i+j) + sgnlt[1]*_tmp1[0][i+j];
	
	a = vloadu(_p);
	v = xsincos(a);
	vstoreu(_p, v.x); // reuse _p for sin
	vstoreu(_c, v.y);
	
	for(j=0; j<VECTLENDP; ++j){
	  exph = _c[j] - I*_p[j];
	  fftw_arr->xa[i+j] = ifo[0].sig.xDatma[i+j]*exph; ///ifo[0].sig.sig2;
	  fftw_arr->xb[i+j] = ifo[0].sig.xDatmb[i+j]*exph; ///ifo[0].sig.sig2;
	}
      } 
#elif defined(YEPPP)
      // use yeppp! library;
      // VLEN is length of vector to be processed
      // for caches L1/L2 64/256kb optimal value is ~2048
      for (j=0; j<bnd; j+=VLEN) {
	for (i=0; i<VLEN; ++i)
	  _p[i] =  het1*(i+j) + sgnlt[1]*_tmp1[0][i+j];
	
	status = yepMath_Sin_V64f_V64f(_p, _s, VLEN);
	assert(status == YepStatusOk);
	status = yepMath_Cos_V64f_V64f(_p, _c, VLEN);
	assert(status == YepStatusOk);
	
	for (i=0; i<VLEN; ++i) {
          exph = _c[i] - I*_s[i];
	  fftw_arr->xa[i+j] = ifo[0].sig.xDatma[i+j]*exph; ///ifo[0].sig.sig2;
	  fftw_arr->xb[i+j] = ifo[0].sig.xDatmb[i+j]*exph; ///ifo[0].sig.sig2;
	}
      }
      // remaining part is shorter than VLEN - no need to vectorize
      for (i=0; i<sett->N-bnd; ++i){
	j = bnd + i;
	_p[i] =  het1*j + sgnlt[1]*_tmp1[0][j];
      }

      status = yepMath_Sin_V64f_V64f(_p, _s, sett->N-bnd);
      assert(status == YepStatusOk);
      status = yepMath_Cos_V64f_V64f(_p, _c, sett->N-bnd);
      assert(status == YepStatusOk);

      for (i=0; i<sett->N-bnd; ++i) {
	j = bnd + i;
	exph = _c[i] - I*_s[i];
	fftw_arr->xa[j] = ifo[0].sig.xDatma[j]*exph; ///ifo[0].sig.sig2;
	fftw_arr->xb[j] = ifo[0].sig.xDatmb[j]*exph; ///ifo[0].sig.sig2;
      }
#elif defined(GNUSINCOS)
      for(i=sett->N-1; i!=-1; --i) {
        phase = het1*i + sgnlt[1]*_tmp1[0][i];
	sincos(phase, &sp, &cp);
	exph = cp - I*sp;
        fftw_arr->xa[i] = ifo[0].sig.xDatma[i]*exph; ///ifo[0].sig.sig2;
        fftw_arr->xb[i] = ifo[0].sig.xDatmb[i]*exph; ///ifo[0].sig.sig2;
      }
#else
      for(i=sett->N-1; i!=-1; --i) {
        phase = het1*i + sgnlt[1]*_tmp1[0][i];
	cp = cos(phase);
      	sp = sin(phase);
	exph = cp - I*sp;
        fftw_arr->xa[i] = ifo[0].sig.xDatma[i]*exph; ///ifo[0].sig.sig2;
        fftw_arr->xb[i] = ifo[0].sig.xDatmb[i]*exph; ///ifo[0].sig.sig2;
      }
#endif

      for(n=1; n<sett->nifo; ++n) {
#if defined(SLEEF)
      // use simd sincos from the SLEEF library;
      // VECTLENDP is a simd vector length defined in the SLEEF library
      // and it depends on selected instruction set e.g. -DENABLE_AVX
	for (i=0; i<sett->N; i+=VECTLENDP) {
	  for(j=0; j<VECTLENDP; j++)
	    _p[j] =  het1*(i+j) + sgnlt[1]*_tmp1[n][i+j];
	
	  a = vloadu(_p);
	  v = xsincos(a);
	  vstoreu(_p, v.x); // reuse _p for sin
	  vstoreu(_c, v.y);
	
	  for(j=0; j<VECTLENDP; ++j){
	    exph = _c[j] - I*_p[j];
	    fftw_arr->xa[i+j] = ifo[n].sig.xDatma[i+j]*exph; //*sig2inv;
	    fftw_arr->xb[i+j] = ifo[n].sig.xDatmb[i+j]*exph; //*sig2inv;
	  }
	} 
#elif defined(YEPPP)
	// use yeppp! library;
	// VLEN is length of vector to be processed
	// for caches L1/L2 64/256kb optimal value is ~2048
	for (j=0; j<bnd; j+=VLEN) {
	  for (i=0; i<VLEN; ++i)
	    _p[i] =  het1*(i+j) + sgnlt[1]*_tmp1[n][i+j];
	
	  status = yepMath_Sin_V64f_V64f(_p, _s, VLEN);
	  assert(status == YepStatusOk);
	  status = yepMath_Cos_V64f_V64f(_p, _c, VLEN);
	  assert(status == YepStatusOk);
	
	  for (i=0; i<VLEN; ++i) {
	    exph = _c[i] - I*_s[i];
	    fftw_arr->xa[i+j] += ifo[n].sig.xDatma[i+j]*exph; //*sig2inv;
	    fftw_arr->xb[i+j] += ifo[n].sig.xDatmb[i+j]*exph; //*sig2inv;
	  }
	}
	// remaining part is shorter than VLEN - no need to vectorize
	for (i=0; i<sett->N-bnd; ++i){
	  j = bnd + i;
	  _p[i] =  het1*j + sgnlt[1]*_tmp1[n][j];
	}

	status = yepMath_Sin_V64f_V64f(_p, _s, sett->N-bnd);
	assert(status == YepStatusOk);
	status = yepMath_Cos_V64f_V64f(_p, _c, sett->N-bnd);
	assert(status == YepStatusOk);

	for (i=0; i<sett->N-bnd; ++i) {
	  j = bnd + i;
	  exph = _c[i] - I*_s[i];
	  fftw_arr->xa[j] += ifo[n].sig.xDatma[j]*exph; //*sig2inv;
	  fftw_arr->xb[j] += ifo[n].sig.xDatmb[j]*exph; //*sig2inv;
	}

#elif defined(GNUSINCOS)
	for(i=sett->N-1; i!=-1; --i) {
	  phase = het1*i + sgnlt[1]*_tmp1[n][i];
	  sincos(phase, &sp, &cp);
	  exph = cp - I*sp;
	  fftw_arr->xa[i] += ifo[n].sig.xDatma[i]*exph; //*sig2inv;
	  fftw_arr->xb[i] += ifo[n].sig.xDatmb[i]*exph; //*sig2inv;
	}
#else
	for(i=sett->N-1; i!=-1; --i) {
	  phase = het1*i + sgnlt[1]*_tmp1[n][i];
	  cp = cos(phase);
	  sp = sin(phase);
	  exph = cp - I*sp;
	  fftw_arr->xa[i] += ifo[n].sig.xDatma[i]*exph; //*sig2inv;
	  fftw_arr->xb[i] += ifo[n].sig.xDatmb[i]*exph; //*sig2inv;
	}
#endif

      } 

      // Zero-padding 
      for(i = sett->fftpad*sett->nfft-1; i != sett->N-1; --i)
	    fftw_arr->xa[i] = fftw_arr->xb[i] = 0.; 

      fftw_execute (plans->plan);
      fftw_execute (plans->plan2);

      (*FNum)++;

    // Computing F-statistic 
      for (i=sett->nmin; i<sett->nmax; i++) {
	F[i] = (sqr(creal(fftw_arr->xa[i])) 
		+ sqr(cimag(fftw_arr->xa[i])))/aa  
	  + (sqr(creal(fftw_arr->xb[i])) 
             + sqr(cimag(fftw_arr->xb[i])))/bb;
      }

#if 0
      FILE *f1 = fopen("fraw-1.dat", "w");
      for(i=sett->nmin; i<sett->nmax; i++)
	fprintf(f1, "%d   %lf   %lf\n", i, F[i], 2.*M_PI*i/((double) sett->fftpad*sett->nfft) + sgnl0);
      fclose(f1);
#endif 

#ifndef NORMTOMAX
#define NAVFSTAT 4096
      // Normalize F-statistics 
      if(!(opts->white_flag))  // if the noise is not white noise
        FStat(F + sett->nmin, sett->nmax - sett->nmin, NAVFSTAT, 0);

      // f1 = fopen("fnorm-4096-1.dat", "w");
      //for(i=sett->nmin; i<sett->nmax; i++)
      //fprintf(f1, "%d   %lf   %lf\n", i, F[i], 2.*M_PI*i/((double) sett->fftpad*sett->nfft) + sgnl0);
      //fclose(f1);
      //      exit(EXIT_SUCCESS);

      for(i=sett->nmin; i<sett->nmax; i++) {
        if ((Fc = F[i]) > opts->trl) { // if F-stat exceeds trl (critical value)
	  printf("i=%d\n",i);
          // Find local maximum for neighboring signals 
          ii = i;

	  while (++i < sett->nmax && F[i] > opts->trl) {
	    if(F[i] >= Fc) {
	      ii = i;
	      Fc = F[i];
	    } // if F[i] 
	  } // while i 
	  // Candidate signal frequency
	  sgnlt[0] = 2.*M_PI*ii/((double) sett->fftpad*sett->nfft) + sgnl0;
	  // Signal-to-noise ratio
	  sgnlt[4] = sqrt(2.*(Fc-sett->nd));
	  
	  (*sgnlc)++; // increase found number
	  
	  // Add new parameters to output array 
	  sgnlv = (FLOAT_TYPE *)realloc(sgnlv, NPAR*(*sgnlc)*sizeof(FLOAT_TYPE));
	  
	  for (j=0; j<NPAR; ++j) // save new parameters
	    sgnlv[NPAR*(*sgnlc-1)+j] = (FLOAT_TYPE)sgnlt[j];

#ifdef VERBOSE
	  printf ("\nSignal %d: %d %d %d %d %d \tsnr=%.2f\n", 
		  *sgnlc, pm, mm, nn, ss, ii, sgnlt[4]);
#endif 
	  
	} // if Fc > trl 
      } // for i
      
#else // new version

      imax = -1;
      // find local maxima first
      //printf("nmin=%d   nmax=%d    nfft=%d   nblocks=%d\n", sett->nmin, sett->nmax, nfft, nfft/blksize);
      for(iblk=0; iblk < nfft/blksize; ++iblk) {
	blkavg = 0.;
	blkstart = sett->nmin + iblk*blksize; // block start index in F 
	// in case the last block is shorter than blksize, include its elements in the previous block
	if(iblk==(nfft/blksize-1)) {blksize = sett->nmax - blkstart; /*printf("blksize modified to %d\n", blksize);*/}
	imax0 = imax+1;// index of first maximum in current block
	//printf("\niblk=%d   blkstart=%d   blksize=%d    imax0=%d\n", iblk, blkstart, blksize, imax0);
	for(i=1; i <= blksize; ++i) { // include first element of the next block,  1..1024
	  ii = blkstart + i;                                                   //  1..1025
	  if(ii < sett->nmax) {ihi=ii+1;} else {ihi = sett->nmax; /*printf("ihi=%d  ii=%d\n", ihi, ii);*/};
	  if(F[ii] > F[ii-1] && F[ii] > F[ihi]) {
	    Fmax[++imax] = ii;
	    blkavg += F[ii];
	    ++i; // next element can't be maximum - skip it
	  }
	} // i
	// now imax points to the last element of Fmax
	// normalize in blocks 
	blkavg /= (double)(imax - imax0 + 1);
	for(i=imax0; i <= imax; ++i)
	  F[Fmax[i]] /= blkavg;

      } // iblk

      //f1 = fopen("fmax.dat", "w");
      //for(i=1; i < imax; i++)
      //fprintf(f1, "%d   %lf \n", Fmax[i], F[Fmax[i]]);
      //fclose(f1);
      //exit(EXIT_SUCCESS);

      // apply threshold limit
      for(i=0; i <= imax; ++i){
	//if(F[Fmax[i]] > opts->trl) {
	if(F[Fmax[i]] > threshold) {
	  sgnlt[0] = 2.*M_PI*i/((double) sett->fftpad*sett->nfft) + sgnl0;
	  // Signal-to-noise ratio
	  sgnlt[4] = sqrt(2.*(F[Fmax[i]] - sett->nd));
	  (*sgnlc)++; // increase found number
	  // Add new parameters to output array 
	  sgnlv = (FLOAT_TYPE *)realloc(sgnlv, NPAR*(*sgnlc)*sizeof(FLOAT_TYPE));
	  for (j=0; j<NPAR; ++j) // save new parameters
	    sgnlv[NPAR*(*sgnlc-1)+j] = (FLOAT_TYPE)sgnlt[j];
#ifdef VERBOSE
	  printf ("\nSignal %d: %d %d %d %d %d \tsnr=%.2f\n", 
		  *sgnlc, pm, mm, nn, ss, Fmax[i], sgnlt[4]);
#endif 
	}
      } // i
      
#endif // old/new version


#if TIMERS>2
      tend = get_current_time();
      spindown_timer += get_time_difference(tstart, tend);
      spindown_counter++;
#endif

    } // if sgnlt[1] 
  } // for ss 
  
#ifndef VERBOSE 
	    printf("Number of signals found: %d\n", *sgnlc); 
#endif 

#if TIMERS>2
  printf("\nTotal spindown loop time: %e s, mean spindown time: %e s (%d runs)\n",
    spindown_timer, spindown_timer/spindown_counter, spindown_counter);
#endif

  return sgnlv;

} // jobcore

