#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <string.h> 
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <malloc.h>
#include <gsl/gsl_vector.h>
#include <complex.h>
#include <fftw3.h>
#include <signal.h>

/* JobCore file */
#include "jobcore.h"
#include "auxi.h"
#include "settings.h"
#include "timer.h"

#include <assert.h>
#if defined(SLEEF)
//#include "sleef-2.80/purec/sleef.h"
//#include <sleefsimd.h>
#include <sleef.h>
#define DORENAME
#ifdef ENABLE_AVX
#define CONFIG 1
#include "helperavx.h"
#include "renameavx.h"
typedef Sleef___m256d_2 vdouble2;
typedef Sleef___m256_2 vfloat2;
#endif
#elif defined(YEPPP)
#include <yepMath.h>
#include <yepLibrary.h>

#endif

#include <omp.h>

extern volatile sig_atomic_t save_state;

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
	    int *FNum ) {

  // struct stat buffer;
  struct flock lck;
  
  int pm, mm, nn;       // hemisphere, sky positions 
  int sgnlc=0;          // number of candidates
  FLOAT_TYPE *sgnlv;    // array with candidates data
  long totsgnl;        // total number of candidates

  char outname[1100];
  int fd, status;
  FILE *state;

#ifdef YEPPP
  status = yepLibrary_Init();
  assert(status == YepStatusOk);
#endif

#ifdef TIMERS
  struct timespec tstart = get_current_time(CLOCK_REALTIME), tend;
#endif
  
  // Allocate buffer for triggers
  sgnlv = (FLOAT_TYPE *)calloc(NPAR*2*sett->nfft, sizeof(FLOAT_TYPE));

  // open mode for trig file
  int tmode = O_WRONLY|O_CREAT|O_APPEND;

  state = NULL;
  if(opts->checkp_flag) state = fopen (opts->qname, "w");
  
  /* Loop over hemispheres */ 
  
  for (pm=s_range->pst; pm<=s_range->pmr[1]; ++pm) {

    sprintf (outname, "%s/triggers_%03d_%04d%s_%d.bin", 
	     opts->prefix, opts->ident, opts->band, opts->label, pm);
    // remove existing trigger file if checkpointing is disabled
    if(! opts->checkp_flag) remove(outname);

    totsgnl = 0;
    /* Two main loops over sky positions */ 
    
    for (mm=s_range->mst; mm<=s_range->mr[1]; ++mm) {	
      for (nn=s_range->nst; nn<=s_range->nr[1]; ++nn) {	
	
	/* Loop over spindowns is inside job_core() */
	status = job_core(
			  pm,           // hemisphere
			  mm,           // grid 'sky position'
			  nn,           // other grid 'sky position'
			  sett,         // search settings
			  opts,         // cmd opts
			  s_range,      // range for searching
			  plans,        // fftw plans 
			  fftw_arr,     // arrays for fftw
			  aux,          // auxiliary arrays
			  &sgnlc,       // current number of candidates
			  sgnlv,        // candidate array
			  FNum);        // candidate signal number
	
	// Get back to regular spin-down range
	s_range->sst = s_range->spndr[0];

	/* Add trigger parameters to a file */
	// if enough signals found (no. of signals > half length of buffer)
	if (sgnlc > sett->nfft || save_state == 1) {
	     if((fd = open (outname, tmode, S_IRUSR|S_IWUSR|S_IRGRP)) < 0) {
		  perror(outname);
		  return;
	     }

#ifdef USE_LOCKING
	     lck.l_type = F_WRLCK;
	     lck.l_whence = 0;
	     lck.l_start = 0L;
	     lck.l_len = 0L;
	     if (fcntl (fd, F_SETLKW, &lck) < 0) perror ("fcntl()");
#endif
	     write(fd, (void *)(sgnlv), sgnlc*NPAR*sizeof(FLOAT_TYPE));
	     totsgnl += sgnlc;
	     if (close(fd) < 0) perror ("close()");
	     sgnlc=0;
	     
	     if(opts->checkp_flag) {
		  ftruncate(fileno(state), 0);  
		  fprintf(state, "%d %d %d %d %d\n", pm, mm, nn+1, s_range->sst, *FNum);
		  fseek(state, 0, SEEK_SET);
		  if (save_state == 1) {
		       //printf("%d %d %d %d %d\n", pm, mm, nn+1, s_range->sst, *FNum);
		       printf("\nState saved after signal\nExiting\n");
		       exit(EXIT_SUCCESS);
		  }
	     }
	     save_state = 0;
	     
	} /* if sgnlc > sett-nfft */
      } // for nn
      s_range->nst = s_range->nr[0];
    } // for mm
    s_range->mst = s_range->mr[0]; 

    // Write the leftover from the last iteration of the buffer 
    if((fd = open(outname, tmode, S_IRUSR|S_IWUSR|S_IRGRP)) < 0) {
	 perror(outname);
	 return; 
    }

#ifdef USE_LOCKING
    lck.l_type = F_WRLCK;
    lck.l_whence = 0;
    lck.l_start = 0L;
    lck.l_len = 0L;
    if (fcntl (fd, F_SETLKW, &lck) < 0) perror ("fcntl()");
#endif
    write(fd, (void *)(sgnlv), sgnlc*NPAR*sizeof(FLOAT_TYPE));
    totsgnl += sgnlc;
    printf("\n### Total number of signals in %s = %ld\n\n", outname, totsgnl);
    if (close(fd) < 0) perror ("close()");
    sgnlc=0; 
    
  } // for pm
  
  //#mb state file has to be modified accordingly to the buffer
  if(opts->checkp_flag) 
    fclose(state); 

  // Free triggers buffer
  free(sgnlv);
  
#ifdef TIMERS
  tend = get_current_time(CLOCK_REALTIME);
  double time_elapsed = get_time_difference(tstart, tend);
  printf("\nwalltime = %e s | ncpus = %d | cputime = %e\n", time_elapsed, omp_get_max_threads(), time_elapsed*omp_get_max_threads());
#endif

  printf("\nEND\n");

}


  /* Main job */ 

int job_core(int pm,                   // Hemisphere
	     int mm,                   // Grid 'sky position'
	     int nn,                   // Second grid 'sky position'
	     Search_settings *sett,    // Search settings
	     Command_line_opts *opts,  // Search options 
	     Search_range *s_range,    // Range for searching
	     FFTW_plans *plans,        // Plans for fftw
	     FFTW_arrays *fftw_arr,    // Arrays for fftw
	     Aux_arrays *aux,          // Auxiliary arrays
	     int *sgnlc,               // Candidate trigger parameters 
	     FLOAT_TYPE *sgnlv,        // Candidate array 
	     int *FNum) {              // Candidate signal number

  int i, j, n;
  int smin = s_range->sst, smax = s_range->spndr[1];
  double al1, al2, sinalt, cosalt, sindelt, cosdelt, sgnlt[NPAR], 
    nSource[3], het0, sgnl0, ft;
  //double _tmp1[sett->nifo][sett->N];
  static double **_tmp1;
  if (!_tmp1) {
    _tmp1 = (double **)malloc(sett->nifo*sizeof(double *));
    for (n=0; n < sett->nifo; n++) _tmp1[n] = (double *)calloc(sett->N, sizeof(double));
  }
#undef NORMTOMAX
#ifdef NORMTOMAX
  double blkavg, threshold = 6.;
  int imax, imax0, iblk, blkstart, ihi;
  int blksize = 1024;
  int nfft = sett->nmax - sett->nmin;
  static int *Fmax;
  if (!Fmax) Fmax = (int *) malloc(nfft*sizeof(int));
#endif

  struct timespec tstart, tend;
  double spindown_timer = 0;
  int spindown_counter  = 0;
  
  //tstart = get_current_time(CLOCK_REALTIME);

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

  // check if the search is in an appropriate region of the grid
  // if not, returns NULL
  if ((sqr(al1)+sqr(al2))/sqr(sett->oms) > 1.) return 0;

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

#define CHUNK 4
#pragma omp parallel default(shared) private(phase,cp,sp,exph)
    {
#pragma omp for schedule(static,CHUNK)
      for(i=0; i<sett->N; ++i) {
	ifo[n].sig.shft[i] = nSource[0]*ifo[n].sig.DetSSB[i*3]
	                   + nSource[1]*ifo[n].sig.DetSSB[i*3+1]
	                   + nSource[2]*ifo[n].sig.DetSSB[i*3+2];
	ifo[n].sig.shftf[i] = ifo[n].sig.shft[i] - shft1;
	_tmp1[n][i] = aux->t2[i] + (double)(2*i)*ifo[n].sig.shft[i];
      }

#pragma omp for schedule(static,CHUNK) 
      for(i=0; i<sett->N; ++i) {
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
	ifo[n].sig.xDatma[i] = ifo[n].sig.xDat[i]*ifo[n].sig.aa[i]*exph;
	ifo[n].sig.xDatmb[i] = ifo[n].sig.xDat[i]*ifo[n].sig.bb[i]*exph;
  
      }

      /* Resampling using spline interpolation:
       * This will double the sampling rate 
       */ 
#pragma omp for schedule(static,CHUNK)  
      for(i=0; i < sett->N; ++i) {
	fftw_arr->xa[i] = ifo[n].sig.xDatma[i];
	fftw_arr->xb[i] = ifo[n].sig.xDatmb[i];
      }

      // Zero-padding (filling with 0s up to sett->nfft, 
      // the nearest power of 2)
#pragma omp for schedule(static,CHUNK)
      for (i=sett->N; i<sett->nfft; ++i) {
	   fftw_arr->xa[i] = 0.;
	   fftw_arr->xb[i] = 0.;
      }
      
    } //omp parallel

    //printf("before xdatma: %f  %f   %f   %f   %f\n", creal(ifo[n].sig.xDatma[2100]), cimag(ifo[n].sig.xDatma[2100]),  
    //	   creal(ifo[n].sig.xDatma[5000]), cimag(ifo[n].sig.xDatma[5000]), ifo[n].sig.shftf[2100] );

    fftw_execute_dft(plans->pl_int,fftw_arr->xa,fftw_arr->xa);  //forward fft (len nfft)
    fftw_execute_dft(plans->pl_int,fftw_arr->xb,fftw_arr->xb);  //forward fft (len nfft)

#if 1
    // move frequencies from second half of spectrum; 
    // and zero frequencies higher than nyquist
    // loop length: nfft - nyqst = nfft - nfft/2 - 1 = nfft/2 - 1

    for(i=nyqst + sett->Ninterp - sett->nfft, j=nyqst; i<sett->Ninterp; ++i, ++j)
	 fftw_arr->xa[i] = fftw_arr->xa[j];

    for(i=nyqst; i<nyqst + sett->Ninterp - sett->nfft; ++i)
	 fftw_arr->xa[i] = 0.;
    
    for(i=nyqst + sett->Ninterp - sett->nfft, j=nyqst; i<sett->Ninterp; ++i, ++j)
	 fftw_arr->xb[i] = fftw_arr->xb[j];

    for(i=nyqst; i<nyqst + sett->Ninterp - sett->nfft; ++i)
	 fftw_arr->xb[i] = 0.;

    // Backward fft (len Ninterp = nfft*interpftpad)
    fftw_execute_dft(plans->pl_inv,fftw_arr->xa,fftw_arr->xa);
    fftw_execute_dft(plans->pl_inv,fftw_arr->xb,fftw_arr->xb);

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

/*
    // alternative linear interpolation
    linterp(fftw_arr->xa, ifo[n].sig.shftf, sett->N, 
	      sett->interpftpad, ifo[n].sig.xDatma);   
    linterp(fftw_arr->xb, ifo[n].sig.shftf, sett->N, 
	      sett->interpftpad, ifo[n].sig.xDatmb);
*/
#endif
#if 0
    // alternative trigonometric interpolation
    triginterp(fftw_arr->xa, fftw_arr->xb, ifo[n].sig.shftf, sett->N, sett->nfft, ifo[n].sig.xDatma, ifo[n].sig.xDatmb);
    printf("after xdatma: %f  %f   %f   %f\n", creal(ifo[n].sig.xDatma[2100]), cimag(ifo[n].sig.xDatma[2100]),  
	   creal(ifo[n].sig.xDatma[5000]), cimag(ifo[n].sig.xDatma[5000] ) );
//    exit(1);
#endif


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

#ifdef YEPPP
#define VLEN 1024
  int bnd = (sett->N/VLEN)*VLEN;
#endif

  // Check if the signal is added to the data 
  // or the range file is given:  
  // if not, proceed with the wide range of spindowns 
  // if yes, use smin = s_range->sst, smax = s_range->spndr[1]  
  if(!strcmp(opts->addsig, "") && !strcmp(opts->range, "")) {

      // Spindown range defined using Smin and Smax (settings.c)  
      smin = trunc((sett->Smin - nn*sett->M[9] - mm*sett->M[13])/sett->M[5]);
      smax = trunc(-(nn*sett->M[9] + mm*sett->M[13] + sett->Smax)/sett->M[5]);

      // swapping smin and smax in case when grid matrix  
      // values are defined with opposite signs than ''usual''
      if(smin > smax) { 
    
        smin = smin + smax ;
        smax = smin - smax ; 
        smin = smin - smax ; 

      }
  } 

  
  if(opts->s0_flag) smin = smax;
  // if spindown parameter is taken into account, smin != smax

  int s_stride = 1;
  printf ("\n>>%d\t%d\t%d\t[%d..%d:%d]\n", *FNum, mm, nn, smin, smax, s_stride);

  static fftw_complex *fxa, *fxb;
  static double *F;
#pragma omp threadprivate(fxa,fxb, F)
#pragma omp threadprivate(F)

  //private loop counter: ss
  //private (declared inside): ii,Fc,het1,k,veto_status,a,v,_p,_c,_s,status
  //shared default: nn,mm,sett,_tmp1,ifo,het0,bnd,plans,opts,aa,bb,
  //                fftw_arr (zostawiamy i robimy nowe), FNum (atomic!)
  //we use shared plans and  fftw_execute with 'new-array' interface
#pragma omp parallel default(shared)				\
  private(i, j, n, sgnl0, exph, phase, cp, sp, tstart, tend)	\
  firstprivate(sgnlt)						\
  reduction(+ : spindown_timer, spindown_counter)

  {
#ifdef YEPPP
    Yep64f _p[VLEN], _s[VLEN], _c[VLEN];
    enum YepStatus status;
#endif
#ifdef SLEEF
    double _p[VECTLENDP], _c[VECTLENDP];
    vdouble2 v;
    vdouble a;
#endif
    
    if (!fxa) fxa = (fftw_complex *)fftw_malloc(fftw_arr->arr_len*sizeof(fftw_complex));
    if (!fxb) fxb = (fftw_complex *)fftw_malloc(fftw_arr->arr_len*sizeof(fftw_complex));
    if (!F) F = (double *)calloc(2*sett->nfft, sizeof(double));
  
    /* Spindown loop  */

    //#pragma omp for schedule(static,4)
#pragma omp for schedule(static,4)
    for(ss=smin; ss<=smax; ss += s_stride) {

#if TIMERS>2
      tstart = get_current_time(CLOCK_PROCESS_CPUTIME_ID);
#endif 

      // Spindown parameter
      sgnlt[1] = ss*sett->M[5] + nn*sett->M[9] + mm*sett->M[13];

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
      for(i=0; i<sett->N; i+=VECTLENDP) {
	for(j=0; j<VECTLENDP; j++)
	  _p[j] =  het1*(i+j) + sgnlt[1]*_tmp1[0][i+j];
	
	a = vloadu_vd_p(_p);
	v = xsincos_u1(a);
	vstoreu_v_p_vd(_p, v.x); // reuse _p for sin
	vstoreu_v_p_vd(_c, v.y);
	
	for(j=0; j<VECTLENDP; ++j){
	  exph = _c[j] - I*_p[j];
	  fxa[i+j] = ifo[0].sig.xDatma[i+j]*exph; //ifo[0].sig.sig2;
	  fxb[i+j] = ifo[0].sig.xDatmb[i+j]*exph; //ifo[0].sig.sig2;
	}
      } 
#elif defined(YEPPP)
      // use yeppp! library;
      // VLEN is length of vector to be processed
      // for caches L1/L2 64/256kb optimal value is ~2048
      for (j=0; j<bnd; j+=VLEN) {
	//double *_tmp2 = &_tmp1[0][j];
	for (i=0; i<VLEN; ++i)
	  //_p[i] =  het1*(i+j) + sgnlt[1]*_tmp2[i];
       	  _p[i] =  het1*(i+j) + sgnlt[1]*_tmp1[0][i+j];
	
	status = yepMath_Sin_V64f_V64f(_p, _s, VLEN);
	assert(status == YepStatusOk);
	status = yepMath_Cos_V64f_V64f(_p, _c, VLEN);
	assert(status == YepStatusOk);

	for (i=0; i<VLEN; ++i) {
	  //	  exph = _c[i] - I*_s[i];
	  fxa[i+j] = ifo[0].sig.xDatma[i+j]*_c[i]-I*ifo[0].sig.xDatma[i+j]*_s[i];
	  fxb[i+j] = ifo[0].sig.xDatmb[i+j]*_c[i]-I*ifo[0].sig.xDatmb[i+j]*_s[i];
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
	//exph = _c[i] - I*_s[i];
	//fxa[j] = ifo[0].sig.xDatma[j]*exph;
	//fxb[j] = ifo[0].sig.xDatmb[j]*exph;
	fxa[j] = ifo[0].sig.xDatma[j]*_c[i]-I*ifo[0].sig.xDatma[j]*_s[i];
	fxb[j] = ifo[0].sig.xDatmb[j]*_c[i]-I*ifo[0].sig.xDatmb[j]*_s[i];
      }
#elif defined(GNUSINCOS)
      for(i=sett->N-1; i!=-1; --i) {
        phase = het1*i + sgnlt[1]*_tmp1[0][i];
	sincos(phase, &sp, &cp);
	exph = cp - I*sp;
        fxa[i] = ifo[0].sig.xDatma[i]*exph; //ifo[0].sig.sig2;
        fxb[i] = ifo[0].sig.xDatmb[i]*exph; //ifo[0].sig.sig2;
      }
#else
      for(i=sett->N-1; i!=-1; --i) {
        phase = het1*i + sgnlt[1]*_tmp1[0][i];
	cp = cos(phase);
      	sp = sin(phase);
	exph = cp - I*sp;
        fxa[i] = ifo[0].sig.xDatma[i]*exph; //ifo[0].sig.sig2;
        fxb[i] = ifo[0].sig.xDatmb[i]*exph; //ifo[0].sig.sig2;
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
	  
	  a = vloadu_vd_p(_p);
	  v = xsincos_u1(a);
	  vstoreu_v_p_vd(_p, v.x); // reuse _p for sin
	  vstoreu_v_p_vd(_c, v.y);
	
	  for(j=0; j<VECTLENDP; ++j){
	    exph = _c[j] - I*_p[j];
	    fxa[i+j] += ifo[n].sig.xDatma[i+j]*exph;
	    fxb[i+j] += ifo[n].sig.xDatmb[i+j]*exph;
	  }
	} 
#elif defined(YEPPP)
	// use yeppp! library;
	// VLEN is length of vector to be processed
	// for caches L1/L2 64/256kb optimal value is ~2048
	for (j=0; j<bnd; j+=VLEN) {
	  //double *_tmp2 = &_tmp1[n][j];
	  for (i=0; i<VLEN; ++i)
	    //  _p[i] =  het1*(i+j) + sgnlt[1]*_tmp2[i];
	    _p[i] =  het1*(j+i) + sgnlt[1]*_tmp1[n][j+i];
	
	  status = yepMath_Sin_V64f_V64f(_p, _s, VLEN);
	  assert(status == YepStatusOk);
	  status = yepMath_Cos_V64f_V64f(_p, _c, VLEN);
	  assert(status == YepStatusOk);
	
	  for (i=0; i<VLEN; ++i) {
	    //exph = _c[i] - I*_s[i];
	    //fxa[j+i] += ifo[n].sig.xDatma[j+i]*exph;
	    //fxb[j+i] += ifo[n].sig.xDatmb[j+i]*exph;
	    fxa[i+j] += ifo[n].sig.xDatma[i+j]*_c[i]-I*ifo[n].sig.xDatma[i+j]*_s[i];
	    fxb[i+j] += ifo[n].sig.xDatmb[i+j]*_c[i]-I*ifo[n].sig.xDatmb[i+j]*_s[i];
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
	  //exph = _c[i] - I*_s[i];
	  //fxa[j] += ifo[n].sig.xDatma[j]*exph;
	  //fxb[j] += ifo[n].sig.xDatmb[j]*exph;
	  fxa[j] += ifo[n].sig.xDatma[j]*_c[i]-I*ifo[n].sig.xDatma[j]*_s[i];
	  fxb[j] += ifo[n].sig.xDatmb[j]*_c[i]-I*ifo[n].sig.xDatmb[j]*_s[i];
	}

#elif defined(GNUSINCOS)
	for(i=sett->N-1; i!=-1; --i) {
	  phase = het1*i + sgnlt[1]*_tmp1[n][i];
	  sincos(phase, &sp, &cp);
	  exph = cp - I*sp;
	  fxa[i] += ifo[n].sig.xDatma[i]*exph;
	  fxb[i] += ifo[n].sig.xDatmb[i]*exph;
	}
#else
	for(i=sett->N-1; i!=-1; --i) {
	  phase = het1*i + sgnlt[1]*_tmp1[n][i];
	  cp = cos(phase);
	  sp = sin(phase);
	  exph = cp - I*sp;
	  fxa[i] += ifo[n].sig.xDatma[i]*exph;
	  fxb[i] += ifo[n].sig.xDatmb[i]*exph;
	}
#endif

      } 

      // Zero-padding 
      for(i = sett->fftpad*sett->nfft-1; i != sett->N-1; --i)
	fxa[i] = fxb[i] = 0.; 

      fftw_execute_dft(plans->plan, fxa, fxa);
      fftw_execute_dft(plans->plan, fxb, fxb);
      
      // Computing F-statistic 
      for (i=sett->nmin; i<sett->nmax; i++) {
	F[i] = (sqr(creal(fxa[i])) + sqr(cimag(fxa[i])))/aa +
	       (sqr(creal(fxb[i])) + sqr(cimag(fxb[i])))/bb;
      }

      //      for (i=sett->nmin; i<sett->nmax; i++) 
      //	F[i] += (sqr(creal(fxb[i])) + sqr(cimag(fxb[i])))/bb;

#pragma omp atomic
      (*FNum)++;

      
#undef FSTATDEB
#ifdef FSTATDEB
      // warning: use with nthreads=1
      static double *fraw;
      static int ifile=0;
      if (!fraw) fraw = (double *) malloc((sett->nmax-sett->nmin)*sizeof(double));
      memcpy(fraw, F+sett->nmin, (sett->nmax-sett->nmin)*sizeof(double));
#endif

#ifndef NORMTOMAX
      double pxout=0.;
      // Normalize F-statistics 
      if(!(opts->white_flag))  // if the noise is not white noise
	   pxout=FStat(F + sett->nmin, sett->nmax - sett->nmin, NAVFSTAT, 0);

      for(i=sett->nmin; i<sett->nmax; i++) {
        if ((Fc = F[i]) > opts->trl) { // if F-stat exceeds trl (critical value)
          // Find local maximum for neighboring signals 
          ii = i;

	  while (++i < sett->nmax && F[i] > opts->trl) {
	    if(F[i] >= Fc) {
	      ii = i;
	      Fc = F[i];
	    } // if F[i] 
	  } // while i 
	  // Candidate signal frequency
	  sgnlt[0] = 2.*M_PI*ii/((FLOAT_TYPE)sett->fftpad*sett->nfft) + sgnl0;
	  // Signal-to-noise ratio
	  sgnlt[4] = sqrt(2.*(Fc-sett->nd));
	  
#ifdef FSTATDEB
	  if (pxout < 0.6) {
	       printf("pxout  mm=%d  nn=%d  ss=%d\n", mm, nn, ss );
	  }
	  //if (sgnlt[4] > 7.1 && sgnlt[0]< 2.984513) {
	  //if ( (mm==-10 && nn==-18 && ss==100) || (mm==-10 && nn==-18 && ss==-175)) {
	  if ( (mm==-61 && nn==-32 && ss==273) ) {
	       char f1name[32];
	       ifile++;
	       sprintf(f1name, "fraw-%d.dat", ifile);
	       FILE *f1 = fopen(f1name, "w");
	       for(i=0; i<(sett->nmax-sett->nmin); i++)
		    fprintf(f1, "%d   %lf   %lf  %lf  %lf  %lf %lf\n", i, fraw[i], 2.*M_PI*i/((double) sett->fftpad*sett->nfft) + sgnl0, 
			    //	    sqr(creal(fxa[i])), sqr(cimag(fxa[i])), sqr(creal(fxb[i])), sqr(cimag(fxb[i])) );
			    sqr(creal( ifo[0].sig.xDatma[2*i] )), sqr(cimag(ifo[0].sig.xDatmb[2*i])), sqr(creal(ifo[1].sig.xDatma[2*i])), sqr(cimag(ifo[1].sig.xDatmb[2*i])) );
	       fclose(f1);

	       sprintf(f1name, "fnorm-%d.dat", ifile);
	       f1 = fopen(f1name, "w");
	       for(i=sett->nmin; i<sett->nmax; i++)
		    fprintf(f1, "%d   %lf   %lf\n", i, F[i], 2.*M_PI*i/((double) sett->fftpad*sett->nfft) + sgnl0);
	       fclose(f1);
	       printf("Dump mm=%d  nn=%d  ss=%d  al1=%.17g   al2=%.17g   oms=%.17g   one=%.17g\n", mm, nn, ss, al1, al2, sett->oms, (sqr(al1)+sqr(al2))/sqr(sett->oms) );
	       //exit(EXIT_SUCCESS);
	  }
#endif
	  // Checking if signal is within a known instrumental line 
	  int k, veto_status = 0; 
	  for(k=0; k<sett->numlines_band; k++){
	    if(sgnlt[0]>=sett->lines[k][0] && sgnlt[0]<=sett->lines[k][1]) {
	      veto_status=1; 
	      break; 
	    }   
	  }
	  int _sgnlc;
	  if(!veto_status) {

	    /* 
#pragma omp critical
	    {
	      (*sgnlc)++; // increase found number
	      // Add new parameters to output array 
	      for (j=0; j<NPAR; ++j)    // save new parameters
		sgnlv[NPAR*(*sgnlc-1)+j] = (FLOAT_TYPE)sgnlt[j];
	    }
	    */

#pragma omp atomic capture
	    {
	      (*sgnlc)++; // increase found number
	      _sgnlc = *sgnlc;
	    }
	    // Add new parameters to output array 
	    for (j=0; j<NPAR; ++j)    // save new parameters
	      sgnlv[NPAR*(_sgnlc-1)+j] = (FLOAT_TYPE)sgnlt[j];
	    
#ifdef VERBOSE
	    printf ("\nSignal %d: %d %d %d %d %d snr=%.2f\n", 
		    *sgnlc, pm, mm, nn, ss, ii, sgnlt[4]);
#endif
	  }
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
	if(iblk==(nfft/blksize-1)) {blksize = sett->nmax - blkstart;}
	imax0 = imax+1; // index of first maximum in current block
	//printf("\niblk=%d   blkstart=%d   blksize=%d    imax0=%d\n", iblk, blkstart, blksize, imax0);
	for(i=1; i <= blksize; ++i) { // include first element of the next block
	  ii = blkstart + i;
	  if(ii < sett->nmax) 
	    {ihi=ii+1;} 
	  else 
	    {ihi = sett->nmax; /*printf("ihi=%d  ii=%d\n", ihi, ii);*/};
	  if(F[ii] > F[ii-1] && F[ii] > F[ihi]) {
	    blkavg += F[ii];
	    Fmax[++imax] = ii;
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
	  sgnlt[0] = 2.*M_PI*i/((FLOAT_TYPE)sett->fftpad*sett->nfft) + sgnl0;
	  // Signal-to-noise ratio
	  sgnlt[4] = sqrt(2.*(F[Fmax[i]] - sett->nd));

	  // Checking if signal is within a known instrumental line 
	  int k, veto_status=0; 
	  for(k=0; k<sett->numlines_band; k++)
	    if(sgnlt[0]>=sett->lines[k][0] && sgnlt[0]<=sett->lines[k][1]) { 
	      veto_status=1; 
	      break; 
	    }   
	  
	  if(!veto_status) { 
	    
	    (*sgnlc)++; // increase number of found candidates
	    // Add new parameters to buffer array 
	    for (j=0; j<NPAR; ++j)
	      sgnlv[NPAR*(*sgnlc-1)+j] = (FLOAT_TYPE)sgnlt[j];
#ifdef VERBOSE
	    printf ("\nSignal %d: %d %d %d %d %d snr=%.2f\n", 
		    *sgnlc, pm, mm, nn, ss, Fmax[i], sgnlt[4]);
#endif 
	  }
	}
      } // i
#endif // old/new version
      
#if TIMERS>2
      tend = get_current_time(CLOCK_PROCESS_CPUTIME_ID);
      spindown_timer += get_time_difference(tstart, tend);
      spindown_counter++;
#endif
      
    } // for ss 
  } // omp parallel
  
#ifndef VERBOSE 
  printf("Number of signals found: %d\n", *sgnlc); 
#endif 

  //  tend = get_current_time(CLOCK_REALTIME);
  //time_elapsed = get_time_difference(tstart, tend);
  //printf("Parallel part: %e  ( per thread %e ) s\n", time_elapsed, time_elapsed/omp_get_max_threads());


#if TIMERS>2
  printf("\nTotal spindown loop time: %e s, mean spindown cpu-time: %e s (%d runs)\n",
	 spindown_timer, spindown_timer/spindown_counter, spindown_counter);
#endif

  return 0;
  
} // jobcore
