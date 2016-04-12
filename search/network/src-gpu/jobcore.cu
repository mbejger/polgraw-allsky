#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <malloc.h>
#include <complex.h>

/* JobCore file */
#include "struct.h"
#include "jobcore.h"
#include "auxi.h"
#include "settings.h"
#include "timer.h"
#include "kernels.h"
#include "spline_z.h"
#include "cuda_error.h"
#include <assert.h>

__constant__ double maa_d, mbb_d;  // mean Amplitude modulation functions

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
  FFT_plans *plans,
  FFT_arrays *fft_arr,
  Aux_arrays *aux,
  int *FNum,
  double *F_d) {

  // struct stat buffer;
  struct flock lck;

  int pm, mm, nn;    // hemisphere, sky positions 
  int sgnlc;         // number of candidates
  FLOAT_TYPE *sgnlv;    // array with candidates data

  char outname[512];
  int fd;
  FILE *state;

  //cublas handle for scaling
  cublasHandle_t scale;
  cublasCreate (&scale);

#if TIMERS>0
  struct timespec tstart = get_current_time(), tend;
#endif

  /* copy ampl. mod. coef. to device constant memory */
  copy_amod_coeff(sett->nifo);
  //printf("after copy_amod\n");

  int cand_buffer_count = 0;

  //allocate vector for FStat_gpu
  int nav_blocks = (sett->nmax - sett->nmin)/NAV;     //number of nav-blocks
  int blocks = (sett->nmax - sett->nmin)/BLOCK_SIZE;  //number of blocks for Fstat-smoothing

  CudaSafeCall ( cudaMalloc((void**)&aux->mu_t_d, sizeof(FLOAT_TYPE)*blocks) );
  CudaSafeCall ( cudaMalloc((void**)&aux->mu_d, sizeof(FLOAT_TYPE)*nav_blocks) );


  state=NULL;
  if(opts->checkp_flag) 
    state = fopen (opts->qname, "w");

  /* Loop over hemispheres  */ 
  for (pm=s_range->pst; pm<=s_range->pmr[1]; ++pm) {
    
    sprintf (outname, "%s/triggers_%03d_%03d%s_%d.bin", 
	     opts->prefix, opts->ident, opts->band, opts->label, pm);
    
    /* Two main loops over sky positions */ 

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
			 fft_arr,      // arrays for fftw
			 aux,          // auxiliary arrays
			 F_d,            // F-statistics array
			 &sgnlc,       // reference to array with the parameters
			 // of the candidate signal
			 // (used below to write to the file)
			 FNum,        // Candidate signal number
			 scale         //handle for scaling
			 );
	
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

  cublasDestroy(scale);

#if TIMERS>0
  tend = get_current_time();
  // printf("tstart = %d . %d\ntend = %d . %d\n", tstart.tv_sec, tstart.tv_usec, tend.tv_sec, tend.tv_usec);
  double time_elapsed = get_time_difference(tstart, tend);
  printf("Time elapsed: %e s\n", time_elapsed);
#endif

}


/* Main job   */ 

FLOAT_TYPE* job_core(
		     int pm,                    // Hemisphere
		     int mm,                    // Grid 'sky position'
		     int nn,                    // Second grid 'sky position'
		     Search_settings *sett,     // Search settings
		     Command_line_opts *opts,   // Search options 
		     Search_range *s_range,     // Range for searching
		     FFT_plans *plans,         // Plans for fftw
		     FFT_arrays *fft_arr,     // Arrays for fftw
		     Aux_arrays *aux,           // Auxiliary arrays
		     double *F_d,                 // F-statistics array
		     int *sgnlc,                // Candidate trigger parameters 
		     int *FNum,                // Candidate signal number
		     cublasHandle_t scale      //handle for scaling
		     ) {

  int i, j, n;
  int smin = s_range->sst, smax = s_range->spndr[1];
  double al1, al2, sinalt, cosalt, sindelt, cosdelt, sgnlt[NPAR], 
    nSource[3], het0, sgnl0, ft;
  double _tmp1[sett->nifo][sett->N];
  FLOAT_TYPE *sgnlv; 

  /// temp for testing  
  static double *F;
  if (F == NULL) F = (double *)malloc(2*sett->nfft*sizeof(double));

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
    //printf("before modvir\n");
    modvir_gpu(sinalt, cosalt, sindelt, cosdelt, 
	       sett->N, &ifo[n], aux, n);

    // Calculate detector positions with respect to baricenter
    nSource[0] = cosalt*cosdelt;
    nSource[1] = sinalt*cosdelt;
    nSource[2] = sindelt;

    shft1 = nSource[0]*ifo[n].sig.DetSSB[0]
          + nSource[1]*ifo[n].sig.DetSSB[1]
          + nSource[2]*ifo[n].sig.DetSSB[2];


    tshift_pmod_kern<<<(sett->nfft + BLOCK_SIZE - 1)/BLOCK_SIZE, BLOCK_SIZE>>>
      ( shft1, het0, nSource[0], nSource[1], nSource[2],
	ifo[n].sig.xDat_d, fft_arr->xa_d, fft_arr->xb_d, 
	ifo[n].sig.shft_d, ifo[n].sig.shftf_d, 
	aux->tshift_d,
	ifo[n].sig.aa_d, ifo[n].sig.bb_d,
	ifo[n].sig.DetSSB_d,
	sett->oms, sett->N, sett->nfft, sett->interpftpad );
    CudaCheckError();
    
    cudaDeviceSynchronize();
    
    //fftw_execute(plans->pl_int);  //forward fft (len nfft)
    //fftw_execute(plans->pl_int2); //forward fft (len nfft)
    cufftExecZ2Z(plans->pl_int, fft_arr->xa_d, fft_arr->xa_d, CUFFT_FORWARD);
    cufftExecZ2Z(plans->pl_int, fft_arr->xb_d, fft_arr->xb_d, CUFFT_FORWARD);


    //shift frequencies and remove those over Nyquist
    resample_postfft<<<(sett->Ninterp + BLOCK_SIZE - 1)/BLOCK_SIZE, BLOCK_SIZE>>>
      ( fft_arr->xa_d, fft_arr->xb_d, sett->nfft, sett->Ninterp, nyqst );
    CudaCheckError();

    // Backward fft (len Ninterp = nfft*interpftpad)
    //fftw_execute (plans->pl_inv);     
    //fftw_execute (plans->pl_inv2); 
    cufftExecZ2Z(plans->pl_inv, fft_arr->xa_d, fft_arr->xa_d, CUFFT_INVERSE);
    cufftExecZ2Z(plans->pl_inv, fft_arr->xb_d, fft_arr->xb_d, CUFFT_INVERSE);

    //scale fft with cublas
    ft = (double)sett->interpftpad / sett->Ninterp;
    cublasZdscal( scale, sett->Ninterp, &ft, fft_arr->xa_d, 1);
    cublasZdscal( scale, sett->Ninterp, &ft, fft_arr->xb_d, 1);
    CudaCheckError();

    // Spline interpolation to xDatma, xDatmb arrays
    gpu_interp(fft_arr->xa_d,       //input data
               sett->Ninterp,       //input data length
	       aux->tshift_d,       //output time domain
	       ifo[n].sig.xDatma_d, //output values
	       sett->N,             //output data length
	       aux->diag_d,         //diagonal
	       aux->ldiag_d,        //lower diagonal
	       aux->udiag_d,        //upper diagonal
	       aux->B_d);           //coefficient matrix

    gpu_interp(fft_arr->xb_d,       //input data
               sett->Ninterp,       //input data length
	       aux->tshift_d,       //output time domain
	       ifo[n].sig.xDatmb_d, //output values
	       sett->N,             //output data length
	       aux->diag_d,         //diagonal
	       aux->ldiag_d,        //lower diagonal
	       aux->udiag_d,        //upper diagonal
	       aux->B_d);           //coefficient matrix

    ft = 1./ifo[n].sig.sig2;
    cublasZdscal( scale, sett->N, &ft, ifo[n].sig.xDatma_d, 1);
    cublasZdscal( scale, sett->N, &ft, ifo[n].sig.xDatmb_d, 1);
    CudaCheckError();

  } // end of detector loop 

  double _maa = 0.;
  double _mbb = 0.;
    
  for(n=0; n<sett->nifo; ++n) {
    double aatemp = 0., bbtemp = 0.;
    // square sums of modulation factors
    cublasDdot (scale, sett->N , 
      (const double *)ifo[n].sig.aa_d, 1,
      (const double *)ifo[n].sig.aa_d, 1,
      &aatemp);
    cublasDdot (scale, sett->N , 
      (const double *)ifo[n].sig.bb_d, 1,
      (const double *)ifo[n].sig.bb_d, 1,
      &bbtemp);
    
    /* or use sqr( cublasSnrm2()) */
    _maa += aatemp/ifo[n].sig.sig2;
    _mbb += bbtemp/ifo[n].sig.sig2;
  }
  CudaCheckError();

  //  printf("maa_d=%f", _maa);
  cudaMemcpyToSymbol(maa_d, &_maa, sizeof(double), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(mbb_d, &_mbb, sizeof(double), 0, cudaMemcpyHostToDevice);
  CudaCheckError();  

  
  /* Spindown loop */
  
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
    sgnlt[1] = ss*sett->M[5] + nn*sett->M[9] + mm*sett->M[13];
    
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
      //      printf("%d  %d\n", BLOCK_SIZE, (sett->N + BLOCK_SIZE - 1)/BLOCK_SIZE );
      
      
      phase_mod_1<<<(sett->N + BLOCK_SIZE - 1)/BLOCK_SIZE, BLOCK_SIZE>>>
        ( fft_arr->xa_d, fft_arr->xb_d,
          ifo[0].sig.xDatma_d, ifo[0].sig.xDatmb_d,
          het1, sgnlt[1], ifo[0].sig.shft_d,
          sett->N );
      
      cudaDeviceSynchronize();

      for(n=1; n<sett->nifo; ++n) {
	phase_mod_2<<<(sett->N + BLOCK_SIZE - 1)/BLOCK_SIZE, BLOCK_SIZE>>>
	  ( fft_arr->xa_d, fft_arr->xb_d,
	    ifo[n].sig.xDatma_d, ifo[n].sig.xDatmb_d,
	    het1, sgnlt[1], ifo[n].sig.shft_d,
	    sett->N );
      }

      // initialize arrays to 0. with integer 0
      // assuming double , remember to change when switching to float
      cuMemsetD32Async((CUdeviceptr) (fft_arr->xa_d + sett->N), 0,
		       (sett->nfftf - sett->N)*2*(sizeof(double)/4), NULL);
      cuMemsetD32Async((CUdeviceptr) (fft_arr->xb_d + sett->N), 0,
		       (sett->nfftf - sett->N)*2*(sizeof(double)/4), NULL);
      CudaCheckError();

      //fft length fftpad*nfft
      cufftExecZ2Z(plans->plan, fft_arr->xa_d, fft_arr->xa_d, CUFFT_FORWARD);
      cufftExecZ2Z(plans->plan, fft_arr->xb_d, fft_arr->xb_d, CUFFT_FORWARD);

      (*FNum)++;
      
      compute_Fstat<<<(sett->nmax-sett->nmin + BLOCK_SIZE - 1)/BLOCK_SIZE, BLOCK_SIZE>>>
        ( fft_arr->xa_d + sett->nmin,
          fft_arr->xb_d + sett->nmin,
          F_d + sett->nmin,
          sett->nmax - sett->nmin );
      CudaCheckError();

#define GPUFSTAT_NO
#ifdef GPUFSTAT
      if(!(opts->white_flag))  // if the noise is not white noise
	FStat_gpu(F_d+sett->nmin, sett->nmax - sett->nmin, NAV, aux->mu_d, aux->mu_t_d);
      
#else
      if(!(opts->white_flag))  // if the noise is not white noise
	FStat_gpu_simple(F_d + sett->nmin, sett->nmax - sett->nmin, NAVFSTAT);
#endif
     
      CudaSafeCall ( cudaMemcpy(F, F_d, 2*sett->nfft*sizeof(FLOAT_TYPE), cudaMemcpyDeviceToHost));
      /*      
      FILE *f1 = fopen("fstat-gpu.dat", "w");
      for(i=sett->nmin; i<sett->nmax; i++)
	fprintf(f1, "%d   %lf\n", i, F[i]);

      fclose(f1);
      printf("wrote fstat-gpu.dat | ss=%d  \n", ss);
      //exit(EXIT_SUCCESS);
      */

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


void modvir_gpu(double sinal, double cosal, double sindel, double cosdel,
		int Np, Detector_settings *ifoi, Aux_arrays *aux, int idet){

  double cosalfr, sinalfr, c2d, c2sd, c, s, c2s, cs;

  cosalfr = cosal*(ifoi->sig.cphir) + sinal*(ifoi->sig.sphir);
  sinalfr = sinal*(ifoi->sig.cphir) - cosal*(ifoi->sig.sphir);
  c2d = sqr(cosdel);
  c2sd = sindel*cosdel;

  modvir_kern<<<BLOCK_DIM(Np, BLOCK_SIZE), BLOCK_SIZE>>>
    ( ifoi->sig.aa_d, ifoi->sig.bb_d, cosalfr, sinalfr, c2d, c2sd, 
      aux->sinmodf_d, aux->cosmodf_d, sindel, cosdel, Np, idet );

  return;

}

/* just simple patch - to be replaced */
void FStat_gpu_simple(FLOAT_TYPE *F_d, int nfft, int nav) {
  
  int blocks = nfft/nav;
  fstat_norm_simple<<<blocks, 1>>>(F_d, nav);
  CudaCheckError();

}


#ifdef GPUFSTAT
/* WARNING
   This won't necessarily work for other values than:
   NAV = 4096
   N = nmax - nmin = 507904
   For those given above it works fine.
   Reason is the "reduction depth", i.e. number of required \
   reduction kernel calls.
*/
void FStat_gpu(FLOAT_TYPE *F_d, int N, int nav, FLOAT_TYPE *mu_d, FLOAT_TYPE *mu_t_d) {

  int nav_blocks = N/nav;           //number of blocks
  int nav_threads = nav/BLOCK_SIZE; //number of blocks computing one nav-block
  int blocks = N/BLOCK_SIZE;

  //    CudaSafeCall ( cudaMalloc((void**)&cu_mu_t, sizeof(float)*blocks) );
  //    CudaSafeCall ( cudaMalloc((void**)&cu_mu, sizeof(float)*nav_blocks) );

  //sum fstat in blocks
  reduction_sum<BLOCK_SIZE_RED><<<blocks, BLOCK_SIZE_RED, BLOCK_SIZE_RED*sizeof(FLOAT_TYPE)>>>(F_d, mu_t_d, N);
  CudaCheckError();

  //sum blocks computed above and return 1/mu (number of divisions: blocks), then fstat_norm doesn't divide (potential number of divisions: N)
  reduction_sum<<<nav_blocks, nav_threads, nav_threads*sizeof(FLOAT_TYPE)>>>(mu_t_d, mu_d, blocks);
  CudaCheckError();
  
  //divide by mu/(2*NAV)
  fstat_norm<<<blocks, BLOCK_SIZE>>>(F_d, mu_d, N, nav);
  CudaCheckError();
  
}
#endif
