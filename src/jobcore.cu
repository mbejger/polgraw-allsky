#include <math.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <malloc.h>
#include <complex.h>

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cufft.h>

/* JobCore file */
#include "jobcore.h"
#include "auxi.h"
#include "settings.h"
#include "timer.h"
#include "kernels.h"
#include "spline_z.h"


//__constant__ Detector_settings *cu_sett;
//__constant__ Command_line_opts *cu_opts;
//__constant__ Ampl_mod_coeff *cu_amod;

__constant__ double cu_c[9];

//main searching function (loops inside)
void search(
	    Detector_settings *sett,
	    Command_line_opts *opts,
	    Search_range *s_range,
	    Arrays *arr,
	    FFT_plans *plans,
	    Ampl_mod_coeff *amod,
	    int *FNum,
	    FLOAT_TYPE *cu_F
	    )
{

  struct timespec tstart = get_current_time(), tend;
  double time_elapsed;

  //copy ampl. mod. coef. to constant memory
  copy_amod_coeff(amod);	
	
  int cand_buffer_count = 0;

  //allocate vector for FStat_gpu
  int nav_blocks = (sett->nmax - sett->nmin)/NAV;     //number of nav-blocks
  int blocks = (sett->nmax - sett->nmin)/BLOCK_SIZE;  //number of blocks for Fstat-smoothing

  CudaSafeCall ( cudaMalloc((void**)&arr->cu_mu_t, sizeof(FLOAT_TYPE)*blocks) );
  CudaSafeCall ( cudaMalloc((void**)&arr->cu_mu, sizeof(FLOAT_TYPE)*nav_blocks) );

  //allocate auxiliary vectors for modvir
  int modvir_blocks = BLOCK_DIM(sett->N, BLOCK_SIZE_RED);
  CudaSafeCall ( cudaMalloc((void**)&arr->cu_o_aa, sizeof(double)*modvir_blocks));
  CudaSafeCall ( cudaMalloc((void**)&arr->cu_o_bb, sizeof(double)*modvir_blocks));
  CudaSafeCall ( cudaMalloc((void**)&arr->cu_o_aa2, sizeof(double)*modvir_blocks));
  CudaSafeCall ( cudaMalloc((void**)&arr->cu_o_bb2, sizeof(double)*modvir_blocks));

  int pm, mm, nn;

  char outname[64], qname[64];
  
  FILE *state;
  if(opts->checkp_flag) {
    if(opts->hemi)
      sprintf (qname, "state_%03d_%03d%s_%d.dat", opts->ident, opts->band, opts->label, opts->hemi);
    else
      sprintf (qname, "state_%03d_%03d%s.dat", opts->ident, opts->band, opts->label);

    state = fopen (qname, "w");
  }
  
  for (pm=s_range->pst; pm<=s_range->pmr[1]; pm++) {	// loop over hemispheres
    sprintf (outname, "%s/triggers_%03d_%03d%s_%d.bin", opts->prefix, opts->ident,
	     opts->band, opts->label, pm);
    
    for (mm=s_range->mst; mm<=s_range->mr[1]; mm++) {	 // 2 loops over
      for (nn=s_range->nst; nn<=s_range->nr[1]; nn++) {	 // sky positions

	if(opts->checkp_flag) {
	  fprintf (state, "%d %d %d %d %d\n", pm, mm, nn, s_range->sst, *FNum);
	  fflush (state);
	}

	/* Loop over Spindows here */
	job_core(
		 pm,	              // hemisphere
		 mm,		      // grid 'sky position'
		 nn,		      // other grid 'sky position'
		 sett,		      // detector settings
		 opts,		      // cmd opts
		 s_range,	      // range for searching
		 arr,
		 plans,	              // plans for fftw
		 amod,	              //amplitude modulation functions coefficients
		 cu_F,
		 // of the candidate signal
		 // (used below to write to the file)
		 FNum,		      // Candidate signal number
		 &cand_buffer_count,
		 outname
		 );

	//get back to regular spin-down range
	s_range->sst = s_range->spndr[0];
				
      } // for nn
			
      //if there's non-saved data
      if (cand_buffer_count > 0) {
	save_candidates(arr->cu_cand_buffer, arr->cand_buffer, &cand_buffer_count, outname);
      }
      
      s_range->nst = s_range->nr[0];
    } // for mm
    s_range->mst = s_range->mr[0];
  } // for pm

  if(opts->checkp_flag) {
    fclose (state);
  }

  //cleanup
  cudaFree(arr->cu_mu);
  cudaFree(arr->cu_mu_t);
  CudaSafeCall ( cudaFree(arr->cu_o_aa) );
  CudaSafeCall ( cudaFree(arr->cu_o_bb) );
  CudaSafeCall ( cudaFree(arr->cu_o_aa2) );
  CudaSafeCall ( cudaFree(arr->cu_o_bb2) );

  tend = get_current_time();
  time_elapsed = get_time_difference(tstart, tend);
  printf("Time elapsed for searching: %e s\n", time_elapsed);

}


double* job_core(
		 int pm,		     // hemisphere
		 int mm,		     // grid 'sky position'
		 int nn,		     // other grid 'sky position'
		 Detector_settings *sett,    // detector settings
		 Command_line_opts *opts,    // cmd opts
		 Search_range *s_range,	     // range for searching
		 Arrays *arr,
		 FFT_plans *plans,	     // plans for fftw
		 Ampl_mod_coeff *amod,	     //amplitude modulation functions coefficients
		 FLOAT_TYPE *cu_F,
		 // of the candidate signal
		 // (used below to write to the file)
		 int *FNum,		     // Candidate signal number
		 int *cand_buffer_count,
		 const char* outname
		 )
{
  struct timespec tstart, tend;

  int j;
  int smin = s_range->sst, smax = s_range->spndr[1];
  double al1, al2, sinalt, cosalt, sindelt, cosdelt, sgnlt[NPAR],
    nSource[3], het0, sgnl0, ft;


  //sgnlt includes current parameters

  /* Matrix	M(.,.) (defined on page 22 of PolGrawCWAllSkyReview1.pdf file)
     defines the transformation form integers (bin, ss, nn, mm) determining
     a grid point to linear coordinates omega, omegadot, alpha_1, alpha_2),
     where bin is the frequency bin number and alpha_1 and alpha_2 are
     defined on p. 22 of PolGrawCWAllSkyReview1.pdf file.

     [omega]
     [bin]
     [omegadot]		    = M(.,.) \times [ss]
     [alpha_1/omega]	    [nn]
     [alpha_2/omega]	    [mm]

     Array M[.] is related to matrix M(.,.) in the following way;

     [ M[0] M[4] M[8]	M[12] ]
     M(.,.) =  [ M[1] M[5] M[9]	 M[13] ]
               [ M[2] M[6] M[10] M[14] ]
               [ M[3] M[7] M[11] M[15] ]

     and

     M[1] = M[2] = M[3] = M[6] = M[7] = 0

  */

  //grid positions: alpha1 and alpha2, as described in paper
  al1 = nn*sett->M[10]+mm*sett->M[14];
  al2 = nn*sett->M[11]+mm*sett->M[15];


  /*
    ############ Checking coordinates ################
  */

  // check if the search is in an appropriate region of the grid
  // if not, returns NULL
  if((sqr(al1)+sqr(al2))/sqr(sett->oms) > 1.) return NULL ;


  int ss;
  double shft1;

  //change linear (grid) coordinates to real coordinates
  lin2ast (al1/sett->oms, al2/sett->oms, pm, sett->sepsm, sett->cepsm,
	   &sinalt, &cosalt, &sindelt, &cosdelt);

  // calculate declination and right ascencion
  // which will be written in file as candidate signal sky positions
  sgnlt[2] = asin (sindelt);
  sgnlt[3] = fmod (atan2 (sinalt, cosalt)+2.*M_PI, 2.*M_PI);


  // amplitude modulation function
  modvir_gpu(sinalt, cosalt, sindelt, cosdelt, sett->sphir,
	     sett->cphir, arr->cu_aa, arr->cu_bb, sett->N, arr);

  //calculate detector positions with respect to baricenter
  nSource[0] = cosalt*cosdelt;
  nSource[1] = sinalt*cosdelt;
  nSource[2] = sindelt;

  shft1 = 0.;
  for (j=0; j<3; j++)
    shft1 += nSource[j]*arr->DetSSB[j];
  het0 = fmod (nn*sett->M[8]+mm*sett->M[12], sett->M[0]);

  //compute shifted time, modulate phase and pad zeros to size 'nfft'
  shift_time_mod<<<(sett->nfft + BLOCK_SIZE - 1)/BLOCK_SIZE, BLOCK_SIZE>>>
    ( shft1, het0, nSource[0], nSource[1], nSource[2],
      arr->cu_xDat, arr->cu_xa, arr->cu_xb, arr->cu_shft,
      arr->cu_shftf, arr->cu_tshift, arr->cu_aa, arr->cu_bb, arr->cu_DetSSB,
      sett->oms, sett->N, sett->nfft, sett->interpftpad );
  CudaCheckError();

  /* Resampling */
  /*
    This part is performed to double the sampling rate (get twice more samples)
  */

  cudaDeviceSynchronize();

  int nyqst = (sett->nfft+2)/2; // Nyquist frequency

  //fft
  cufftExecZ2Z(plans->pl_int, arr->cu_xa, arr->cu_xa, CUFFT_FORWARD);
  cufftExecZ2Z(plans->pl_int, arr->cu_xb, arr->cu_xb, CUFFT_FORWARD);

  //shift frequencies and remove those over Nyquist
  resample_postfft<<<(sett->Ninterp + BLOCK_SIZE - 1)/BLOCK_SIZE, BLOCK_SIZE>>>
    ( arr->cu_xa, arr->cu_xb, sett->nfft, sett->Ninterp, nyqst );
  CudaCheckError();

  //fft back
  cufftExecZ2Z(plans->pl_inv, arr->cu_xa, arr->cu_xa, CUFFT_INVERSE);
  cufftExecZ2Z(plans->pl_inv, arr->cu_xb, arr->cu_xb, CUFFT_INVERSE);

  //scale fft
  ft = (double)sett->interpftpad / sett->Ninterp; //scale FFT
  scale_fft<<<(sett->Ninterp + BLOCK_SIZE - 1)/BLOCK_SIZE, BLOCK_SIZE>>>
    ( arr->cu_xa, arr->cu_xb, ft, sett->Ninterp );
  CudaCheckError();

  /* spline interpolation here */

  //interpolation of xa
  gpu_interp(arr->cu_xa,	//input data
	     sett->Ninterp,	//input data length
	     arr->cu_tshift,	//output time domain
	     arr->cu_xar,	//output values
	     sett->N,		//output data length
	     arr->cu_d,		//diagonal
	     arr->cu_dl,	//upper diagonal
	     arr->cu_du,	//lower diagonal
	     arr->cu_B);	//coefficient matrix

  //interpolation of xb
  gpu_interp(arr->cu_xb,	//input data
	     sett->Ninterp,	//input data length
	     arr->cu_tshift,	//output time domain
	     arr->cu_xbr,	//output values
	     sett->N,		//output data length
	     arr->cu_d,		//diagonal
	     arr->cu_dl,	//upper diagonal
	     arr->cu_du,	//lower diagonal
	     arr->cu_B);	//coefficient matrix


  /* here we have arrays cu_xar and cu_xbr with spline-interpolated signal */

  //copy arrays to float (or to double, again...)
  double_to_float<<<(sett->N + BLOCK_SIZE - 1) / BLOCK_SIZE, BLOCK_SIZE>>>
    (arr->cu_xar, arr->cu_xbr,
     arr->cu_xar_f, arr->cu_xbr_f, sett->N);
  
  double spindown_timer = 0;
  int spindown_counter  = 0;

  /*
    ############ Spindown loop ################
  */

  printf ("\n>>%d\t%d\t%d\t[%d..%d]\n", *FNum, mm, nn, smin, smax);

  // if no-spindown
  if (opts->s0_flag) smin = smax;
  // if spindown parameter is taken into account,
  // smin != smax
  for (ss=smin; ss<=smax; ss++) {
		
    // Spindown parameter
    sgnlt[1] = (opts->s0_flag ? 0. : ss * sett->M[5] + nn * sett->M[9] + mm * sett->M[13]);

    if (sgnlt[1] >= -sett->Smax && sgnlt[1] <= 0.) { //look only for negative-valued spindowns
      tstart = get_current_time();

      double  het1;

      // print a 'dot' every new spindown
      // (now it slows down execution)
      putchar('.');
      fflush (stdout);

      het1 = fmod(ss*sett->M[4], sett->M[0]);
      if (het1<0) het1 += sett->M[0];

      sgnl0 = het0+het1;

      // this kernel modulates phase
      // xar -> xa
      phase_mod_2<<<(sett->N + BLOCK_SIZE - 1)/BLOCK_SIZE, BLOCK_SIZE>>>
	( arr->cu_xa_f, arr->cu_xb_f,
	  arr->cu_xar_f, arr->cu_xbr_f,
	  het1, sgnlt[1], arr->cu_shft,
	  sett->N );
      
      // arrays used in different interpolation methods aren't the same
      COMPLEX_TYPE *cu_xa_final, *cu_xb_final;
			
      if (opts->fftinterp == INT) { //interpolation by interbinning
	//pad zeros from N to nfft
		pad_zeros<<<(sett->nfft - sett->N + BLOCK_SIZE - 1)/BLOCK_SIZE, BLOCK_SIZE>>>
	  ( arr->cu_xa_f+sett->N, arr->cu_xb_f+sett->N, sett->nfft - sett->N );

	CudaCheckError();
				
	// FFT length nfft
	CUFFT_EXEC_FFT(plans->plan, arr->cu_xa_f, arr->cu_xa_f, CUFFT_FORWARD);
	CUFFT_EXEC_FFT(plans->plan, arr->cu_xb_f, arr->cu_xb_f, CUFFT_FORWARD);

	// make a gap between every bin
	interbinning_gap<<<(sett->nfft/2 + BLOCK_SIZE - 1)/BLOCK_SIZE, BLOCK_SIZE>>>
	  ( arr->cu_xa_f, arr->cu_xb_f, arr->cu_xa2_f, arr->cu_xb2_f, sett->nfft );
	CudaCheckError();
				
	//interpolate values between every bin
	interbinning_interp<<<( (sett->nfft-2)/2 + BLOCK_SIZE - 1)/BLOCK_SIZE, BLOCK_SIZE>>>
	  ( arr->cu_xa2_f, arr->cu_xb2_f, sett->nfft );
	CudaCheckError();
				
	cu_xa_final = arr->cu_xa2_f;
	cu_xb_final = arr->cu_xb2_f;
				
      } else { //interpolation by FFT (fftpad-times-longer FFT)
	//pad zeros from N to nfft*fftpad
	//pad_zeros<<<(sett->nfft*sett->fftpad - sett->N + BLOCK_SIZE - 1)/BLOCK_SIZE, BLOCK_SIZE>>>
	// ( arr->cu_xa_f+sett->N, arr->cu_xb_f+sett->N, sett->nfft*sett->fftpad - sett->N );
	/*
	pad_zeros1<<<(sett->nfft*sett->fftpad - sett->N + BLOCK_SIZE - 1)/BLOCK_SIZE, BLOCK_SIZE>>>
	  ( arr->cu_xa_f+sett->N, sett->nfft*sett->fftpad - sett->N );
	pad_zeros1<<<(sett->nfft*sett->fftpad - sett->N + BLOCK_SIZE - 1)/BLOCK_SIZE, BLOCK_SIZE>>>
	  ( arr->cu_xb_f+sett->N, sett->nfft*sett->fftpad - sett->N );
	*/

	cuMemsetD32Async((CUdeviceptr) (arr->cu_xa_f+sett->N), 0, 
			 (sett->nfft*sett->fftpad - sett->N)*2, NULL);
	cuMemsetD32Async((CUdeviceptr) (arr->cu_xb_f+sett->N), 0, 
			 (sett->nfft*sett->fftpad - sett->N)*2, NULL);

	CudaCheckError();
				
	//fft length fftpad*nfft
	CUFFT_EXEC_FFT(plans->plan, arr->cu_xa_f, arr->cu_xa_f, CUFFT_FORWARD);
	CUFFT_EXEC_FFT(plans->plan, arr->cu_xb_f, arr->cu_xb_f, CUFFT_FORWARD);
				
	cu_xa_final = arr->cu_xa_f;
	cu_xb_final = arr->cu_xb_f;
      }
			
      (*FNum)++;

      //compute F-stat
      compute_Fstat<<<(sett->nmax-sett->nmin + BLOCK_SIZE - 1)/BLOCK_SIZE, BLOCK_SIZE>>>
	( cu_xa_final + sett->nmin,
	  cu_xb_final + sett->nmin,
	  cu_F + sett->nmin,
	  1/(sett->crf0),
	  sett->nmax - sett->nmin );
      CudaCheckError();

      // normalize F-statistics
      if (sett->sig2 < 0.0) { // when noise is not white-noise 
	FStat_gpu(cu_F+sett->nmin, sett->nmax - sett->nmin, NAV, arr->cu_mu, arr->cu_mu_t);
      } else { // when noise is white-noise
	kernel_norm_Fstat_wn<<<(sett->nmax - sett->nmin + BLOCK_SIZE - 1)/BLOCK_SIZE, BLOCK_SIZE>>>
	  ( cu_F + sett->nmin,
	    1/(sett->sig2),
	    sett->nmax - sett->nmin );
      }
      
      //zero candidates count
      cudaMemset(arr->cu_cand_count, 0, sizeof(int));

      /* FINDING AND SAVING CANDIDATES */
      find_candidates<<<(sett->nmax - sett->nmin + BLOCK_SIZE - 1)/BLOCK_SIZE, BLOCK_SIZE>>>
	( cu_F,				// F-statistic array
	  arr->cu_cand_params,	        // local parameters buffer
	  arr->cu_cand_count,	        // local found candidates count
	  opts->trl,			// threshold value
	  sett->nmin,			// minimum F-statistic array index
	  sett->nmax,			// maximum -||-
	  1/((double)sett->fftpad),			// FFT padding for calculation of freq.
	  1/((double)sett->nfft),			// FFT length
	  sgnl0,			// base frequency parameter
	  sett->nd,			// number of degrees of freedom
	  sgnlt[1],			// spindown
	  sgnlt[2],			// sky 
	  sgnlt[3] );			// positions


/////////////////////////
      
      int cand_count;
      cudaMemcpy(&cand_count, arr->cu_cand_count, sizeof(int), cudaMemcpyDeviceToHost);
			
      if (cand_count>0) { //if something was found
	printf ("\nSome signals found in %d %d %d %d\n", pm, mm, nn, ss);
	cudaMemcpy(arr->cu_cand_buffer + *cand_buffer_count * NPAR,
		   arr->cu_cand_params,
		   sizeof(FLOAT_TYPE) * cand_count*NPAR,
		   cudaMemcpyDeviceToDevice);
	*cand_buffer_count += cand_count;
				
	// check if it's time to save results
	// (place left in buffer is less than one-spindown buffer length)
	if (*cand_buffer_count >= arr->cand_buffer_size - arr->cand_params_size) { 
	  save_candidates(arr->cu_cand_buffer, arr->cand_buffer, cand_buffer_count, outname);
	}
      }

      tend = get_current_time();
      spindown_timer += get_time_difference(tstart, tend);
      spindown_counter++;

    } /* if sgnlt[1] */
  } /* for ss */

  printf("\nTotal spindown loop time: %e s, mean spindown time: %e s (%d runs)\n",
	 spindown_timer, spindown_timer/spindown_counter, spindown_counter);

  return 0;
} /* JobCore() */



void save_candidates(FLOAT_TYPE* cu_cand_buffer, FLOAT_TYPE* cand_buffer, 
		     int *cand_buffer_count, const char* outname) 
{
  CudaSafeCall
    ( cudaMemcpy(cand_buffer, cu_cand_buffer, sizeof(FLOAT_TYPE)*NPAR * (*cand_buffer_count), 
		 cudaMemcpyDeviceToHost)
      );

  struct flock lck;
  int fd;

  if ((fd = open (outname, O_WRONLY|O_CREAT|O_APPEND,
		  S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH)) < 0) {
    perror (outname);
    return;
  }
  lck.l_type = F_WRLCK;
  lck.l_whence = 0;
  lck.l_start = 0L;
  lck.l_len = 0L;
  if (fcntl (fd, F_SETLKW, &lck) < 0) perror ("fcntl()");

  //#mb 
  // For debbuging
  /*
  int jj, kk;
  printf("\nParameters to save:\n");
  for (jj=0; jj<*cand_buffer_count; jj++) {
    for(kk=0; kk<NPAR; kk++) {
      printf("%lf ", cand_buffer[jj*NPAR+kk]);
    }
    printf("\n");
  }
  */ 

  write (fd, (void *)(cand_buffer), *cand_buffer_count*NPAR*sizeof(FLOAT_TYPE));

  if (close (fd) < 0) perror ("close()");
  *cand_buffer_count = 0;
}



//not used
//left for reference
void
modvir (double sinal, double cosal, double sindel, double cosdel, double sphir, double cphir, 
	double *a, double *b, int Np, Ampl_mod_coeff *amod, Arrays *arr)
{
  /* Amplitude modulation functions */
  int t;
  double cosalfr, sinalfr, c2d, c2sd, c, s, c2s, cs, as, bs;

  double c1 = amod->c1,
         c2 = amod->c2,
         c3 = amod->c3,
         c4 = amod->c4,
         c5 = amod->c5,
         c6 = amod->c6,
         c7 = amod->c7,
         c8 = amod->c8,
         c9 = amod->c9;

  cosalfr = cosal*cphir+sinal*sphir;
  sinalfr = sinal*cphir-cosal*sphir;
  c2d = sqr(cosdel);
  c2sd = sindel*cosdel;

  as = bs = 0.;
  for (t=0; t<Np; t++) {
    c = cosalfr * arr->cosmodf[t]+sinalfr * arr->sinmodf[t];
    s = sinalfr * arr->cosmodf[t]-cosalfr * arr->sinmodf[t];
    c2s = 2.*sqr(c);
    cs = c*s;
    /* modulation factors */
    a[t] = c1*(2.-c2d)*c2s+c2*(2.-c2d)*2.*cs+c3*c2sd*c+c4*c2sd*s-c1*(2.-c2d)+c5*c2d;
    b[t] = c6*sindel*c2s+c7*sindel*2.*cs+c8*cosdel*c+c9*cosdel*s-c6*sindel;
    as += sqr(a[t]);
    bs += sqr(b[t]);
  } /* for t */
  as /= Np;
  bs /= Np;
  as = sqrt (as);
  bs = sqrt (bs);
  printf("Sa, Sb: %e %e (CPU)\n", as, bs);
  for (t=0; t<Np; t++) {
    a[t] /= as;
    b[t] /= bs;
  }
	
} /* modvir() */



void modvir_gpu (double sinal, double cosal, double sindel, double cosdel,	
		 double sphir, double cphir, double *cu_a, double *cu_b, int N, Arrays *arr) 
{

  double cosalfr = cosal*cphir+sinal*sphir;
  double sinalfr = sinal*cphir-cosal*sphir;
  double c2d = sqr(cosdel);
  double c2sd = sindel*cosdel;
	
  //compute aa, bb
  compute_modvir<<<BLOCK_DIM(N, BLOCK_SIZE), BLOCK_SIZE>>>
    ( cu_a, cu_b, cosalfr, sinalfr, c2d, c2sd, arr->cu_sinmodf,
      arr->cu_cosmodf, sindel, cosdel, N );
	

  int n = N;
  int blocks = BLOCK_DIM(n,BLOCK_SIZE_RED);

  bool turn = false;
  //normalization
  //first, compute sum of squares and sum in every block
  reduction_sumsq<<<blocks, BLOCK_SIZE_RED, sizeof(double)*2*BLOCK_SIZE_RED>>>
    (cu_a, arr->cu_o_aa, n);
  reduction_sumsq<<<blocks, BLOCK_SIZE_RED, sizeof(double)*2*BLOCK_SIZE_RED>>>
    (cu_b, arr->cu_o_bb, n);
  CudaCheckError();	

  double *ina, *inb, *outa, *outb;
  while(blocks > 1) {
    n = blocks;
    blocks = BLOCK_DIM(blocks,BLOCK_SIZE_RED);
    //reduce partial sums
    if (turn == false) {
      ina = arr->cu_o_aa;
      inb = arr->cu_o_bb;
      outa = arr->cu_o_aa2;
      outb = arr->cu_o_bb2;
    } else { 
      ina = arr->cu_o_aa2;
      inb = arr->cu_o_bb2;
      outa = arr->cu_o_aa;
      outb = arr->cu_o_bb;
    }

    reduction_sum<<<blocks, BLOCK_SIZE_RED, sizeof(double)*2*BLOCK_SIZE_RED>>>
      (ina, outa, n);
    reduction_sum<<<blocks, BLOCK_SIZE_RED, sizeof(double)*2*BLOCK_SIZE_RED>>>
      (inb, outb, n);
    CudaCheckError();

    turn=!turn;
  } 
	
  double s_a, s_b;
  //copy reduction results
  CudaSafeCall ( cudaMemcpy(&s_a, outa, sizeof(double), cudaMemcpyDeviceToHost));
  CudaSafeCall ( cudaMemcpy(&s_b, outb, sizeof(double), cudaMemcpyDeviceToHost));
	
  s_a = sqrt(s_a/N);
  s_b = sqrt(s_b/N);

  //	printf("Sa, Sb: %e %e (GPU)\n", s_a, s_b);
	
  //normalize
  modvir_normalize<<<BLOCK_DIM(N, BLOCK_SIZE), BLOCK_SIZE>>>(cu_a, cu_b, 1/(s_a), 1/(s_b), N);
  CudaCheckError();
	
}


/* WARNING
   This won't necessarily work for other values than:
   NAV = 4096
   N = nmax - nmin = 507904
   For those given above it works fine.
   Reason is the "reduction depth", i.e. number of required \
   reduction kernel calls.
*/
void FStat_gpu(FLOAT_TYPE *cu_F, int N, int nav, float *cu_mu, float *cu_mu_t) {
	
  int nav_blocks = N/nav;           //number of blocks
  int nav_threads = nav/BLOCK_SIZE; //number of blocks computing one nav-block
  int blocks = N/BLOCK_SIZE;
	
  //	CudaSafeCall ( cudaMalloc((void**)&cu_mu_t, sizeof(float)*blocks) );
  //	CudaSafeCall ( cudaMalloc((void**)&cu_mu, sizeof(float)*nav_blocks) );
	
  //sum fstat in blocks
  reduction_sum<<<blocks, BLOCK_SIZE, BLOCK_SIZE*sizeof(float)>>>
    (cu_F, cu_mu_t, N);
  CudaCheckError();
		
  //sum blocks computed above
  reduction_sum<<<nav_blocks, nav_threads, nav_threads*sizeof(float)>>>
    (cu_mu_t, cu_mu, blocks);
  CudaCheckError();
	
  //divide by mu/(2*NAV)
  fstat_norm<<<blocks, BLOCK_SIZE>>>(cu_F, cu_mu, N, nav);
  CudaCheckError();
	
}
