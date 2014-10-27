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


//main searching function (loops inside)
void search(
	    Detector_settings *sett,
	    Command_line_opts *opts,
	    Search_range *s_range,
	    FFTW_arrays *fftw_arr,
	    Signals *sig,
	    FFTW_plans *plans,
	    Aux_arrays *aux,
	    Ampl_mod_coeff *amod,
	    int *FNum,
	    double *F
	    )
{
  struct flock lck;
  //struct stat buffer;

  //	clock_t cstart, cend; //clock
  //	double time_elapsed; //for measuring time

  int pm, mm, nn;
  int sgnlc; //number of canditates
  double *sgnlv; //array with candidates data

  char outname[64], qname[64];
  int fd;
	
  FILE *state;
  if(opts->checkp_flag) {
    if(opts->hemi)
      sprintf (qname, "state_%03d_%03d%s_%d.dat", opts->ident, opts->band, opts->label, opts->hemi);
    else
      sprintf (qname, "state_%03d_%03d%s.dat", opts->ident, opts->band, opts->label);

    state = fopen (qname, "w");
  }

	
  struct timespec tstart = get_current_time(), tend;
		
  for (pm=s_range->pst; pm<=s_range->pmr[1]; pm++) {	// loop over hemispheres
    sprintf (outname, "%s/triggers_%03d_%03d%s_%d.bin", opts->prefix, opts->ident,
	     opts->band, opts->label, pm);

    for (mm=s_range->mst; mm<=s_range->mr[1]; mm++) {	// 2 loops over
      for (nn=s_range->nst; nn<=s_range->nr[1]; nn++) {	// sky positions

	if(opts->checkp_flag) {
	  fprintf (state, "%d %d %d %d %d\n", pm, mm, nn, s_range->sst, *FNum);
	  fflush (state);
	}

	/* Loop over Spindows here */
	sgnlv = job_core(
			 pm,          // hemisphere
			 mm,	      // grid 'sky position'
			 nn,	      // other grid 'sky position'
			 sett,        // detector settings
			 opts,        // cmd opts
			 s_range,     // range for searching
			 sig,         // signals
			 fftw_arr,    // arrays for fftw
			 plans,       // plans for fftw
			 aux, 	      // auxiliary arrays
			 amod,        // amplitude modulation functions coefficients 
			 F,	      // F-statistics array
			 &sgnlc,      // reference to array with the parameters
			              // of the candidate signal
			              // (used below to write to the file)
			 FNum	      // Candidate signal number
			 );
				
	//get back to regular spin-down range
	s_range->sst = s_range->spndr[0];
	/* Add trigger parameters to a file */

	//if any signals found (Fstat>Fc)
	if (sgnlc) {
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

	  /* For debbuging */
	  int jj, kk;
	  printf("\nParameters to save:\n");
	  for (jj=0; jj<sgnlc; jj++) {
	    for(kk=0; kk<NPAR; kk++) {
	      printf("%lf ", sgnlv[jj*NPAR+kk]);
	    }
	    printf("\n");
	  }
	  /* */										
	  write (fd, (void *)(sgnlv), sgnlc*NPAR*sizeof (double));

	  if (close (fd) < 0) perror ("close()");

	} /* if sgnlc */
	free (sgnlv);
      } // for nn
      s_range->nst = s_range->nr[0];
    } // for mm
    s_range->mst = s_range->mr[0];
  } // for pm

  if(opts->checkp_flag) {
    fclose (state);
  }

  tend = get_current_time();
	
  // printf("tstart = %d . %d\ntend = %d . %d\n", tstart.tv_sec, tstart.tv_usec, tend.tv_sec, tend.tv_usec);
	
  double time_elapsed = get_time_difference(tstart, tend);
  printf("Time elapsed: %e s\n", time_elapsed);

}


double* job_core(
		 int pm,		  // hemisphere
		 int mm,		  // grid 'sky position'
		 int nn,		  // other grid 'sky position'
		 Detector_settings *sett, // detector settings
		 Command_line_opts *opts, // cmd opts
		 Search_range *s_range,	  // range for searching
		 Signals *sig,		  // signals
		 FFTW_arrays *fftw_arr,   // arrays for fftw
		 FFTW_plans *plans,       // plans for fftw
		 Aux_arrays *aux, 	  // auxiliary arrays
		 Ampl_mod_coeff *amod,    // amplitude modulation functions coefficients 
		 double *F,		  // F-statistics array
		 int *sgnlc,		  // reference to array with the parameters
		                          // of the candidate signal
		                          // (used below to write to the file)
		 int *FNum		  // Candidate signal number
		 )
{

  int j, i;
  int smin = s_range->sst, smax = s_range->spndr[1];
  double al1, al2, sinalt, cosalt, sindelt, cosdelt, sgnlt[NPAR], \
    nSource[3], het0, sgnl0, *sgnlv, ft;


  //sgnlt includes current parameters
	
  /* Matrix	M(.,.) (defined on page 22 of PolGrawCWAllSkyReview1.pdf file)
     defines the transformation form integers (bin, ss, nn, mm) determining
     a grid point to linear coordinates omega, omegadot, alpha_1, alpha_2),
     where bin is the frequency bin number and alpha_1 and alpha_2 are
     defined on p. 22 of PolGrawCWAllSkyReview1.pdf file.

     [omega]													 [bin]
     [omegadot]				= M(.,.) \times [ss]
     [alpha_1/omega]									 [nn]
     [alpha_2/omega]									 [mm]

     Array M[.] is related to matrix M(.,.) in the following way;

     [ M[0] M[4] M[8]	M[12] ]
     M(.,.) =	    [ M[1] M[5] M[9]	M[13] ]
     [ M[2] M[6] M[10] M[14] ]
     [ M[3] M[7] M[11] M[15] ]

     and

     M[1] = M[2] = M[3] = M[6] = M[7] = 0

  */

  //grid positions
  al1 = nn*sett->M[10]+mm*sett->M[14];
  al2 = nn*sett->M[11]+mm*sett->M[15];

  sgnlv = NULL;
  *sgnlc = 0;

  /*
    ############ Checking coordinates ################
  */

  // check if the search is in an appropriate region of the grid
  // if not, returns NULL
  if((sqr(al1)+sqr(al2))/sqr(sett->oms) > 1.) return NULL ;

  //	printf("%lf %lf\n", al1, al2);

  int ss;
  double shft1, phase, cp, sp;
  complex double exph;


  //	fftw_arr->xDb = fftw_arr->xDa + sett->nfft; //position of 'b' array
  //	fftw_arr->rDb = fftw_arr->rDa + sett->Ninterp; //same
  //	fftw_arr->xb = fftw_arr->xbo =


  //!
  // Switch for different interpolation methods
  // (set while calling JobCore() from JobNAllSky.c)
  //	switch (opts->fftinterp) {
  //	case INT:
  //		fftw_arr->xb = fftw_arr->xa + sett->nfft;
  //		fftw_arr->xbo = fftw_arr->xao + sett->nfft;
  //		break;
  //	case FFT:
  fftw_arr->xb = fftw_arr->xa + fftw_arr->arr_len;// + sett->fftpad * sett->nfft;
  //	} /* case fftinterp */

  //change linear (grid) coordinates to real coordinates
  lin2ast (al1/sett->oms, al2/sett->oms, pm, sett->sepsm, sett->cepsm,
	   &sinalt, &cosalt, &sindelt, &cosdelt);


  // calculate declination and right ascention
  // written in file as candidate signal sky positions
  sgnlt[2] = asin (sindelt);
  sgnlt[3] = fmod (atan2 (sinalt, cosalt)+2.*M_PI, 2.*M_PI);


  /* amplitude modulation functions */
  modvir (sinalt, cosalt, sindelt, cosdelt,
	  sett->sphir, sett->cphir, aux->aa, aux->bb, sett->N, amod, aux);
	
  //calculate detector positions with respect to baricenter
  nSource[0] = cosalt*cosdelt;
  nSource[1] = sinalt*cosdelt;
  nSource[2] = sindelt;
  //some magic
  shft1 = 0.;
  for (j=0; j<3; j++)
    shft1 += nSource[j]*aux->DetSSB[j];
  het0 = fmod (nn*sett->M[8]+mm*sett->M[12], sett->M[0]);

  for (i=0; i<sett->N; i++) {
    aux->shft[i] = 0.;
    for (j=0; j<3; j++)
      aux->shft[i] += nSource[j] * aux->DetSSB[i*3+j];
    aux->shftf[i] = aux->shft[i]-shft1;

    /* Phase modulation */
    phase = het0*i+sett->oms*aux->shft[i];
    cp = cos (phase);
    sp = sin (phase);
    exph = cp - I*sp;
    /* Matched filter */
    sig->xDatma[i] = sig->xDat[i]*aux->aa[i]*exph;
    sig->xDatmb[i] = sig->xDat[i]*aux->bb[i]*exph;
  } /* for i */



  /* Resampling */
  /*
    This part is performed to double the sampling rate (get twice more samples)
  */
  for (i=0; i < sett->N; i++) { //rewrite data
    fftw_arr->xa[i] = sig->xDatma[i];
    fftw_arr->xb[i] = sig->xDatmb[i];
  }

  for (i=sett->N; i<sett->nfft; i++) { //zero padding (filling) to size of nearest power of 2
    fftw_arr->xa[i] = 0.;
    fftw_arr->xb[i] = 0.;
  }
	
  //	save_array(fftw_arr->xa, sett->nfft, "xa0.dat");
  //	save_array(fftw_arr->xb, sett->nfft, "xb0.dat");

  fftw_execute (plans->pl_int);  //forward fft (len nfft)
  fftw_execute (plans->pl_int2); //forward fft (len nfft)

  //	save_array(fftw_arr->xa, sett->nfft, "xa1.dat");
  //	save_array(fftw_arr->xb, sett->nfft, "xb1.dat");



  int nyqst = (sett->nfft+2)/2; // Nyquist frequency
	
  //move frequencies from second half of spectrum; loop length: nfft - nyqst =
  // = nfft - nfft/2 - 1 = nfft/2 - 1
  for (i=nyqst + sett->Ninterp - sett->nfft, j=nyqst; i < sett->Ninterp; i++, j++) {
    fftw_arr->xa[i] = fftw_arr->xa[j];
    fftw_arr->xb[i] = fftw_arr->xb[j];
  }
	
  //zero frequencies higher than nyquist
  for (i=nyqst; i< nyqst + sett->Ninterp - sett->nfft; i++) {
    fftw_arr->xa[i] = 0;
    fftw_arr->xb[i] = 0;
  }
	
  //	save_array(fftw_arr->xa, sett->nfft, "xa2.dat");
  //	save_array(fftw_arr->xb, sett->nfft, "xb2.dat");


  fftw_execute (plans->pl_inv); //backward fft (len Ninterp = nfft  *interpftpad
  fftw_execute (plans->pl_inv2); //backward fft (len Ninterp = nfft  *interpftpad


  ft = (double)sett->interpftpad / sett->Ninterp; //scale FFT
  for (i=0; i < sett->Ninterp; i++) {
    fftw_arr->xa[i] *= ft;
    fftw_arr->xb[i] *= ft;
  }

	

  //	struct timeval tstart = get_current_time(), tend;

  splintpad (fftw_arr->xa, aux->shftf, sett->N, sett->interpftpad, sig->xDatma); //spline interpolation to xDatm(a/b) arrays
  splintpad (fftw_arr->xb, aux->shftf, sett->N, sett->interpftpad, sig->xDatmb);

  //	tend = get_current_time();
	
  //	double time_elapsed = get_time_difference(tstart, tend);
  //	printf("Time elapsed, 2 splines: %e s\n", time_elapsed);
	
  //	save_array(sig->xDatma, sett->N, "xa.dat");


  /*
    ############ Spindown loop ################
  */

  struct timespec tstart, tend;
  double spindown_timer=0;
  int spindown_counter = 0;


  printf ("\n>>%d\t%d\t%d\t[%d..%d]\n", *FNum, mm, nn, smin, smax);

  // if no-spindown
  if (opts->s0_flag) smin = smax;
  // if spindown parameter is taken into account,
  // smin != smax
  for (ss=smin; ss<=smax; ss++) {
    tstart = get_current_time();
    // Spindown parameter
    sgnlt[1] = (opts->s0_flag ? 0. : ss * sett->M[5] + nn * sett->M[9] + mm * sett->M[13]);

    if (sgnlt[1] >= -sett->Smax && sgnlt[1] <= 0.) { //look only for negative-valued spindowns
      int ii;
      double phase2, cosPH, sinPH, Fc, het1;

      //print a 'dot' every new spindown
      printf (".");
      fflush (stdout);

      // ?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
      het1 = fmod(ss*sett->M[4], sett->M[0]);
      if (het1<0) het1 += sett->M[0];

      sgnl0 = het0+het1;

      //			FILE *fphase= fopen("phase.dat", "w");
      // phase modulation before fft
      for (i=0; i < sett->N; i++) {
	phase2 = het1*i + sgnlt[1]*(aux->t2[i]+2.*i*aux->shft[i]);
				
	// fprintf(fphase, "%d %lf %lf %lf %lf %lf\n",i, -phase2, het1, sgnlt[1], aux->shft[i], aux->t2[i]);
	cosPH = cos (phase2);
	sinPH = sin (phase2);
	exph = cosPH - I*sinPH;
	fftw_arr->xa[i] = sig->xDatma[i]*exph;
	fftw_arr->xb[i] = sig->xDatmb[i]*exph;
      } /* for i */
      // fclose(fphase);


      //!
      /*			switch (opts->fftinterp) {*/
      for (i = sett->N; i < sett->fftpad * sett->nfft; i++)
	fftw_arr->xa[i] = fftw_arr->xb[i] = 0.; //pad zeros
			
      //	save_array(fftw_arr->xa, sett->nfft, "xa-prefft.dat");
      //	save_array(fftw_arr->xb, sett->nfft, "xb-prefft.dat");

      //			fftw_execute (plans->plan);
      fftw_execute (plans->plan);
      fftw_execute (plans->plan2);

      //	save_array(fftw_arr->xa, sett->nfft, "xa-postfft.dat");
      //	save_array(fftw_arr->xb, sett->nfft, "xb-postfft.dat");


      (*FNum)++;

      /* Computing F-STATISTICS from Fa, Fb */
      for (i=sett->nmin; i<sett->nmax; i++) {
	F[i] = (sqr (creal(fftw_arr->xa[i])) + sqr (cimag(fftw_arr->xa[i])) +
		sqr (creal(fftw_arr->xb[i])) + sqr (cimag(fftw_arr->xb[i])))/sett->crf0;
      }

      /* Normalize F-statistics */
      if ( sett->sig2 < 0.)	// if the noise is not white noise
	FStat (F + sett->nmin, sett->nmax - sett->nmin, NAV, 0);
      else
	for (i = sett->nmin; i < sett->nmax; i++)
	  F[i] /= sett->sig2; //normalize by variance

      //			save_array_double(F, sett->nfft, "F.dat");



      for ( i = sett->nmin;  i < sett->nmax; i++) {
	if ((Fc = F[i]) > opts->trl) { //if F-stat exceeds trl (critical value)
	  /* Find local maximum for neighboring signals */
	  ii = i;
	  while (++i < sett->nmax && F[i] > opts->trl) {
	    if (F[i] >= Fc) {
	      ii = i;
	      Fc = F[i];
	    } /* if F[i] */
	  } /* while i */
	  // Candidate signal frequency
					
	  sgnlt[0] = 2.*M_PI*ii/((double) sett->fftpad * sett->nfft)+sgnl0;
	  // Signal-to-noise ratio
	  sgnlt[4] = sqrt (2.*(Fc-sett->nd));

	  printf("\n%lf %lf %lf %lf %lf\n", 
		 sgnlt[0], sgnlt[1], sgnlt[2], sgnlt[3], sgnlt[4]);
	  (*sgnlc)++; //increase found number
	  /* Add new parameters to output array */
	  sgnlv = (double *) realloc (sgnlv,			\
				      NPAR*(*sgnlc)*sizeof (double));
	  for (j=0; j<NPAR; j++) // save new parameters
	    sgnlv[NPAR*(*sgnlc-1)+j] = sgnlt[j];

	  printf ("\nSignal %d: %d %d %d %d %d \tsnr=%.2f\n", \
		  *sgnlc, pm, mm, nn, ss, ii, sgnlt[4]);
	} /* if Fc > trl */
      } /* for i */
      tend = get_current_time();
      spindown_timer += get_time_difference(tstart, tend);
      spindown_counter++;

      //if (spindown_counter %10 == 0 && spindown_counter > 0 ) {
      //	printf("\nTotal spindown loop time: %e s, mean spindown time: %e s (%d runs)\n",
      //		spindown_timer, spindown_timer/spindown_counter, spindown_counter);
      //	}

    } /* if sgnlt[1] */
  } /* for ss */

  printf("\nTotal spindown loop time: %e s, mean spindown time: %e s (%d runs)\n",
	 spindown_timer, spindown_timer/spindown_counter, spindown_counter);

  return sgnlv;
} /* JobCore() */








void modvir (double sinal, double cosal, double sindel, double cosdel, 
	     double sphir, double cphir, double *a, double *b, int Np, 
	     Ampl_mod_coeff *amod, Aux_arrays *aux) {
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
  for (t=0; t<Np; t++) { //for every time step
    c = cosalfr * aux->cosmodf[t] + sinalfr * aux->sinmodf[t];
    s = sinalfr * aux->cosmodf[t] - cosalfr * aux->sinmodf[t];
    c2s = 2.*sqr(c);
    cs = c*s;
    /* modulation factors */ // compute values in time t
    a[t] = c1*(2.-c2d)*c2s + c2*(2.-c2d)*2.*cs +
           c3*c2sd*c + c4*c2sd*s - c1*(2.-c2d) + c5*c2d;
    b[t] = c6*sindel*c2s + c7*sindel*2.*cs + 
           c8*cosdel*c + c9*cosdel*s - c6*sindel;
    as += sqr(a[t]); // sum square of values
    bs += sqr(b[t]);
  } /* for t */
  as /= Np;       //scale sum of squares by length
  bs /= Np;
  as = sqrt (as); // take square root of sum squares
  bs = sqrt (bs);
  for (t=0; t<Np; t++) { // normalize a[] and b[]
    a[t] /= as;
    b[t] /= bs;
  }
} /* modvir() */
