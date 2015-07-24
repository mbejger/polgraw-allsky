#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <getopt.h>
#include <gsl/gsl_linalg.h>
#include <time.h>

#include "init.h"
#include "struct.h"
#include "settings.h"
#include "auxi.h"



	/*	Command line options handling 
	 */ 
	
void handle_opts(
  Search_settings *sett, 
  Command_line_opts *opts,
	int argc, 
	char* argv[]) {
	
  opts->hemi=0;
  opts->wd=NULL;
  opts->trl=20;
  opts->fftinterp=INT;
	
  strcpy (opts->prefix, TOSTR(PREFIX));
  strcpy (opts->dtaprefix, TOSTR(DTAPREFIX));
  opts->label[0]  = '\0';
  opts->addsig[0] = '\0';
	
  // Initial value of starting frequency
  // set to a negative quantity. If this is not
  // changed by the command line value
  // fpo is calculated from the band number b.
  sett->fpo = -1;

  opts->help_flag=0;
  opts->white_flag=0;
  opts->s0_flag=0;
  opts->checkp_flag=0;

  static int help_flag=0, white_flag=0, s0_flag=0, checkp_flag=1;

  // Reading arguments 

  while (1) {
    static struct option long_options[] = {
      {"help", no_argument, 			&help_flag, 1},
      {"whitenoise", no_argument, 	&white_flag, 1},
      {"nospindown", no_argument, 	&s0_flag, 1},
      {"nocheckpoint", no_argument, &checkp_flag, 0},
      // frame number
      {"ident", required_argument, 0, 'i'},
      // frequency band number
      {"band", required_argument, 0, 'b'},
      // output directory
      {"output", required_argument, 0, 'o'},
      // input data directory
      {"data", required_argument, 0, 'd'},
      // non-standard label for naming files
      {"label", required_argument, 0, 'l'},
      // Spotlight grid range parameter file
      {"range", required_argument, 0, 'r'},
      // change directory parameter
      {"cwd", required_argument, 0, 'c'},
      // interpolation method
      {"threshold", required_argument, 0, 't'},
      // hemisphere
      {"hemisphere", required_argument, 0, 'h'},
      // fpo value
      {"fpo value", required_argument, 0, 'p'},
      // add signal parameters
      {"addsig", required_argument, 0, 'x'},
      {0, 0, 0, 0}
    };

    if (help_flag) {

      printf("*** Continuous GW search code using the F-statistic ***\n");
      printf("Usage: ./search -[switch1] <value1> -[switch2] <value2> ...\n") ;
      printf("Switches are:\n\n");
      printf("-d	Data directory (default is .)\n");
      printf("-o	Output directory (default is ./candidates)\n");
      printf("-i	Frame number\n");
      printf("-b	Band number\n");
      printf("-l	Custom label for the input and output files\n");
      printf("-r	Spotlight search file with sky and spindown grid points\n");
      printf("-c	Change to directory <dir>\n");
      printf("-t	Threshold for the F-statistic (default is 20)\n");
      printf("-h	Hemisphere (default is 0 - does both)\n");
      printf("-p	fpo (starting frequency) value\n");
      printf("-x	Add signal with parameters from <file>\n\n");

      printf("Also:\n\n");
      printf("--whitenoise	white Gaussian noise assumed\n");
      printf("--nospindown	spindowns neglected\n");
      printf("--nocheckpoint	state file won't be created (no checkpointing)\n");
      printf("--help		This help\n");

      exit (0);
    }

    int option_index = 0;
    int c = getopt_long (argc, argv, "i:b:o:d:l:r:c:t:h:p:x:", long_options, &option_index);
    if (c == -1)
      break;

    switch (c) {
    case 'i':
      opts->ident = atoi (optarg);
      break;
    case 't':
      opts->trl = atof(optarg);
      break;
    case 'h':
      opts->hemi = atof(optarg);
      break;
    case 'b':
      opts->band = atoi (optarg);
      break;
    case 'o':
      strcpy (opts->prefix, optarg);
      break;
    case 'd':
      strcpy (opts->dtaprefix, optarg);
      break;
    case 'l':
      opts->label[0] = '_';
      strcpy (1+opts->label, optarg);
      break;
    case 'r':
      strcpy (opts->spotlight, optarg);
      break;
    case 'c':
      opts->wd = (char *) malloc (1+strlen(optarg));
      strcpy (opts->wd, optarg);
      break;
    case 'p':
      sett->fpo = atof(optarg);
      break;
    case 'x':
      strcpy (opts->addsig, optarg);
      break;
    case '?':
      break;
    default:
      break ;
    } /* switch c */
  } /* while 1 */

  opts->white_flag = white_flag;
  opts->s0_flag = s0_flag;
  opts->checkp_flag = checkp_flag;	
	
  printf ("Data directory is %s\n", opts->dtaprefix);
  printf ("Output directory is %s\n", opts->prefix);
  printf ("Frame number is %d\n", opts->ident);
  printf ("Band number is %d\n", opts->band);

  if (opts->white_flag)
    printf ("Assuming white Gaussian noise\n");  
  printf ("Using fftinterp=FFT (FFT interpolation by zero-padding)\n");
  if(opts->trl!=20)
    printf ("Threshold for the F-statistic is %lf\n", opts->trl);
  if(opts->hemi)
    printf ("Search for hemisphere %d\n", opts->hemi);
  if (opts->s0_flag)
    printf ("Assuming s_1 = 0.\n");
  if (strlen(opts->label))
    printf ("Using '%s' as data label\n", opts->label);

  if (strlen(opts->spotlight))
    printf("Obtaining spotlight grid points from '%s'\n", opts->spotlight);
  else { 
    printf("No spotlight grid range file provided! Exiting...\n");   
    abort(); 
  }

  if (strlen(opts->addsig))
    printf ("Adding signal from '%s'\n", opts->addsig);
  if (opts->wd) {
    printf ("Changing working directory to %s\n", opts->wd);
    if (chdir(opts->wd)) {
      perror (opts->wd);
      abort ();
    }
  }

  // Starting band frequency:
  // fpo_val is optionally read from the command line
  // Its initial value is set to -1
  if(!(sett->fpo >= 0))
    // The usual definition:
    sett->fpo = 100. + 0.96875 * opts->band;

  printf("The reference frequency (fpo) is %f\n", sett->fpo);

} // end of command line options handling 


	/* Generate grid from the M matrix (grid.bin)
	 */ 

void read_grid(
	Search_settings *sett, 
	Command_line_opts *opts) {

  sett->M = (double *) calloc (16, sizeof (double));

  FILE *data;
  char filename[512];
  sprintf (filename, "%s/%03d/grid.bin", opts->dtaprefix, opts->ident);
	if ((data=fopen (filename, "r")) != NULL) {
    fread ((void *)&sett->fftpad, sizeof (int), 1, data);

	printf("Using fftpad from the grid file: %d\n", sett->fftpad); 
	
    // M: vector of 16 components consisting of 4 rows
    // of 4x4 grid-generating matrix
    fread ((void *)sett->M, sizeof (double), 16, data);
    fclose (data);
  } else {
	  perror (filename);
      exit(EXIT_FAILURE);
  }

} // end of read grid 


  /* Array initialization 
	 */ 

void init_arrays(
  Search_settings *sett, 
  Command_line_opts *opts,
  Aux_arrays *aux_arr,
  double** F) {

  int i, status; 
  // Allocates and initializes to zero the data, detector ephemeris
  // and the F-statistic arrays

  char filename[512];
  FILE *data;

  for(i=0; i<sett->nifo; i++) { 

    ifo[i].sig.xDat = (double *) calloc(sett->N, sizeof(double));

    // Input time-domain data handling
    sprintf (filename, "%s/%03d/%s/xdatc_%03d_%03d%s.bin", 
	      opts->dtaprefix, opts->ident, ifo[i].name, 
	      opts->ident, opts->band, opts->label);
 
    if((data = fopen(filename, "r")) != NULL) {
      status = fread((void *)(ifo[i].sig.xDat), 
               sizeof(double), sett->N, data);
      fclose (data);

    } else {
      perror (filename);
      exit(EXIT_FAILURE); 
    }

    int j, Nzeros=0;
    // Checking for null values in the data
    for(j=0; j < sett->N; j++)
      if(!ifo[i].sig.xDat[j]) Nzeros++;

    ifo[i].sig.Nzeros = Nzeros; 

    // factor N/(N - Nzeros) to account for null values in the data
    ifo[i].sig.crf0 = (double)sett->N/(sett->N - ifo[i].sig.Nzeros);

    // Estimation of the variance for each detector 
    ifo[i].sig.sig2 = (ifo[i].sig.crf0)*var(ifo[i].sig.xDat, sett->N);

    ifo[i].sig.DetSSB = (double *) calloc(3*sett->N, sizeof(double));

    // Ephemeris file handling
    sprintf (filename, "%s/%03d/%s/DetSSB.bin", 
        opts->dtaprefix, opts->ident, ifo[i].name);

    if((data = fopen(filename, "r")) != NULL) {
      // Detector position w.r.t Solar System Baricenter
      // for every datapoint
      status = fread((void *)(ifo[i].sig.DetSSB), 
               sizeof(double), 3*sett->N, data);

      // Deterministic phase defining the position of the Earth
      // in its diurnal motion at t=0 
      status = fread((void *)(&ifo[i].sig.phir), 
               sizeof(double), 1, data);

      // Earth's axis inclination to the ecliptic at t=0
      status = fread((void *)(&ifo[i].sig.epsm), 
               sizeof(double), 1, data);
      fclose (data);

    } else {
      perror (filename);
      return ;
    }

    // sincos 
    ifo[i].sig.sphir = sin(ifo[i].sig.phir);
    ifo[i].sig.cphir = cos(ifo[i].sig.phir);
    ifo[i].sig.sepsm = sin(ifo[i].sig.epsm);
    ifo[i].sig.cepsm = cos(ifo[i].sig.epsm);

    sett->sepsm = ifo[i].sig.sepsm; 
    sett->cepsm = ifo[i].sig.cepsm; 

    ifo[i].sig.xDatma = 
      (complex double *) calloc(sett->N, sizeof(complex double));
    ifo[i].sig.xDatmb = 
      (complex double *) calloc(sett->N, sizeof(complex double));

    ifo[i].sig.aa = (double *) calloc(sett->N, sizeof(double));
    ifo[i].sig.bb = (double *) calloc(sett->N, sizeof(double));

    ifo[i].sig.shft = (double *) calloc(sett->N, sizeof(double));
    ifo[i].sig.shftf = (double *) calloc(sett->N, sizeof(double));
 
  } // end loop for detectors 

  // Check if the ephemerids have the same epsm parameter
  for(i=1; i<sett->nifo; i++) {  
    if(!(ifo[i-1].sig.sepsm == ifo[i].sig.sepsm)) { 
      printf("The parameter epsm (DetSSB.bin) differs for detectors %s and %s. Aborting...\n", ifo[i-1].name, ifo[i].name); 
      exit(EXIT_FAILURE);

    } 

  } 

  // if all is well with epsm, take the first value 
  sett->sepsm = ifo[0].sig.sepsm;
  sett->cepsm = ifo[0].sig.cepsm;

  *F = (double *) calloc(2*sett->nfft, sizeof(double));

  // Auxiliary arrays, Earth's rotation
  aux_arr->t2 = (double *) calloc(sett->N, sizeof (double));
  aux_arr->cosmodf = (double *) calloc(sett->N, sizeof (double));
  aux_arr->sinmodf = (double *) calloc(sett->N, sizeof (double));
  double omrt;

  for (i=0; i<sett->N; i++) {
    omrt = (sett->omr)*i;     // Earth angular velocity * dt * i
    aux_arr->t2[i] = sqr((double)i);
    aux_arr->cosmodf[i] = cos(omrt);
    aux_arr->sinmodf[i] = sin(omrt);

  }

} // end of init arrays 


  /* Add signal to data
   */ 

//#mb not used right now 
/*
void add_signal(
  Search_settings *sett,
  Command_line_opts *opts,
  Aux_arrays *aux_arr,
  Search_range *s_range) {

  int i, j, n, gsize; 
  double h0, cof; 
  double sinaadd, cosaadd, sindadd, cosdadd, phaseadd, shiftadd, signadd; 
  double nSource[3], sgnlo[10], sgnlol[4];

  FILE *data;

  // Signal parameters are read
  if ((data=fopen (opts->addsig, "r")) != NULL) {
	
    fscanf (data, "%le %d %d", &h0, &gsize, s_range->pmr);     
	  for(i=0; i<10; i++)
			fscanf(data, "%le",i+sgnlo); 	
		fclose (data);
                    
    } else {
      perror (opts->addsig);
    }
  
 		// VSR1 search-specific parametrization of freq. 
		// for the software injection
		// snglo[0]: frequency, sgnlo[1]: frequency. derivative  
	  //#mb fixme
	  //sgnlo[0] += - 2.*sgnlo[1]*(sett->N)*(68 - opts->ident); 

    cof = sett->oms + sgnlo[0]; 
        			  
    for(i=0; i<2; i++) sgnlol[i] = sgnlo[i]; 
	  
    sgnlol[2] = sgnlo[8]*cof; 
	  sgnlol[3] = sgnlo[9]*cof;  
		 	
		// solving a linear system in order to translate 
		// sky position, frequency and spindown (sgnlo parameters) 
		// into the position in the grid
		 
		double *MM ; 
		MM = (double *) calloc (16, sizeof (double));
		for(i=0; i<16; i++) MM[i] = sett->M[i] ;
		
		gsl_vector *x = gsl_vector_alloc (4);     
		int s;
		
		gsl_matrix_view m = gsl_matrix_view_array (MM, 4, 4);
		gsl_matrix_transpose (&m.matrix) ; 
		gsl_vector_view b = gsl_vector_view_array (sgnlol, 4);
		gsl_permutation *p = gsl_permutation_alloc (4);
     
		gsl_linalg_LU_decomp (&m.matrix, p, &s);
		gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
     
		s_range->spndr[0] = round(gsl_vector_get(x,1)); 
		s_range->nr[0] 	= round(gsl_vector_get(x,2));
		s_range->mr[0] 	= round(gsl_vector_get(x,3));
       
		gsl_permutation_free (p);
		gsl_vector_free (x);
		free (MM);
       
		// Define the grid range in which the signal will be looked for
		s_range->spndr[1] = s_range->spndr[0] + gsize; 
    s_range->spndr[0] -= gsize;
		s_range->nr[1] = s_range->nr[0] + gsize; 
    s_range->nr[0] -= gsize;
		s_range->mr[1] = s_range->mr[0] + gsize; 
    s_range->mr[0] -= gsize;
		s_range->pmr[1] = s_range->pmr[0]; 
       	  
		// sgnlo[2]: declination, snglo[3]: right ascension 
		sindadd = sin(sgnlo[2]); 
		cosdadd = cos(sgnlo[2]); 
		sinaadd = sin(sgnlo[3]);  
		cosaadd = cos(sgnlo[3]); 
		
    // Loop for each detector 
    for(n=0; n<sett->nifo; n++) {

      modvir(sinaadd, cosaadd, sindadd, cosdadd,
             sett->N, &ifo[n], aux_arr);

      // Normalization of the modulation amplitudes 
      double as = 0, bs = 0;
      for (i=0; i<sett->N; i++) {
        as += sqr(ifo[n].sig.aa[i]); 
        bs += sqr(ifo[n].sig.bb[i]);
      }

      as /= sett->N; bs /= sett->N;
      as = sqrt (as); bs = sqrt (bs);
 
      for (i=0; i<sett->N; i++) {
        ifo[n].sig.aa[i] /= as;
        ifo[n].sig.bb[i] /= bs; 
      }

		  nSource[0] = cosaadd*cosdadd;
		  nSource[1] = sinaadd*cosdadd;
		  nSource[2] = sindadd;
								
      // adding signal to data (point by point)  								
		  for (i=0; i<sett->N; i++) {

			  shiftadd = 0.; 					 
			  for (j=0; j<3; j++)
				  shiftadd += nSource[j]*ifo[n].sig.DetSSB[i*3+j];		 
					 
			  phaseadd = sgnlo[0]*i + sgnlo[1]*aux_arr->t2[i]  
				  	 + (sett->oms + sgnlo[0] + 2.*sgnlo[1]*i)*shiftadd;

			  signadd = sgnlo[4]*(ifo[n].sig.aa[i])*cos(phaseadd) 
				        + sgnlo[6]*(ifo[n].sig.aa[i])*sin(phaseadd) 
				        + sgnlo[5]*(ifo[n].sig.bb[i])*cos(phaseadd) 
				        + sgnlo[7]*(ifo[n].sig.bb[i])*sin(phaseadd);
						
			  if(ifo[n].sig.xDat[i]) { 
				  ifo[n].sig.xDat[i] += h0*signadd;
//				  thsnr   += pow(signadd, 2.);
			  }	 
		  }
  
    }



} 
*/

	/* Search range 
	 */ 

void set_search_range(
	Search_settings *sett, 
	Command_line_opts *opts, 
	Search_range *s_range) { 

  FILE *data;

  if(strlen(opts->spotlight)) {

    if ((data=fopen (opts->spotlight, "r")) != NULL) {

      int skypos, i; 
      int aqq = fscanf(data, "%d %d", 
                  &s_range->spotlight_pm, &s_range->spotlight_skypos); 
  
      if(aqq != 2) { 
        printf("Problem with the spotlight range file %s. Exiting...\n", opts->spotlight); 
        exit(EXIT_FAILURE);
      }  

      s_range->spotlight_mm   = (int *)calloc(s_range->spotlight_skypos, sizeof(int));
      s_range->spotlight_nn   = (int *)calloc(s_range->spotlight_skypos, sizeof(int));
      s_range->spotlight_noss = (int *)calloc(s_range->spotlight_skypos+1, sizeof(int));

      s_range->spotlight_ss   = (int *)calloc(s_range->spotlight_skypos*MAX_SPOTLIGHT, sizeof(int));

      skypos=0;       
      while(aqq != EOF) { 

        aqq = fscanf(data, "%d %d %d", 
          s_range->spotlight_mm+skypos, s_range->spotlight_nn+skypos, s_range->spotlight_noss+skypos);  
 
        for(i=0; i<s_range->spotlight_noss[skypos]; i++) 
          aqq = fscanf(data, "%d", s_range->spotlight_ss+MAX_SPOTLIGHT*skypos+i); 

        skypos++; 

      }  

//#mb testing printout 
/* 
      printf("Hemisphere: %d\n", s_range->spotlight_pm); 
      printf("No. of sky positions: %d\n", s_range->spotlight_skypos); 

      for(skypos=0; skypos<s_range->spotlight_skypos; skypos++) { 

        printf("%d %d %d\n", 
          s_range->spotlight_mm[skypos], s_range->spotlight_nn[skypos], s_range->spotlight_noss[skypos]); 

        for(i=0; i<s_range->spotlight_noss[skypos]; i++) 
          printf("%d ", s_range->spotlight_ss[MAX_SPOTLIGHT*skypos+i]); 

        printf("\n"); 

      }  
*/

      fclose (data);

    } else { 
      perror (opts->spotlight);
      exit(EXIT_FAILURE);
    }  

  } 

} // end of set search range



  /* FFT Plans 
	 */

void plan_fftw(
  Search_settings *sett, 
	Command_line_opts *opts,
	FFTW_plans *plans, 
	FFTW_arrays *fftw_arr, 
	Aux_arrays *aux_arr) {

  char hostname[512], wfilename[512];
  FILE *wisdom;

  /* Imports a "wisdom file" containing information 
   * (previous tests) about how to optimally compute Fourier 
   * transforms on a given machine. If wisdom file is not present, 
   * it will be created after the test (measure) runs 
   * of the fft_plans are performed below 
   * (see http://www.fftw.org/fftw3_doc/Wisdom.html)
   */ 

  gethostname(hostname, 512);
  sprintf (wfilename, "wisdom-%s.dat", hostname);
  if((wisdom = fopen (wfilename, "r")) != NULL) {
    fftw_import_wisdom_from_file(wisdom);
    fclose (wisdom);
  }

  sett->Ninterp = sett->interpftpad*sett->nfft; 

  // array length (xa, xb) is max{fftpad*nfft, Ninterp}
  fftw_arr->arr_len = (sett->fftpad*sett->nfft > sett->Ninterp 
                    ? sett->fftpad*sett->nfft : sett->Ninterp);

  fftw_arr->xa = fftw_malloc(2*fftw_arr->arr_len*sizeof(fftw_complex));
  fftw_arr->xb = fftw_arr->xa + fftw_arr->arr_len;

  sett->nfftf = sett->fftpad*sett->nfft;

  plans->plan = fftw_plan_dft_1d(sett->nfftf, fftw_arr->xa, fftw_arr->xa, FFTW_FORWARD, FFTW_MEASURE);
  plans->plan2 = fftw_plan_dft_1d(sett->nfftf, fftw_arr->xb, fftw_arr->xb, FFTW_FORWARD, FFTW_MEASURE);
	                             
  plans->pl_int = fftw_plan_dft_1d(sett->nfft, fftw_arr->xa, fftw_arr->xa, FFTW_FORWARD, FFTW_MEASURE);
  plans->pl_int2 = fftw_plan_dft_1d(sett->nfft, fftw_arr->xb, fftw_arr->xb, FFTW_FORWARD, FFTW_MEASURE);
	                             
  plans->pl_inv = fftw_plan_dft_1d(sett->Ninterp, fftw_arr->xa, fftw_arr->xa, FFTW_BACKWARD, FFTW_MEASURE);
  plans->pl_inv2 = fftw_plan_dft_1d(sett->Ninterp, fftw_arr->xb, fftw_arr->xb, FFTW_BACKWARD, FFTW_MEASURE);
	                             
  // Generates a wisdom FFT file if there is none
  if((wisdom = fopen(wfilename, "r")) == NULL) {
    wisdom = fopen(wfilename, "w");
    fftw_export_wisdom_to_file(wisdom);
  }

  fclose (wisdom);

} // end of FFT plans 


  /* Checkpointing
	 */

/*
void read_checkpoints(
	Command_line_opts *opts, 
  Search_range *s_range, 
	int *FNum) {

  if(opts->checkp_flag) {
		
    // filename of checkpoint state file, depending on the hemisphere
    if(opts->hemi)
      sprintf(opts->qname, "state_%03d_%03d%s_%d.dat",  
	            opts->ident, opts->band, opts->label, opts->hemi);
    else
      sprintf(opts->qname, "state_%03d_%03d%s.dat", 
	            opts->ident, opts->band, opts->label);

    FILE *state;
    if((state = fopen(opts->qname, "r")) != NULL) {

      // Scan the state file to get last recorded parameters
      if((fscanf(state, "%d %d %d %d %d", &s_range->pst, &s_range->mst,
		      &s_range->nst, &s_range->sst, FNum)) == EOF) {

        // This means that state file is empty (=end of the calculations)
		    fprintf (stderr, "State file empty: nothing to do...\n");
		    fclose (state);
		    return;

      }

      fclose (state);

    // No state file - start from the beginning
    } else {
      s_range->pst = s_range->pmr[0];
      s_range->mst = s_range->mr[0];
      s_range->nst = s_range->nr[0];
      s_range->sst = s_range->spndr[0];
      *FNum = 0;
    } // if state

  } else {
    s_range->pst = s_range->pmr[0];
    s_range->mst = s_range->mr[0];
    s_range->nst = s_range->nr[0];
    s_range->sst = s_range->spndr[0];
    *FNum = 0;
  } // if checkp_flag

} // end reading checkpoints
*/

  /* Cleanup & memory free 
	 */

void cleanup(
	Search_settings *sett,
	Command_line_opts *opts,
	Search_range *s_range,
	FFTW_plans *plans,
	FFTW_arrays *fftw_arr,
	Aux_arrays *aux,
	double *F) {

  int i; 

  for(i=0; i<sett->nifo; i++) {
    free(ifo[i].sig.xDat);
    free(ifo[i].sig.xDatma);
    free(ifo[i].sig.DetSSB);
    free(ifo[i].sig.aa);
    free(ifo[i].sig.bb);
    free(ifo[i].sig.shftf);
    free(ifo[i].sig.shft);
  } 
	
  free(aux->sinmodf);
  free(aux->cosmodf);
  free(aux->t2);
  free(F);
	
  fftw_free(fftw_arr->xa);
	
  free(sett->M);
	
  fftw_destroy_plan(plans->plan);
  fftw_destroy_plan(plans->pl_int);
  fftw_destroy_plan(plans->pl_inv);

  //spotlight free 
  free(s_range->spotlight_mm); 
  free(s_range->spotlight_nn);
  free(s_range->spotlight_noss); 
  free(s_range->spotlight_ss); 

} // end of cleanup & memory free 
