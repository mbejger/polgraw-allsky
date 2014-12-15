#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <getopt.h>
#include <time.h>

#include "init.h"
#include "struct.h"
#include "settings.h"
#include "auxi.h"
#include "spline_z.h"

#include "cuda_error.h"


/*
 * Command line options
 */
void handle_opts(int argc, char* argv[], Command_line_opts *opts, Detector_settings *sett) {
	
  opts->hemi=0;
  opts->wd=NULL;
  opts->trl=20;
  opts->fftinterp=INT;
	
  strcpy (opts->prefix, TOSTR(PREFIX));
  strcpy (opts->dtaprefix, TOSTR(DTAPREFIX));
  opts->label[0] = '\0';
  opts->range[0] = '\0';
	
  // Initial value of starting frequency
  // set to a negative quantity. If this is not
  // changed by the command line value
  // fpo is calculated from the band number b.
  sett->fpo = -1;

  // Default choice of the detector is Virgo:
  strcpy(opts->ifo_choice, "V1");



  opts->help_flag=0;
  opts->white_flag=0;
  opts->s0_flag=0;
  opts->checkp_flag=0;

  /*
    ############ Reading arguments ################
  */
  static int help_flag=0, white_flag=0, s0_flag=0, checkp_flag=1;


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
      // narrower grid range parameter file
      {"range", required_argument, 0, 'r'},
      // change directory parameter
      {"cwd", required_argument, 0, 'c'},
      // interpolation method
      {"int/fft", required_argument, 0, 'f'},
      // interpolation method
      {"threshold", required_argument, 0, 't'},
      // hemisphere
      {"hemisphere", required_argument, 0, 'h'},
      // detector
      {"detector", required_argument, 0, 'a'},
      // fpo value
      {"fpo value", required_argument, 0, 'p'},
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
      printf("-r	File with grid range or pulsar position\n");
      printf("-c	Change to directory <dir>\n");
      printf("-f	Intepolation method (INT [default] or FFT)\n");
      printf("-t	Threshold for the F-statistic (default is 20)\n");
      printf("-h	Hemisphere (default is 0 - does both)\n");
      printf("-a	Detector (L1, H1 or V1); default is V1\n");
      printf("-p	fpo (starting frequency) value\n\n");
      printf("Also:\n\n");
      printf("--whitenoise	white Gaussian noise assumed\n");
      printf("--nospindown	spindowns neglected\n");
      printf("--nocheckpoint	state file won't be created (no checkpointing)\n");
      printf("--help		This help\n");

      exit (0);
    }

    int option_index = 0;
    int c = getopt_long (argc, argv, "i:b:o:d:l:r:c:t:f:h:a:p:", long_options, &option_index);
    if (c == -1)
      break;

    switch (c) {
    case 'i':
      opts->ident = atoi (optarg);
      break;
    case 'f':
      if(!strcmp(optarg, "FFT")) opts->fftinterp=FFT;
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
      strcpy (opts->range, optarg);
      break;
    case 'c':
      opts->wd = (char *) malloc (1+strlen(optarg));
      strcpy (opts->wd, optarg);
      break;
    case 'p':
      sett->fpo = atof(optarg);
      break;
    case 'a':
      strcpy (opts->ifo_choice, optarg);
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
  if (opts->fftinterp==INT)
    printf ("Using fftinterp=INT (FFT interpolation by interbinning)\n");
  else
    printf ("Using fftinterp=FFT (FFT interpolation by zero-padding)\n");
  if(opts->trl!=20)
    printf ("Threshold for the F-statistic is %f\n", opts->trl);
  if(opts->hemi)
    printf ("Search for hemisphere %d\n", opts->hemi);
  if (opts->s0_flag)
    printf ("Assuming s_1 = 0.\n");
  if (strlen(opts->label))
    printf ("Using '%s' as data label\n", opts->label);
  if (strlen(opts->range))
    printf ("Obtaining grid range from '%s'\n", opts->range);
  if (opts->wd) {
    printf ("Changing working directory to %s\n", opts->wd);
    if (chdir(opts->wd)) {
      perror (opts->wd);
      abort ();
    }
  }

  /*
    ############ End reading arguments ################
  */


  // Starting band frequency:
  // fpo_val is optionally read from the command line
  // Its initial value is set to -1
  if(!(sett->fpo >= 0))
    // The usual definition:
    sett->fpo = 100. + 0.96875 * opts->band;

  printf("The reference frequency (fpo) is %f\n", sett->fpo);

}





void read_grid(Detector_settings *sett, Command_line_opts *opts)
{

  /*
    ############ Grid-generating matrix ################
  */

  sett->M = (double *) calloc (16, sizeof (double));

  FILE *data;
  char filename[CHAR_BUFFER_SIZE];
  sprintf (filename, "%s/%03d/grid.bin", opts->dtaprefix, opts->ident);

  if ((data=fopen (filename, "r")) != NULL) {
    // fftpad: used to zero padding to fftpad*nfft data points
    fread ((void *)&sett->fftpad, sizeof (int), 1, data);
    // M: vector of 16 components consisting of 4 rows
    // of 4x4 grid-generating matrix
    fread ((void *)sett->M, sizeof (double), 16, data);
    fclose (data);
  } else {
    perror (filename);
    printf("Problem with %s... Exiting...\n", filename);
    exit(1); 
  }

}



void init_arrays(Arrays *arr, FLOAT_TYPE** cu_F, 
		 Command_line_opts *opts, Detector_settings *sett) 
{

  // Allocates and initializes to zero the data, detector ephemeris
  // and the F-statistic arrays
  arr->xDat = (double *) calloc (sett->N, sizeof (double));
  CudaSafeCall ( cudaMalloc((void**)&arr->cu_xDat, sizeof(double)*sett->N));

  arr->DetSSB = (double *) calloc (3*sett->N, sizeof (double));
  CudaSafeCall ( cudaMalloc((void**)&arr->cu_DetSSB, sizeof(double)*3*sett->N));

  CudaSafeCall ( cudaMalloc((void**)cu_F, sizeof(FLOAT_TYPE)*sett->fftpad*sett->nfft));
  CudaSafeCall ( cudaMemset(*cu_F, 0, sizeof(FLOAT_TYPE)*sett->fftpad*sett->nfft));

  char filename[CHAR_BUFFER_SIZE];
  FILE *data;
  // Input time-domain data handling
  sprintf (filename, "%s/%03d/xdat_%03d_%03d%s.bin", opts->dtaprefix, opts->ident, \
	   opts->ident, opts->band, opts->label);
  if ((data = fopen (filename, "r")) != NULL) {
    fread ((void *)(arr->xDat), sizeof (double), sett->N, data); // !!! wczytanie danych
    fclose (data);
  } else {
    perror (filename);
    printf("Problem with %s... Exiting...\n", filename);
    exit(1); 
  }
  //copy to device
  CudaSafeCall ( cudaMemcpy(arr->cu_xDat, arr->xDat, sizeof(double)*sett->N, cudaMemcpyHostToDevice));


  int Nzeros=0;
  int i;
  // Checking for null values in the data
  for(i=0; i < sett->N; i++)
    if(!arr->xDat[i]) Nzeros++;

  // factor N/(N - Nzeros) to account for null values in the data
  sett->crf0 = (double)sett->N/(sett->N-Nzeros);


  //if white noise...
  if (opts->white_flag)
    sett->sig2 = sett->N*var (arr->xDat, sett->N);
  else
    sett->sig2 = -1.;

  double epsm, phir;

  /*
    ############ Efemerydy ################
  */

  // Ephemeris file handling
  sprintf (filename, "%s/%03d/DetSSB.bin", opts->dtaprefix, opts->ident);
  if ((data = fopen (filename, "r")) != NULL) {
    // Detector position w.r.t solar system baricenter
    // for every datapoint
    fread ((void *)(arr->DetSSB), sizeof (double), 3*sett->N, data);
    // Deterministic phase defining the position of the Earth
    // in its diurnal motion at t=0
    fread ((void *)(&phir), sizeof (double), 1, data);
    // Earth's axis inclination to the ecliptic at t=0
    fread ((void *)(&epsm), sizeof (double), 1, data);
    fclose (data);
  } else {
    perror (filename);
    printf("Problem with %s... Exiting...\n", filename);
    exit(1);
  }

  //copy DetSSB to device
  CudaSafeCall ( cudaMemcpy(arr->cu_DetSSB, arr->DetSSB, sizeof(double)*sett->N*3, cudaMemcpyHostToDevice));


  /*
    ############ Sincos ################
  */


  sett->sphir = sin (phir);
  sett->cphir = cos (phir);
  sett->sepsm = sin (epsm);
  sett->cepsm = cos (epsm);

  //misc. arrays
  arr->aa = (double*) malloc(sizeof(double)*sett->N);
  arr->bb = (double*) malloc(sizeof(double)*sett->N);
  CudaSafeCall ( cudaMalloc((void**)&arr->cu_aa, sizeof(double)*sett->nfft));
  CudaSafeCall ( cudaMalloc((void**)&arr->cu_bb, sizeof(double)*sett->nfft));

  CudaSafeCall ( cudaMalloc((void**)&arr->cu_shft, sizeof(double)*sett->N));
  CudaSafeCall ( cudaMalloc((void**)&arr->cu_shftf, sizeof(double)*sett->N));
  CudaSafeCall ( cudaMalloc((void**)&arr->cu_tshift, sizeof(double)*sett->N));
	
  //for splines
  init_spline_matrices(&arr->cu_d,
		       &arr->cu_dl,
		       &arr->cu_du,
		       &arr->cu_B,
		       sett->Ninterp);

  arr->cand_params_size = (sett->nmax - sett->nmin);
  arr->cand_buffer_size = (sett->nmax - sett->nmin)*CANDIDATE_BUFFER_SCALE;

  //parameters of found signal

  CudaSafeCall (cudaMalloc((void**)&arr->cu_cand_params, sizeof(FLOAT_TYPE)*arr->cand_params_size));
  CudaSafeCall (cudaMalloc((void**)&arr->cu_cand_buffer, sizeof(FLOAT_TYPE)*arr->cand_buffer_size));
  CudaSafeCall (cudaMalloc((void**)&arr->cu_cand_count, sizeof(int)));
	
  arr->cand_buffer = (FLOAT_TYPE*)malloc(sizeof(FLOAT_TYPE)*arr->cand_buffer_size);

}



void set_search_range(Search_range *s_range, Command_line_opts *opts, Detector_settings *sett) 
{

  // Hemispheres (with respect to the ecliptic)
  if(opts->hemi) {
    s_range->pmr[0] = opts->hemi;
    s_range->pmr[1] = opts->hemi;

  } else {
    s_range->pmr[0] = 1;
    s_range->pmr[1] = 2;
  }

  /*
    ############ Setting search range ################
  */

  FILE *data;

  // If the parameter range is invoked, the search is performed
  // within the range of grid parameters from an ascii file
  // ("-r range_file" from the command line)
  if (strlen (opts->range)) {

    if ((data=fopen (opts->range, "r")) != NULL) {

      int foo = fscanf (data, "%d %d %d %d %d %d %d %d",
			s_range->spndr, 1+s_range->spndr, s_range->nr,
			1+s_range->nr, s_range->mr, 1+s_range->mr,
			s_range->pmr, 1+s_range->pmr);

      // this supresses warnings... :)
      if (foo) {
			
      }
      /*
      // the case when range file does not contain 8 integers
      // describing the grid ranges, but other values:
      // the pulsar position, frequency, and spindowns.
      if(range_status!=8) {

      rewind(data);
      range_status = fscanf (data, "%le %le %le %le %le %le %d",
      &pepoch, &alpha, &delta, &f0, &f1, &f2, &gsize);

      // GPS time of the first sample
      double gps1;
      sprintf (filename, "%s/%03d/starting_date", dtaprefix, ident);
      if ((data2 = fopen (filename, "r")) != NULL) {
      fscanf (data2, "%le", &gps1);
      fclose(data2);
      } else {
      perror (filename);
      return 1;
      }

      // Conversion of mjd to gps time
      double pepoch_gps = (pepoch - 44244)*86400 - 51.184;
      //			 gps1 = (gps1 - 44244)*86400 - 51.184;

      // Interpolation of ephemeris parameters to the starting time
      double *sgnlo;
      sgnlo = (double *) calloc (4, sizeof (double));

      double *be;
      be = (double *) calloc (2, sizeof (double));

      // ast2lin (auxi.c) returns the hemisphere number
      // and the vector be (used for sky position in linear coords.)
      pmr[0] = ast2lin(alpha, delta, epsm, be);

      sgnlo[0] = f0 + f1*(gps1 - pepoch_gps) + f2*pow(gps1 - pepoch_gps, 2)/2.;
      sgnlo[0] = 2*M_PI*2*sgnlo[0]*dt - oms;

      sgnlo[1] = f1 + f2*(gps1 - pepoch_gps);
      sgnlo[1] = M_PI*2*sgnlo[1]*dt*dt;

      sgnlo[2] = be[0]*(oms + sgnlo[0]);
      sgnlo[3] = be[1]*(oms + sgnlo[0]);

      // solving a linear system in order to translate
      // sky position, frequency and spindown (sgnlo parameters)
      // into the position in the grid

      gsl_vector *x = gsl_vector_alloc (4);
      int s;

      gsl_matrix_view m = gsl_matrix_view_array (M, 4, 4);
      gsl_matrix_transpose (&m.matrix) ;
      gsl_vector_view b = gsl_vector_view_array (sgnlo, 4);
      gsl_permutation *p = gsl_permutation_alloc (4);

      gsl_linalg_LU_decomp (&m.matrix, p, &s);
      gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);

      spndr[0] = round(gsl_vector_get(x, 1));
      nr[0]		= round(gsl_vector_get(x, 2));
      mr[0]		= round(gsl_vector_get(x, 3));

      gsl_permutation_free (p);
      gsl_vector_free (x);
      free (be);
      free (sgnlo);


      // Warnings and infos
      if(hemi)

      printf("Warning: -h switch hemisphere choice (%d) may be altered\nby the choice of -r grid range...\n", hemi);

      // Define the grid range in which the signal will be looked for
      spndr[1] = spndr[0] + gsize ;
      spndr[0] -= gsize;
      nr[1] = nr[0] + gsize ;
      nr[0] -= gsize;
      mr[1] = mr[0] + gsize ;
      mr[0] -= gsize;
      pmr[1] = pmr[0];

      }
      */
      fclose (data);

    } else {
      perror (opts->range);
      return;
    }

    /*
      ############ End setting range ################
    */

    printf("The following grid range is used\n");
    printf("(spndr, nr, mr, pmr pairs): %d %d %d %d %d %d %d %d\n", \
	   s_range->spndr[0], s_range->spndr[1], s_range->nr[0], s_range->nr[1],
	   s_range->mr[0], s_range->mr[1], s_range->pmr[0], s_range->pmr[1]);

  } else {

    // Establish the grid range in which the search will be performed
    // with the use of the M matrix from grid.bin
    gridr (sett->M, s_range->spndr, s_range->nr, s_range->mr, sett->oms, sett->Smax);
  }

}




void plan_fft(FFT_plans *plans, Arrays *arr, 
	      Detector_settings *sett, Command_line_opts *opts)
{
  /*
    ############ FFT Plans ################
  */

  //arrlen is maximum of Ninterp and fftpad*nfft
  arr->arr_len = (sett->fftpad * sett->nfft > sett->Ninterp ? sett->fftpad * sett->nfft : sett->Ninterp);

  CudaSafeCall ( cudaMalloc((void**)&arr->cu_xa, arr->arr_len*sizeof(cufftDoubleComplex)) );
  CudaSafeCall ( cudaMalloc((void**)&arr->cu_xb, arr->arr_len*sizeof(cufftDoubleComplex)) );
	
  CudaSafeCall ( cudaMalloc((void**)&arr->cu_xar, arr->arr_len*sizeof(cufftDoubleComplex)) );
  CudaSafeCall ( cudaMalloc((void**)&arr->cu_xbr, arr->arr_len*sizeof(cufftDoubleComplex)) );


  CudaSafeCall ( cudaMalloc((void**)&arr->cu_xa_f, arr->arr_len*sizeof(COMPLEX_TYPE)) );
  CudaSafeCall ( cudaMalloc((void**)&arr->cu_xb_f, arr->arr_len*sizeof(COMPLEX_TYPE)) );
	
  CudaSafeCall ( cudaMalloc((void**)&arr->cu_xar_f, arr->arr_len*sizeof(COMPLEX_TYPE)) );
  CudaSafeCall ( cudaMalloc((void**)&arr->cu_xbr_f, arr->arr_len*sizeof(COMPLEX_TYPE)) );

  if (opts->fftinterp == INT) { //interbinning
    CudaSafeCall ( cudaMalloc((void**)&arr->cu_xa2_f, sett->nfft*sizeof(COMPLEX_TYPE)) );
    CudaSafeCall ( cudaMalloc((void**)&arr->cu_xb2_f, sett->nfft*sizeof(COMPLEX_TYPE)) );
  }

  sett->nfftf = sett->fftpad*sett->nfft;

  if (opts->fftinterp == INT) { //interbinning
    cufftPlan1d( &(plans->plan),
		 sett->nfft,
		 CUFFT_TRANSFORM_TYPE, 1);
  } else { //fft & zero padding
    cufftPlan1d( &(plans->plan),
		 sett->nfftf,
		 CUFFT_TRANSFORM_TYPE, 1);

  }
	
  //plans for interpolation with splines
	
  cufftPlan1d(&(plans->pl_int),
	      sett->nfft,
	      CUFFT_Z2Z, 1);
		
  cufftPlan1d(&(plans->pl_inv),
	      sett->Ninterp,
	      CUFFT_Z2Z, 1);

  /*
    ############ FFT plans end ################
  */

}




void read_checkpoints(Search_range *s_range, int *FNum, Command_line_opts *opts) 
{

  /*
    ############ Checkpoints ################
  */

  // Checkpointing
  if(opts->checkp_flag) {
		
    // filename of checkpoint state file, depending on the hemisphere
    if(opts->hemi)
      sprintf (opts->qname, "state_%03d_%03d%s_%d.dat", opts->ident, opts->band, opts->label, opts->hemi);
    else
      sprintf (opts->qname, "state_%03d_%03d%s.dat", opts->ident, opts->band, opts->label);

    FILE *state;
    if ((state = fopen (opts->qname, "r")) != NULL) {

      // Scan the state file to get last recorded parameters
      if ((fscanf (state, "%d %d %d %d %d", &s_range->pst, &s_range->mst,
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


  /*
    ############ Checkpoints end ################
  */

}


void cleanup(Detector_settings *sett,
	     Command_line_opts *opts,
	     Search_range *s_range,
	     Arrays *arr,
	     FFT_plans *plans,
	     Ampl_mod_coeff *amod,
	     FLOAT_TYPE *cu_F) 
{

  free(arr->xDat);
	
  free(arr->sinmodf);
  free(arr->cosmodf);
  free(arr->aa);
  free(arr->bb);
  free(arr->DetSSB);
	
  free(arr->cand_buffer);
	
  cudaFree(arr->cu_xa);
  cudaFree(arr->cu_xb);
  cudaFree(arr->cu_xar);
  cudaFree(arr->cu_xbr);
  cudaFree(arr->cu_xa_f);
  cudaFree(arr->cu_xb_f);
  cudaFree(arr->cu_xar_f);
  cudaFree(arr->cu_xbr_f);
  cudaFree(arr->cu_xDat);
	
  cudaFree(arr->cu_aa);
  cudaFree(arr->cu_bb);
	
  cudaFree(arr->cu_shft);
  cudaFree(arr->cu_shftf);
  cudaFree(arr->cu_tshift);
  cudaFree(arr->cu_DetSSB);
	
  cudaFree(arr->cu_d);
  cudaFree(arr->cu_dl);
  cudaFree(arr->cu_du);
  cudaFree(arr->cu_B);

  cudaFree(cu_F);

  cudaFree(arr->cu_sinmodf);
  cudaFree(arr->cu_cosmodf);

  cudaFree(arr->cu_cand_buffer);
  cudaFree(arr->cu_cand_params);
  cudaFree(arr->cu_cand_count);
	
  free(sett->M);

  if (opts->fftinterp == INT ) {//interbinning
    cudaFree(arr->cu_xa2_f);
    cudaFree(arr->cu_xb2_f);
  }
	
}
