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



void read_checkpoints(Search_range *s_range, int *FNum, Command_line_opts *opts) {

  /*
    ############ Checkpointy ################
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
    ############ Koniec checkpointow ################
  */

}




void plan_fftw(FFTW_plans *plans, FFTW_arrays *fftw_arr, Signals *sig, 
	       Aux_arrays *aux_arr, Detector_settings *sett, Command_line_opts *opts)
{
  /*
    ############ FFT Plans ################
  */

  char hostname[32], wfilename[64];
  FILE *wisdom;

  // Imports a "wisdom file" containing information about how to optimally
  // compute Fourier transforms on a given machine. If such file is not
  // present, it will be created after the measure runs of the fft_plans
  // are performed below
  // (more info at http://www.fftw.org/fftw3_doc/Wisdom.html)
  gethostname (hostname, 32);
  sprintf (wfilename, "wisdom-%s.dat", hostname);
  if ((wisdom = fopen (wfilename, "r")) != NULL) {
    fftw_import_wisdom_from_file (wisdom);
    fclose (wisdom);
  }

  int N = sett->N;


  //array length (xa, xb) is max{fftpad*nfft , Ninterp}

  sett->Ninterp = sett->interpftpad*sett->nfft;
  fftw_arr->arr_len = (sett->fftpad * sett->nfft > sett->Ninterp ? sett->fftpad * sett->nfft : sett->Ninterp);

  fftw_arr->xa = fftw_malloc( 2 * fftw_arr->arr_len*sizeof(fftw_complex));
  fftw_arr->xb = fftw_arr->xa + fftw_arr->arr_len;

  sett->nfftf = sett->fftpad*sett->nfft;



   plans->plan = fftw_plan_dft_1d(sett->nfftf, fftw_arr->xa, fftw_arr->xa, FFTW_FORWARD, FFTW_MEASURE);
  plans->plan2 = fftw_plan_dft_1d(sett->nfftf, fftw_arr->xb, fftw_arr->xb, FFTW_FORWARD, FFTW_MEASURE);
	                             
  plans->pl_int = fftw_plan_dft_1d(sett->nfft, fftw_arr->xa, fftw_arr->xa, FFTW_FORWARD, FFTW_MEASURE);
  plans->pl_int2 = fftw_plan_dft_1d(sett->nfft, fftw_arr->xb, fftw_arr->xb, FFTW_FORWARD, FFTW_MEASURE);
	                             
  plans->pl_inv = fftw_plan_dft_1d(sett->Ninterp, fftw_arr->xa, fftw_arr->xa, FFTW_BACKWARD, FFTW_MEASURE);
  plans->pl_inv2 = fftw_plan_dft_1d(sett->Ninterp, fftw_arr->xb, fftw_arr->xb, FFTW_BACKWARD, FFTW_MEASURE);
	                             

  // Generates a 'wisdom' FFT file if there is none
  if ((wisdom = fopen (wfilename, "r")) == NULL) {

    wisdom = fopen (wfilename, "w");
    fftw_export_wisdom_to_file (wisdom);
  }
  fclose (wisdom);

  /*
    ############ Koniec planow FFT ################
  */

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
    ############ Zmiana zakresu wyszukiwania ################
  */

  FILE *data;

  // If the parameter range is invoked, the search is performed
  // within the range of grid parameters from an ascii file
  // ("-r range_file" from the command line)
  if (strlen (opts->range)) {

    if ((data=fopen (opts->range, "r")) != NULL) {

      int aqq = fscanf (data, "%d %d %d %d %d %d %d %d",
			s_range->spndr, 1+s_range->spndr, s_range->nr,
			1+s_range->nr, s_range->mr, 1+s_range->mr,
			s_range->pmr, 1+s_range->pmr);

      if (aqq) {
			
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
      ############ Koniec zmiany zakresu wyszukiwania ################
    */

    /*
      ############ Generacja grida ################
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




void read_grid(Detector_settings *sett, Command_line_opts *opts)
{

  /*
    ############ Macierz generująca grid ################
  */


  sett->M = (double *) calloc (16, sizeof (double));

  FILE *data;
  char filename[64];
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
    return;
  }
}



void init_arrays(Signals *sig, Aux_arrays *aux_arr, double** F, 
		 Command_line_opts *opts, Detector_settings *sett) 
{

  // Allocates and initializes to zero the data, detector ephemeris
  // and the F-statistic arrays
  sig->xDat = (double *) calloc (sett->N, sizeof (double));
  aux_arr->DetSSB = (double *) calloc (3*sett->N, sizeof (double));
  *F = (double *) calloc (2*sett->nfft, sizeof (double));

  aux_arr->aa = (double *) calloc (sett->N, sizeof (double));
  aux_arr->bb = (double *) calloc (sett->N, sizeof (double));

  aux_arr->shft = (double *) calloc (sett->N, sizeof (double));
  aux_arr->shftf = (double *) calloc (sett->N, sizeof (double));
  sig->xDatma = (complex double *) calloc (sett->N, sizeof (complex double));
  sig->xDatmb = (complex double *) calloc (sett->N, sizeof (complex double));



  char filename[64];
  FILE *data;
  // Input time-domain data handling
  sprintf (filename, "%s/%03d/xdat_%03d_%03d%s.bin", opts->dtaprefix, opts->ident, \
	   opts->ident, opts->band, opts->label);
  if ((data = fopen (filename, "r")) != NULL) {
    fread ((void *)(sig->xDat), sizeof (double), sett->N, data); // !!! wczytanie danych
    fclose (data);
  } else {
    perror (filename);
    return ;
  }


  int Nzeros=0;
  int i;
  // Checking for null values in the data
  for(i=0; i < sett->N; i++)
    if(!sig->xDat[i]) Nzeros++;

  // factor N/(N - Nzeros) to account for null values in the data
  sett->crf0 = (double)sett->N/(sett->N-Nzeros);


  //if white noise...
  if (opts->white_flag)
    sett->sig2 = sett->N*var (sig->xDat, sett->N);
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
    fread ((void *)(aux_arr->DetSSB), sizeof (double), 3*sett->N, data);
    // Deterministic phase defining the position of the Earth
    // in its diurnal motion at t=0
    fread ((void *)(&phir), sizeof (double), 1, data);
    // Earth's axis inclination to the ecliptic at t=0
    fread ((void *)(&epsm), sizeof (double), 1, data);
    fclose (data);
  } else {
    perror (filename);
    return ;
  }


  /*
    ############ Sincos ################
  */


  sett->sphir = sin (phir);
  sett->cphir = cos (phir);
  sett->sepsm = sin (epsm);
  sett->cepsm = cos (epsm);
}




void handle_opts(int argc, char* argv[], Command_line_opts *opts, Detector_settings *sett) 
{
	
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
    printf ("Threshold for the F-statistic is %lf\n", opts->trl);
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



void cleanup(Detector_settings *sett,
	     Command_line_opts *opts,
	     Search_range *s_range,
	     FFTW_arrays *fftw_arr,
	     Signals *sig,
	     FFTW_plans *plans,
	     Aux_arrays *aux,
	     Ampl_mod_coeff *amod,
	     double *F) 
{

  free(sig->xDat);
  free(sig->xDatma);
	
  free(aux->sinmodf);
  free(aux->cosmodf);
  free(aux->t2);
  free(aux->aa);
  free(aux->bb);
  free(aux->shftf);
  free(aux->shft);
  free(aux->DetSSB);
  free(F);
	
  fftw_free(fftw_arr->xa);
  //	if (opts->fftinterp ==) {
  //		fftw_free(fftw_arr->xao);
  //	}
	
  free(sett->M);
	
  fftw_destroy_plan(plans->plan);
  fftw_destroy_plan(plans->pl_int);
  fftw_destroy_plan(plans->pl_inv);
}
