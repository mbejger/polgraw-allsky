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
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <time.h>

#include "init.h"
#include "struct.h"
#include "settings.h"
#include "auxi.h"



	/*  Command line options handling: search 
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
  opts->range[0]  = '\0';
  opts->addsig[0] = '\0';
	
  // Initial value of starting frequency set to a negative quantity. 
  // If this is not changed by the command line value, fpo is calculated 
  // from the band number b (fpo = fpo = 100. + 0.96875*b)
  sett->fpo = -1;

  // Default initial value of the data sampling time 
  sett->dt = 0.5; 

  opts->help_flag=0;
  opts->white_flag=0;
  opts->s0_flag=0;
  opts->checkp_flag=0;

  static int help_flag=0, white_flag=0, s0_flag=0, checkp_flag=1;

  // Reading arguments 

  while (1) {
    static struct option long_options[] = {
      {"help", no_argument, &help_flag, 1},
      {"whitenoise", no_argument, &white_flag, 1},
      {"nospindown", no_argument, &s0_flag, 1},
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
      {"threshold", required_argument, 0, 't'},
      // hemisphere
      {"hemisphere", required_argument, 0, 'h'},
      // fpo value
      {"fpo", required_argument, 0, 'p'},
      // add signal parameters
      {"addsig", required_argument, 0, 'x'},
      // data sampling time 
      {"dt", required_argument, 0, 's'},
      {0, 0, 0, 0}
    };

    if (help_flag) {

      printf("polgraw-allsky CGW search code using the F-statistic\n");
      printf("Usage: ./search -[switch1] <value1> -[switch2] <value2> ...\n") ;
      printf("Switches are:\n\n");
      printf("-d, -data         Data directory (default is .)\n");
      printf("-o, -output       Output directory (default is ./candidates)\n");
      printf("-i, -ident        Frame number\n");
      printf("-b, -band         Band number\n");
      printf("-l, -label        Custom label for the input and output files\n");
      printf("-r, -range        File with grid range or pulsar position\n");
      printf("-c, -cwd          Change to directory <dir>\n");
      printf("-t, -threshold    Threshold for the F-statistic (default is 20)\n");
      printf("-h, -hemisphere   Hemisphere (default is 0 - does both)\n");
      printf("-p, -fpo          Reference band frequency fpo value\n");
      printf("-s, -dt           Data sampling time dt (default value: 0.5)\n");
      printf("-x, -addsig       Add signal with parameters from <file>\n\n");

      printf("Also:\n\n");
      printf("--whitenoise      White Gaussian noise assumed\n");
      printf("--nospindown      Spindowns neglected\n");
      printf("--nocheckpoint    State file won't be created (no checkpointing)\n");
      printf("--help            This help\n");

      exit (0);
    }

    int option_index = 0;
    int c = getopt_long_only (argc, argv, "i:b:o:d:l:r:c:t:h:p:x:s:", long_options, &option_index);
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
      opts->band = atoi(optarg);
      break;
    case 'o':
      strcpy(opts->prefix, optarg);
      break;
    case 'd':
      strcpy(opts->dtaprefix, optarg);
      break;
    case 'l':
      opts->label[0] = '_';
      strcpy(1+opts->label, optarg);
      break;
    case 'r':
      strcpy(opts->range, optarg);
      break;
    case 'c':
      opts->wd = (char *) malloc (1+strlen(optarg));
      strcpy(opts->wd, optarg);
      break;
    case 'p':
      sett->fpo = atof(optarg);
      break;
    case 'x':
      strcpy(opts->addsig, optarg);
      break;
    case 's':
      sett->dt = atof(optarg);
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
	
  printf("Input data directory is %s\n", opts->dtaprefix);
  printf("Output directory is %s\n", opts->prefix);
  printf("Frame and band numbers are %d and %d\n", opts->ident, opts->band);

  // Starting band frequency:
  // fpo_val is optionally read from the command line
  // Its initial value is set to -1
  if(!(sett->fpo >= 0))
    // The usual definition:
    sett->fpo = 100. + 0.96875 * opts->band;

  printf("The reference frequency fpo is %f\n", sett->fpo);
  printf("The data sampling time dt is  %f\n", sett->dt); 

  if (opts->white_flag)
    printf ("Assuming white Gaussian noise\n");

  // For legacy: FFT is now the only option 
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

  if (strlen(opts->addsig))
    printf ("Adding signal from '%s'\n", opts->addsig);
  if (opts->wd) {
    printf ("Changing working directory to %s\n", opts->wd);
    if (chdir(opts->wd)) { perror (opts->wd); abort (); }
  }

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

  FILE *data;

  for(i=0; i<sett->nifo; i++) { 

    ifo[i].sig.xDat = (double *) calloc(sett->N, sizeof(double));

    // Input time-domain data handling
    // 
    // The file name ifo[i].xdatname is constructed 
    // in settings.c, while looking for the detector 
    // subdirectories

    if((data = fopen(ifo[i].xdatname, "r")) != NULL) {
      status = fread((void *)(ifo[i].sig.xDat), 
               sizeof(double), sett->N, data);
      fclose (data);

    } else {
      perror (ifo[i].xdatname);
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
    char filename[512];
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

      printf("Using %s as detector %s ephemerids...\n", filename, ifo[i].name);

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


	/* Search range 
	 */ 

void set_search_range(
	Search_settings *sett, 
	Command_line_opts *opts, 
	Search_range *s_range) { 

  // Hemispheres (with respect to the ecliptic)
  if(opts->hemi) {
    s_range->pmr[0] = opts->hemi;
    s_range->pmr[1] = opts->hemi;

  } else {
    s_range->pmr[0] = 1;
    s_range->pmr[1] = 2;
  }

  // If the parameter range is invoked, the search is performed
  // within the range of grid parameters from an ascii file
  // ("-r range_file" from the command line)
  FILE *data;
  if (strlen (opts->range)) {

    if ((data=fopen (opts->range, "r")) != NULL) {

      int aqq = fscanf (data, "%d %d %d %d %d %d %d %d",
			s_range->spndr, 1+s_range->spndr, s_range->nr,
			1+s_range->nr, s_range->mr, 1+s_range->mr,
			s_range->pmr, 1+s_range->pmr);

      if (aqq) {
			
      }
      /*
      //#mb commented-out for now - useful for tests 

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
      exit(EXIT_FAILURE);
    }

    // Grid range is established from above
    printf("The following grid range is used\n");
    printf("(spndr, nr, mr, pmr pairs): %d %d %d %d %d %d %d %d\n", \
	   s_range->spndr[0], s_range->spndr[1], s_range->nr[0], s_range->nr[1],
	   s_range->mr[0], s_range->mr[1], s_range->pmr[0], s_range->pmr[1]);

  } else {

    // Establish the grid range in which the search will be performed
    // with the use of the M matrix from grid.bin
    gridr(
	    sett->M, 
	    s_range->spndr,
	    s_range->nr,
	    s_range->mr,
	    sett->oms,
	    sett->Smax);

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

  // Change FFTW_MEASURE to FFTW_PATIENT for more optimized plan
  // (takes more time to generate the wisdom file)
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
    free(ifo[i].sig.xDatmb);
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

} // end of cleanup & memory free 


	/*	Command line options handling: coincidences  
	 */ 
	
void handle_opts_coinc(
    Search_settings *sett, 
    Command_line_opts_coinc *opts,
	int argc, 
	char* argv[]) {
	
  opts->wd=NULL;
	
  strcpy (opts->prefix, TOSTR(PREFIX));
  strcpy (opts->dtaprefix, TOSTR(DTAPREFIX));

  // Default initial value of the data sampling time 
  sett->dt = 0.5;

  opts->help_flag=0;
  static int help_flag=0;

  // Reading arguments 

  while (1) {
    static struct option long_options[] = {
      {"help", no_argument, &help_flag, 1},
      // shift of cells 
      {"shift", required_argument, 0, 's'},
      // Cell scalling in frequency
      {"scale_f", required_argument, 0, 'f'},
      // Cell scalling in spindown 
      {"scale_s", required_argument, 0, 'z'},
      // Cell scalling in right ascension
      {"scale_a", required_argument, 0, 'a'},
      // Cell scalling in declination
      {"scale_d", required_argument, 0, 'b'},
      // Reference frame number 
      {"refr", required_argument, 0, 'r'},
      // output directory
      {"output", required_argument, 0, 'o'},
      // input data directory
      {"data", required_argument, 0, 'd'},
      // fpo value
      {"fpo", required_argument, 0, 'p'},
      // data sampling time 
      {"dt", required_argument, 0, 't'},
      // triggers' name prefactor 
      {"trigname", required_argument, 0, 'e'},
      // Location of the reference grid.bin and starting_date files  
      {"refloc", required_argument, 0, 'g'},
      {0, 0, 0, 0}
    };

    if (help_flag) {

      printf("polgraw-allsky CGW code for concidences among candidates\n");
      printf("Usage: ./search -[switch1] <value1> -[switch2] <value2> ...\n") ;
      printf("Switches are:\n\n");
      printf("-data         Data directory (default is ./candidates)\n");
      printf("-output       Output directory (default is ./coinc-results)\n");
      printf("-shift        Shift of cells\n");
      printf("-scale_f      Cell scalling in frequency\n");
      printf("-scale_s      Cell scalling in spindown\n");
      printf("-scale_a      Cell scalling in right ascenscion\n");
      printf("-scale_d      Cell scalling in declination\n");
      printf("-refr         Reference frame number\n");
      printf("-fpo          Reference band frequency fpo value\n");
      printf("-dt           Data sampling time dt (default value: 0.5)\n");
      printf("-trigname     Triggers' name prefactor\n");
      printf("-refloc       Location of the reference grid.bin and starting_date files\n\n");

      printf("Also:\n\n");
      printf("--help		This help\n");

      exit (0);
    }

    int option_index = 0;
    int c = getopt_long_only (argc, argv, "f:p:o:d:b:s:a:z:r:t:e:g:", long_options, &option_index);
    if (c == -1)
      break;

    switch (c) {
    case 'p':
      sett->fpo = atof(optarg);
      break;
    case 's':
      opts->shift = atof(optarg);
      break;
    case 'f':   // -scale_f 
      opts->scale[0] = atoi(optarg);
      break;
    case 'z':   // -scale_s 
      opts->scale[1] = atoi(optarg);
      break;
    case 'b':  // -scale_d 
      opts->scale[2] = atoi(optarg);
      break;
    case 'a':  // -scale_a 
      opts->scale[3] = atoi(optarg);
      break;
    case 'r':
      opts->refr = atoi(optarg);
      break;
    case 'o':
      strcpy(opts->prefix, optarg);
      break;
    case 'd':
      strcpy(opts->dtaprefix, optarg);
      break;
    case 't':
      sett->dt = atof(optarg);
      break;
    case 'e':
      strcpy(opts->trigname, optarg);
      break;
    case 'g':
      strcpy(opts->refloc, optarg);
      break;
    case '?':
      break;
    default:
      break ;
    } /* switch c */
  } /* while 1 */

  printf("#mb add info at the beginning...\n"); 

} // end of command line options handling: coincidences  



	/* Manage grid matrix (read from grid.bin, find eigenvalues 
	 * and eigenvectors) and reference GPS time from starting_time
     * (expected to be in the same directory)    
	 */ 

void manage_grid_matrix(
	Search_settings *sett, 
	Command_line_opts_coinc *opts) {

  sett->M = (double *)calloc(16, sizeof (double));

  FILE *data;
  char filename[512];
  sprintf (filename, "%s/grid.bin", opts->refloc);

  if ((data=fopen (filename, "r")) != NULL) {

    printf("Reading the reference grid.bin at %s\n", opts->refloc);

    fread ((void *)&sett->fftpad, sizeof (int), 1, data);

	printf("fftpad from the grid file: %d\n", sett->fftpad); 
	
    fread ((void *)sett->M, sizeof(double), 16, data);
    // We actually need the second (Fisher) matrix from grid.bin, 
    // hence the second fread: 
    fread ((void *)sett->M, sizeof(double), 16, data);
    fclose (data);
  } else {
	  perror (filename);
      exit(EXIT_FAILURE);
  }

/* //#mb seems not needed at the moment 
  sprintf (filename, "%s/starting_date", opts->refloc);
  
  if ((data=fopen (filename, "r")) != NULL) {
    fscanf(data, "%le", &opts->refgps);

    printf("Reading the reference starting_date file at %s The GPS time is %12f\n", opts->refloc, opts->refgps);
    fclose (data);
  } else {
      perror (filename);
      exit(EXIT_FAILURE);
  }
*/ 

  // Calculating the eigenvectors and eigenvalues 
  gsl_matrix_view m = gsl_matrix_view_array(sett->M, 4, 4);

  gsl_vector *eval = gsl_vector_alloc(4);
  gsl_matrix *evec = gsl_matrix_alloc(4, 4);

  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(4); 
  gsl_eigen_symmv(&m.matrix, eval, evec, w);
  gsl_eigen_symmv_free(w);

  double eigval[4], eigvec[4][4]; 
  // Saving the results to the settings struct sett->vedva[][]
  { int i, j;
    for(i=0; i<4; i++) { 
      eigval[i] = gsl_vector_get(eval, i); 
      gsl_vector_view evec_i = gsl_matrix_column(evec, i);

      for(j=0; j<4; j++)   
        eigvec[j][i] = gsl_vector_get(&evec_i.vector, j);               
    } 

    // This is an auxiliary matrix composed of the eigenvector 
    // columns multiplied by a matrix with sqrt(eigenvalues) on diagonal  
    for(i=0; i<4; i++)
      for(j=0; j<4; j++)
        sett->vedva[i][j]  = eigvec[i][j]*sqrt(eigval[j]);  

  } 

  gsl_vector_free (eval);
  gsl_matrix_free (evec);

} // end of manage grid matrix  



static int colnum = 0; // columns to compare when searching coincidences
int way2compare(const void *a, const void *b){

    int** x = (int**) a;
	int** y = (int**) b;

	int xval, yval;

	xval = *(*(x) + colnum);
	yval = *(*(y) + colnum);

	if(xval < yval) return 1;
    else if(xval > yval) return -1;
	else return 0;
		
}

// Allocation of memory for martix with given number of rows and columns
int** matrix(int rows, int cols) {

	int k, **m;
	m = (int **)malloc(rows*sizeof(int *));

	for (k=0; k<rows; k++){
		m[k] = (int *)calloc(cols, sizeof(int));
	}

	return m;

}


  /* Convert triggers to linear integer coordinates 
   * with the use of eigenvectors and eigenvalues 
   * obtained in manage_grid_matrix()
   */ 

void convert_to_linear(
  Search_settings *sett,
  Command_line_opts_coinc *opts,
  Candidate_triggers *trig) {

  int i, j, k, numtr, shift[4], **MTR;
  double sqrN, tmp[4], v[4][4], be[2]; 

  sqrN = pow(sett->N, 2);

  // Total number of triggers 
  numtr = trig->num_of_trig; 

  // Calculating the shifts from opts->shift 
  int val = opts->shift;
  for(i=0; i<4; i++) shift[i] = 0; 
  i=3; 
  while (val > 0) { 
    if(val%10) shift[i] = val%10; 
    i--; val /= 10;
  }

  // Allocation of the matrix for the coincidence sorting and search 
  MTR = matrix(numtr, trig->size);

  // Dividing the transformation matrix elements by the corresponding 
  // scale factors outside the loop over candidates to save computational time  
  for(j=0; j<4; j++)  
    for(i=0; i<4; i++)
        v[i][j] = sett->vedva[i][j]/(opts->scale[j]);   
  
  // Loop over all candidates
  for(i=0; i<numtr; i++) { 

    //#mb account for a possibility  
    //#mb that the trigger may exit the band   

    // Conversion to linear parameters 
    tmp[0] = trig->f[i]*sett->N; 
    tmp[1] = trig->s[i]*sqrN; 

    // Transformation of astronomical to linear coordinates 
    // C_EPSMA, an average value of epsm, is defined in settings.h  
    ast2lin(trig->a[i], trig->d[i], C_EPSMA, be);

    // tmp[2] corresponds to trig->d[i]
    // tmp[3] corresponds to trig->a[i]
    tmp[2] = sett->oms*sett->N*be[0]; 
    tmp[3] = sett->oms*sett->N*be[1]; 

    // Saving the integer candidate values (0=f, 1=s, 2=d, 3=a)  
    for(j=0; j<4; j++)
      MTR[i][j] = round(tmp[0]*v[0][j] + tmp[1]*v[1][j] + tmp[2]*v[2][j] 
                      + tmp[3]*v[3][j] + 0.5*shift[j]);  

    MTR[i][4] = trig->fr[i]; 

    // printf("%d %d %d %d\n", MTR[i][0], MTR[i][1], MTR[i][2], MTR[i][3]); 

  }

  // Sorting 
  for (i=0; i<trig->size; i++) {
    colnum = i;
	qsort(MTR, numtr, sizeof(int *), way2compare);
  }
 
  // Finding the duplicate rows in the MTR matrix
  int **IMTR, diff=0, weight=0, coindx=0, numcoincidences, mincoincidences=2;

  numcoincidences = numtr/mincoincidences; // maximal possible amount of coincidents is the half or correlation with weight of coincidences

  IMTR = matrix(numcoincidences, 2); 

  for (i=0; i<numtr-1; i++) {
    diff = (MTR[i][0] - MTR[i+1][0]);
	
	if(!diff) { // if first column is not equal then rest could be ignored
	  for(j=0; j<trig->size; j++) 
	    // distance measure to compare rows  
        diff += (MTR[i][j] - MTR[i+1][j])*(MTR[i][j] - MTR[i+1][j]); 
			
	}

  if(diff==0) {
	  weight++;
	  if(weight==mincoincidences - 1) { //only coincidences with weight above 1 will be stored
	    coindx++;		
	  }
	  IMTR[coindx][0] = i;
	  IMTR[coindx][1] = weight;

	} else weight = 0;
		
	}

	colnum = 1;
	qsort(IMTR, numcoincidences, sizeof(int *), way2compare);

	for (i = 0; i < numcoincidences; i++) {
        int idxmtr = IMTR[i][0];
	printf("i %d weight %d coincidence %d %d %d %d floatings %8.10e %8.10e %8.10e %8.10e\n", idxmtr, IMTR[i][1],MTR[idxmtr][0],MTR[idxmtr][1],MTR[idxmtr][2],MTR[idxmtr][3], trig->f[idxmtr], trig->s[idxmtr], trig->d[idxmtr], trig->a[idxmtr]);
	}


}
