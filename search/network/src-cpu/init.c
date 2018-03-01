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

void handle_opts( Search_settings *sett, 
		  Command_line_opts *opts,
		  int argc, 
		  char* argv[]) {
	
  opts->hemi=0;
  opts->wd=NULL;

  // Default F-statistic threshold 
  opts->trl=20;

  strcpy (opts->prefix, TOSTR(PREFIX));
  strcpy (opts->dtaprefix, TOSTR(DTAPREFIX));

  opts->label[0]    = '\0';
  opts->range[0]    = '\0';
  opts->getrange[0] = '\0';
  opts->usedet[0]   = '\0';
  opts->addsig[0]   = '\0';
  opts->addline[0]  = '\0';

  	
  // Initial value of starting frequency set to a negative quantity. 
  // If this is not changed by the command line value, fpo is calculated 
  // from the band number b (fpo = fpo = fstart + 0.96875*b/(2dt))
  sett->fpo = -1;

  // Default initial value of the data sampling time 
  sett->dt = 0.5; 

  // Default value of the narrow-down parameter 
  opts->narrowdown=0.5; 

  // Initial value of the number of days is set to 0
  sett->nod = 0; 

  opts->help_flag=0;
  opts->white_flag=0;
  opts->s0_flag=0;
  opts->checkp_flag=0;
  opts->veto_flag=0; 

  static int help_flag=0, white_flag=0, s0_flag=0, 
             checkp_flag=1, veto_flag=0;

  // Reading arguments 

  while (1) {
    static struct option long_options[] = {
      {"help", no_argument, &help_flag, 1},
      {"whitenoise", no_argument, &white_flag, 1},
      {"nospindown", no_argument, &s0_flag, 1},
      {"nocheckpoint", no_argument, &checkp_flag, 0},
      {"vetolines", no_argument, &veto_flag, 1}, 
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
      // write full grid range to file
      {"getrange", required_argument, 0, 'g'},
      // change directory parameter
      {"cwd", required_argument, 0, 'c'},
      // threshold for the F-statistic
      {"threshold", required_argument, 0, 't'},
      // hemisphere
      {"hemisphere", required_argument, 0, 'h'},
      // fpo value
      {"fpo", required_argument, 0, 'p'},
      // add signal parameters
      {"addsig", required_argument, 0, 'x'},
      // add stationary line parameters
      {"addline", required_argument, 0, 'z'},
      // number of days in the time-domain segment 
      {"nod", required_argument, 0, 'y'},
      // which detectors to use
      {"usedet", required_argument, 0, 'u'}, 
      // data sampling time 
      {"dt", required_argument, 0, 's'},
      // Narrow down the frequency band (+- the center of band) 
      {"narrowdown", required_argument, 0, 'n'},
      {0, 0, 0, 0}
    };

    if (help_flag) {

      printf("polgraw-allsky periodic GWs: search for candidate signals with the F-statistic\n");
      printf("Usage: ./search -[switch1] <value1> -[switch2] <value2> ...\n") ;
      printf("Switches are:\n\n");
      printf("-data         Data directory (default is .)\n");
      printf("-output       Output directory (default is ./candidates)\n");
      printf("-ident        Frame number\n");
      printf("-band         Band number\n");
      printf("-label        Custom label for the input and output files\n");
      printf("-range        Use file with grid range or pulsar position\n");
      printf("-getrange     Write grid ranges & save fft wisdom & exit (ignore -r)\n");
      printf("-cwd          Change to directory <dir>\n");
      printf("-threshold    Threshold for the F-statistic (default is 20)\n");
      printf("-hemisphere   Hemisphere (default is 0 - does both)\n");
      printf("-fpo          Reference band frequency fpo value\n");
      printf("-dt           Data sampling time dt (default value: 0.5)\n");
      printf("-usedet       Use only detectors from string (default is use all available)\n");
      printf("-addsig       Add signal with parameters from <file>\n");
      printf("-addline      Add stationary line with parameters from <file>\n");
      printf("-nod          Number of days\n");
      printf("-narrowdown   Narrow-down the frequency band (range [0, 0.5] +- around center)\n\n");


      printf("Also:\n\n");
      printf("--whitenoise      White Gaussian noise assumed\n");
      printf("--nospindown      Spindowns neglected\n");
      printf("--nocheckpoint    State file won't be created (no checkpointing)\n");
      printf("--vetolines       Veto known lines from files in data directory\n");
      printf("--help            This help\n");

      exit(EXIT_SUCCESS);
    }

    int option_index = 0;
    int c = getopt_long_only(argc, argv, "i:b:o:d:l:r:g:c:t:h:p:x:z:y:s:u:n:", 
			     long_options, &option_index);
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
    case 'g':
      strcpy(opts->getrange, optarg);
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
    case 'z':
      strcpy(opts->addline, optarg);
      break;
    case 'y':
      sett->nod = atoi(optarg);
      break;
    case 's':
      sett->dt = atof(optarg);
      break;
    case 'u':
      strcpy(opts->usedet, optarg);
      break;
    case 'n':
      opts->narrowdown = atof(optarg);
      break;
    case '?':
      break;
    default:
      break ;
    } /* switch c */
  } /* while 1 */

  // Check if sett->nod was set up, if not, exit
  if(!(sett->nod)) { 
    printf("Number of days not set... Exiting\n"); 
    exit(EXIT_FAILURE); 
  } 

  printf("Number of days is %d\n", sett->nod); 

  // Putting the parameter in triggers' frequency range [0, pi] 
  opts->narrowdown *= M_PI; 

  opts->white_flag = white_flag;
  opts->s0_flag = s0_flag;
  opts->checkp_flag = checkp_flag;
  opts->veto_flag = veto_flag; 

  printf("Input data directory is %s\n", opts->dtaprefix);
  printf("Output directory is %s\n", opts->prefix);
  printf("Frame and band numbers are %d and %d\n", opts->ident, opts->band);

  // Starting band frequency:
  // fpo_val is optionally read from the command line
  // Its initial value is set to -1
  if(!(sett->fpo >= 0))

    // The usual definition (multiplying the offset by B=1/(2dt))
    // !!! in RDC_O1 the fstart equals 10, not 100 like in VSR1 !!! 
    // 
    sett->fpo = 10. + 0.96875*opts->band*(0.5/sett->dt);

  printf("The reference frequency fpo is %f\n", sett->fpo);
  printf("The data sampling time dt is %f\n", sett->dt); 

  if (opts->white_flag)
    printf ("Assuming white Gaussian noise\n");

  // For legacy: FFT is now the only option 
  printf ("Using fftinterp=FFT (FFT interpolation by zero-padding)\n");

  if(opts->trl!=20)
    printf ("Threshold for the F-statistic is %lf\n", opts->trl);
  if(opts->hemi)
    printf ("Search for hemisphere %d\n", opts->hemi);
  if(opts->s0_flag)
    printf ("Assuming s_1 = 0.\n");
  if(strlen(opts->label))
    printf ("Using '%s' as data label\n", opts->label);

  if(strlen(opts->getrange)){
    printf ("Writing full grid ranges to '%s'\n", opts->getrange);
    if(strlen(opts->range)) {
      opts->range[0] = '\0';
      printf ("     WARNING! -r option will be ignored...\n");
    }
  }

  if (strlen(opts->range))
    printf ("Obtaining grid range from '%s'\n", opts->range);

  if (strlen(opts->addsig))
    printf ("Adding signal from '%s'\n", opts->addsig);

  if (strlen(opts->addline))
    printf ("Adding stationary line from '%s'\n", opts->addline);

  if (opts->wd) {
    printf ("Changing working directory to %s\n", opts->wd);
    if (chdir(opts->wd)) { perror (opts->wd); abort (); }
  }

  if(opts->veto_flag) 
    printf("Known lines will be vetoed (reading from files in the data directory)\n");

} // end of command line options handling 


	/* Generate grid from the M matrix (grid.bin)
	 */ 

void read_grid(
	       Search_settings *sett, 
	       Command_line_opts *opts) {

  sett->M = (double *) calloc (16, sizeof (double));

  FILE *data;
  char filename[512];

  // In case when -usedet option is used for one detector
  // i.e. opts->usedet has a length of 2 (e.g. H1 or V1), 
  // read grid.bin from this detector subdirectory 
  // (see detectors_settings() in settings.c for details) 
  if(strlen(opts->usedet)==2)
    sprintf (filename, "%s/%03d/%s/grid.bin", opts->dtaprefix, opts->ident, opts->usedet);
  else 
    sprintf (filename, "%s/%03d/grid.bin", opts->dtaprefix, opts->ident);


  if ((data=fopen (filename, "r")) != NULL) {
    printf("Using grid file from %s\n", filename);
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


  /* Array initialization */ 

void init_arrays(
		 Search_settings *sett, 
		 Command_line_opts *opts,
		 Aux_arrays *aux_arr,
		 double** F) {

  int i; 
  size_t status; 
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
    /* 
    const size_t array_bytes = 3*sett->N*sizeof(double);
    ifo[i].sig.DetSSB = NULL;
    if ( posix_memalign((void**)&ifo[i].sig.DetSSB, 32, array_bytes) ) exit (1);
    */

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


  /* Add a signal to the data */ 

void add_signal(
		Search_settings *sett,
		Command_line_opts *opts,
		Aux_arrays *aux_arr,
		Search_range *s_range) {

  int i, j, n, gsize, reffr; 
  double snr=0, sum = 0., h0=0, cof, d1; 
  double sigma_noise = 1.0;
  double be[2];
  double sinaadd, cosaadd, sindadd, cosdadd, phaseadd, shiftadd; 
  double nSource[3], sgnlo[8], sgnlol[4];
 
  char amporsnr[3];  
 
  FILE *data;
  
  // Signal parameters are read
  if ((data=fopen (opts->addsig, "r")) != NULL) {
	
    // Fscanning for the GW amplitude h0 or signal-to-noise,  
    // the grid size and the reference frame 
    // (for which the signal freq. is not spun-down/up)

    fscanf (data, "%s", amporsnr);    

    if(!strcmp(amporsnr, "amp")) { 
      fscanf (data, "%le %d %d", &h0, &gsize, &reffr); 
      printf("add_signal(): GW amplitude h0 is %le\n", h0); 
    } else if(!strcmp(amporsnr, "snr")) { 
      fscanf (data, "%le %d %d", &snr, &gsize, &reffr); 
      printf("add_signal(): GW (network) signal-to-noise ratio is %le\n", snr); 
    } else { 
      printf("Problem with the signal file. Exiting...\n"); 
      exit(0); 
    } 

    // Fscanning signal parameters: f, fdot, delta, alpha (sgnlo[0], ..., sgnlo[3])
    // four amplitudes sgnlo[4], ..., sgnlo[7] 
    // (see sigen.c and Phys. Rev. D 82, 022005 2010, Eqs. 2.13a-d) 

    for(i=0; i<8; i++)
      fscanf(data, "%le",i+sgnlo); 
    
    fclose (data);
                 
  } else {
    perror (opts->addsig);
  }
  
  // Search-specific parametrization of freq. 
  // for the software injections
  // sgnlo[0]: frequency, sgnlo[1]: frequency. derivative  
 
  sgnlo[0] += -2.*sgnlo[1]*(sett->N)*(reffr - opts->ident); 
 
  // Check if the signal is in band 
  if(sgnlo[0]<0) exit(171);          // &laquo;  
  else if (sgnlo[0]>M_PI) exit(187); // &raquo;

  cof = sett->oms + sgnlo[0]; 
  
  for(i=0; i<2; i++) sgnlol[i] = sgnlo[i]; 
  
  // Calculate the hemisphere and be vector 
  s_range->pmr[0] = ast2lin(sgnlo[3], sgnlo[2], C_EPSMA, be);

  sgnlol[2] = be[0]*cof;  
  sgnlol[3] = be[1]*cof; 

 		 	
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
  
  printf("add_signal(): following grid range is used (spndr, nr, mr, pmr pairs)\n");
  printf("%d %d %d %d %d %d %d %d\n", \
   s_range->spndr[0], s_range->spndr[1], s_range->nr[0], s_range->nr[1],
   s_range->mr[0], s_range->mr[1], s_range->pmr[0], s_range->pmr[1]);

  // sgnlo[2]: declination, sgnlo[3]: right ascension 
  sindadd = sin(sgnlo[2]); 
  cosdadd = cos(sgnlo[2]); 
  sinaadd = sin(sgnlo[3]);  
  cosaadd = cos(sgnlo[3]); 
	
  // To keep coherent phase between time segments  
  double phaseshift = sgnlo[0]*sett->N*(reffr - opts->ident)   
    + sgnlo[1]*pow(sett->N*(reffr - opts->ident), 2); 


  // Allocate arrays for added signal, for each detector 
  double **signadd = malloc((sett->nifo)*sizeof(double *));
  for(n=0; n<sett->nifo; n++)
    signadd[n] = malloc((sett->N)*sizeof(double));

  // Loop for each detector - sum calculations
  for(n=0; n<sett->nifo; n++) {
    
    modvir(sinaadd, cosaadd, sindadd, cosdadd,
	   sett->N, &ifo[n], aux_arr);

    nSource[0] = cosaadd*cosdadd;
    nSource[1] = sinaadd*cosdadd;
    nSource[2] = sindadd;
					
    for (i=0; i<sett->N; i++) {
      
      shiftadd = 0.; 					 
      for (j=0; j<3; j++)
      	shiftadd += nSource[j]*ifo[n].sig.DetSSB[i*3+j];		 
      
      // Phase 
      phaseadd = sgnlo[0]*i + sgnlo[1]*aux_arr->t2[i] 
        + (cof + 2.*sgnlo[1]*i)*shiftadd
        - phaseshift; 

      // The whole signal with 4 amplitudes and modulations 
      signadd[n][i] = sgnlo[4]*(ifo[n].sig.aa[i])*cos(phaseadd) 
                    + sgnlo[6]*(ifo[n].sig.aa[i])*sin(phaseadd) 
                    + sgnlo[5]*(ifo[n].sig.bb[i])*cos(phaseadd) 
                    + sgnlo[7]*(ifo[n].sig.bb[i])*sin(phaseadd);

      // Sum over signals
      sum += pow(signadd[n][i], 2.);
    
    } // data loop
   
  } // detector loop


  // Signal amplitude h0 from the snr 
  // (currently only makes sense for Gaussian noise with fixed sigma)
  if(snr)
    h0 = (snr*sigma_noise)/(sqrt(sum));

  // Loop for each detector - adding signal to data (point by point)  								
  for(n=0; n<sett->nifo; n++) {
    for (i=0; i<sett->N; i++) {

      // Adding the signal to the data vector 
      if(ifo[n].sig.xDat[i]) { 
        ifo[n].sig.xDat[i] += h0*signadd[n][i];

      } 

    } // data loop

  } // detector loop

  // printf("snr=%le h0=%le\n", snr, h0);

  // Free auxiliary 2d array 
  for(n=0; n<sett->nifo; n++) 
    free(signadd[n]);
  free(signadd);
 
} // add_signal()


  /* Add a stationary line to the data */ 

void add_stationary_line(
		Search_settings *sett,
		Command_line_opts *opts,
		Aux_arrays *aux_arr,
		Search_range *s_range) {

  int i, j, n, gsize, reffr; 
  double snr=0, sum = 0., h0=0, cof, d1; 
  double sigma_noise = 1.0;
  double be[2];
  double sinaadd, cosaadd, sindadd, cosdadd, phaseadd, shiftadd; 
  double nSource[3], sgnlo[8], sgnlol[4];
 
  char amporsnr[3];  
 
  FILE *data;
  
  // Line parameters: this function uses files 
  // generated by sigen.c 
  // Some parameters are omitted 
  if ((data=fopen (opts->addline, "r")) != NULL) {
	
    // Fscanning for the GW amplitude h0 or signal-to-noise,  
    // the grid size and the reference frame 
    // (for which the signal freq. is not spun-down/up)

    fscanf (data, "%s", amporsnr);    

    if(!strcmp(amporsnr, "amp")) { 
      fscanf (data, "%le %d %d", &h0, &gsize, &reffr); 
      printf("add_stationary_line(): Amplitude h0 is %le\n", h0); 
    } else if(!strcmp(amporsnr, "snr")) { 
      fscanf (data, "%le %d %d", &snr, &gsize, &reffr); 
      printf("add_stationary_line(): Network signal-to-noise ratio is %le\n", snr); 
    } else { 
      printf("Problem with the line file. Exiting...\n"); 
      exit(0); 
    } 

    // Fscanning line frequency f and other parameters: 
    // fdot, delta, alpha (sgnlo[0], ..., sgnlo[3])
    // four amplitudes sgnlo[4], ..., sgnlo[7] 
    // (see sigen.c and Phys. Rev. D 82, 022005 2010, Eqs. 2.13a-d) 
    // Here f is used, fdot == 0

    for(i=0; i<8; i++)
      fscanf(data, "%le",i+sgnlo); 
    
    fclose (data);
                 
  } else {
    perror (opts->addline);
  }
  
  // Line is stationary, does not change frequency 
  // sgnlo[0]: frequency, sgnlo[1]: frequency derivative  
  sgnlo[1] = 0; 

  cof = sett->oms + sgnlo[0]; 
  
  for(i=0; i<2; i++) sgnlol[i] = sgnlo[i]; 

  // We are looking at a specific part of the sky  
  // Calculate the hemisphere and be vector 
  s_range->pmr[0] = ast2lin(sgnlo[3], sgnlo[2], C_EPSMA, be);

  sgnlol[2] = be[0]*cof;  
  sgnlol[3] = be[1]*cof; 


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
  
  printf("add_stationary_file(): following grid range is used (spndr, nr, mr, pmr pairs)\n");
  printf("%d %d %d %d %d %d %d %d\n", \
   s_range->spndr[0], s_range->spndr[1], s_range->nr[0], s_range->nr[1],
   s_range->mr[0], s_range->mr[1], s_range->pmr[0], s_range->pmr[1]);

  // sgnlo[2]: declination, sgnlo[3]: right ascension 
  sindadd = sin(sgnlo[2]); 
  cosdadd = cos(sgnlo[2]); 
  sinaadd = sin(sgnlo[3]);  
  cosaadd = cos(sgnlo[3]); 
	
  //#mb not sure if needed 
  // To keep coherent phase between time segments  
  double phaseshift = sgnlo[0]*sett->N*(reffr - opts->ident);    

  // Allocate arrays for added signal, for each detector 
  double **signadd = malloc((sett->nifo)*sizeof(double *));
  for(n=0; n<sett->nifo; n++)
    signadd[n] = malloc((sett->N)*sizeof(double));

  // Loop for each detector - sum calculations
  for(n=0; n<sett->nifo; n++) {
					
    for (i=0; i<sett->N; i++) {
      
      // Phase 
      phaseadd = sgnlo[0]*i - phaseshift; 

      // The whole signal with 4 amplitudes and modulations
      // Let's just use one amplitude e.g. sgnlo[4]  
      signadd[n][i] = sgnlo[4]*cos(phaseadd);  
/*                    + sgnlo[6]*sin(phaseadd) 
                    + sgnlo[5]*cos(phaseadd) 
                    + sgnlo[7]*sin(phaseadd);
*/ 

      // Sum over signals
      sum += pow(signadd[n][i], 2.);
    
    } // data loop
   
  } // detector loop


  // Signal amplitude h0 from the snr 
  // (currently only makes sense for Gaussian noise with fixed sigma)
  if(snr)
    h0 = (snr*sigma_noise)/(sqrt(sum));

  // Loop for each detector - adding signal to data (point by point)  								
  for(n=0; n<sett->nifo; n++) {
    for (i=0; i<sett->N; i++) {

      // Adding the signal to the data vector 
      if(ifo[n].sig.xDat[i]) { 
        ifo[n].sig.xDat[i] += h0*signadd[n][i];

      } 

    } // data loop

  } // detector loop

  // printf("snr=%le h0=%le\n", snr, h0);

  // Free auxiliary 2d array 
  for(n=0; n<sett->nifo; n++) 
    free(signadd[n]);
  free(signadd);
 
} // add_stationary_line()



/* Search range */ 

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

      int aqq = fscanf(data, "%d %d %d %d %d %d %d %d",
		       s_range->spndr, 1+s_range->spndr, 
		       s_range->nr, 1+s_range->nr, 
		       s_range->mr, 1+s_range->mr,
		       s_range->pmr, 1+s_range->pmr);
      
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

    if (strlen(opts->getrange)) {

      FILE *data;
      if ((data=fopen (opts->getrange, "w")) != NULL) {
	fprintf(data, "%d %d\n%d %d\n%d %d\n%d %d\n",
		s_range->spndr[0], s_range->spndr[1],
		s_range->nr[0], s_range->nr[1],
		s_range->mr[0], s_range->mr[1],
		s_range->pmr[0], s_range->pmr[1] );
	
	printf("Wrote input data grid ranges to %s\n", opts->getrange);
	fclose (data);
	//	exit(EXIT_SUCCESS);
	
      } else {
	
	printf("Can't open %s file for writing\n", opts->getrange);
       	exit(EXIT_FAILURE);
	
      }
    }

  }

  printf("set_search_range() - the grid ranges are maximally this:\n");
  printf("(spndr, nr, mr, pmr pairs): %d %d %d %d %d %d %d %d\n",	\
	 s_range->spndr[0], s_range->spndr[1], s_range->nr[0], s_range->nr[1],
	 s_range->mr[0], s_range->mr[1], s_range->pmr[0], s_range->pmr[1]);
  
  printf("Smin: %le, -Smax: %le\n", sett->Smin, sett->Smax); 

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
	                             
  plans->pl_int = fftw_plan_dft_1d(sett->nfft, fftw_arr->xa, fftw_arr->xa, FFTW_FORWARD, FFTW_MEASURE);
	                             
  plans->pl_inv = fftw_plan_dft_1d(sett->Ninterp, fftw_arr->xa, fftw_arr->xa, FFTW_BACKWARD, FFTW_MEASURE);
	                             
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
      sprintf(opts->qname, "state_%03d_%04d%s_%d.dat",  
	            opts->ident, opts->band, opts->label, opts->hemi);
    else
      sprintf(opts->qname, "state_%03d_%04d%s.dat", 
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
  fftw_destroy_plan(plans->plan2);
  fftw_destroy_plan(plans->pl_int);
  fftw_destroy_plan(plans->pl_int2);
  fftw_destroy_plan(plans->pl_inv);
  fftw_destroy_plan(plans->pl_inv2);

  fftw_forget_wisdom();
  fftw_cleanup();


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

  // Initial value of the number of days is set to 0
  sett->nod = 0;

  // Default initial value of the data sampling time 
  sett->dt = 0.5;

  opts->help_flag=0;
  static int help_flag=0;

  // Default value of the minimal number of coincidences 
  opts->mincoin=3; 

  // Default value of the narrow-down parameter 
  opts->narrowdown=0.5; 

  // Default value of the cell shift: 0000 (no shifts)
  opts->shift=0;

  // Default value of the cell scaling: 1111 (no scaling)
  opts->scale=1111;

  // Default signal-to-noise threshold cutoff
  opts->snrcutoff=6;

  // Reading arguments 

  while (1) {
    static struct option long_options[] = {
      {"help", no_argument, &help_flag, 1},
      // Cell shifts  
      {"shift", required_argument, 0, 's'},
      // Cell scaling 
      {"scale", required_argument, 0, 'z'},
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
      // Minimal number of coincidences recorded in the output  
      {"mincoin", required_argument, 0, 'm'},
      // Narrow down the frequency band (+- the center of band) 
      {"narrowdown", required_argument, 0, 'n'},
      // Signal-to-noise threshold cutoff  
      {"snrcutoff", required_argument, 0, 'c'},
      // number of days in the time-domain segment 
      {"nod", required_argument, 0, 'y'},
      {0, 0, 0, 0}
    };

    if (help_flag) {

      printf("polgraw-allsky periodic GWs: search for concidences among candidates\n");
      printf("Usage: ./coincidences -[switch1] <value1> -[switch2] <value2> ...\n") ;
      printf("Switches are:\n\n");
      printf("-data         Data directory (default is ./candidates)\n");
      printf("-output       Output directory (default is ./coinc-results)\n");
      printf("-shift        Cell shifts in fsda directions (4 digit number, e.g. 0101, default 0000)\n");
      printf("-scale        Cell scaling in fsda directions (4 digit number, e.g. 4824, default 1111)\n");
      printf("-refr         Reference frame number\n");
      printf("-fpo          Reference band frequency fpo value\n");
      printf("-dt           Data sampling time dt (default value: 0.5)\n");
      printf("-trigname     Part of triggers' name (for identifying files)\n");
      printf("-refloc       Location of the reference grid.bin and starting_date files\n");
      printf("-mincoin      Minimal number of coincidences recorded\n");
      printf("-narrowdown   Narrow-down the frequency band (range [0, 0.5] +- around center)\n");
      printf("-nod          Number of days\n");
      printf("-snrcutoff    Signal-to-noise threshold cutoff (default value: 6)\n\n");

      printf("Also:\n\n");
      printf("--help		This help\n");

      exit (0);
    }

    int option_index = 0;
    int c = getopt_long_only (argc, argv, "p:o:d:s:z:r:t:e:g:m:n:c:y:", long_options, &option_index);
    if (c == -1)
      break;

    switch (c) {
    case 'p':
      sett->fpo = atof(optarg);
      break;
    case 's': // Cell shifts 
      opts->shift = atof(optarg);
      break;
    case 'z': // Cell scaling   
      opts->scale = atoi(optarg);
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
    case 'm':
      opts->mincoin = atoi(optarg);
      break;
    case 'n':
      opts->narrowdown = atof(optarg);
      break;
    case 'c':
      opts->snrcutoff = atof(optarg);
      break;
    case 'y':
      sett->nod = atoi(optarg);
      break;
    case '?':
      break;
    default:
      break ;
    } /* switch c */
  } /* while 1 */

  // Putting the parameter in triggers' frequency range [0, pi] 
  opts->narrowdown *= M_PI; 

  // Check if sett->nod was set up, if not, exit
  if(!(sett->nod)) { 
    printf("Number of days not set... Exiting\n"); 
    exit(EXIT_FAILURE); 
  } 

  printf("Number of days is %d\n", sett->nod); 

  printf("The SNR threshold cutoff is %.12f, ", opts->snrcutoff); 
  printf("corresponding to F-statistic value of %.12f\n", 
    pow(opts->snrcutoff, 2)/2. + 2); 

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
    for(i=0; i<4; i++) { 
      for(j=0; j<4; j++) { 
        sett->vedva[i][j]  = eigvec[i][j]*sqrt(eigval[j]); 
//        printf("%.12le ", sett->vedva[i][j]); 
      } 
//      printf("\n"); 
    }
      
  } 

/* 
  //#mb matrix generated in matlab, for tests 
  double _tmp[4][4] = { 
    {-2.8622034614137332e-001, -3.7566564762376159e-002, -4.4001551065376701e-012, -3.4516253934827171e-012}, 
    {-2.9591999145463371e-001, 3.6335210834374479e-002, 8.1252443441098394e-014, -6.8170555119669981e-014}, 
    {1.5497867603229576e-005, 1.9167007413107127e-006, 1.0599051611325639e-008, -5.0379548388381567e-008}, 
    {2.4410008440913992e-005, 3.2886518554938671e-006, -5.7338464150027107e-008, -9.3126913365595100e-009},
  };
  { int i,j; 
  for(i=0; i<4; i++) 
    for(j=0; j<4; j++) 
      sett->vedva[i][j]  = _tmp[i][j]; 
  }
  printf("\n"); 
  { int i, j; 
  for(i=0; i<4; i++) { 
    for(j=0; j<4; j++) {
        printf("%.12le ", sett->vedva[i][j]);
      }
      printf("\n"); 
  } 
 } 
*/ 

  gsl_vector_free (eval);
  gsl_matrix_free (evec);

} // end of manage grid matrix  
