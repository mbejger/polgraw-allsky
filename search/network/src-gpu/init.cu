#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <getopt.h>
#include <time.h>

#include "cuda_error.h"
#include "struct.h"
#include "init.h"
#include "settings.h"
#include "auxi.h"
#include "kernels.h"
#include "spline_z.h"

#include <cuda_runtime_api.h>
#include <cufft.h>

/*  Command line options handling: search  */ 

void handle_opts( Search_settings *sett, 
		  Command_line_opts *opts,
		  int argc, char* argv[]) {
  
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

      printf("polgraw-allsky periodic GWs: search for candidate signals with the F-statistic\n");
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
    // The usual definition (multiplying the offset by B=1/(2dt) ):
    sett->fpo = 100. + 0.96875*opts->band*(0.5/sett->dt);

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


/* Generate grid from the M matrix (grid.bin) */ 

void read_grid( Search_settings *sett, 
	        Command_line_opts *opts) {

  sett->M = (double *) calloc (16, sizeof (double));

  FILE *data;
  char filename[512];
  sprintf (filename, "%s/%03d/grid.bin", opts->dtaprefix, opts->ident);
  if ((data=fopen (filename, "r")) != NULL) {
    fread ((void *)&sett->fftpad, sizeof(int), 1, data);

    printf("Using fftpad from the grid file: %d\n", sett->fftpad); 
	
    // M: vector of 16 components consisting of 4 rows
    // of 4x4 grid-generating matrix
    fread ((void *)sett->M, sizeof(double), 16, data);
    fclose (data);
  } else {
    perror (filename);
    exit(EXIT_FAILURE);
  }

} // end of read grid 


 /* Array initialization */ 

void init_arrays( Search_settings *sett, 
		  Command_line_opts *opts,
		  Aux_arrays *aux_arr,
		  double **F_d) {

  int i, status; 
  // Allocates and initializes to zero the data, detector ephemeris
  // and the F-statistic arrays

  FILE *data;

  sett->Ninterp = sett->interpftpad*sett->nfft;
  sett->nfftf = sett->fftpad*sett->nfft;


  for(i=0; i<sett->nifo; i++) { 

    /// ifo[i].sig.xDat = (double *) calloc(sett->N, sizeof(double));
    /// mapped memory works for CUDART_VERSION >= 2020
    /// we should test if it's available, if not copy data explicitly to device
    CudaSafeCall( cudaHostAlloc((void **)&(ifo[i].sig.xDat), sett->N*sizeof(double), 
				cudaHostAllocMapped) );
    CudaSafeCall( cudaHostGetDevicePointer((void **)&(ifo[i].sig.xDat_d), 
					   (void *)ifo[i].sig.xDat, 0) );

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

    //ifo[i].sig.DetSSB = (double *) calloc(3*sett->N, sizeof(double));
    CudaSafeCall( cudaHostAlloc((void **)&(ifo[i].sig.DetSSB), 3*sett->N*sizeof(double), 
				cudaHostAllocMapped) );
    CudaSafeCall( cudaHostGetDevicePointer((void **)&(ifo[i].sig.DetSSB_d), 
					   (void *)ifo[i].sig.DetSSB, 0) );


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
      printf("phir=%f\n", ifo[i].sig.phir );
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


    /*    ifo[i].sig.xDatma = 
      (complex double *) calloc(sett->N, sizeof(complex double));
    ifo[i].sig.xDatmb = 
      (complex double *) calloc(sett->N, sizeof(complex double));

    ifo[i].sig.aa = (double *) calloc(sett->N, sizeof(double));
    ifo[i].sig.bb = (double *) calloc(sett->N, sizeof(double));

    ifo[i].sig.shft = (double *) calloc(sett->N, sizeof(double));
    ifo[i].sig.shftf = (double *) calloc(sett->N, sizeof(double));
    */
    CudaSafeCall( cudaMalloc((void**)&ifo[i].sig.xDatma_d,
				 sizeof(cufftDoubleComplex)*sett->N) );
    CudaSafeCall( cudaMalloc((void**)&ifo[i].sig.xDatmb_d, 
				 sizeof(cufftDoubleComplex)*sett->N) );

    CudaSafeCall( cudaMalloc((void**)&(ifo[i].sig.aa_d), 
				 sizeof(double)*sett->N) );
    CudaSafeCall( cudaMalloc((void**)&(ifo[i].sig.bb_d), 
				 sizeof(double)*sett->N) );

    CudaSafeCall( cudaMalloc((void**)&(ifo[i].sig.shft_d), 
				 sizeof(double)*sett->N) );
    CudaSafeCall( cudaMalloc((void**)&(ifo[i].sig.shftf_d), 
				 sizeof(double)*sett->N) );

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

  //  *F = (double *) calloc(2*sett->nfft, sizeof(double));
  CudaSafeCall ( cudaMalloc((void **)F_d, 2*sett->nfft*sizeof(double)));

  // Auxiliary arrays, Earth's rotation
  //aux_arr->t2 = (double *) calloc(sett->N, sizeof (double));
  //aux_arr->cosmodf = (double *) calloc(sett->N, sizeof (double));
  //aux_arr->sinmodf = (double *) calloc(sett->N, sizeof (double));

  CudaSafeCall( cudaMalloc((void**)&(aux_arr->t2_d),
			   sizeof(double)*sett->N) );
  CudaSafeCall( cudaMalloc((void**)&(aux_arr->cosmodf_d), 
			   sizeof(double)*sett->N) );
  CudaSafeCall( cudaMalloc((void**)&(aux_arr->sinmodf_d), 
			   sizeof(double)*sett->N) );

  CudaSafeCall( cudaMalloc((void**)&(aux_arr->tshift_d),
			   sizeof(double)*sett->N) );

  init_spline_matrices(&aux_arr->diag_d, &aux_arr->ldiag_d, &aux_arr->udiag_d, 
		       &aux_arr->B_d, sett->Ninterp);

  compute_sincosmodf<<<sett->N/256+1,256>>>(aux_arr->sinmodf_d, aux_arr->cosmodf_d, 
					    sett->omr, sett->N);
  /*
  double omrt;

  for (i=0; i<sett->N; i++) {
    omrt = (sett->omr)*(double)i;     // Earth angular velocity * dt * i
    printf(">>> omrt=%f\n", omrt);
    (aux_arr->t2_d)[i] = sqr((double)i);
    aux_arr->cosmodf_d[i] = cos(omrt);
    aux_arr->sinmodf_d[i] = sin(omrt);
  }
  */
} // end of init arrays 


  /* Add signal to data */ 

/// disabled in gpu version - would need to rewrite for gpu 
/// due to sig.aa, sig.bb dependency
#if 0
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
	//  thsnr   += pow(signadd, 2.);
      }	 
    }

  }

} 
#endif

   /* Search range */ 

void set_search_range( Search_settings *sett, 
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



void plan_fft (Search_settings *sett, 
	       // Command_line_opts *opts,
	       FFT_plans *plans, 
	       FFT_arrays *fft_arr
	       // Aux_arrays *aux_arr
	       ) {

  //  sett->Ninterp = sett->interpftpad*sett->nfft; //moved to init_arrays

  fft_arr->arr_len = (sett->fftpad*sett->nfft > sett->Ninterp 
		       ? sett->fftpad*sett->nfft : sett->Ninterp);

  CudaSafeCall ( cudaMalloc((void **)&fft_arr->xa_d, 2*fft_arr->arr_len*sizeof(cufftDoubleComplex)) );
  fft_arr->xb_d = fft_arr->xa_d + fft_arr->arr_len;

  //  sett->nfftf = sett->fftpad*sett->nfft; // moved to init_arrays

  // no need for plans '2' - dimaensions are the same
  cufftPlan1d( &(plans->plan), sett->nfftf, CUFFT_Z2Z, 1);
  cufftPlan1d( &(plans->pl_int), sett->nfft, CUFFT_Z2Z, 1);
  cufftPlan1d( &(plans->pl_inv), sett->Ninterp, CUFFT_Z2Z, 1);

  CudaSafeCall ( cudaMalloc((void **)&fft_arr->xa_d, 2*fft_arr->arr_len*sizeof(cufftDoubleComplex)) );

}


/* Checkpointing */

void read_checkpoints(Command_line_opts *opts, 
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


   /* Cleanup & memory free */

void cleanup(
	     Search_settings *sett,
	     Command_line_opts *opts,
	     Search_range *s_range,
	     FFT_plans *plans,
	     FFT_arrays *fft_arr,
	     Aux_arrays *aux,
	     double *F_d) {

  int i; 
  
  for(i=0; i<sett->nifo; i++) {
    CudaSafeCall( cudaFreeHost(ifo[i].sig.xDat) );
    CudaSafeCall( cudaFreeHost(ifo[i].sig.DetSSB) );

    CudaSafeCall( cudaFree(ifo[i].sig.xDatma_d) );
    CudaSafeCall( cudaFree(ifo[i].sig.xDatmb_d) );

    CudaSafeCall( cudaFree(ifo[i].sig.aa_d) );
    CudaSafeCall( cudaFree(ifo[i].sig.bb_d) );

    CudaSafeCall( cudaFree(ifo[i].sig.shft_d) );
    CudaSafeCall( cudaFree(ifo[i].sig.shftf_d) );
  } 

  CudaSafeCall( cudaFree(aux->cosmodf_d) );
  CudaSafeCall( cudaFree(aux->sinmodf_d) );
  CudaSafeCall( cudaFree(aux->t2_d) );

  CudaSafeCall( cudaFree(F_d) );

  CudaSafeCall( cudaFree(fft_arr->xa_d) );

  free(sett->M);

  cufftDestroy(plans->plan);
  cufftDestroy(plans->pl_int);
  cufftDestroy(plans->pl_inv);

} // end of cleanup & memory free 



 /* Command line options handling: coincidences  */ 

void handle_opts_coinc( Search_settings *sett, 
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

  // Default value of the minimal number of coincidences 
  opts->mincoin=3; 

  // Default value of the narrow-down parameter 
  opts->narrowdown=0.5; 

  // Default value of the cell shift: 0000 (no shifts)
  opts->shift=0;

  // Default value of the cell scaling: 1111 (no scaling)
  opts->scale=1111;

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
      printf("-narrowdown   Narrow-down the frequency band (range [0,1], +- around center)\n\n");

      printf("Also:\n\n");
      printf("--help		This help\n");

      exit (0);
    }

    int option_index = 0;
    int c = getopt_long_only (argc, argv, "p:o:d:s:z:r:t:e:g:m:n:", long_options, &option_index);
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
    case '?':
      break;
    default:
      break ;
    } /* switch c */
  } /* while 1 */

  // Putting the parameter in triggers' frequency range [0, pi] 
  opts->narrowdown *= M_PI; 

  printf("#mb add info at the beginning...\n"); 

} // end of command line options handling: coincidences  



#if 0
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

#endif

/*---------------------------------------------------------------------------*/

/*
  Initialize CUDA: cuinit
  - sets cuda device to (in priority order): cdev, 0 
  - returns: device id or -1 on error
*/
int cuinit(int cdev)
{
  int dev, deviceCount = 0;
  cudaDeviceProp deviceProp;
  
  if (cudaGetDeviceCount(&deviceCount) != cudaSuccess) {
    printf("ERROR: cudaGetDeviceCount FAILED CUDA Driver and Runtime version may be mismatched.\n");
    return(-1);
  }
  if (deviceCount == 0) {
    printf("ERROR: There is no device supporting CUDA\n");
    return(-1);
  }
  if (cdev < 0 && cdev >= deviceCount) {
    printf("\nWARNING: Device %d is not available! Trying device 0\n", cdev);
    cdev = 0;
  }

  printf("__________________________________CUDA devices___________________________________\n");
  printf("Set | ID |        Name        |   Gmem(B)   | Smem(B) | Cmem(B) | C.Cap. | Thr/bl |\n");
  
  for (dev = 0; dev < deviceCount; ++dev) {
    cudaGetDeviceProperties(&deviceProp, dev);
    if (deviceProp.major == 9999 && deviceProp.minor == 9999) {
      printf("- | %1d | %16s | Error | Error | Error | Error | Error |\n", dev, deviceProp.name );
      if ( dev==cdev ) {
	printf("ERROR: Can't set device %d\n", cdev);
	return(-1);
      }
    }
    if (dev==cdev) {
      printf(" *  |");
      cudaSetDevice(cdev);
    } else {
      printf("    |");
    }
    printf(" %1d  | %18.18s | %11Zu | %7Zu | %7Zu |   %d.%d  | %6d |\n", 
	   dev, deviceProp.name, deviceProp.totalGlobalMem, deviceProp.sharedMemPerBlock, 
	   deviceProp.totalConstMem, deviceProp.major, deviceProp.minor, deviceProp.maxThreadsPerBlock );
  }
  printf("---------------------------------------------------------------------------------\n");
  
  /* enable mapped memory */
  cudaSetDeviceFlags(cudaDeviceMapHost);

  /* force initialization */
  cudaThreadSynchronize();
  return(cdev);
}
