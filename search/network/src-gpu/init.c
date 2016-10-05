// MSVC macro to include constants, such as M_PI (include before math.h)
#define _USE_MATH_DEFINES

// Standard C includes
#include <stdio.h>
#include <stdlib.h>
//#include <unistd.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <getopt.h>
#include <time.h>

#include <direct.h>


//#include "cuda_error.h"
#include "init.h"
#include "struct.h"
#include "settings.h"
#include "auxi.h"
//#include "kernels.h"
#include "spline_z.h"

//#include <cuda_runtime_api.h>
//#include <cufft.h>

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

  opts->label[0]  = '\0';
  opts->range[0]  = '\0';
  opts->getrange[0] = '\0';
  opts->usedet[0]   = '\0';
  opts->addsig[0] = '\0';
	
  // Initial value of starting frequency set to a negative quantity. 
  // If this is not changed by the command line value, fpo is calculated 
  // from the band number b (fpo = fpo = fstart + 0.96875*b/(2dt))
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
      // write full grid range to file
      {"getrange", required_argument, 0, 'g'},
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
      // which detectors to use
      {"usedet", required_argument, 0, 'u'}, 
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
      printf("-r, -range        Use file with grid range or pulsar position\n");
      printf("-g, -getrange     Write grid ranges & exit (ignore -r)\n");
      printf("-c, -cwd          Change to directory <dir>\n");
      printf("-t, -threshold    Threshold for the F-statistic (default is 20)\n");
      printf("-h, -hemisphere   Hemisphere (default is 0 - does both)\n");
      printf("-p, -fpo          Reference band frequency fpo value\n");
      printf("-s, -dt           data sampling time dt (default value: 0.5)\n");
      printf("-u, -usedet       Use only detectors from string (default is use all available)\n");
      printf("-x, -addsig       Add signal with parameters from <file>\n\n");


      printf("Also:\n\n");
      printf("--whitenoise      White Gaussian noise assumed\n");
      printf("--nospindown      Spindowns neglected\n");
      printf("--nocheckpoint    State file won't be created (no checkpointing)\n");
      printf("--help            This help\n");

      exit(EXIT_SUCCESS);
    }

    int option_index = 0;
    int c = getopt_long_only(argc, argv, "i:b:o:d:l:r:g:c:t:h:p:x:s:u:", 
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
    case 's':
      sett->dt = atof(optarg);
      break;
    case 'u':
      strcpy(opts->usedet, optarg);
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

    // The usual definition (multiplying the offset by B=1/(2dt))
    // !!! in RDC_O1 the fstart equals 10, not 100 like in VSR1 !!! 
    // 
    sett->fpo = 10. + 0.96875*opts->band*(0.5/sett->dt);

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
  if (opts->wd) {
    printf ("Changing working directory to %s\n", opts->wd);
#ifdef WIN32
    if (_chdir(opts->wd)) { perror(opts->wd); abort(); }
#else
    if (chdir(opts->wd)) { perror (opts->wd); abort (); }
#endif
  }

} // end of command line options handling 


/* Generate grid from the M matrix (grid.bin) 
 */ 

void read_grid( Search_settings *sett, 
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
    //CudaSafeCall( cudaHostAlloc((void **)&(ifo[i].sig.xDat), sett->N*sizeof(double), 
	//			cudaHostAllocMapped) );
    //CudaSafeCall( cudaHostGetDevicePointer((void **)&(ifo[i].sig.xDat_d), 
	//				   (void *)ifo[i].sig.xDat, 0) );

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

    //CudaSafeCall( cudaHostAlloc((void **)&(ifo[i].sig.DetSSB), 3*sett->N*sizeof(double), 
	//			cudaHostAllocMapped) );
    //CudaSafeCall( cudaHostGetDevicePointer((void **)&(ifo[i].sig.DetSSB_d), 
	//				   (void *)ifo[i].sig.DetSSB, 0) );


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

    //CudaSafeCall( cudaMalloc((void**)&ifo[i].sig.xDatma_d,
	//			 sizeof(cufftDoubleComplex)*sett->N) );
    //CudaSafeCall( cudaMalloc((void**)&ifo[i].sig.xDatmb_d, 
	//			 sizeof(cufftDoubleComplex)*sett->N) );
    //
    //CudaSafeCall( cudaMalloc((void**)&(ifo[i].sig.aa_d), 
	//			 sizeof(double)*sett->N) );
    //CudaSafeCall( cudaMalloc((void**)&(ifo[i].sig.bb_d), 
	//			 sizeof(double)*sett->N) );
    //
    //CudaSafeCall( cudaMalloc((void**)&(ifo[i].sig.shft_d), 
	//			 sizeof(double)*sett->N) );
    //CudaSafeCall( cudaMalloc((void**)&(ifo[i].sig.shftf_d), 
	//			 sizeof(double)*sett->N) );

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
  //CudaSafeCall ( cudaMalloc((void **)F_d, 2*sett->nfft*sizeof(double)));

  // Auxiliary arrays, Earth's rotation
 
//  CudaSafeCall( cudaMalloc((void**)&(aux_arr->t2_d),
//			   sizeof(double)*sett->N) );
//  CudaSafeCall( cudaMalloc((void**)&(aux_arr->cosmodf_d), 
//			   sizeof(double)*sett->N) );
//  CudaSafeCall( cudaMalloc((void**)&(aux_arr->sinmodf_d), 
//			   sizeof(double)*sett->N) );
//
//  CudaSafeCall( cudaMalloc((void**)&(aux_arr->tshift_d),
//			   sizeof(double)*sett->N) );
//
//  init_spline_matrices(&aux_arr->diag_d, &aux_arr->ldiag_d, &aux_arr->udiag_d, 
//		       &aux_arr->B_d, sett->Ninterp);
//
//  compute_sincosmodf<<<sett->N/256+1,256>>>(aux_arr->sinmodf_d, aux_arr->cosmodf_d, 
//					    sett->omr, sett->N);

} // end of init arrays 


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
      
      if (aqq != 8) {
	printf("Error when reading range file!\n");
	exit(EXIT_FAILURE);
      }
      
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

void plan_fft (Search_settings *sett, 
	       // Command_line_opts *opts,
	       FFT_plans *plans, 
	       FFT_arrays *fft_arr
	       // Aux_arrays *aux_arr
	       ) {

  //  sett->Ninterp = sett->interpftpad*sett->nfft; //moved to init_arrays

  fft_arr->arr_len = (sett->fftpad*sett->nfft > sett->Ninterp 
		       ? sett->fftpad*sett->nfft : sett->Ninterp);

  //CudaSafeCall ( cudaMalloc((void **)&fft_arr->xa_d, 2*fft_arr->arr_len*sizeof(cufftDoubleComplex)) );
  fft_arr->xb_d = fft_arr->xa_d + fft_arr->arr_len;

  //  sett->nfftf = sett->fftpad*sett->nfft; // moved to init_arrays

  // no need for plans '2' - dimaensions are the same
  //cufftPlan1d( &(plans->plan), sett->nfftf, CUFFT_Z2Z, 1);
  //cufftPlan1d( &(plans->pl_int), sett->nfft, CUFFT_Z2Z, 1);
  //cufftPlan1d( &(plans->pl_inv), sett->Ninterp, CUFFT_Z2Z, 1);
  //
  //CudaSafeCall ( cudaMalloc((void **)&fft_arr->xa_d, 2*fft_arr->arr_len*sizeof(cufftDoubleComplex)) );

}


/* Checkpointing */

void read_checkpoints(Command_line_opts *opts, 
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
    //CudaSafeCall( cudaFreeHost(ifo[i].sig.xDat) );
    //CudaSafeCall( cudaFreeHost(ifo[i].sig.DetSSB) );
    //
    //CudaSafeCall( cudaFree(ifo[i].sig.xDatma_d) );
    //CudaSafeCall( cudaFree(ifo[i].sig.xDatmb_d) );
    //
    //CudaSafeCall( cudaFree(ifo[i].sig.aa_d) );
    //CudaSafeCall( cudaFree(ifo[i].sig.bb_d) );
    //
    //CudaSafeCall( cudaFree(ifo[i].sig.shft_d) );
    //CudaSafeCall( cudaFree(ifo[i].sig.shftf_d) );
  } 

  //CudaSafeCall( cudaFree(aux->cosmodf_d) );
  //CudaSafeCall( cudaFree(aux->sinmodf_d) );
  //CudaSafeCall( cudaFree(aux->t2_d) );
  //
  //CudaSafeCall( cudaFree(F_d) );
  //
  //CudaSafeCall( cudaFree(fft_arr->xa_d) );

  free(sett->M);

  //cufftDestroy(plans->plan);
  //cufftDestroy(plans->pl_int);
  //cufftDestroy(plans->pl_inv);

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
      printf("-snrcutoff    Signal-to-noise threshold cutoff (default value: 6)\n\n");

      printf("Also:\n\n");
      printf("--help		This help\n");

      exit (0);
    }

    int option_index = 0;
    int c = getopt_long_only (argc, argv, "p:o:d:s:z:r:t:e:g:m:n:c:", long_options, &option_index);
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
    case '?':
      break;
    default:
      break ;
    } /* switch c */
  } /* while 1 */

  // Putting the parameter in triggers' frequency range [0, pi] 
  opts->narrowdown *= M_PI; 

  printf("#mb add info at the beginning...\n"); 
  printf("The SNR threshold cutoff is %.12f, ", opts->snrcutoff); 
  printf("corresponding to F-statistic value of %.12f\n", 
    pow(opts->snrcutoff, 2)/2. + 2); 

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
  //int dev, deviceCount = 0;
  //cudaDeviceProp deviceProp;
  
//  if (cudaGetDeviceCount(&deviceCount) != cudaSuccess) {
//    printf("ERROR: cudaGetDeviceCount FAILED CUDA Driver and Runtime version may be mismatched.\n");
//    return(-1);
//  }
//  if (deviceCount == 0) {
//    printf("ERROR: There is no device supporting CUDA\n");
//    return(-1);
//  }
//  if (cdev < 0 && cdev >= deviceCount) {
//    printf("\nWARNING: Device %d is not available! Trying device 0\n", cdev);
//    cdev = 0;
//  }
//
//  printf("__________________________________CUDA devices___________________________________\n");
//  printf("Set | ID |        Name        |   Gmem(B)   | Smem(B) | Cmem(B) | C.Cap. | Thr/bl |\n");
//  
//  for (dev = 0; dev < deviceCount; ++dev) {
//    cudaGetDeviceProperties(&deviceProp, dev);
//    if (deviceProp.major == 9999 && deviceProp.minor == 9999) {
//      printf("- | %1d | %16s | Error | Error | Error | Error | Error |\n", dev, deviceProp.name );
//      if ( dev==cdev ) {
//	printf("ERROR: Can't set device %d\n", cdev);
//	return(-1);
//      }
//    }
//    if (dev==cdev) {
//      printf(" *  |");
//      cudaSetDevice(cdev);
//    } else {
//      printf("    |");
//    }
//    printf(" %1d  | %18.18s | %11Zu | %7Zu | %7Zu |   %d.%d  | %6d |\n", 
//	   dev, deviceProp.name, deviceProp.totalGlobalMem, deviceProp.sharedMemPerBlock, 
//	   deviceProp.totalConstMem, deviceProp.major, deviceProp.minor, deviceProp.maxThreadsPerBlock );
//  }
//  printf("---------------------------------------------------------------------------------\n");
//  
//  /* enable mapped memory */
//  cudaSetDeviceFlags(cudaDeviceMapHost);
//
//  /* force initialization */
//  cudaThreadSynchronize();
  return(cdev);
}
