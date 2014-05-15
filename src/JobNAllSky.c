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

#include "auxi.h"
#include "settings.h"

static int white_flag=0; 							//  white noise flag
static int s0_flag=0;								// no spin-down flag
static int help_flag=0; 

/* Default output and data directories */

#ifndef PREFIX
#define PREFIX .
#endif

#ifndef DTAPREFIX
#define DTAPREFIX ../data
#endif

double *cosmodf, *sinmodf, *t2, *aa, *bb, \
	   *shftf, *shft, *xDat, *DetSSB, *F;
complex double *xDatma, *xDatmb;

int
JobNAllSky (int argc, char *argv[]) {
  int i, pm, mm, nn, pst, mst, nst, sst, sgnlc, fd, hemi=0, Nzeros=0,	\
    Ninterp, FNum, nmin, nmax, spndr[2], nr[2], mr[2], pmr[2], c,	\
    ident=0, band=0, nfftf, 
	fftinterp=INT; // default value
  char hostname[32], wfilename[96], filename[64], outname[64], qname[64], 
    prefix[64], dtaprefix[64], label[64], range[64], ifo_choice[2], *wd=NULL;
  double *sgnlv, omrt, coft, epsm, sepsm, cepsm, phir, sphir, cphir, *M, 
	trl=20., // default value for the F-statistic threshold
	fpo, fpo_val, sig2, crf0;
  fftw_complex *xa, *xao, *xDa, *rDa;
  fftw_plan plan, pl_int, pl_inv;
  FILE *wisdom, *state, *data;
  struct flock lck;
  struct stat buffer;

  strcpy (prefix, TOSTR(PREFIX));
  strcpy (dtaprefix, TOSTR(DTAPREFIX));
  label[0] = '\0';
  range[0] = '\0';

  // Initial value of starting frequency 
  // set to a negative quantity. If this is not 
  // changed by the command line value 
  // fpo is calculated from the band number b.
  fpo_val = -1; 

  // Default choice of the detector is Virgo:
  strcpy(ifo_choice, "V1"); 

  while (1) {
    static struct option long_options[] = {
      {"help", no_argument, &help_flag, 1},	
      {"whitenoise", no_argument, &white_flag, 1},
      {"nospindown", no_argument, &s0_flag, 1},
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

  if(help_flag) {

     printf("*** Continuous GW search code using the F-statistic ***\n"); 
     printf("Usage: ./search -[switch1] <value1> -[switch2] <value2> ...\n") ;
     printf("Switches are:\n\n"); 
     printf("-d	Data directory (default is .)\n"); 
     printf("-o	Output directory (default is ./candidates)\n"); 
     printf("-i	Frame number\n"); 
     printf("-b	Band number\n"); 
     printf("-l	Custom label for the input and output files\n"); 
     printf("-r	Grid range filename\n"); 
     printf("-c	Change to directory <dir>\n"); 
     printf("-f	Intepolation method (INT [default] or FFT)\n"); 
     printf("-t	Threshold for the F-statistic (default is 20)\n");
     printf("-h	Hemisphere (default is 0 - does both)\n");
	 printf("-a	Detector (L1, H1 or V1); default is V1\n");     
     printf("-p	fpo (starting frequency) value\n\n");
     printf("Also:\n\n"); 
     printf("--whitenoise	Will assume white Gaussian noise\n"); 
     printf("--nospindown	Will neglect spindowns\n"); 
     printf("--help	This help\n"); 		

     exit (0);

  }

    int option_index = 0;
    c = getopt_long (argc, argv, "i:b:o:d:l:r:c:t:f:h:a:p:", long_options, &option_index);
    if (c == -1)
      break;

    switch (c) {
    case 'i':
      ident = atoi (optarg);
      break;
    case 'f': 
      if(!strcmp(optarg, "FFT")) fftinterp=FFT;
      break;
    case 't':
      trl = atof(optarg);
      break;
    case 'h':
      hemi = atof(optarg);
      break;
    case 'b':
      band = atoi (optarg);
      break;
    case 'o':
      strcpy (prefix, optarg);
      break;
    case 'd':
      strcpy (dtaprefix, optarg);
      break;
    case 'l':
      label[0] = '_';
      strcpy (1+label, optarg);
      break;
    case 'r':
      strcpy (range, optarg);
      break;
    case 'c':
      wd = (char *) malloc (1+strlen(optarg));
      strcpy (wd, optarg);
      break;
    case 'p':
      fpo_val = atof(optarg); 
      break; 	
	case 'a':
	  strcpy (ifo_choice, optarg); 
	  break; 
    case '?':
      break;
    default: break ; 
    } /* switch c */
  } /* while 1 */

  printf ("Data directory is %s\n", dtaprefix);
  printf ("Output directory is %s\n", prefix);
  printf ("Frame number is %d\n", ident);
  printf ("Band number is %d\n", band);    

  if (white_flag)
    printf ("Will assume white Gaussian noise\n");
  if (fftinterp==INT)
    printf ("Will use fftinterp=INT (FFT interpolation by interbinning)\n");
  else 
    printf ("Will use fftinterp=FFT (FFT interpolation by zero-padding)\n");
  if(trl!=20)
    printf ("Threshold for the F-statistic is %lf\n", trl);
  if(hemi)
    printf ("Search for hemisphere %d\n", hemi);
  if (s0_flag)
    printf ("Will assume s_1 = 0.\n");
  if (strlen(label))
    printf ("Will use '%s' as data label\n", label);
  if (strlen(range))
    printf ("Will read grid range from '%s'\n", range);
  if (wd) {
    printf ("Changing working directory to %s\n", wd);
    if (chdir(wd)) {
      perror (wd);
      abort ();
    }
  }

  // Output data handling 
  if (stat(prefix, &buffer) == -1) {
    if (errno == ENOENT) {
      /* Output directory apparently does not exist, try to create one */
      if (mkdir(prefix, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH	\
		| S_IXOTH) == -1) {
	perror (prefix);
	return 1;
      }
    } else { /* can't access output directory */
      perror (prefix);
      return 1;
    }
  }

  M = (double *) calloc (16, sizeof (double));

  sprintf (filename, "%s/%03d/grid.bin", dtaprefix, ident);
  if ((data=fopen (filename, "r")) != NULL) {
    // fftpad: used to zero padding to fftpad*nfft data points 
    fread ((void *)&fftpad, sizeof (int), 1, data);
    // M: vector of 16 components consisting of 4 rows 
    // of 4x4 grid-generating matrix  
    fread ((void *)M, sizeof (double), 16, data);
    fclose (data);
  } else {
    perror (filename);
    return 1;
  }

  // Starting band frequency: 
  // fpo_val is optionally read from the command line
  // Its initial value is set to -1 
  if(fpo_val>=0)
      fpo = fpo_val; 
  else 
  // The usual definition: 
      fpo = 100. + 0.96875 * band;

  // Detector, ephemerides, constants 
  settings (fpo, ifo_choice);
  // Amplitude modulation functions coefficients
  rogcvir ();
  // Establish grid range in which the search will be performed 
  // with the use of the M matrix from grid.bin  
  gridr (M, spndr, nr, mr);

  // Hemispheres (with respect to the ecliptic)
  if(hemi) { 

	pmr[0] = hemi; pmr[1] = hemi; 

  } else { 

  	pmr[0] = 1;
  	pmr[1] = 2;

  }

  // If the parameter range is invoked, the search is performed 
  // within the range of grid parameters from an ascii file 
  // ("-r range_file" from the command line) 
  if (strlen (range)) {
    if ((data=fopen (range, "r")) != NULL) {
      fscanf (data, "%d %d", spndr, 1+spndr);
      fscanf (data, "%d %d", nr, 1+nr);
      fscanf (data, "%d %d", mr, 1+mr);
      fscanf (data, "%d %d", pmr, 1+pmr);
      fclose (data);
    } else {
      perror (range);
      return 1;
    }
  }

  // Allocates and initializes to zero the data, detector ephemeris 
  // and the F-statistic arrays
  xDat = (double *) calloc (N, sizeof (double));
  DetSSB = (double *) calloc (3*N, sizeof (double));
  F = (double *) calloc (2*nfft, sizeof (double));

  // Input time-domain data handling 
  sprintf (filename, "%s/%03d/xdatc_%03d_%03d%s.bin", dtaprefix, ident,	\
	   ident, band, label);
  if ((data = fopen (filename, "r")) != NULL) {
    fread ((void *)(xDat), sizeof (double), N, data);
    fclose (data);
  } else {
    perror (filename);
    return 1;
  }

  // Checking for null values in the data 
  for(i=0; i<N; i++)
	if(!xDat[i]) Nzeros++;

  // factor N/(N - Nzeros) to account for null values in the data 
  crf0 = (double)N/(N-Nzeros); 
  
  if (white_flag)
    sig2 = N*var (xDat, N);
  else
    sig2 = -1.;

  // Ephemeris file handling 
  sprintf (filename, "%s/%03d/DetSSB.bin", dtaprefix, ident);
  if ((data = fopen (filename, "r")) != NULL) {
    // Detector position w.r.t solar system baricenter
    // for every datapoint
    fread ((void *)(DetSSB), sizeof (double), 3*N, data);
    // Deterministic phase defining the position of the Earth 
    // in its diurnal motion at t=0  
    fread ((void *)(&phir), sizeof (double), 1, data);
    // Earth's axis inclination to the ecliptic at t=0  
    fread ((void *)(&epsm), sizeof (double), 1, data);
    fclose (data);
  } else {
    perror (filename);
    return 1;
  }

#ifdef HAVE_SINCOS
  sincos (phir, &sphir, &cphir);
  sincos (epsm, &sepsm, &cepsm);
#else
  sphir = sin (phir);
  cphir = cos (phir);
  sepsm = sin (epsm);
  cepsm = cos (epsm);
#endif
   
  t2 = (double *) calloc (Nv, sizeof (double));
  cosmodf = (double *) calloc (Nv, sizeof (double));
  sinmodf = (double *) calloc (Nv, sizeof (double));
  for (i=0; i<Nv; i++) {
    omrt = omr*i;
    t2[i] = sqr ((double)i);
    cosmodf[i] = cos (omrt);
    sinmodf[i] = sin (omrt);
  }

  // Because of frequency-domain filters, we search 
  // F-statistic in range (nmin+1, nmax) of data points 
  nmin = fftpad*NAV;
  nmax = (nfft/2 - NAV)*fftpad;

  coft = oms; 

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

  aa = (double *) calloc (Nv, sizeof (double));
  bb = (double *) calloc (Nv, sizeof (double));

    shft = (double *) calloc (N, sizeof (double));
    shftf = (double *) calloc (N, sizeof (double));
    xDatma = (complex double *) calloc (N, sizeof (complex double));
    xDatmb = (complex double *) calloc (N, sizeof (complex double));

  xa = fftw_malloc (4 * fftpad * nfft * sizeof (fftw_complex));

  // FFT plans vary w.r.t the method of calculation: 
  // case INT (simple interpolation [interbinning] of a shorter Fourier transform)
  if(fftinterp==INT) { 

	xao = fftw_malloc (2 * nfft * sizeof (fftw_complex));
 
	// Plans a multidimensional DFT, where the input variables are:
	plan = fftw_plan_many_dft 
			(1,	// dimension of the transform (rank) 
			&nfft,	// size of the transform 
			2,	// how many transforms (the k-th transform 
				// is of the array starting at xa + k*idist
				// and xao + k*odist) 
			xa,	// input array
			NULL,	// inembed (defines auxilary array 
				// of size nfft)
			1,	// input stride (the j-th element of 
				// the input array is located at j*istride)
			nfft,	// idist (needed if many transforms) 
			xao,	// output array
			NULL,	// onembed (defines auxilary array 
                                // of size nfft)
			1,	// output stride (the j-th element of 
                                // the output array is located at j*ostride)
			nfft,	// odist (needed if many transforms)
			FFTW_FORWARD,	// sign of the transform 
					// (-1 in the exponent)
			FFTW_MEASURE);	// to get the optimal plan
	// (more info at http://www.fftw.org/fftw3_doc/Advanced-Complex-DFTs.html)

  // Case FFT (longer Fourier transforms, interpolation by zero padding)
  } else { 

	nfftf = fftpad*nfft ;

	// Plans a multidimensional DFT, where the input variables are:
	plan = fftw_plan_many_dft 
                        (1,     // dimension of the transform (rank) 
                        &nfftf, // size of the transform 
                        2,      // how many transforms (the k-th transform 
                                // is of the array starting at xa + k*idist
                                // and xao + k*odist) 
                        xa,     // input array
                        NULL,   // inembed (defines auxilary array 
                                // of size nfft)
                        1,      // input stride (the j-th element of 
                                // the input array is located at j*istride)
                        nfftf,  // idist (needed if many transforms) 
                        xa,     // output array
                        NULL,   // onembed (defines auxilary array 
                                // of size nfft)
                        1,      // output stride (the j-th element of 
                                // the output array is located at j*ostride)
                        nfftf,  // odist (needed if many transforms)
                        FFTW_FORWARD,   // sign of the transform 
                                        // (-1 in the exponent)
                        FFTW_MEASURE);  // to get the optimal plan
  	// (more info at http://www.fftw.org/fftw3_doc/Advanced-Complex-DFTs.html)

 } 


  // These two plans below are used in the resampling 
  // procedure in JobCore() 
  xDa = xa;
  pl_int = fftw_plan_many_dft (1, &nfft, 2, xDa, NULL, 1, nfft, xDa,	\
			       NULL, 1, nfft, FFTW_FORWARD,		\
			       FFTW_MEASURE);

  rDa = xa;
  Ninterp = interpftpad*nfft;
  pl_inv = fftw_plan_many_dft (1, &Ninterp, 2, rDa, NULL, 1, Ninterp,	\
			       rDa, NULL, 1, Ninterp, FFTW_BACKWARD,	\
			       FFTW_MEASURE);

  // Generates a 'wisdom' FFT file if there is none  
  if ((wisdom = fopen (wfilename, "r")) == NULL) {

	wisdom = fopen (wfilename, "w");
	fftw_export_wisdom_to_file (wisdom);
  }
  fclose (wisdom);

  if(hemi)
  	sprintf (qname, "state_%03d_%03d%s_%d.dat", ident, band, label, hemi);
  else 
	sprintf (qname, "state_%03d_%03d%s.dat", ident, band, label);

  // Checkpointing 
  if ((state = fopen (qname, "r")) != NULL) {

    // Scan the state file to get last recorded parameters
    if ((fscanf (state, "%d %d %d %d %d", &pst, &mst, &nst, &sst, &FNum) \
	 ) == EOF) {

	// This means that state file is empty (=end of the calculations)
	fprintf (stderr, "state file empty: ending...\n") ;
	fclose (state);
	return 0;

    } fclose (state);

  // No state file - start from the beginning
  } else {
    pst = pmr[0];
    mst = mr[0];
    nst = nr[0];
    sst = spndr[0];
    FNum = 0;
  } /* if state */

  // Main loops 

  for (pm=pst; pm<=pmr[1]; pm++) {	// loop over hemispheres 
    sprintf (outname, "%s/triggers_%03d_%03d%s_%d.bin", prefix, ident,	\
	     band, label, pm);

    for (mm=mst; mm<=mr[1]; mm++) {	// 2 loops over 
      for (nn=nst; nn<=nr[1]; nn++) {	// sky positions 

	  state = fopen (qname, "w");
	  fprintf (state, "%d %d %d %d %d\n", pm, mm, nn, sst, FNum);
	  fclose (state);

	sgnlv = JobCore(pm,		// hemisphere 
			mm,		// grid 'sky position' 		
			nn,		// other grid 'sky position' 
			sst,		// grid spindown coordinate
			spndr[1],	// spindown range limit 	
			M,		// grid generating matrix 
			DetSSB,		// ephemerides array 
			xDat,		// time-domain input data array 
			N,		// Number of data points
			Ninterp,	// interpftpad*nfft (for resampling)
			nfft,		// size of the FFT transform 
			xDa,		// Array for resampling
			rDa,		// Array for resampling
			xa,		// input array for plan
			xao,		// output array for plan
			pl_int,		// fftw_plan needed for resampling
			pl_inv,		// fftw_plan needed for resampling
					// (inverse transformation)
			plan,		// main fftw_plan 
			nmin,		// nmin+1: first point of 
					// the F-statistic to search for a signal
			nmax,		// nmax: last point of 
                                        // the F-statistic to search for a signal
			sepsm,		// sin(epsm)
			cepsm,		// cos(epsm)
			sphir,		// sin(phi_r)
			cphir,		// cos(phi_r)
			&sgnlc,		// reference to array with the parameters 
					// of the candidate signal 
					// (used below to write to the file) 
			1,		// std output writing flag
			fftinterp,	// interpolation flag 
					// (INT of FFT, see settings.h) 
			&FNum,		// Candidate signal number 
			coft,		// = oms 
			trl,		// F-statistic threshold
			sig2,		// N*(variance of xDat) if white_flag
					// else sig2=-1 
			crf0, 		// factor N/(N - Nzeros)
			s0_flag		// No-spindown flag
			);
	sst = spndr[0];

	/* Add trigger parameters to a file */

	if (sgnlc) {
	  if ((fd = open (outname, O_WRONLY|O_CREAT|O_APPEND,
			  S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH)) < 0) {
	    perror (outname);
	    return 1;
	  }
	  lck.l_type = F_WRLCK;
	  lck.l_whence = 0;
	  lck.l_start = 0L;
	  lck.l_len = 0L;
	  if (fcntl (fd, F_SETLKW, &lck) < 0) perror ("fcntl()");

	  write (fd, (void *)(sgnlv), sgnlc*NPAR*sizeof (double));

	  if (close (fd) < 0) perror ("close()");

	} /* if sgnlc */
	free (sgnlv);
      } /* for nn */
      nst = nr[0];
    } /* for mm */
    mst = mr[0];
  } /* for pm */

    state = fopen (qname, "w");
    fclose (state);

  free (aa);
  free (bb);
  free (shft);
  free (shftf);
  free (xDatma);
  free (xDatmb);
  free (xDat);
  free (DetSSB);
  free (F);
  
  free (t2);
  free (cosmodf);
  free (sinmodf);
  fftw_free (xa);
  // This array is defined only in the case of fftinterp==INT
  if(fftinterp==INT) fftw_free (xao);
  free (M);
  
  return 0;
} /* JobNAllSky() */
