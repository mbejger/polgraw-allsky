#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include "auxi.h"
#include "settings.h" 

// Dummy device code
__global__ void Multiply(double* A, double* B, int N) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < N)
        B[i] = A[i]*A[i];
}

static int white_flag=0;                            //  white noise flag
static int s0_flag=0;                               // no spin-down flag
static int help_flag=0;

#ifndef PREFIX
#define PREFIX .
#endif

#ifndef DTAPREFIX
#define DTAPREFIX .
#endif

// Host code
int main(int argc, char *argv[]) {

	int i, pm, mm, nn, pst, mst, nst, sst, sgnlc, fd, hemi=0, Nzeros=0,   \
    	Ninterp, FNum, nmin, nmax, c,       \
    	ident=0, band=0, nfftf,
    	fftinterp=INT; // default value

	double *sgnlv, omrt, coft, epsm, phir, *M, 
	trl=20.,        // default value for the F-statistic threshold
	fpo, sig2;

	char prefix[64], dtaprefix[64], filename[64], label[64], range[64], *ifo = "V", *wd = NULL; 
	FILE *data;

    strcpy (prefix, TOSTR(PREFIX));
    strcpy (dtaprefix, TOSTR(DTAPREFIX));
	label[0] = '\0';
  	range[0] = '\0';

	while (1) {
		static struct option long_options[] = {

			{"help", no_argument, &help_flag, 1},
		    {"whitenoise", no_argument, &white_flag, 1},
      	 	{"nospindown", no_argument, &s0_flag, 1},
	        // frequency band number 
      	 	{"band", required_argument, 0, 'b'},
	        // change directory parameter
      	 	{"cwd", required_argument, 0, 'c'},
	 	    // input data directory 
	        {"data", required_argument, 0, 'd'},
	        // interpolation method
      		{"int/fft", required_argument, 0, 'f'},
	        // hemisphere
      		{"hemisphere", required_argument, 0, 'h'},
	        // frame number
      		{"ident", required_argument, 0, 'i'},
	        // non-standard label for naming files 
      		{"label", required_argument, 0, 'l'},
	        // output directory
      		{"output", required_argument, 0, 'o'},
	        // narrower grid range parameter file
      		{"range", required_argument, 0, 'r'},
	        // interpolation method
      		{"threshold", required_argument, 0, 't'},
	        // the detector
      		{"detector", required_argument, 0, 'x'},
      		{0, 0, 0, 0}
    	};

		if(help_flag) {

			printf("*** Continuous GW search code using the F-statistic ***\n");
		    printf("Usage: ./search -[switch1] <value1> -[switch2] <value2> ...\n") ;
    	printf("Switches are:\n\n");
		    printf("-b     Band number\n");
		    printf("-c     Change to directory <dir>\n");
		    printf("-d     Data directory (default is .)\n");
    	printf("-f     Intepolation method (INT [default] or FFT)\n");
		    printf("-h     Hemisphere (default is 0 - does both)\n");
		    printf("-i     Frame number\n");
		    printf("-l     Custom label for the input and output files\n");
		    printf("-o     Output directory (default is ./candidates)\n");
		    printf("-r     Grid range filename\n");
		    printf("-t     Threshold for the F-statistic (default is 20)\n\n");
		    printf("-x     Detector (H, L or V; default is V)\n");
		    printf("Also:\n");
		    printf("--whitenoise   Will assume white Gaussian noise\n");
		    printf("--nospindown   Will assume that spindown is not important\n");
		    printf("--help     This help\n");

    	exit (0);

		}

	    int option_index = 0;
	    c = getopt_long (argc, argv, "b:c:d:f:h:i:l:o:r:t:x:", \
			long_options, &option_index);
	
	    if (c == -1) break;
	
		switch (c) {
			case 'b':
				band = atoi (optarg);
				break;
			case 'c':
				wd = (char *) malloc (1+strlen(optarg));
				strcpy (wd, optarg);
				break;
			case 'd':
				strcpy (dtaprefix, optarg);
				break;
			case 'f':
				if(!strcmp(optarg, "FFT")) fftinterp=FFT;
				break;
			case 'h':
				hemi = atof(optarg);
				break;
			case 'i':
				ident = atoi (optarg);
				break;
			case 'l':
				label[0] = '_';
				strcpy (1+label, optarg);
				break;
			case 'o':
				strcpy (prefix, optarg);
				break;
			case 'r':
				strcpy (range, optarg);
				break;
			case 't':
				trl = atof(optarg);
				break;
			case 'x':
				ifo = (char *) malloc (1+strlen(optarg));
				strcpy (ifo, optarg);
				break;
			case '?':
				break;
			default: 
				break;
	
		} // switch c

	} // while 1

	// Screen info printouts
	//----------------------

	// Detector choice 
	switch(ifo[0]) {
    	case 'V':
        	printf ("The detector is Virgo\n");
	        break; 
    	case 'H': 
        	printf ("The detector is LIGO Hanford\n");
	        break;  
    	case 'L':
        	printf ("The detector is LIGO Livingston\n");
	        break; 
    	default: 
        	printf ("Unknown detector %s! Cannot continue...\n", ifo);
        	exit(0); 
  	}

	printf ("Data directory is %s\n", dtaprefix);
	printf ("Output directory is %s\n", prefix);
	printf ("Frame number: %d, band number: %d\n", ident, band);

	if (strlen(label))
		printf ("Using '%s' as data label; datafile is xdat_%02d_%03d%s.bin\n",
		label, ident, band, label);
	else
		printf ("Datafile is xdat_%02d_%03d%s.bin\n", ident, band, label);

	if (white_flag)
		printf ("Assuming white Gaussian noise\n");
	if (fftinterp==INT)
		printf ("fftinterp=INT (FFT interpolation by interbinning)\n");
	else
		printf ("fftinterp=FFT (FFT interpolation by zero-padding)\n");
	if(trl!=20)
		printf ("Threshold for the F-statistic is %lf\n", trl);
	if(hemi)
		printf ("Search for hemisphere %d\n", hemi);
	if (s0_flag)
		printf ("Will assume s_1 = 0.\n");
	if (strlen(label))
		printf ("Will use '%s' as datat label\n", label);
	if (strlen(range))
		printf ("Reading grid range from '%s'\n", range);

	if (wd) {
		printf ("Changing working directory to %s\n", wd);
		if (chdir(wd)) {
			perror (wd);
			abort ();
		}
	}


	// Grid-generating matrix
	//-----------------------

	M = (double *) calloc (16, sizeof (double));
	sprintf (filename, "%s/%02d/grid.bin", dtaprefix, ident);

	if ((data=fopen (filename, "r")) != NULL) {

		// fftpad: used to zero padding to fftpad*nfft data points 
		fread ((void *)&pars.fftpad, sizeof (int), 1, data);
		// M: vector of 16 components consisting of 4 rows 
		// of 4x4 grid-generating matrix 
		//#mb add to pars  
		fread ((void *)M, sizeof (double), 16, data);
		fclose (data);
		} else {
		perror (filename);
		return 1;
	}
	
	// Starting band frequency
	fpo = 100. + 0.96875 * band;
	// Detector, ephemerides, constants, amplitude modulation functions coefficients
	set_search_parameters (fpo, ifo);
	// Establish grid range in which the search will be performed 
	// with the use of the M matrix from grid.bin  
	gridr (M, pars.spndr, pars.nr, pars.mr);
	
	// Hemispheres (with respect to the ecliptic)
	if(hemi) {
	
		pars.pmr[0] = hemi; pars.pmr[1] = hemi;
	
	} else {
	
		pars.pmr[0] = 1;
		pars.pmr[1] = 2;
	
	}
	
	// If the parameter range is invoked, the search is performed 
	// within the range of grid parameters from an ascii file 
	// ("-r range_file" from the command line) 
	if (strlen (range)) {

		if ((data=fopen (range, "r")) != NULL) {

			fscanf (data, "%d %d", pars.spndr, 1+pars.spndr);
			fscanf (data, "%d %d", pars.nr, 1+pars.nr);
			fscanf (data, "%d %d", pars.mr, 1+pars.mr);
			fscanf (data, "%d %d", pars.pmr, 1+pars.pmr);
			fclose (data);
		} else {

		perror (range);
		return 1;
		}
	}

	set_search_parameters(100, ifo); 

/* //#mb just for tests
	printf("%d %lf %le\n", pars.N, pars.c1, pars.Smax) ; 
	printf("%lf %d %d %d %d %d %d %d %d\n", fpo, pars.spndr[0], pars.spndr[1], 
		pars.nr[0], pars.nr[1], pars.mr[0], pars.mr[1], pars.pmr[0], pars.pmr[1]) ; 
*/ 


	int N = pars.N; 
    size_t size = N * sizeof(double);
	double *xDat; 

	// Allocate memory for input data array on host 
	xDat = (double *) malloc(size);

	// Input time-domain data handling 
  	sprintf (filename, "%s/%02d/xdat_%02d_%03d%s.bin", dtaprefix, ident,	\
	   ident, band, label);

	if ((data = fopen (filename, "r")) != NULL) {
    	fread ((void *)(xDat), sizeof (double), N, data);
    	fclose (data);
  	} else {
    	perror (filename);
    	return 1;
  	}

/*
  // Checking for null values in the data 
  //#mb to be done on the device? 
  for(i=0; i<N; i++) 
	if(!xDat[i])
		Nzeros++; 

  // In case of white noise, factor N/(N - Nzeros) 
  // accounts for null values in the data  
  //#mb to be done on the device?
  if (white_flag) sig2 = N*var (xDat, N)*N/(N - Nzeros);
  else sig2 = -1.;
*/

	double *DetSSB; 
      
	// Allocate memory for input data array on host 
	DetSSB = (double *) malloc(size);

	// Ephemeris file handling 
	sprintf (filename, "%s/%02d/DetSSB.bin", dtaprefix, ident);
	if ((data = fopen (filename, "r")) != NULL) {
		// Detector position w.r.t solar system baricenter for every datapoint
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

//#mb check if sincos is really worth using 
#ifdef HAVE_SINCOS
	sincos (phir, &pars.sphir, &pars.cphir);
	sincos (epsm, &pars.sepsm, &pars.cepsm);
#else
	pars.sphir = sin (phir);
	pars.cphir = cos (phir);
	pars.sepsm = sin (epsm);
	pars.cepsm = cos (epsm);
#endif


	//#mb maybe it will be better to do this 
	// on the device - to be checked
	// (is done only once) 
	double *t2, *cosmodf, *sinmodf; 

	t2 = (double *) calloc (N, sizeof (double));
	cosmodf = (double *) calloc (N, sizeof (double));
	sinmodf = (double *) calloc (N, sizeof (double));
	
	for (i=0; i<N; i++) {
		omrt = pars.omr*i;
		t2[i] = sqr ((double)i);
		cosmodf[i] = cos (omrt);
		sinmodf[i] = sin (omrt);
	}

	for(int i=0; i<10; i++) printf("%lf\n", xDat[i]); 
	
	// Allocate memory for input data array on device 
	double *xDat_d;
    cudaMalloc(&xDat_d, size);

	double *out_d;
    cudaMalloc(&out_d, size);

	// Copy vectors from host memory to device memory
    cudaMemcpy(xDat_d, xDat, size, cudaMemcpyHostToDevice);

	// Invoke kernel
    int threadsPerBlock = 256;
    int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
    Multiply<<<blocksPerGrid, threadsPerBlock>>>(xDat_d, out_d, N);

	// Copy result from device memory to host memory
    cudaMemcpy(xDat, out_d, size, cudaMemcpyDeviceToHost);

	for(int i=0; i<10; i++) printf("%lf\n", xDat[i]) ;

	// Free device memory
    cudaFree(xDat_d);

	// Free host memory 
	free (xDat);

	return 0; 

} 
