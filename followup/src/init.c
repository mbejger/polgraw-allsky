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
#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h>
#include <time.h>
#include <sys/time.h>

#include "init.h"
#include "struct.h"
#include "settings.h"
#include "auxi.h"


/***************************************************************
Command line options handling
***************************************************************/
void handle_opts( Search_settings *sett, Command_line_opts *opts, int argc, char* argv[]) {
	opts->wd=NULL;
// Default F-statistic threshold 
  	opts->trl = 20;

  	strcpy (opts->prefix, TOSTR(PREFIX));
  	strcpy (opts->dtaprefix, TOSTR(DTAPREFIX));

  	opts->label[0]    = '\0';
  	opts->usedet[0]   = '\0';
  	opts->addsig[0]   = '\0';
  	opts->candidates[0]   = '\0';
// Initial value of starting frequency set to a negative quantity. 
// If this is not changed by the command line value, fpo is calculated 
// from the band number b (fpo = fpo = fstart + 0.96875*b/(2dt))
  	sett->fpo = -1;
// Default initial value of the data sampling time 
  	sett->dt = 0.5; 

// Initial number of reference frame
  	opts->refr = -1;

  	opts->help_flag=0;
  	opts->simplex_flag=0;
  	opts->onepoint_flag=0;
  	opts->mads_flag=0;
  	opts->skymads_flag=0;
  	opts->gauss_flag=0;
  	opts->neigh_flag=0;
  	opts->naive_flag=0;

  	static int help_flag=0, simplex_flag=0, onepoint_flag=0, mads_flag=0, skymads_flag=0, gauss_flag=0, neigh_flag=0, naive_flag=0;

// Reading arguments 
  	while (1) {
		static struct option long_options[] = {
			{"help", no_argument, &help_flag, 1},
//Simplex maximum search
			{"simplex", no_argument, &simplex_flag, 1},
//invertedMADS/MADS maximum search
			{"mads", no_argument, &mads_flag, 1},
//skyMADS search - only in frequency and spindown parameters
			{"skymads", no_argument, &skymads_flag, 1},
//Calculate F-statistics only for initial point (don't generate grid)
			{"onepoint", no_argument, &onepoint_flag, 1},
//Generate Gaussian noise instead of reading data
			{"gauss", no_argument, &gauss_flag, 1},
//Uniform grid around initial point
			{"neigh", no_argument, &neigh_flag, 1},
//Uniform grid between grid from grid.bin closest points
			{"naive", no_argument, &naive_flag, 1},
//Frame number
			{"ident", required_argument, 0, 'i'},
// Frequency band number
			{"band", required_argument, 0, 'b'},
// Output directory
	      		{"output", required_argument, 0, 'o'},
// Input data directory
	      		{"data", required_argument, 0, 'd'},
// Non-standard label for naming files
	      		{"label", required_argument, 0, 'l'},
// Change directory parameter
	      		{"cwd", required_argument, 0, 'c'},
// Threshold value
	      		{"threshold", required_argument, 0, 't'},
// fpo value
	      		{"fpo", required_argument, 0, 'p'},
// Number of days in the time-domain segment 
	      		{"nod", required_argument, 0, 'y'},
// Add signal parameters
	      		{"addsig", required_argument, 0, 'x'},
// Which detectors to use
	      		{"usedet", required_argument, 0, 'u'}, 
// Data sampling time 
	      		{"dt", required_argument, 0, 's'},
// Use candidates from file 
	      		{"candidates", required_argument, 0, 'a'},
// Reference frame number 
	      		{"refr", required_argument, 0, 'r'},
// File with precise Fisher matrix 
	      		{"fisher", required_argument, 0, 'f'},
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
      		printf("-c, -cwd          Change to directory <dir>\n");
      		printf("-t, -threshold    Threshold for the F-statistic (default is 20)\n");
      		printf("-p, -fpo          Reference band frequency fpo value\n");
      		printf("-s, -dt           data sampling time dt (default value: 0.5)\n");
      		printf("-u, -usedet       Use only detectors from string (default is use all available)\n");
      		printf("-y, -nod          Number of days\n");
      		printf("-x, -addsig       Add signal with parameters from <file>\n");
      		printf("-a, -candidates   As a starting point in followup use parameters from <file>\n");
     		printf("-r, -refr         Reference frame number\n");
      		printf("-f, -fisher       Read precise, optimal Fisher matrix from <file> (put just name; <file> should be in the same directory as data)\n\n");

      		printf("Also:\n\n");
      		printf("--simplex       Direct search of maximum using Nelder-Mead (simplex) algorithm\n");
      		printf("--mads       	Direct search of maximum using invertedMADS/MADS algorithm\n");
      		printf("--skymads      	Direct search of maximum using 2D MADS algorithm (search only in frequency and spindown; keep sky position fixed).\n");
      		printf("--gauss		Generate Gaussian noise instead of reading data. Amplitude and sigma of the noise declared in init.c\n");
      		printf("--neigh		Function neigh() generate area as %% from initial value instead of taking it from grid.bin\n");
      		printf("--naive		Function naive() generate area as +/- points taking it from grid.bin and divide it into smaller grid.\n");
      		printf("--onepoint	Calculate Fstatistic only in one point taken from file with candidates (without generating any grid).\n");
      		printf("--help            This help\n");

      		exit(EXIT_SUCCESS);
    		}

    		int option_index = 0;
    		int c = getopt_long_only(argc, argv, "i:b:o:d:l:c:t:p:x:a:r:s:e:u:", long_options, &option_index);
    		if (c == -1){
     			break;
		}
    		switch (c) {
    			case 'i':
	      			opts->ident = atoi (optarg);
	      			break;
    			case 't':
		      		opts->trl = atof(optarg);
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
    			case 'c':
      				opts->wd = (char *) malloc (1+strlen(optarg));
      				strcpy(opts->wd, optarg);
      				break;
    			case 'p':
      				sett->fpo = atof(optarg);
      				break;
    			case 'y':
      				sett->nod = atoi(optarg);
      				break;
    			case 'x':
      				strcpy(opts->addsig, optarg);
      				break;
    			case 'a':
      				strcpy(opts->candidates, optarg);
      				break;
    			case 'r':
      				opts->refr = atoi(optarg);
      				break;
    			case 's':
     	 			sett->dt = atof(optarg);
      				break;
    			case 'u':
      				strcpy(opts->usedet, optarg);
      				break;
    			case 'f':
      				strcpy(opts->fisher, optarg);
      				break;
			case '?':
      				break;
    			default:
      				break ;
    		} /* switch c */
  	} /* while 1 */

	opts->simplex_flag = simplex_flag;
  	opts->onepoint_flag = onepoint_flag;
  	opts->mads_flag = mads_flag;
  	opts->skymads_flag = skymads_flag;
  	opts->gauss_flag = gauss_flag;
  	opts->neigh_flag = neigh_flag;
  	opts->naive_flag = naive_flag;
// Check if sett->nod was set up, if not, exit
  	if(!(sett->nod)) { 
    		printf("Number of days not set... Exiting\n"); 
    		exit(EXIT_FAILURE); 
  	} 

  	printf("Number of days is %d\n", sett->nod); 
	printf("Input data directory is %s\n", opts->dtaprefix);
  	printf("Output directory is %s\n", opts->prefix);
  	printf("Frame and band numbers are %d and %d\n", opts->ident, opts->band);

// Starting band frequency:
// fpo_val is optionally read from the command line
// Its initial value is set to -1
  	if(!(sett->fpo >= 0))
// The usual definition (multiplying the offset by B=1/(2dt))
// !!! in RDC_O1 the fstart equals 10, not 100 like in VSR1 !!! 
    	sett->fpo = 10. + 0.96875*opts->band*(0.5/sett->dt);
//Print information on the screen
  	printf("The reference frequency fpo is %f\n", sett->fpo);
  	printf("The data sampling time dt is %f\n", sett->dt); 

  	if(opts->trl!=20){
    		printf ("Threshold for the F-statistic is %lf\n", opts->trl);
	}
  	if(opts->refr > 0){
    		printf ("Reference frame is %d\n", opts->refr);
	}
  	if(strlen(opts->label)){
    		printf ("Using '%s' as data label\n", opts->label);
	}
  	if (strlen(opts->addsig)){
    		printf ("Adding signal from '%s'\n", opts->addsig);
	}
  	if (strlen(opts->candidates)){
    		printf ("Starting point for followup taken from '%s'\n", opts->candidates);
	}
  	if (opts->wd) {
    		printf ("Changing working directory to %s\n", opts->wd);
    		if (chdir(opts->wd)) { 
			perror (opts->wd); abort (); 
		}
  	}
  	if(opts->gauss_flag){ 
    		printf("Gaussian noise will be generated instead of reading data from file\n");
	}
  	if(opts->neigh_flag){ 
		printf("Area of calculation will be defined as +/- from initial value instead of taking it from grid.bin\n");
	}
	if(opts->naive_flag){ 
    		printf("Area of calculation will be defined as +/- points from grid.bin\n Then area will be divided into bins\n (number of points and bins defined in followup.c)\n");
	}
	if(opts->mads_flag){ 
    		printf("(inverted)MADS direct maximum search\n");
	}
	if(opts->skymads_flag){ 
    		printf("skyMADS direct maximum search - only in frequency and spindown! Sky position is fixed!\n");
	}
  	if(opts->simplex_flag){ 
    		printf("Simplex direct maximum search\n");
	}
} // end of command line options handling 

/*************************************************************** 
Generate grid from the M matrix (grid.bin)
***************************************************************/ 
void read_grid(Search_settings *sett, Command_line_opts *opts) {
  	sett->M = (double *) calloc (16, sizeof (double));
  	sett->gamrn = (double *) calloc (16, sizeof (double));
  	FILE *data;
  	char filename[512];
  	int i;

// In case when -usedet option is used for one detector
// i.e. opts->usedet has a length of 2 (e.g. H1 or V1), 
// read grid.bin from this detector subdirectory 
// (see detectors_settings() in settings.c for details) 
  	if(strlen(opts->usedet)==2){
    		sprintf (filename, "%s/%03d/%s/grid.bin", opts->dtaprefix, opts->ident, opts->usedet);
	}
  	else{ 
    		sprintf (filename, "%s/%03d/grid.bin", opts->dtaprefix, opts->ident);
	}
  	if ((data=fopen (filename, "r")) != NULL) {
    		printf("Using grid file from %s\n", filename);
    		fread ((void *)&sett->fftpad, sizeof (int), 1, data);
    		printf("Using fftpad from the grid file: %d\n", sett->fftpad); 
// M: vector of 16 components consisting of 4 rows
// of 4x4 grid-generating matrix
    		fread ((void *)sett->M, sizeof (double), 16, data);
    		fread ((void *)sett->gamrn, sizeof (double), 16, data);
    		fclose (data);
  	} 
	else {
    		perror (filename);
    		exit(EXIT_FAILURE);
  	}
} // end of read grid 

/*************************************************************** 
Array initialization 
***************************************************************/
void init_arrays(Search_settings *sett, Command_line_opts *opts, Aux_arrays *aux_arr) {
	int i, status;
//Noise parameters 
  	double amplitude = 1.0;
  	double sigma = 1.0;

// Allocates and initializes to zero the data, detector ephemeris
// and the F-statistic arrays

  	FILE *data;

  	for(i=0; i<sett->nifo; i++) { 
    		ifo[i].sig.xDat = (double *) calloc(sett->N, sizeof(double));
    		if(!opts->gauss_flag){	    
// Input time-domain data handling
// The file name ifo[i].xdatname is constructed in settings.c,
// while looking for the detector subdirectories
			if((data = fopen(ifo[i].xdatname, "r")) != NULL) {
         			status = fread((void *)(ifo[i].sig.xDat), 
 	  	       		sizeof(double), sett->N, data);
         			fclose (data);
       
       			} 
			else {
         			perror (ifo[i].xdatname);
         			exit(EXIT_FAILURE); 
       			}
    		}
    		else {
//If gauss_flag: generating Gaussian noise instead of reading data
      			gauss_xdat(sett, amplitude, sigma, i);	
    		}
    		int j, Nzeros=0;
// Checking for null values in the data
    		for(j=0; j < sett->N; j++){
      			if(!ifo[i].sig.xDat[j]){ 
				Nzeros++;
			}
		}
		ifo[i].sig.Nzeros = Nzeros; 
// factor N/(N - Nzeros) to account for null values in the data
    		ifo[i].sig.crf0 = (double)sett->N/(sett->N - ifo[i].sig.Nzeros);
// Estimation of the variance for each detector 
    		ifo[i].sig.sig2 = (ifo[i].sig.crf0)*var(ifo[i].sig.xDat, sett->N);
		ifo[i].sig.DetSSB = (double *) calloc(3*sett->N, sizeof(double));
// Ephemeris file handling
    		char filename[512];
    		sprintf (filename, "%s/%03d/%s/DetSSB.bin", opts->dtaprefix, opts->ident, ifo[i].name);
    		if((data = fopen(filename, "r")) != NULL) {
// Detector position w.r.t Solar System Baricenter for every datapoint
      			status = fread((void *)(ifo[i].sig.DetSSB), 
               		sizeof(double), 3*sett->N, data);
// Deterministic phase defining the position of the Earth in its diurnal motion at t=0 
			status = fread((void *)(&ifo[i].sig.phir), sizeof(double), 1, data);
// Earth's axis inclination to the ecliptic at t=0
      			status = fread((void *)(&ifo[i].sig.epsm), sizeof(double), 1, data);
      			fclose (data);
			printf("Using %s as detector %s ephemerids...\n", filename, ifo[i].name);
		} 
		else {
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

    		ifo[i].sig.xDatma = (complex double *) calloc(sett->N, sizeof(complex double));
    		ifo[i].sig.xDatmb = (complex double *) calloc(sett->N, sizeof(complex double));

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

// If all is well with epsm, take the first value 
  	sett->sepsm = ifo[0].sig.sepsm;
  	sett->cepsm = ifo[0].sig.cepsm;
  
// Auxiliary arrays, Earth's rotation
  	aux_arr->t2 = (double *) calloc(sett->N, sizeof (double));
  	aux_arr->cosmodf = (double *) calloc(sett->N, sizeof (double));
  	aux_arr->sinmodf = (double *) calloc(sett->N, sizeof (double));
  	double omrt;

  	for (i=0; i<sett->N; i++) {
// Earth angular velocity * dt * i
    		omrt = (sett->omr)*i;     
    	aux_arr->t2[i] = sqr((double)i);
    	aux_arr->cosmodf[i] = cos(omrt);
    	aux_arr->sinmodf[i] = sin(omrt);
  	}
} // end of init arrays 

/*************************************************************** 
Get random seed (needed for Gaussian noise generator)
***************************************************************/
unsigned long int random_seed(){
	unsigned int seed;
 	struct timeval tv;
 	FILE *devrandom;

 	if ((devrandom = fopen("/dev/random","r")) == NULL) {
   		gettimeofday(&tv,0);
   		seed = tv.tv_sec + tv.tv_usec;
 	} 
	else {
   		fread(&seed,sizeof(seed),1,devrandom);
   		fclose(devrandom);
 	}
	return(seed);
}//end random_seed()

/*************************************************************** 
Generate Gaussian noise with GSL function (keep it only in memory; do not save) 
***************************************************************/
void gauss_xdat(Search_settings *sett, double amplitude, double sigma, int i){
  	gsl_rng * r;
  	int j;
  	unsigned long mySeed;
  	mySeed = random_seed();
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	 T = gsl_rng_default;
  	r = gsl_rng_alloc (T);
  	gsl_rng_set(r, mySeed);
// Generate normal distribution (around 0, with amplitude and sigma as parameters)
  	for(j = 0; j < sett->N; j++){ 
    		ifo[i].sig.xDat[j] = amplitude*gsl_ran_gaussian_ziggurat(r, sigma);
	}
	gsl_rng_free(r); 
} //end of Gaussian noise generator

/***************************************************************  
Add signal to data (signal parameters are read from file)
***************************************************************/ 
void add_signal(Search_settings *sett, Command_line_opts *opts, Aux_arrays *aux_arr) {
	puts("Adding signal from file");
	int i, j, n, gsize, reffr, k, hem; 
  	double snr=0., sum = 0., h0=0., cof, thsnr = 0., d1; 
  	double sigma_noise = 1.0;
  	double be[2];
  	double sinaadd, cosaadd, sindadd, cosdadd, phaseadd, shiftadd; 
  	double nSource[3], sgnlo[8];
  	char amporsnr[3];  
//Amplitude modulation factors
  	double **sigaa, **sigbb; 
  	sigaa = (double **)malloc(sett->nifo*sizeof(double *));
  	for (k=0; k < sett->nifo; k++){ 
		sigaa[k] = (double *)calloc(sett->N, sizeof(double));
	}
  	sigbb = (double **)malloc(sett->nifo*sizeof(double *));
  	for (k=0; k < sett->nifo; k++){ 
		sigbb[k] = (double *)calloc(sett->N, sizeof(double));
	}
	FILE *data;
// Signal parameters are read
  	if ((data=fopen (opts->addsig, "r")) != NULL) {
// Fscanning for the GW amplitude h0 or signal-to-noise, the grid size and
// the reference frame (for which the signal freq. is not spun-down/up)
		fscanf (data, "%s", amporsnr);  
// Check if first value is an amplitude or SNR (file should have flag 'amp' or 'snr')  
		if(!strcmp(amporsnr, "amp")) { 
      			fscanf (data, "%le %d %d", &h0, &gsize, &reffr); 
      			printf("add_signal(): GW amplitude h0 is %le\n", h0); 
    		} 
		else if(!strcmp(amporsnr, "snr")) { 
      			fscanf (data, "%le %d %d", &snr, &gsize, &reffr); 
      			printf("add_signal(): GW (network) signal-to-noise ratio is %le\n", snr); 
    		} 
		else { 
      			printf("Problem with the signal file. Exiting...\n"); 
      			exit(0); 
    		} 
// This is checkpoint for tests
		puts("-1000 -1000 -1000 -1000 -1000 -1000");
// Fscanning signal parameters: f, fdot, delta, alpha (sgnlo[0], ..., sgnlo[3])
// four amplitudes sgnlo[4], ..., sgnlo[7] 
// (see sigen.c and Phys. Rev. D 82, 022005 2010, Eqs. 2.13a-d) 

    		for(i=0; i<8; i++){
      			fscanf(data, "%le",i+sgnlo); 
		}   
    		fclose (data);
	} 
	else {
    		perror (opts->addsig);
  	}
// Search-specific parametrization of freq. for the software
// injections sgnlo[0]: frequency, sgnlo[1]: frequency. derivative  
	sgnlo[0] += -2.*sgnlo[1]*(sett->N)*(reffr - opts->ident);  
// Check if the signal is in band 
  	if(sgnlo[0]<0){			// &laquo; 
		exit(171); 
	} 
  	else if(sgnlo[0]>M_PI){  	// &raquo;
		 exit(187);
	}
	cof = sett->oms + sgnlo[0];   
//Hemisphere and be vector 
//(previously was fscanned from sigfile, now calculated here)
	hem = ast2lin(sgnlo[3], sgnlo[2], C_EPSMA, be);
// sgnlo[2]: declination, sgnlo[3]: right ascension 
  	sindadd = sin(sgnlo[2]); 
  	cosdadd = cos(sgnlo[2]); 
  	sinaadd = sin(sgnlo[3]);  
  	cosaadd = cos(sgnlo[3]); 	
// To keep coherent phase between time segments  
  	double phaseshift = sgnlo[0]*sett->N*(reffr - opts->ident) + sgnlo[1]*pow(sett->N*(reffr - opts->ident), 2); 
// Allocate arrays for added signal, for each detector 
  	double **signadd = malloc((sett->nifo)*sizeof(double *));
  	for(n=0; n<sett->nifo; n++){
    		signadd[n] = malloc((sett->N)*sizeof(double));
	}
// Loop for each detector - sum calculations
  	for(n=0; n<sett->nifo; n++) {
		modvir(sinaadd, cosaadd, sindadd, cosdadd, sett->N, &ifo[n], aux_arr, sigaa[n], sigbb[n]);

    		nSource[0] = cosaadd*cosdadd;
    		nSource[1] = sinaadd*cosdadd;
    		nSource[2] = sindadd;
					
    		for (i=0; i<sett->N; i++) {
			shiftadd = 0.; 					 
      			for (j=0; j<3; j++){
      				shiftadd += nSource[j]*ifo[n].sig.DetSSB[i*3+j];	
			}	 
// Phase 
      			phaseadd = sgnlo[0]*i + sgnlo[1]*aux_arr->t2[i] 
        		+ (cof + 2.*sgnlo[1]*i)*shiftadd
        		- phaseshift; 
// The whole signal with 4 amplitudes and modulations 
			signadd[n][i] = sgnlo[4]*(sigaa[n][i])*cos(phaseadd) 
        		+ sgnlo[6]*(sigaa[n][i])*sin(phaseadd) 
        		+ sgnlo[5]*(sigbb[n][i])*cos(phaseadd) 
        		+ sgnlo[7]*(sigbb[n][i])*sin(phaseadd);
// Sum over signals
      			sum += pow(signadd[n][i], 2.);
		} // data loop
	} // detector loop
// Signal amplitude h0 from the snr 
// (currently only makes sense for Gaussian noise with fixed sigma)
  	if(snr){
    		h0 = (snr*sigma_noise)/(sqrt(sum));
	}
// Loop for each detector - adding signal to data (point by point)  								
  	for(n=0; n<sett->nifo; n++) {
    		for (i=0; i<sett->N; i++) {
// Adding the signal to the data vector 
		      	if(ifo[n].sig.xDat[i]) { 
				ifo[n].sig.xDat[i] += h0*signadd[n][i];
			} 
		} // data loop
	} // detector loop
// Free auxiliary 2d array 
  	for(n=0; n<sett->nifo; n++){ 
    		free(signadd[n]);
  		free(signadd);
	}
	for (i = 0; i < sett->nifo; i++){ 
		free(sigaa[i]);
	}
	free(sigaa);
	for (i = 0; i < sett->nifo; i++){ 
		free(sigbb[i]);
	}
	free(sigbb); 
} // add_signal()

/*************************************************************** 
Cleanup & memory free 
***************************************************************/
void cleanup_followup(Search_settings *sett, Command_line_opts *opts, Aux_arrays *aux) {
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
  	free(sett->M);
  	free(sett->gamrn);
  	free(aux->sinmodf);
  	free(aux->cosmodf);
  	free(aux->t2);
}//cleanup()


