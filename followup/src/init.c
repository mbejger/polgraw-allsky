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



/*  Command line options handling: search 
 */ 

void handle_opts( Search_settings *sett, 
		  Command_line_opts *opts,
		  int argc, 
		  char* argv[]) {
	
  opts->wd=NULL;

  // Default F-statistic threshold 
  opts->trl=20;

  strcpy (opts->prefix, TOSTR(PREFIX));
  strcpy (opts->dtaprefix, TOSTR(DTAPREFIX));

  opts->label[0]    = '\0';
  opts->usedet[0]   = '\0';
  opts->addsig[0]   = '\0';
//  opts->glue[0]   = '\0';
	
  // Initial value of starting frequency set to a negative quantity. 
  // If this is not changed by the command line value, fpo is calculated 
  // from the band number b (fpo = fpo = fstart + 0.96875*b/(2dt))
  sett->fpo = -1;

  // Default initial value of the data sampling time 
  sett->dt = 0.5; 

  opts->help_flag=0;
  opts->veto_flag=0; 
  opts->simplex_flag=0;
  opts->mads_flag=0;
  opts->gauss_flag=0;

  static int help_flag=0, veto_flag=0, simplex_flag=0, mads_flag=0, gauss_flag=0;

  // Reading arguments 

  while (1) {
    static struct option long_options[] = {
      {"help", no_argument, &help_flag, 1},
      {"vetolines", no_argument, &veto_flag, 1}, 
      {"simplex", no_argument, &simplex_flag, 1},
      {"mads", no_argument, &mads_flag, 1},
      {"gauss", no_argument, &gauss_flag, 1},
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
      // change directory parameter
      {"cwd", required_argument, 0, 'c'},
      // interpolation method
      {"threshold", required_argument, 0, 't'},
      // fpo value
      {"fpo", required_argument, 0, 'p'},
      // add signal parameters
      {"addsig", required_argument, 0, 'x'},
      // glue frames together
//      {"glue", required_argument, 0, 'e'},
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
      printf("-c, -cwd          Change to directory <dir>\n");
      printf("-t, -threshold    Threshold for the F-statistic (default is 20)\n");
      printf("-p, -fpo          Reference band frequency fpo value\n");
      printf("-s, -dt           data sampling time dt (default value: 0.5)\n");
      printf("-u, -usedet       Use only detectors from string (default is use all available)\n");
//      printf("-e, -glue		Glue chosen frames together. Names of frames from <file>\n");
      printf("-x, -addsig       Add signal with parameters from <file>\n\n");

      printf("Also:\n\n");
      printf("--vetolines       Veto known lines from files in data directory\n");
      printf("--simplex         Direct search of maximum using Nelder-Mead (simplex) algorithm\n");
      printf("--mads       	Direct search of maximum using MADS algorithm\n");
      printf("--gauss		Generate Gaussian noise instead of reading data. Amplitude and sigma of the noise declared in init.c\n");
      printf("--help            This help\n");

      exit(EXIT_SUCCESS);
    }

    int option_index = 0;
    int c = getopt_long_only(argc, argv, "i:b:o:d:l:c:t:p:x:s:e:u:", 
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
    case 'x':
      strcpy(opts->addsig, optarg);
      break;
    case 's':
      sett->dt = atof(optarg);
      break;
/*    case 'e':
      strcpy(opts->glue, optarg);
      break; */
    case 'u':
      strcpy(opts->usedet, optarg);
      break;


    case '?':
      break;
    default:
      break ;
    } /* switch c */
  } /* while 1 */

  opts->veto_flag = veto_flag; 
  opts->simplex_flag = simplex_flag;
  opts->mads_flag = mads_flag;
  opts->gauss_flag = gauss_flag;

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

  if(opts->trl!=20)
    printf ("Threshold for the F-statistic is %lf\n", opts->trl);
  if(strlen(opts->label))
    printf ("Using '%s' as data label\n", opts->label);
  if (strlen(opts->addsig))
    printf ("Adding signal from '%s'\n", opts->addsig);
//  if(strlen(opts->glue))
//    printf ("Gluing together frames from '%s'\n", opts->glue);
  if (opts->wd) {
    printf ("Changing working directory to %s\n", opts->wd);
    if (chdir(opts->wd)) { perror (opts->wd); abort (); }
  }

  if(opts->veto_flag) 
    printf("Known lines will be vetoed (reading from files in the data directory)\n");

  if(opts->gauss_flag) 
    printf("Gaussian noise will be generated instead of reading data from file\n");

  if(opts->mads_flag) 
    printf("MADS direct maximum search\n");

  if(opts->simplex_flag) 
    printf("Simplex direct maximum search\n");

} // end of command line options handling 

  /* Array initialization */ 

void init_arrays(
		 Search_settings *sett, 
		 Command_line_opts *opts,
		 Aux_arrays *aux_arr) {

  int i, status;
//Noise parameters 
  double amplitude = 1.0;
  double sigma = 1.0;

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
  
//  *F = (double *) calloc(2*sett->nfft, sizeof(double));

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

  /* Gluing frames */

/* void glue(Command_line_opts *opts){
//	char prfx[512]= opts->prefix;  
//	char dataprfx[512]= opts->dtaprefix;
//	char band[10]= opts->label;
	FILE *list; 
	FILE *fp;
	FILE *output;
	int i, j, k, l, flag;
	int data = 1;
	int line;
	int base_size = 2067952; //size of DetSSB.bin (in bytes)
	char xdat[200];
	char listpath[200];
	char path[200];
	char part[100];
	char out[200];
	char output_tot[512];
	char output_L[512];
	char output_H[512];
	char output_0[512];
	double *cand;
	double *ssb1_l, *ssb2_l, *ssb1_h, *ssb2_h;
	cand = (double *) calloc (data, sizeof(double));
	ssb1_l = (double *) calloc (data, sizeof(double));
	ssb2_l = (double *) calloc (data, sizeof(double));
	ssb1_h = (double *) calloc (data, sizeof(double));
	ssb2_h = (double *) calloc (data, sizeof(double));
//	sprintf(listpath, "%s/list.txt", dataprfx); //list.txt should be in data dir

	
	if ((list = fopen (opts->glue, "r")) != NULL) {
		flag = 0;
		l = 0;
		while (fscanf(list, "%d\n", &line) == 1){
			if (l != 0) flag = 1;
			l++;
;
			sprintf(xdat, "xdatc_%03d_%04d%s.bin", line, opts->band, opts->label);
			mkdir(opts->prefix, 0777);
			sprintf(output_tot, "%s/followup_total_data/", opts->prefix);
			mkdir(output_tot, 0777);
			sprintf(output_0, "%s/followup_total_data/000", opts->prefix);
			mkdir(output_0, 0777);
			sprintf(output_H, "%s/H1/", output_0);
			mkdir(output_H, 0777);
			sprintf(output_L, "%s/L1/", output_0);
			mkdir(output_L, 0777);
			for(i = 0; i < 4; i++){
				if(i == 0){
					sprintf(part, "/L1/DetSSB.bin");
					sprintf(out, "%s%s", output_0, part);
				}
				if(i == 1){
					sprintf(part, "/H1/DetSSB.bin");
					sprintf(out, "%s%s", output_0, part);
				}
				if(i == 2){
					sprintf(part, "/L1/%s", xdat);
					sprintf(out, "%s/L1/xdatc_000_%04d%s.bin", output_0, opts->band, opts->label);
				}
				if(i == 3){
					sprintf(part, "/H1/%s", xdat);
					sprintf(out, "%s/H1/xdatc_000_%04d%s.bin", output_0, opts->band, opts->label);
				}
				sprintf(path, "%s/%03d%s", opts->dtaprefix, line, part);
				k = 1;
				if ((fp = fopen (path, "rb")) != NULL) {
					if ((output = fopen (out, "ab")) != NULL) {			
						if((i == 2)||(i == 3)) {
							while (!feof(fp)) {
								if((fread ((void *)(cand), sizeof (double), data, fp))==data){
									for(j = 0; j < data; j++) fwrite((void *)(cand), sizeof (double), data, output);
								}
							}
						}
						else {
							while (!feof(fp)) {
								if((fread ((void *)(cand), sizeof (double), data, fp))==data){
									if(k < ((base_size-8)/8)){ 
										fwrite((void *)(cand), sizeof (double), data, output);
									}
									if((k == ((base_size-8)/8))&&(flag == 0)&&(i == 0)){ 
										ssb1_l[0] = cand[0]; 									
									}
									if((k == (base_size/8))&&(flag == 0)&&(i == 0)){ 
										ssb2_l[0] = cand[0]; 
									}	
									if((k == ((base_size-8)/8))&&(flag == 0)&&(i == 1)){ 
										ssb1_h[0] = cand[0]; 
									}
									if((k == (base_size/8))&&(flag == 0)&&(i == 1)){ 
										ssb2_h[0] = cand[0]; 
									}
									k++;
								}
							}
							
						}
					}
					else {		
						printf("Problem with %s file!\n", out);
					}
				}
				else {		
					printf("Problem with %s file!\n", path);
				}
				fclose(fp);
				fclose(output);

			}
		}
	}
	else {
		
		perror (listpath);
//		return 1;
	}
	sprintf(out, "%s/DetSSB.bin",output_L);
	if ((output = fopen (out, "ab")) != NULL) {
		fwrite((void *)(ssb1_l), sizeof (double), data, output);
		fwrite((void *)(ssb2_l), sizeof (double), data, output);
	}
	else {		
		printf("Problem with %s file - at the end!\n", out);
	}
	sprintf(out, "%s/DetSSB.bin",output_H);
	if ((output = fopen (out, "ab")) != NULL) {
		fwrite((void *)(ssb1_h), sizeof (double), data, output);
		fwrite((void *)(ssb2_h), sizeof (double), data, output);
	}
	else {		
		printf("Problem with %s file - at the end!\n", out);
	}
//	puts("END OF GLUE FUNCTION");
	fclose(list);
//	return 0;

} */

unsigned long int random_seed()
{

 unsigned int seed;
 struct timeval tv;
 FILE *devrandom;

 if ((devrandom = fopen("/dev/random","r")) == NULL) {
   gettimeofday(&tv,0);
   seed = tv.tv_sec + tv.tv_usec;
 } else {
   fread(&seed,sizeof(seed),1,devrandom);
   fclose(devrandom);
 }

 return(seed);

}

  /* Generate Gaussian noise (keep it only in memory; do not save) */

void gauss_xdat(Search_settings *sett, double amplitude, double sigma, int i){
puts("Generate Gaussian noise");
  gsl_rng * r;
  int j;
  unsigned long mySeed;
  mySeed = random_seed();

  const gsl_rng_type * T;
  
  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_set(r, mySeed);
  // Generate normal distribution (around 0, 
  // with amplitude and sigma as parameters)
  for(j = 0; j < sett->N; j++) 
    ifo[i].sig.xDat[j] = amplitude*gsl_ran_gaussian_ziggurat(r, sigma);
 
  gsl_rng_free(r); 


}

  /* Add signal to data   */ 


void add_signal(
		Search_settings *sett,
		Command_line_opts *opts,
		Aux_arrays *aux_arr) {
puts("Adding signal from file");
  int i, j, n, gsize, reffr, k, hem; 
  double snr, sum = 0., h0, cof, thsnr = 0., d1; 
  double sigma_noise = 1.0;
  double be[2];
  double sinaadd, cosaadd, sindadd, cosdadd, phaseadd, shiftadd, signadd; 
  double nSource[3], sgnlo[8];
  double **sigaa, **sigbb;   // aa[nifo][N]
  sigaa = (double **)malloc(sett->nifo*sizeof(double *));
  for (k=0; k < sett->nifo; k++) sigaa[k] = (double *)calloc(sett->N, sizeof(double));
  sigbb = (double **)malloc(sett->nifo*sizeof(double *));
  for (k=0; k < sett->nifo; k++) sigbb[k] = (double *)calloc(sett->N, sizeof(double));
  
  FILE *data;
  
  // Signal parameters are read
  if ((data=fopen (opts->addsig, "r")) != NULL) {
	
    // Fscanning for the GW snr, grid size and the reference
    // frame (for which the signal freq. is not spun-down/up)
    fscanf (data, "%le %d %d", &snr, &gsize, &reffr);    

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
  //#mb For VSR1 reffr=67
 
  sgnlo[0] += -2.*sgnlo[1]*(sett->N)*(reffr - opts->ident); 
 
  // Check if the signal is in band 
  if(sgnlo[0]<0) exit(171);          // &laquo;  
  else if (sgnlo[0]>M_PI) exit(187); // &raquo;

  cof = sett->oms + sgnlo[0]; 

//Hemisphere an be vector 
//(previously was fscanned from sigfile, now calculated here)

  hem = ast2lin(sgnlo[3], sgnlo[2], C_EPSMA, be);

  printf("%le %le %d\n", be[0], be[1], hem);


  // sgnlo[2]: declination, sgnlo[3]: right ascension 
  sindadd = sin(sgnlo[2]); 
  cosdadd = cos(sgnlo[2]); 
  sinaadd = sin(sgnlo[3]);  
  cosaadd = cos(sgnlo[3]); 
	
  // To keep coherent phase between time segments  
  double phaseshift = sgnlo[0]*sett->N*(reffr - opts->ident)   
    + sgnlo[1]*pow(sett->N*(reffr - opts->ident), 2); 

  // Loop for each detector - sum calculations
  for(n=0; n<sett->nifo; n++) {
    
    modvir(sinaadd, cosaadd, sindadd, cosdadd,
	   sett->N, &ifo[n], aux_arr, sigaa[n], sigbb[n]);

    nSource[0] = cosaadd*cosdadd;
    nSource[1] = sinaadd*cosdadd;
    nSource[2] = sindadd;
					
    // adding signal to data (point by point)  								
    for (i=0; i<sett->N; i++) {
      shiftadd = 0.; 					 
      for (j=0; j<3; j++)
      	shiftadd += nSource[j]*ifo[n].sig.DetSSB[i*3+j];		 
      
      // Phase 
      phaseadd = sgnlo[0]*i + sgnlo[1]*aux_arr->t2[i] 
        + (cof + 2.*sgnlo[1]*i)*shiftadd
        - phaseshift; 

      // The whole signal with 4 amplitudes and modulations 
      signadd = sgnlo[4]*(sigaa[n][i])*cos(phaseadd) 
        + sgnlo[6]*(sigaa[n][i])*sin(phaseadd) 
        + sgnlo[5]*(sigbb[n][i])*cos(phaseadd) 
        + sgnlo[7]*(sigbb[n][i])*sin(phaseadd);

// Sum over signals
      sum += pow(signadd, 2.);
	 
    } //data loop

    // Write the data+signal to file   
/*    FILE *dataout;
    char xxx[512]; 
    sprintf(xxx, "%s/%03d/%s/xdatc_%03d_%04d%s.bin",
          opts->dtaprefix, opts->ident, ifo[n].name, 
          opts->ident, opts->band, opts->label);
    printf("Flag 1: try to write xDat to fole %s\n", xxx);
    if((dataout = fopen(xxx, "wb")) != NULL) {
    	fwrite(ifo[n].sig.xDat, sizeof(*ifo[n].sig.xDat), sett->N, dataout);
    }
    else{
    	printf("Problem with %s file!\n", xxx);
    }
    fclose(dataout); */
    
  } //detector loop

//Signal amplitude

  h0 = (snr*sigma_noise)/(sqrt(sum));

// Loop for each detector - adding signal to data (point by point)  								
  for(n=0; n<sett->nifo; n++) {
    
    modvir(sinaadd, cosaadd, sindadd, cosdadd,
	   sett->N, &ifo[n], aux_arr, sigaa[n], sigbb[n]);

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
      signadd = sgnlo[4]*(sigaa[n][i])*cos(phaseadd) 
        + sgnlo[6]*(sigaa[n][i])*sin(phaseadd) 
        + sgnlo[5]*(sigbb[n][i])*cos(phaseadd) 
        + sgnlo[7]*(sigbb[n][i])*sin(phaseadd);

      // Adding the signal to the data vector 
      if(ifo[n].sig.xDat[i]) { 
        ifo[n].sig.xDat[i] += h0*signadd;
      } // if xDat
 
    } //data loop
  } //detector loop

//printf("snr=%le h0=%le\n", snr, h0);

// Free memory
for (i = 0; i < sett->nifo; i++) free(sigaa[i]);
free(sigaa);
for (i = 0; i < sett->nifo; i++) free(sigbb[i]);
free(sigbb);
//exit(0);
} //add_signal

  /* Cleanup & memory free 
	 */
void cleanup_followup(
	Search_settings *sett,
	Command_line_opts *opts,
	Aux_arrays *aux) {

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

}


