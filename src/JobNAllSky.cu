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
#include <cuda.h>
#include <cufft.h>
#include "auxi.h"
#include "settings.h"

static int white_flag=0; 	//  white noise flag
static int s0_flag=0;		// no spin-down flag
static int help_flag=0; 

/* Default output and data directories */

#ifndef PREFIX
#define PREFIX .
#endif

#ifndef DTAPREFIX
#define DTAPREFIX ../data
#endif

double *cosmodf, *sinmodf, *aa, *bb, \
  *shftf, *shft, *xDat, *DetSSB, *F;

//#mb not needed
//double *t2; 

complex double *xDatma, *xDatmb;

__global__ void 
zerotrig(int *trig) {
  int i=threadIdx.x+blockIdx.x*blockDim.x;
  if(i==0) {
    trig[0]=0;
  }
}

__global__ void 
zapisz(int *trig, cufftComplex *du,cufftComplex *d,cufftComplex *dl, int N){
  int i=threadIdx.x+blockIdx.x*blockDim.x;
  if(i==0) {
    trig[0]=0;
    du[0].x=1.0;
    du[0].y=0.0;
    dl[0].x=0.0;
    dl[0].y=0.0;
    d[0].x=4.0;
    d[0].y=0.0;
  }
  if(i==(N-1)) {
    du[i].x=0.0;
    du[i].y=0.0;
    dl[i].x=1.0;
    dl[i].y=0.0;
    d[i].x=4.0;
    d[i].y=0.0;
  }
  if(i>0 && i<(N-1)) {
    du[i].x=1.0;
    du[i].y=0.0;
    dl[i].x=1.0;
    dl[i].y=0.0;
    d[i].x=4.0;
    d[i].y=0.0;
  }
}

int
JobNAllSky (int argc, char *argv[]) {
  int i, pm, mm, nn, pst, mst, nst, sst, sgnlc, hemi=0, Nzeros=0,
    Ninterp, FNum, nmin, nmax, spndr[2], nr[2], mr[2], pmr[2], c,
    ident=0, band=0, nfftf, 
    fftinterp=INT; // default value
  char filename[64], outname[64], qname[64], 
    prefix[64], dtaprefix[64], label[64], range[64], *ifo="V", *wd=NULL;
  double *sgnlv, omrt, coft, epsm, sepsm, cepsm, phir, sphir, cphir, *M, 
    trl=20.,        // default value for the F-statistic threshold
    fpo, sig2;
  complex float *xa, *xao;

#if TIMERS==1
  struct timespec t0, t1;
#define NANO_INV 1000000000L
#endif

  cufftHandle plan, pl_int, pl_inv;
  FILE *state, *data;
  FILE *candid;

  struct stat buffer;

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
      printf("-b		Band number\n"); 
      printf("-c		Change to directory <dir>\n"); 
      printf("-d		Data directory (default is .)\n"); 
      printf("-f		Intepolation method (INT [default] or FFT)\n"); 
      printf("-h		Hemisphere (default is 0 - does both)\n");
      printf("-i		Frame number\n"); 
      printf("-l		Custom label for the input and output files\n"); 
      printf("-o		Output directory (default is ./candidates)\n"); 
      printf("-r		Grid range filename\n"); 
      printf("-t		Threshold for the F-statistic (default is 20)\n\n");
      printf("-x		Detector (H, L or V; default is V)\n");  
      printf("Also:\n"); 
      printf("--whitenoise	Will assume white Gaussian noise\n"); 
      printf("--nospindown	Will assume that spindown is not important\n"); 
      printf("--help		This help\n"); 		

      exit (0);
    }

    int option_index = 0;
    c = getopt_long (argc, argv, "b:c:d:f:h:i:l:o:r:t:x:", long_options, \
		     &option_index);
    if (c == -1)
      break;

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
    default: break ; 
    } /* switch c */
  } /* while 1 */


  // Detector choice, by first letter 
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
  sprintf (filename, "%s/%02d/grid.bin", dtaprefix, ident);
  if ((data=fopen (filename, "r")) != NULL) {
    // fftpad: used to zero padding to fftpad*nfft data points 
    if(fread ((void *)&fftpad, sizeof (int), 1, data)==0) {
      printf("\n Failed to read a grid file \n");
      abort();
    }
    // M: vector of 16 components consisting of 4 rows 
    // of 4x4 grid-generating matrix  
    if(fread ((void *)M, sizeof (double), 16, data)==0) {
      printf("\n Failed to read a grid file \n");
      abort();
    }
    fclose (data);
  } else {
    perror (filename);
    return 1;
  }

  // Starting band frequency
  fpo = 100. + 0.96875 * band;
  // Detector, ephemerides, constants 
  settings (fpo, ifo);
  // Amplitude modulation functions coefficients
  //rogcvir ();
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
      if(fscanf (data, "%d %d", spndr, 1+spndr)==0 ||
	 fscanf (data, "%d %d", nr, 1+nr)==0||
	 fscanf (data, "%d %d", mr, 1+mr)==0||
	 fscanf (data, "%d %d", pmr, 1+pmr)==0)
	{
	  printf("\n Failed reading the range file \n");
	  abort();
	}
      fclose (data);
    } else {
      perror (range);
      return 1;
    }
  }

  // hard work start here
  if (cuinit(CUDA_DEV) == -1) {
    printf("\nGPU device initialization error!\n");
    exit(EXIT_FAILURE);
  }
#if TIMERS==1
  cudaThreadSynchronize();
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);  
#endif

  // Allocates and initializes to zero the data, detector ephemeris 
  // and the F-statistic arrays
  xDat = (double *) calloc (N, sizeof (double));
  DetSSB = (double *) calloc (3*N, sizeof (double));
  //F = (double *) calloc (2*nfft, sizeof (double));
  cudaMallocHost((void**)&F, 2*nfft*sizeof(double));
  double *cuF;
  cudaMalloc((void**)&cuF,2*nfft*sizeof(double));
  // pci wtf? po co kopiować śmieci? 
  // F na cpu jest globalna - mogłaby być lokalna w tej funkcji
  // na gpu wystarczy tylko w device bo nie będzie używana poza kernelem
  // i nie będzie potrzeby jej kopiować do i z device
  cudaMemcpy(cuF,F, 2*nfft*sizeof(double),cudaMemcpyHostToDevice);
  // Input time-domain data handling 
  sprintf (filename, "%s/%02d/xdat_%02d_%03d%s.bin", dtaprefix, ident,	\
	   ident, band, label);
  if ((data = fopen (filename, "r")) != NULL) {
    if(fread ((void *)(xDat), sizeof (double), N, data)==0) {
      printf("\n Failed reading data file \n");
      abort();
    }
    fclose (data);
  } else {
    perror (filename);
    return 1;
  }
  double *cuxDat;
  cudaMalloc((void**)&cuxDat,N*sizeof(double));
  // pci wtf? po co kopiować śmieci? 
  cudaMemcpy(cuxDat,xDat, N*sizeof(double),cudaMemcpyHostToDevice);
  // Checking for null values in the data 
  for(i=0; i<N; i++) 
    if(!xDat[i])
      Nzeros++; 

  // In case of white noise, factor N/(N - Nzeros) 
  // accounts for null values in the data 
  if (white_flag) sig2 = N*var (xDat, N)*N/(N - Nzeros);
  else sig2 = -1.;

  // Ephemeris file handling 
  sprintf (filename, "%s/%02d/DetSSB.bin", dtaprefix, ident);
  if ((data = fopen (filename, "r")) != NULL) {
    // Detector position w.r.t solar system baricenter for every datapoint
    if(fread ((void *)(DetSSB), sizeof (double), 3*N, data)==0) {
      printf("\n Failed reading detector ephemerides file 1 \n");
      abort();
    }
    // Deterministic phase defining the position of the Earth 
    // in its diurnal motion at t=0  
    if(fread ((void *)(&phir), sizeof (double), 1, data)==0) {
      printf("\n Failed reading detector ephemerides file 2  \n");
      abort();
    }
    // Earth's axis inclination to the ecliptic at t=0  
    if(fread ((void *)(&epsm), sizeof (double), 1, data)==0) {
      printf("\n Failed reading detector ephemerides file 3 \n");
      abort();
    }
    fclose (data);
  } else {
    perror (filename);
    return 1;
  }

  double *cuDetSSB;
  cudaMalloc((void**)&cuDetSSB, 3*N*sizeof(double));
  // pci wtf? po co kopiować śmieci? 
  cudaMemcpy(cuDetSSB,DetSSB,3*N*sizeof(double), cudaMemcpyHostToDevice);

#ifdef HAVE_SINCOS
  sincos (phir, &sphir, &cphir);
  sincos (epsm, &sepsm, &cepsm);
#else
  sphir = sin (phir);
  cphir = cos (phir);
  sepsm = sin (epsm);
  cepsm = cos (epsm);
#endif
   
  //#mb not needed
  //t2 = (double *) calloc (Nv, sizeof (double));

  cosmodf = (double *) calloc (Nv, sizeof (double));
  sinmodf = (double *) calloc (Nv, sizeof (double));

  for (i=0; i<Nv; i++) {
    omrt = omr*i;

    //#mb not needed
    //t2[i] = sqr ((double)i);

    cosmodf[i] = cos (omrt);
    sinmodf[i] = sin (omrt);
  }

  // Because of frequency-domain filters, we search 
  // F-statistic in range (nmin+1, nmax) of data points 
  nmin = 2*NAV;
  nmax = nfft-2*NAV;

  coft = oms; 

  aa = (double *) calloc (Nv, sizeof (double));
  bb = (double *) calloc (Nv, sizeof (double));

  shft = (double *) calloc (N, sizeof (double));
  shftf = (double *) calloc (N, sizeof (double));
  xDatma = (complex double *) calloc (N, sizeof (complex double));
  xDatmb = (complex double *) calloc (N, sizeof (complex double));

  xa =(complex float *) malloc (4 * fftpad * nfft * sizeof (complex float));

  // FFT plans vary w.r.t the method of calculation: 
  // case INT (simple interpolation [interbinning] of a shorter Fourier transform)
  if(fftinterp==INT) { 
    xao = (complex float *) malloc (2 * nfft * sizeof (complex float));
    cufftPlanMany(&plan,1,&nfft,NULL,1, nfft,NULL,1,nfft, CUFFT_C2C, 2);
    // Case FFT (longer Fourier transforms, interpolation by zero padding)
  } else { 
    nfftf = fftpad*nfft ; 
    cufftPlanMany(&plan,1,&nfftf,NULL,1, nfftf,NULL,1,nfftf, CUFFT_C2C, 2);
  } 


  // These two plans below are used in the resampling 
  // procedure in JobCore() 
  // xDa = xa;
  
  cufftPlanMany(&pl_int,1,&nfft,NULL,1, nfft,NULL,1,nfft, CUFFT_C2C, 2);
  
  //rDa = xa;
  Ninterp = interpftpad*nfft;
  
  cufftPlanMany(&pl_inv,1,&Ninterp,NULL,1, Ninterp,NULL,1,Ninterp, CUFFT_C2C, 2);

  //#mb disabling state file for tests 
  /*
  if(hemi)
    sprintf (qname, "state_%02d_%03d%s_%d.dat", ident, band, label, hemi);
  else 
    sprintf (qname, "state_%02d_%03d%s.dat", ident, band, label);
  
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

  */ 

    pst = pmr[0];
    mst = mr[0];
    nst = nr[0];
    sst = spndr[0];
    FNum = 0;

  //#mb disabling state file for tests
  //} /* if state */
  
  cufftComplex *nakarcie1, *curDa;
  cudaMalloc((void**)&nakarcie1,4 * fftpad * nfft*sizeof(cufftComplex));
  cudaMalloc((void**)&curDa,4 * fftpad * nfft*sizeof(cufftComplex));
  cufftComplex *du, *d, *dl, *splvec;
  cudaMalloc((void**)&du,(Ninterp-2)*sizeof(cufftComplex));
  cudaMalloc((void**)&d,(Ninterp-2)*sizeof(cufftComplex));
  cudaMalloc((void**)&dl,(Ninterp-2)*sizeof(cufftComplex));
  cudaMalloc((void**)&splvec,(Ninterp-2)*sizeof(cufftComplex));
  int *trigcount;
  cudaMalloc((void**)&trigcount,sizeof(int));
  zapisz<<<(int)((float)(Ninterp-2)/256.0)+1,256>>>(trigcount,du,d,dl,Ninterp-2);

  float *humem, *humem_host;
  cudaMalloc((void**)&humem,2000*5*sizeof(float));
  humem_host=(float*)malloc(2000*5*sizeof(float));
  int triggers_cpu;

  double *cuaa, *cubb;
  cudaMalloc((void**)&cuaa,N*sizeof(double)); 
  cudaMalloc((void**)&cubb,N*sizeof(double)); 

  cufftComplex *cuxao;
  cudaMalloc((void**)&cuxao,2*fftpad*nfft*sizeof(cufftComplex));

  // Main loops 

  for (pm=pst; pm<=pmr[1]; pm++) {	// loop over hemispheres 
   
    sprintf (outname, "%s/triggers_%02d_%03d%s_%d.bin", prefix, ident,	\
	     band, label, pm);
    
    for (mm=mst; mm<=mr[1]; mm++) {	// 2 loops over 
      for (nn=nst; nn<=nr[1]; nn++) {	// sky positions 

    //#mb disabling state file for tests
    /* 	
	state = fopen (qname, "w");
	fprintf (state, "%d %d %d %d %d\n", pm, mm, nn, sst, FNum);
	fclose (state);
	*/ 

	sgnlv = JobCore(pm,		// hemisphere 
			mm,		// grid 'sky position' 		
			nn,		// other grid 'sky position' 
			sst,		// grid spindown coordinate
			spndr[1],	// spindown range limit 	
			M,		// grid generating matrix 
			cuDetSSB,		// ephemerides array 
			cuxDat,		// time-domain input data array 
			N,		// Number of data points
			Ninterp,	// interpftpad*nfft (for resampling)
			nfft,		// size of the FFT transform 
			//xDa,		// Array for resampling
			//rDa,		// Array for resampling
			xa,		// input array for plan
			cuxao,		// output array for plan
			pl_int,		// cufftHandle needed for resampling
			pl_inv,		// cufftHandle needed for resampling
					// (inverse transformation)
			plan,		// main cufftHandle 
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
			s0_flag,	// No-spindown flag
			cuF,
			nakarcie1, curDa, du, d, dl, splvec, cuaa, cubb, trigcount, humem 	//F stat on GPU
			);
	sst = spndr[0];

	

	/* if sgnlc */
	free (sgnlv);
      } /* for nn */
      nst = nr[0];
    } /* for mm */
    mst = mr[0];

    /* Add trigger parameters to a file */
    cudaMemcpy(&triggers_cpu, trigcount,sizeof(int),cudaMemcpyDeviceToHost);
    cudaMemcpy(humem_host,humem,triggers_cpu*5*sizeof(float),cudaMemcpyDeviceToHost);
    candid=fopen(outname,"a");
    fwrite(humem_host, sizeof(float),5*triggers_cpu,candid);
    fclose(candid);
    zerotrig<<<1,1>>>(trigcount);
    
    
  } /* for pm */
  
  
  //#mb disabling state file for tests
  /*   
  state = fopen (qname, "w");
  fclose (state);
  */ 

  free (aa);
  free (bb);
  free (shft);
  free (shftf);
  free (xDatma);
  free (xDatmb);
  free (xDat);
  free (DetSSB);


  //#mb things that are not needed
  // free (F);
  // free (t2);

  free (cosmodf);
  free (sinmodf);
  free (xa);
  free(humem_host);
  // This array is defined only in the case of fftinterp==INT
  if(fftinterp==INT) free (xao);
  free (M);
  cudaFree(cuxDat);
  cudaFree(cuDetSSB);
  cudaFree(nakarcie1);
  cudaFree(curDa);
  cudaFree(du);
  cudaFree(dl);
  cudaFree(d);
  cudaFree(splvec);
  cudaFree(cuaa);
  cudaFree(cubb);
  cudaFreeHost(F);
  cudaFree(cuF);
  cudaFree(trigcount);
  cudaFree(humem);
  cudaFree(cuxao);

#if TIMERS==1
  cudaThreadSynchronize();
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
  printf("\nJobNAllSky execution time: %lf sec.\n", 
	 (t1.tv_sec-t0.tv_sec)+(double)(t1.tv_nsec-t0.tv_nsec)/(double)NANO_INV
	 );
#endif

  return 0;
} /* JobNAllSky() */



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
    printf("ERROR: CUDA enabled device not found!\n");
    return(-1);
  }
  if (cdev < 0 || cdev >= deviceCount) {
    printf("\nWARNING: Device %d is not available! Trying device 0\n", cdev);
    cdev = 0;
  }

  printf("__________________________________CUDA devices___________________________________\n");
  printf("Set | ID |       Name       |   Gmem(B)   | Smem(B) | Cmem(B) | C.Cap. | Thr/bl |\n");
  
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
    printf(" %1d  | %16s | %11Zu | %7Zu | %7Zu |   %d.%d  | %6d |\n", 
           dev, deviceProp.name, deviceProp.totalGlobalMem, deviceProp.sharedMemPerBlock, 
           deviceProp.totalConstMem, deviceProp.major, deviceProp.minor, deviceProp.maxThreadsPerBlock );
  }
  printf("---------------------------------------------------------------------------------\n");
  
  /* force initialization */
  cudaThreadSynchronize();
  return(cdev);
}
