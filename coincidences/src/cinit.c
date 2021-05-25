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


/*    Command line options handling for coincidences   */ 
	
void handle_opts_coinc( Search_settings *sett,
			Command_line_opts_coinc *opts,
			int argc,
			char* argv[]) {
	
  opts->wd=NULL;
	
  strcpy (opts->prefix, TOSTR(PREFIX));
  //strcpy (opts->dtaprefix, TOSTR(DTAPREFIX));

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
  opts->scalef=4;
  opts->scales=4;
  opts->scaled=4;
  opts->scalea=4;

  sett->fpo = -1;
  
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
      // frequency band number
      {"band", required_argument, 0, 'b'},
      // fpo value
      {"fpo", required_argument, 0, 'p'},
      // data sampling time 
      {"dt", required_argument, 0, 't'},
      // a file with input files (triggers + grids)
      {"infiles", required_argument, 0, 'i'},
      // Location of the reference frame grid
      {"refgrid", required_argument, 0, 'g'},
      // Minimal number of coincidences recorded in the output  
      {"mincoin", required_argument, 0, 'm'},
      // Narrow down the frequency band (+- the center of band) 
      {"narrowdown", required_argument, 0, 'n'},
      // Signal-to-noise threshold cutoff  
      {"snrcutoff", required_argument, 0, 'c'},
      // number of days in the time-domain segment 
      {"nod", required_argument, 0, 'y'},
      // band overlap
      {"overlap", required_argument, 0, 'v'},
      {0, 0, 0, 0}
    };

    if (help_flag) {

      printf("polgraw-allsky periodic GWs: search for concidences among candidates\n");
      printf("Usage: ./coincidences -[switch1] <value1> -[switch2] <value2> ...\n") ;
      printf("Switches are:\n\n");
      printf("-output       Output directory (default is ./coinc-results)\n");
      printf("-shift        Cell shifts in fsda directions (4 digit number, e.g. 0101, default 0000)\n");
      printf("-scale        Cell scaling in fsda directions (coma separated, e.g. 32,8,4,4, default 4,4,4,4)\n");
      printf("-refr         Reference frame number\n");
      printf("-band         Band number\n");
      printf("-fpo          Reference band frequency fpo value\n");
      printf("-dt           Data sampling time dt (default value: 0.5)\n");
      printf("-infile       File containing the list of trigger and grid files\n");      
      printf("-refgrid      Location of the reference frame grid\n");
      printf("-mincoin      Minimal number of coincidences recorded\n");
      printf("-narrowdown   Narrow-down the frequency band (range [0, 0.5] +- around center)\n");
      printf("-nod          Number of days\n");
      printf("-snrcutoff    Signal-to-noise threshold cutoff (default value: 6)\n");
      printf("-overlap      Band overlap, fpo=10+(1-overlap)*band/(dt*2) ; obligatory if band is used\n\n");
      
      printf("Also:\n\n");
      printf("--help		This help\n");

      exit (0);
    }

    int option_index = 0;
    int c = getopt_long_only (argc, argv, "p:o:s:z:r:t:g:m:n:c:y:b:v:i:", long_options, &option_index);
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
	 opts->scalef = atoi(strtok(optarg,","));
	 opts->scales = atoi(strtok(NULL,","));
	 opts->scaled = atoi(strtok(NULL,","));
	 opts->scalea = atoi(strtok(NULL,","));      
	 break;
    case 'r':
	 opts->refr = atoi(optarg);
	 break;
    case 'o':
	 strcpy(opts->prefix, optarg);
	 break;
    case 't':
	 sett->dt = atof(optarg);
	 break;
    case 'i':
	 strcpy(opts->infile, optarg);
	 break;
    case 'g':
	 strcpy(opts->refgrid, optarg);
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
    case 'b':
	 opts->band = atoi(optarg);
	 break;
    case 'v':
	 opts->overlap = atof(optarg);
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

  if(!(opts->band)) { 
       printf("Band is not set... Exiting\n"); 
       exit(EXIT_FAILURE);
  }
  if(!(opts->overlap)) { 
       printf("Band overlap is not set... Exiting\n"); 
       exit(EXIT_FAILURE);
  }
  
  printf("Band=%04d  Overlap=%f\n", opts->band, opts->overlap);

  // hemi must be decoded from filename
  opts->hemi = -1;
  
  // Starting band frequency:
  // fpo_val is optionally read from the command line
  // Its initial value is set to -1
  if (!(sett->fpo >= 0)) {
       if (opts->band > 0 && opts->overlap >=0.) {
	    sett->fpo = 10. + (1. - opts->overlap)*opts->band*(0.5/sett->dt);
       } else {
	    printf("Band AND overlap or fpo must be specified!\n");
	    exit(EXIT_FAILURE);
       }
  }
  printf("The reference frequency fpo is %f\n", sett->fpo);

  printf("Cell scaling factors are: %d %d %d %d\n", opts->scalef, opts->scales,
	 opts->scaled, opts->scalea);
  
} // end of command line options handling: coincidences  



/* Manage grid matrix (read from grid.bin, find eigenvalues 
 * and eigenvectors) and reference GPS time from starting_time
 * (expected to be in the same directory)    
 */ 

void manage_grid_matrix( Search_settings *sett, char *gridfile ) {

    FILE *data;
    int gdim=4;
    int i, j;
#define DEBINV

    sett->M = (double *)calloc(gdim*gdim, sizeof (double));
    //sett->invM = (double *)calloc(gdim*gdim, sizeof (double));

    if ((data=fopen(gridfile, "r")) != NULL) {
	printf("Reading grid file: %s\n", gridfile);
	fread ((void *)&sett->fftpad, sizeof (int), 1, data);
	printf("fftpad from the grid file: %d\n", sett->fftpad); 
	fread ((void *)sett->M, sizeof(double), gdim*gdim, data);
	fclose (data);
    } else {
	perror(gridfile);
	exit(EXIT_FAILURE);
    }

#ifdef DEBINV
    printf("M=\n");
    for(i=0; i<gdim; i++) {
	for(j=0; j<gdim; j++) {
	    printf("%.12le ", *(sett->M + i*gdim + j) );
	}
	printf("\n");
    }
#endif    
    // Calculating the eigenvectors and eigenvalues 
    gsl_matrix_view m = gsl_matrix_view_array(sett->M, gdim, gdim);

    // for testing
    gsl_matrix *mcopy = gsl_matrix_alloc(gdim, gdim);
    gsl_matrix_memcpy(mcopy, &m.matrix);
    
    // Inverting the matrix 
    gsl_permutation *p = gsl_permutation_alloc(gdim);
    int s;
  
    // Compute the LU decomposition of matrix m
    gsl_linalg_LU_decomp(&m.matrix, p, &s);
  
    // Compute the inverse of the LU decomposition
    gsl_matrix *inv_m = gsl_matrix_alloc(gdim, gdim);
    gsl_linalg_LU_invert(&m.matrix, p, inv_m);
    gsl_permutation_free(p);

    //printf("invM=\n");
    for(i=0; i<gdim; i++) {
	for(j=0; j<gdim; j++) {
	    //*(sett->invM + i*gdim + j)  = gsl_matrix_get(inv_m, i, j);
	    sett->invM[i][j] = gsl_matrix_get(inv_m, i, j);
	    //printf("%.12le ", gsl_matrix_get(inv_m, i, j));
	}
	//printf("\n");
    }

#ifdef DEBINV
    // test 
    printf("one=\n");
    gsl_matrix *ident = gsl_matrix_calloc(gdim, gdim);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mcopy, inv_m, 0.0, ident);
    for(i=0; i<gdim; i++) {
	for(j=0; j<gdim; j++) {
	    printf("%.12le ", gsl_matrix_get(ident, i, j));
	}
	printf("\n");
    }
#endif

} // end of manage grid matrix  
