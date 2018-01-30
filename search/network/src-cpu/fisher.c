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
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <time.h>
#include <dirent.h>

#include "auxi.h"
#include "settings.h"
#include "struct.h"
#include "init.h"

#ifndef CODEVER
#define CODEVER unknown
#endif

#ifndef DTAPREFIX
#define DTAPREFIX .
#endif


int fisher (
  Search_settings *sett,
  Command_line_opts *opts,
  Aux_arrays *aux);  

void cleanup_fisher(
	Search_settings *sett,
	Command_line_opts *opts,
	Aux_arrays *aux);


void handle_opts_fisher (
  Search_settings *sett,
  Command_line_opts *opts,
  int argc,
  char* argv[]); 


int main (int argc, char* argv[]) {

  Command_line_opts opts;
  Search_settings sett;
  Search_range s_range; 
  Aux_arrays aux_arr;
  double *F; 			  // F-statistic array
  int i; 

#define QUOTE(x) #x
#define STR(macro) QUOTE(macro)
#define CVSTR STR(CODEVER)

  printf("Code version : " CVSTR "\n");

  // Command line options 
  handle_opts_fisher(&sett, &opts, argc, argv);  
	
  // Search settings
  search_settings(&sett); 

  // Detector network settings
  detectors_settings(&sett, &opts); 

  // Array initialization and reading the ephemerids 
  init_arrays(&sett, &opts, &aux_arr, &F);

  // Amplitude modulation functions for each detector  
  for(i=0; i<sett.nifo; i++)   
    rogcvir(&ifo[i]); 


  // Main job - Fisher matrix calculations   
  for(i=0; i<sett.nifo; i++) 
    fisher(&sett, &opts, &aux_arr);


  // Cleanup & memory free 
  cleanup_fisher(&sett, &opts, &aux_arr);

  return 0; 
	
}

int fisher(Search_settings *sett,
           Command_line_opts *opts, 
           Aux_arrays *aux) { 


  int i, j, k, n; 

  // Signal parameters: f, fdot, delta, alpha, a1, a2, a3, a4
  // (see Phys. Rev. D 82, 022005 2010, Eqs. 2.13a-d) 
  double sgnlo[8]; 

  FILE *data; 

  // Reading signal parameters 
  if ((data=fopen (opts->addsig, "r")) != NULL) {

    for(i=0; i<8; i++)
      fscanf(data, "%le",i+sgnlo); 
    
    fclose (data);
                 
  } else {
    perror (opts->addsig);
  }

  
  double ma[8][8], mFl[8][8]; 
  double omega0, omega1, domega, a1, a2, a3, a4; 
  double sindelt, cosdelt, sinalt, cosalt; 

  for(k=0; k<8; k++) 
    for(j=0; j<8; j++) { 
      ma[k][j] = 0; 
      mFl[k][j] = 0; 
    }

  omega0 = sgnlo[0]/sett->dt; 
  omega1 = sgnlo[1]; 

  domega = 2*M_PI*sett->fpo*sett->dt; 

  sindelt = sin(sgnlo[2]); 
  cosdelt = cos(sgnlo[2]); 
  sinalt  = sin(sgnlo[3]); 
  cosalt  = cos(sgnlo[3]); 
  
  a1 = sgnlo[4]; 
  a2 = sgnlo[5]; 
  a3 = sgnlo[6]; 
  a4 = sgnlo[7]; 

  // Loop for each detector 
  for(n=0; n<sett->nifo; ++n) { 

    /* Amplitude modulation functions aa and bb 
     * for each detector (in signal sub-struct 
     * of _detector, ifo[n].sig.aa, ifo[n].sig.bb) 
     */

    modvir(sinalt, cosalt, sindelt, cosdelt,
	   sett->N, &ifo[n], aux);


   for(i=0; i<sett->N; ++i) {

      double dpdf, dpds, dpdd, dpda, xet, yet, zet, a, b, t; 

      xet = ifo[n].sig.DetSSB[i*3]; 
      yet = ifo[n].sig.DetSSB[i*3+1]; 
      zet = ifo[n].sig.DetSSB[i*3+2]; 

      a = ifo[n].sig.aa[i]; 
      b = ifo[n].sig.bb[i]; 

      t = i*sett->dt; 
      
      dpdf = t + xet*cosalt*cosdelt + yet*cosdelt*sinalt + zet*sindelt;  
      dpds = t*t;   
      dpdd = (domega + omega0)*(zet*cosdelt - (xet*cosalt + yet*sinalt)*sindelt); 
      dpda = (domega + omega0)*cosdelt*(yet*cosalt - xet*sinalt); 


      // Matrix element calculation 
      ma[0][0] = (a*a*a1*a1*dpdf*dpdf)/2 + (a*a*a3*a3*dpdf*dpdf)/2. + a*a1*a2*b*dpdf*dpdf + a*a3*a4*b*dpdf*dpdf + (a2*a2*b*b*dpdf*dpdf)/2. + (a4*a4*b*b*dpdf*dpdf)/2.; 
 
      ma[0][1] = (a*a*a1*a1*dpds*dpdf)/2 + (a*a*a3*a3*dpds*dpdf)/2. + a*a1*a2*b*dpds*dpdf + a*a3*a4*b*dpds*dpdf + (a2*a2*b*b*dpds*dpdf)/2. + (a4*a4*b*b*dpds*dpdf)/2.;
 
      ma[0][2] = (a*a*a1*a1*dpdf*dpdd)/2 + (a*a*a3*a3*dpdf*dpdd)/2. + a*a1*a2*b*dpdf*dpdd + a*a3*a4*b*dpdf*dpdd + (a2*a2*b*b*dpdf*dpdd)/2. + (a4*a4*b*b*dpdf*dpdd)/2.; 

      ma[0][3] = (a*a*a1*a1*dpdf*dpda)/2. + (a*a*a3*a3*dpdf*dpda)/2. + a*a1*a2*b*dpdf*dpda + a*a3*a4*b*dpdf*dpda + (a2*a2*b*b*dpdf*dpda)/2. + (a4*a4*b*b*dpdf*dpda)/2.;  

      ma[0][4] = (a*a*a3*dpdf)/2. + (a*a4*b*dpdf)/2.;  

      ma[0][5] = (a*a3*b*dpdf)/2. + (a4*b*b*dpdf)/2.;
 
      ma[0][6] = -(a*a*a1*dpdf)/2. - (a*a2*b*dpdf)/2.;
 
      ma[0][7] = -(a*a1*b*dpdf)/2. - (a2*b*b*dpdf)/2.; 

      ma[1][1] = (a*a*a1*a1*dpds*dpds)/2. + (a*a*a3*a3*dpds*dpds)/2. + a*a1*a2*b*dpds*dpds + a*a3*a4*b*dpds*dpds + (a2*a2*b*b*dpds*dpds)/2. + (a4*a4*b*b*dpds*dpds)/2.; 

      ma[1][2] = (a*a*a1*a1*dpds*dpdd)/2. + (a*a*a3*a3*dpds*dpdd)/2. + a*a1*a2*b*dpds*dpdd + a*a3*a4*b*dpds*dpdd + (a2*a2*b*b*dpds*dpdd)/2. + (a4*a4*b*b*dpds*dpdd)/2.;
 
      ma[1][3] = (a*a*a1*a1*dpds*dpda)/2. + (a*a*a3*a3*dpds*dpda)/2. + a*a1*a2*b*dpds*dpda + a*a3*a4*b*dpds*dpda + (a2*a2*b*b*dpds*dpda)/2. + (a4*a4*b*b*dpds*dpda)/2.;

      ma[1][4] = (a*a*a3*dpds)/2. + (a*a4*b*dpds)/2.;
 
      ma[1][5] = (a*a3*b*dpds)/2. + (a4*b*b*dpds)/2.; 
 
      ma[1][6] = -(a*a*a1*dpds)/2. - (a*a2*b*dpds)/2.;
 
      ma[1][7] = -(a*a1*b*dpds)/2. - (a2*b*b*dpds)/2.;
 
      ma[2][2] = (a*a*a1*a1*dpdd*dpdd)/2. + (a*a*a3*a3*dpdd*dpdd)/2. + a*a1*a2*b*dpdd*dpdd + a*a3*a4*b*dpdd*dpdd + (a2*a2*b*b*dpdd*dpdd)/2. + (a4*a4*b*b*dpdd*dpdd)/2.;
 
      ma[2][3] = (a*a*a1*a1*dpdd*dpda)/2. + (a*a*a3*a3*dpdd*dpda)/2. + a*a1*a2*b*dpdd*dpda + a*a3*a4*b*dpdd*dpda + (a2*a2*b*b*dpdd*dpda)/2. + (a4*a4*b*b*dpdd*dpda)/2.;
 
      ma[2][4] = (a*a*a3*dpdd)/2. + (a*a4*b*dpdd)/2.;
 
      ma[2][5] = (a*a3*b*dpdd)/2. + (a4*b*b*dpdd)/2.;
 
      ma[2][6] = -(a*a*a1*dpdd)/2. - (a*a2*b*dpdd)/2.;
 
      ma[2][7] = -(a*a1*b*dpdd)/2. - (a2*b*b*dpdd)/2.;
 
      ma[3][3] = (a*a*a1*a1*dpda*dpda)/2. + (a*a*a3*a3*dpda*dpda)/2. + a*a1*a2*b*dpda*dpda + a*a3*a4*b*dpda*dpda + (a2*a2*b*b*dpda*dpda)/2. + (a4*a4*b*b*dpda*dpda)/2.;
 
      ma[3][4] = (a*a*a3*dpda)/2. + (a*a4*b*dpda)/2.;
 
      ma[3][5] = (a*a3*b*dpda)/2. + (a4*b*b*dpda)/2.;
 
      ma[3][6] = -(a*a*a1*dpda)/2. - (a*a2*b*dpda)/2.;
 
      ma[3][7] = -(a*a1*b*dpda)/2. - (a2*b*b*dpda)/2.;
 
      ma[4][4] = a*a/2.;
 
      ma[4][5] = (a*b)/2.;
 
      ma[4][6] = 0;
 
      ma[4][7] = 0;
 
      ma[5][5] = b*b/2.;
 
      ma[5][6] = 0;
 
      ma[5][7] = 0;
 
      ma[6][6] = a*a/2.;
 
      ma[6][7] = (a*b)/2.;
 
      ma[7][7] = b*b/2.;


      // mFl 
      for(k=0; k<8; k++) 
        for(j=0; j<8; j++) 
          mFl[k][j] += ma[k][j];

    } 

    // Symmetrize mFl 
    for(k=1; k<8; k++) 
      for(j=0; j<k; j++)
        mFl[k][j] = mFl[j][k]; 

    for(k=0; k<8; k++) { 
      printf("[");
      for(j=0; j<8; j++)
        printf("%.16f, ", 0.5*mFl[k][j]);
      printf("],\n");
    }


    printf("\nSubmatrices:\n"); 

    // Calculate the inverse 
    int s; 
    int x = 0, y = 4, nn = 4; 
    gsl_matrix *m = gsl_matrix_alloc (nn, nn);
    gsl_matrix *m2 = gsl_matrix_alloc (nn, nn);
    gsl_matrix *mtest = gsl_matrix_alloc (nn, nn);


    gsl_matrix *inverse = gsl_matrix_alloc (nn, nn);
    gsl_permutation *perm = gsl_permutation_alloc (nn);
  
    // Fill the matrix
    printf("\n1-4 4x4 part\n"); 
    for(k=x; k<y; k++) {  
      printf("["); 
      for(j=x; j<y; j++) { 
        gsl_matrix_set (m, k, j, 0.5*mFl[k][j]);
        gsl_matrix_set (m2, k, j, 0.5*mFl[k][j]);
        printf("%.16f, ", gsl_matrix_get (m, k, j));   
      } 
      printf("],\n"); 
    } 

    // Make LU decomposition of matrix m
    gsl_linalg_LU_decomp (m, perm, &s);
  
    // Invert the matrix m
    gsl_linalg_LU_invert (m, perm, inverse);
 
    printf("\nIts inverse:\n"); 

    for(k=x; k<y; k++) { 
      for(j=x; j<y; j++)  
        printf("%.6e ", gsl_matrix_get (inverse, k, j));   
      printf("\n"); 
    } 

    printf("\nDiagonal values:\n"); 

    for(k=x; k<y; k++) { 
        printf("%.6e ", gsl_matrix_get (inverse, k, k));   
    } 

    printf("\n\nTest:\n"); 

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, m2, inverse, 0.0, mtest);

    for(k=x; k<y; k++) { 
      for(j=x; j<y; j++)  
        printf("%.6e ", gsl_matrix_get (mtest, k, j));   
      printf("\n"); 
    } 

    // Calculate the inverse 
    x = 4; 
    y = 8; 

    // Fill the matrix
    printf("\n5-8 4x4 part\n");
    for(k=x; k<y; k++) {  
      for(j=x; j<y; j++) { 
        gsl_matrix_set (m, k-nn, j-nn, 0.5*mFl[k][j]);
        printf("%.6f ", gsl_matrix_get (m, k-nn, j-nn));   
      } 
      printf("\n"); 
    } 

    // Make LU decomposition of matrix m
    gsl_linalg_LU_decomp (m, perm, &s);
  
    // Invert the matrix m
    gsl_linalg_LU_invert (m, perm, inverse);
 
    printf("\nIts inverse:\n");

    for(k=x; k<y; k++) { 
      for(j=x; j<y; j++)  
        printf("%.6e ", gsl_matrix_get (inverse, k-nn, j-nn));   
      printf("\n"); 
    } 

    printf("\nDiagonal values:\n");

    for(k=x; k<y; k++) { 
        printf("%.6e ", gsl_matrix_get (inverse, k-nn, k-nn));   
    } 

    printf("\n"); 

    gsl_matrix_free(m); 
    gsl_matrix_free(m2); 
    gsl_matrix_free(mtest); 
    gsl_matrix_free(inverse); 
    gsl_permutation_free(perm); 



  } 

  return 0; 

} 


  /* Cleanup & memory free (calculation 
   * of the Fisher matrix) 
	 */

void cleanup_fisher(
	Search_settings *sett,
	Command_line_opts *opts,
	Aux_arrays *aux) {

  int i; 

  for(i=0; i<sett->nifo; i++) {
    free(ifo[i].sig.xDat);
    free(ifo[i].sig.DetSSB);
    free(ifo[i].sig.aa);
    free(ifo[i].sig.bb);
  } 
	
  free(aux->sinmodf);
  free(aux->cosmodf);
  free(aux->t2);


} // end of cleanup & memory free 



/*  Command line options handling: fisher 
 */ 

void handle_opts_fisher( Search_settings *sett, 
		  Command_line_opts *opts,
		  int argc, 
		  char* argv[]) {

  strcpy (opts->prefix, TOSTR(PREFIX));
  strcpy (opts->dtaprefix, TOSTR(DTAPREFIX));

  opts->usedet[0]   = '\0';
  opts->addsig[0]   = '\0';
	
  // Initial value of starting frequency set to a negative quantity. 
  // If this is not changed by the command line value, fpo is calculated 
  // from the band number b (fpo = fpo = fstart + 0.96875*b/(2dt))
  sett->fpo = -1;

  // Default initial value of the data sampling time 
  sett->dt = 0.5; 

  // Initial value of the number of days is set to 0
  sett->nod = 0; 

  opts->help_flag=0;

  static int help_flag=0;

  // Reading arguments 

  while (1) {
    static struct option long_options[] = {
      {"help", no_argument, &help_flag, 1},
      // frame number
      {"ident", required_argument, 0, 'i'},
      // frequency band number
      {"band", required_argument, 0, 'b'},
      // input data directory
      {"data", required_argument, 0, 'd'},
      // fpo value
      {"fpo", required_argument, 0, 'p'},
      // add signal parameters
      {"addsig", required_argument, 0, 'x'},
      // number of days in the time-domain segment 
      {"nod", required_argument, 0, 'y'},
      // which detectors to use
      {"usedet", required_argument, 0, 'u'}, 
      // data sampling time 
      {"dt", required_argument, 0, 's'},
      {0, 0, 0, 0}
    };

    if (help_flag) {

      printf("polgraw-allsky periodic GWs: Fisher matrix calculation\n");
      printf("Usage: ./fisher -[switch1] <value1> -[switch2] <value2> ...\n") ;
      printf("Switches are:\n\n");
      printf("-data         Data directory (default is .)\n");
      printf("-ident        Frame number\n");
      printf("-band         Band number\n");
      printf("-fpo          Reference band frequency fpo value\n");
      printf("-dt           Data sampling time dt (default value: 0.5)\n");
      printf("-usedet       Use only detectors from string (default is use all available)\n");
      printf("-addsig       Add signal with parameters from <file>\n");
      printf("-nod          Number of days\n\n");


      printf("Also:\n\n");
      printf("--help            This help\n");

      exit(EXIT_SUCCESS);
    }

    int option_index = 0;
    int c = getopt_long_only(argc, argv, "i:b:dp:x:y:s:u:", 
			     long_options, &option_index);
    if (c == -1)
      break;

    switch (c) {
    case 'i':
      opts->ident = atoi (optarg);
      break;
    case 'b':
      opts->band = atoi(optarg);
      break;
    case 'd':
      strcpy(opts->dtaprefix, optarg);
      break;
    case 'p':
      sett->fpo = atof(optarg);
      break;
    case 'x':
      strcpy(opts->addsig, optarg);
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

  printf("Input data directory is %s\n", opts->dtaprefix);
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

  if (strlen(opts->addsig))
    printf ("Adding signal from '%s'\n", opts->addsig);

} // end of command line options handling: fisher  

