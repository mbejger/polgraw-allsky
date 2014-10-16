#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gsl/gsl_linalg.h>
#include <getopt.h>

#include "settings.h"

static int help_flag=0;

int main (int argc, char *argv[]) {

  size_t i; 
  int c, fftpad, freadstatus, ident=0, band=0, fact_sp1=1; 
  char filename[512], dtaprefix[512];
  double *M, *gamrn, *Mn; 
  double cputime=0, fpo;
  FILE *data;

  // Initial value of starting frequency 
  // set to a negative quantity. If this is not 
  // changed by the command line value 
  // fpo is calculated from the band number b.
  double fpo_val = -1; 

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
      {"fpo value", required_argument, 0, 'p'},
      // CPU time 
      {"cpu time", required_argument, 0, 't'},
      {0, 0, 0, 0}
    };

  if(help_flag) {

     printf("*** Continuous GW benchmark CPU time estimator (polgraw-allsky) ***\n"); 
     printf("Usage: ./bench -[switch1] <value1> -[switch2] <value2> ...\n") ;
     printf("Switches are:\n\n");
     printf("-d	Data directory (default is .)\n"); 
     printf("-i	Frame number\n"); 
     printf("-b	Band number\n"); 
     printf("-t CPU time\n");
     printf("Also:\n\n"); 
     printf("--help		This help\n"); 		

     exit (0);

  }

    int option_index = 0;
    c = getopt_long (argc, argv, "i:b:d:t:", long_options, &option_index);
    if (c == -1)
      break;

    switch (c) {
    case 'i':
      ident = atoi (optarg);
      break;
    case 'b':
      band = atoi (optarg);
      break;
    case 'd':
      strcpy (dtaprefix, optarg);
      break;
    case 't':
      cputime = atof (optarg);
      break;
    case '?':
      break;
    default: break ; 
    } /* switch c */
  } /* while 1 */

  printf ("Data directory is %s\n", dtaprefix);
  printf ("Frame number is %d\n", ident);
  printf ("Band number is %d\n", band);    

  M = (double *) calloc (16, sizeof (double));
  gamrn = (double *) calloc (16, sizeof (double));
  Mn = (double *) calloc (16, sizeof (double));

  sprintf (filename, "%s/%03d/grid.bin", dtaprefix, ident);
  if ((data=fopen (filename, "r")) != NULL) {
    // fftpad: used to zero padding to fftpad*nfft data points 
    freadstatus = fread ((void *)&fftpad, sizeof (int), 1, data);
    // M: vector of 16 components consisting of 4 rows 
    // of 4x4 grid-generating matrix  
    freadstatus = fread ((void *)M, sizeof (double), 16, data);
    // gamrn: normalised Fisher matrix 
    freadstatus = fread ((void *)gamrn, sizeof (double), 16, data);
    // Mn: normalised grid-generating matrix 
    freadstatus = fread ((void *)Mn, sizeof (double), 16, data);
    fclose (data);

	if(freadstatus) {}
    
    //#mb for tests: 
    printf("Grid read from the '%s' file\n", filename);

    //#mb warning: this fftpad is just for check (the value 
	// with which the grid was generated)
	// It is overwritten as soon as settings() is called 
	printf("fftpad from grid.bin: %d\n", fftpad); 
    printf("Matrix M from grid.bin:\n"); 

	printf("%e %e %e %e\n", M[0], M[1], M[2], M[3]);
	printf("%e %e %e %e\n", M[4], M[5], M[6], M[7]);
	printf("%e %e %e %e\n", M[8], M[9], M[10], M[11]);
	printf("%e %e %e %e\n", M[12], M[13], M[14], M[15]);   

    printf("Matrix Mn from grid.bin:\n");

	printf("%e %e %e %e\n", Mn[0], Mn[1], Mn[2], Mn[3]);
	printf("%e %e %e %e\n", Mn[4], Mn[5], Mn[6], Mn[7]);
	printf("%e %e %e %e\n", Mn[8], Mn[9], Mn[10], Mn[11]);
	printf("%e %e %e %e\n", Mn[12], Mn[13], Mn[14], Mn[15]);   

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

  printf("The reference frequency (fpo) is %f\n", fpo); 

  // Detector, ephemerides, constants (second variable 
  // is the detector - not used here, Virgo as default choice) 	
  settings (fpo, "V1");

  // Observation time [s]  
  double To = nod*SIDday;

  double f_min = fpo; 
  double f_max = fpo + B; 

  // Factorial of s+1
  for (i = 1; i <= s+1; i++)
  	fact_sp1 *= i;

  // Paper III Eq.(71) x (2*pi*To)^2
  double cof = 4*M_PI*M_PI*pow(2, 2*s + 1)*pow(M_PI, s+2)/((s+3.)*fact_sp1); 
  double tof = pow(To, s+3)*pow(To/tau_min, s*(s+1)/2.);

  // Volume of the intrinsic parameter space
  double vol = cof*tof*(pow(f_max, s+3) - pow(f_min, s+3));

  // Paper III Eq.(92) x (2*pi*To)^2
  double coff = 4*M_PI*M_PI*pow(2, 2*s)*pow(M_PI, s+1)/fact_sp1;
  double toff = pow(To, s+2)*pow(To/tau_min, s*(s+1)/2.);

  // Volume of the intrinsic parameter space
  // with frequency parameter taken out
  double volf = coff*toff*pow(f_max, s+2);

  // Area of the sky
  double vols = 4*M_PI*pow(M_PI*To*f_max, 2);

  int s;

  gsl_matrix_view mn = gsl_matrix_view_array (Mn, 4, 4);
  gsl_matrix_view mnr = gsl_matrix_submatrix (&mn.matrix, 1, 1, 3, 3); 

  // gsl_matrix_fprintf(stdout, &mnr.matrix, "%e");

  gsl_permutation *p = gsl_permutation_alloc (3);
  gsl_linalg_LU_decomp (&mnr.matrix, p, &s);

  double vf = fabs(gsl_linalg_LU_det(&mnr.matrix, s));

  // No. of filters
  double Nf = round(volf/vf); 

  printf("No. of templates: %f\nTime per template: %f\n", Nf, cputime/Nf); 

  return 0;
 
} 
