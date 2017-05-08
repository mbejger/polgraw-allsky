#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_linalg.h>
#include <getopt.h>

#include "auxi.h"
#include "settings.h"
#include "struct.h"
#include "init.h"

static int help_flag=0;

int main (int argc, char *argv[]) {

  Command_line_opts opts;
  Search_settings sett;

  size_t i; 
  int c, fftpad, freadstatus, ident=0, band=0, fact_sp1=1; 
  char filename[512], dtaprefix[512];
  double *M, *gamrn, *Mn; 
  double cputime=0, fpo;
  FILE *data;

  // Initial value of starting frequency set to a negative quantity. 
  // If this is not changed by the command line value, fpo is calculated 
  // from the band number b (fpo = fpo = fstart + 0.96875*b/(2dt))
  sett.fpo = -1;

  // Default initial value of the data sampling time 
  sett.dt = 2 ; 

  // Initial value of the number of days is set to 0
  sett.nod = 0; 

  opts.help_flag=0;
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
      // cpu time 
      {"time", required_argument, 0, 't'},
      // fpo value
      {"fpo", required_argument, 0, 'p'},
      // number of days in the time-domain segment 
      {"nod", required_argument, 0, 'y'},
      // data sampling time 
      {"dt", required_argument, 0, 's'},
      {0, 0, 0, 0}
    };

    if (help_flag) {

      printf("polgraw-allsky periodic GWs benchmark CPU time estimator\n");
      printf("Usage: ./bench -[switch1] <value1> -[switch2] <value2> ...\n") ;
      printf("Switches are:\n\n");
      printf("-d, -data         Data directory (default is .)\n");
      printf("-i, -ident        Frame number\n");
      printf("-b, -band         Band number\n");
      printf("-t, -time         CPU time\n");
      printf("-p, -fpo          Reference band frequency fpo value\n");
      printf("-s, -dt           data sampling time dt (default value: 0.5)\n");
      printf("-y, -nod          Number of days\n\n");


      printf("Also:\n\n");
      printf("--help            This help\n");

      exit(EXIT_SUCCESS);
    }

    int option_index = 0;
    int c = getopt_long_only(argc, argv, "i:b:d:t:p:y:s:", 
			     long_options, &option_index);
    if (c == -1)
      break;

    switch (c) {
    case 'i':
      opts.ident = atoi (optarg);
      break;
    case 't':
      cputime = atof(optarg);
      break;
    case 'b':
      opts.band = atoi(optarg);
      break;
    case 'd':
      strcpy(opts.dtaprefix, optarg);
      break;
    case 'p':
      sett.fpo = atof(optarg);
      break;
    case 'y':
      sett.nod = atoi(optarg);
      break;
    case 's':
      sett.dt = atof(optarg);
      break;

    case '?':
      break;
    default:
      break ;
    } /* switch c */
  } /* while 1 */


  // Check if sett.nod was set up, if not, exit
  if(!(sett.nod)) { 
    printf("Number of days not set... Exiting\n"); 
    exit(EXIT_FAILURE); 
  } 

  printf("Number of days is %d\n", sett.nod); 

  printf("Input data directory is %s\n", opts.dtaprefix);
  printf("Frame and band numbers are %d and %d\n", opts.ident, opts.band);

  // Starting band frequency:
  // fpo_val is optionally read from the command line
  // Its initial value is set to -1
  if(!(sett.fpo >= 0))
    sett.fpo = 10. + 0.96875*opts.band*(0.5/sett.dt);

  printf("The reference frequency fpo is %f\n", sett.fpo);
  printf("The data sampling time dt is %f\n", sett.dt); 

  if(cputime)
    printf ("CPU time is %lf [s]\n", cputime);
	
  // Search settings
  search_settings(&sett); 
  Mn = (double *) calloc (16, sizeof (double));

  sprintf (filename, "%s/%03d/grid.bin", opts.dtaprefix, opts.ident);

  if ((data=fopen (filename, "rb")) != NULL) {

    // skipping fftpad (1 int) and M and gamrn matrices (2x16 doubles) 
    fseek(data, sizeof(int) + 32*sizeof(double), SEEK_SET);
    // Mn: normalised grid-generating matrix 
    fread ((void *)Mn, sizeof (double), 16, data);

    fclose (data);
  } else {
    perror (filename);
    exit(EXIT_FAILURE);
  }

  printf("Matrix Mn from grid.bin:\n");

	printf("%e %e %e %e\n", Mn[0], Mn[1], Mn[2], Mn[3]);
	printf("%e %e %e %e\n", Mn[4], Mn[5], Mn[6], Mn[7]);
	printf("%e %e %e %e\n", Mn[8], Mn[9], Mn[10], Mn[11]);
	printf("%e %e %e %e\n", Mn[12], Mn[13], Mn[14], Mn[15]);   

  // Observation time [s]  
  double To = sett.nod*C_SIDDAY;

  double f_min = sett.fpo; 
  double f_max = sett.fpo + sett.B; 

  double tau_min = 1000.*C_YEARSEC ; 

  // number of spindowns 
  int s = 1; 

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

  int signum;
  gsl_matrix_view mn = gsl_matrix_view_array (Mn, 4, 4);
  gsl_matrix_view mnr = gsl_matrix_submatrix (&mn.matrix, 1, 1, 3, 3); 

  gsl_permutation *p = gsl_permutation_alloc (3);
  gsl_linalg_LU_decomp (&mnr.matrix, p, &signum);

  double vf = fabs(gsl_linalg_LU_det(&mnr.matrix, signum));

  // No. of filters
  double Nf = round(volf/vf); 

  printf("\nNo. of templates: %f\nTime per template: %f\n", Nf, cputime/Nf); 

  // free 
  free(Mn); 

  return 0;
 
} 
