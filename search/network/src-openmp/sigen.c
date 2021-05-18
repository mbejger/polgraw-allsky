#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>

#include "auxi.h"
#include "settings.h"
#include "struct.h"

static int help_flag=0;

double get_rand() ;

int main(int argc, char *argv[]) {
  
  int i, numl=0, freq_line_check, c, pm, gsize=2, band=0, reffr, nfrinband; 
  char filename[512], dtaprefix[512], *wd=NULL ; 
  double amp=0, snr=0;
  double freql[32768], linew[32768], sgnlo[8], rrn[2], 
      dvr, fpo_val, be1, be2, fr, lw, freqlp, freqlm,
      sepsm, cepsm, sinalt, cosalt, sindelt, cosdelt,
      iota, ph_o, psik, hop, hoc, overlap;
  
  FILE *data ;
  
  Search_settings sett;

  //#mb fftpad value=1 (default for some time now) 
  sett.fftpad = 1; 

  // Default data sampling time (s)  
  sett.dt = 0.5; 

  // Default reference time frame (frame in which the frequency 
  // of the signal is as selected below, not spun-down/up) is 1
  reffr = 1; 
  
  // Initial value of starting frequency 
  // set to a negative quantity. If this is not 
  // changed by the command line value 
  // fpo is calculated from the band number b.
  fpo_val = -1;

  // Initial value of the number of days is set to 0
  sett.nod = 0; 
  
  overlap = -1.;
  nfrinband = 0;
  
  while (1) {
    static struct option long_options[] = {
      {"help", no_argument, &help_flag, 1},			
      // GW amplitude 
      {"amp", required_argument, 0, 'a'},
      // SNR
      {"snr", required_argument, 0, 'n'},
      // frequency band number 
      {"band", required_argument, 0, 'b'},
      // band overlap 
      {"overlap", required_argument, 0, 'v'},
      // change directory parameter
      {"cwd", required_argument, 0, 'c'},
      // fpo value
      {"fpo", required_argument, 0, 'p'},
      // grid search range 
      {"gsize", required_argument, 0, 'g'},
      // data sampling time 
      {"dt", required_argument, 0, 's'},
      // reference frame () 
      {"reffr", required_argument, 0, 'r'},
      // +/- frames from reference to keep the signal inband
      {"nfrinband", required_argument, 0, 'i'},
      // number of days in the time-domain segment 
      {"nod", required_argument, 0, 'y'}, 
      {0, 0, 0, 0}
    };

    if(help_flag) {
	  
      printf("*** Software injection parameters - cont. GW signal ***\n"); 
      printf("Usage: ./sigen -[switch1] <value1> -[switch2] <value2> ...\n") ;
      printf("Switches are:\n\n"); 
      printf("-amp      GW amplitude of the injected (mutually exclusive with -snr)\n"); 
      printf("-snr      SNR of the injected signal (mutually exclusive with -amp)\n"); 
      printf("-band     Band number\n"); 
      printf("-overlap  Band overlap\n"); 
      printf("-cwd      Change to directory <dir>\n");     
      printf("-fpo      fpo (starting frequency) value\n");
      printf("-gsize    Grid search range (default value: 2)\n");
      printf("-dt       Data sampling time dt (default value: 0.5)\n");
      printf("-nod      Number of days\n");
      printf("-reffr    Reference frame (default value: 1)\n");
      printf("-nfrinband  Force the signal to stay inband for +/- frames \n\n");
      
      printf("--help    This help\n"); 		
      exit (0);
    }

    int option_index = 0;
    c = getopt_long_only (argc, argv, "a:n:b:v:c:p:g:s:r:i:y:", long_options, &option_index);
    
    if (c == -1)
      break;
    switch (c) {
     case 'a':
      amp  = atof (optarg);
      break; 
    case 'n':
      snr  = atof (optarg);
      break; 
    case 'b':
      band = atoi (optarg);
      break;
    case 'v':
      overlap = atof (optarg);
      break;
    case 'c':
      wd = (char *) malloc (1+strlen(optarg));
      strcpy (wd, optarg);
      break;      
    case 'g':
      gsize = atoi (optarg);
      break;
    case 'p':
      fpo_val = atof(optarg);
      break;
    case 's':
      sett.dt = atof(optarg);
      break;
    case 'r':
      reffr = atof(optarg);
      break;
    case 'i':
      nfrinband = atoi(optarg);
      break;
    case 'y':
      sett.nod = atoi(optarg);
      break;
    case '?':
      break;
    default: break ; 

    } /* switch c */
  } /* while 1 */

  // Check if sett->nod was set up, if not, exit
  if(!(sett.nod)) { 
    printf("Number of days not set... Exiting\n"); 
    exit(EXIT_FAILURE); 
  } 
  

  if (wd) {
    printf ("Changing working directory to %s\n", wd);
    if (chdir(wd)) {
      perror (wd);
      abort ();
    }
  }

  // Check if the options are consistent with each other 
  if(!(amp || snr)) { 
    printf("Options -amp or -snr not given. Exiting...\n"); 
    exit(0); 
  } else if(amp && snr) { 
    printf("Options -amp and -snr are mutually exclusive. Exiting...\n"); 
    exit(0); 
  } 

  // Starting band frequency: 
  // fpo_val is optionally read from the command line
  // Its initial value is set to -1 
  if(fpo_val>=0)
      sett.fpo = fpo_val;
  else
      if (band > 0 && overlap >=0.) {
	  sett.fpo = 10. + (1.-overlap)*band*(0.5/sett.dt);
      } else {
	  printf("Band AND overlap or fpo must be specified!\n");
	  exit(EXIT_FAILURE);
      }

  // Search settings
  search_settings(&sett);
  
  fprintf(stderr, "Band number is %04d\n", band);
  fprintf(stderr, "The reference frequency fpo is %f\n", sett.fpo);
  fprintf(stderr, "The data sampling time dt is %f\n", sett.dt);
  fprintf(stderr, "The reference time frame is %d\n", reffr);

  //#mb Changing the limits near the bands border 
  //#mb For O3 UL simulation run 
  rrn[0] = M_PI/20.; 
  rrn[1] = M_PI - rrn[0]; 

  // Random signal parameters 
  //------------------------- 
  double f1,f2;
  int inband=0;
  do {
      // Frequency derivative 
      sgnlo[1] = sett.Smin - (sett.Smin + sett.Smax)*get_rand();
  
      // Frequency in (0, \pi) range  
      sgnlo[0] =  get_rand()*(rrn[1] - rrn[0]) + rrn[0];

      f1 = sgnlo[0] + 2.*sgnlo[1]*sett.N*nfrinband;
      f2 = sgnlo[0] - 2.*sgnlo[1]*sett.N*nfrinband;
 
      // Check if the signal is in band 
      inband = ((f1>0) && (f1<M_PI) && (f2<M_PI) && (f2>0));
      //printf("bad: f1=%f   f2=%f  spnd=%e  N=%d   nfrinband=%d\n", f1, f2, sgnlo[1], sett.N, nfrinband);

  } while (inband != 1);
  
  //printf("f1=%f    f2=%f   spnd=%e  N=%d\n", f1, f2, sgnlo[1], sett.N);

  // Hemisphere 
  pm = round(get_rand()) + 1 ;
  
  // epsma - average inclination of Earth's axis to the ecliptic
  // (see definitions of constants in settings.h) 
  sepsm = sin(C_EPSMA);
  cepsm = cos(C_EPSMA);

  // Random sky position such that the signal is on the grid 
  // (inner loop), but is far enough from the poles (outer loop) 
  //#mv outer loop not used now (O3) 

  // Uniform sphere sampling algorithm 
  double x1, x2, X, Y, Z; 
  //  do {

  do { 
      x1 = 2*get_rand() - 1;
      x2 = 2*get_rand() - 1;
  } while(x1*x1 + x2*x2 >= 1); 
       
  X = 2*x1*sqrt(1 - x1*x1 - x2*x2);
  Y = 2*x2*sqrt(1 - x1*x1 - x2*x2);
  Z = 1 - 2*(x1*x1 + x2*x2);
  
  // Sky position: declination
  sgnlo[2] = M_PI_2 - acos(Z); 
  
  // Right ascension
  sgnlo[3] = atan2(Y, X) + M_PI; 

  // Random phase and polarization of the signal
  ph_o = 2.*M_PI*get_rand();    
  psik = 2.*M_PI*get_rand();    
  hoc =  2.*get_rand() - 1.;   
  hop = (1. + hoc*hoc)/2.; 
  iota = acos(hoc);

  sgnlo[4] =  cos(2.*psik)*hop*cos(ph_o) - sin(2.*psik)*hoc*sin(ph_o) ;
  sgnlo[5] =  sin(2.*psik)*hop*cos(ph_o) + cos(2.*psik)*hoc*sin(ph_o) ;
  sgnlo[6] = -cos(2.*psik)*hop*sin(ph_o) - sin(2.*psik)*hoc*cos(ph_o) ;
  sgnlo[7] = -sin(2.*psik)*hop*sin(ph_o) + cos(2.*psik)*hoc*cos(ph_o) ;

  // Output (GW amplitude or signal-to-noise ratio)  
  if(amp) 
    printf("amp %le\n%d\n%d\n", amp, gsize, reffr);   		 
  else if(snr) 
    printf("snr %le\n%d\n%d\n", snr, gsize, reffr);

  printf("%.16le\n%.16le\n%.16le\n%.16le\n%.16le\n%.16le\n%.16le\n%.16le\n", 
	 sgnlo[0], sgnlo[1], sgnlo[2], sgnlo[3], 
	 sgnlo[4], sgnlo[5], sgnlo[6], sgnlo[7]);

  
  //printf("%.16le %.16le %.16le %.16le\n", sgnlo[0], sgnlo[1], sgnlo[2], sgnlo[3]);

 
  // Testing printouts	
/*
  printf("Random number between 0 and 1: %lf\n", rand1) ; 
  printf("pm, mm, nn, spnd: %d %d %d %d\n", pm, mm, nn, spnd) ;
  printf("%le\n%d\n%d\n", amp, gsize, pm) ;
  printf("rrn[0], rrn[1]: %.16le %.16le\n", rrn[0], rrn[1]) ; 
  
  printf("sgnlo[0]: %.16le\nsgnlo[1]: %.16le\nsgnlo[2]: %.16le\nsgnlo[3]: %.16le\n",
		sgnlo[0], sgnlo[1], sgnlo[2], sgnlo[3]) ;
  printf("sgnlo[4]: %.16le\nsgnlo[5]: %.16le\nsgnlo[6]: %.16le\nsgnlo[7]: %.16le\n",
		sgnlo[4], sgnlo[5], sgnlo[6], sgnlo[7]) ;		
  printf("be1: %.16le\nbe2: %.16le\n", be1, be2) ; 
  
  printf("iota, ph_o, psik: %.8lf %.8lf %.8lf\nhop, hoc: %.8lf %.8lf\n", 
		 iota, ph_o, psik, hop, hoc) ;
 
*/
 		 
  return 0;

} // sigen()


// Random number generator with seed from /dev/urandom
// range: [0,1] 
double get_rand() { 
	
  FILE *urandom;
  unsigned int seed;

  urandom = fopen ("/dev/urandom", "r");
  if (urandom == NULL) {
    fprintf (stderr, "Cannot open /dev/urandom!\n");
    exit(EXIT_FAILURE);
  }

  fread (&seed, sizeof (seed), 1, urandom);
  srand (seed);
  return ((double)rand()/(double)(RAND_MAX));

}	
