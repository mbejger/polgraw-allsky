
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>

//#include <time.h>

#include "auxi.h"
#include "settings.h"
#include "struct.h"

static int help_flag=0;
static int evade_flag=0; 

// Default output and data directories

double get_rand() ;

int main(int argc, char *argv[]) {
	
  int i, numl=0, freq_line_check, c, pm, gsize=2, band=0, reffr; 
  char filename[512], dtaprefix[512], *wd=NULL ; 
  double amp=0, snr=0;
  double freql[32768], linew[32768], sgnlo[8], rrn[2], 
	dvr, fpo_val, be1, be2, fr, lw, freqlp, freqlm, f1, f2,  
	sepsm, cepsm, sinalt, cosalt, sindelt, cosdelt, 
	iota, ph_o, psik, hop, hoc ; 

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
  
  while (1) {
    static struct option long_options[] = {
      {"help", no_argument, &help_flag, 1},			
      // GW amplitude 
      {"amp", required_argument, 0, 'a'},
      // SNR
      {"snr", required_argument, 0, 'n'},
      // frequency band number 
      {"band", required_argument, 0, 'b'},
      // change directory parameter
      {"cwd", required_argument, 0, 'c'},
      // input data directory 
      {"data", required_argument, 0, 'd'},
      // fpo value
      {"fpo", required_argument, 0, 'p'},
      // evade lines 
      {"evade-lines", no_argument, &evade_flag, 1},
      // grid search range 
      {"gsize", required_argument, 0, 'g'},
      // data sampling time 
      {"dt", required_argument, 0, 's'},
      // reference frame () 
      {"reffr", required_argument, 0, 'r'},
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
     printf("-cwd      Change to directory <dir>\n");     
     printf("-data     Data directory in case of lines (default is .)\n"); 
     printf("-fpo      fpo (starting frequency) value\n");
     printf("-gsize    Grid search range (default value: 2)\n");
     printf("-dt       Data sampling time dt (default value: 0.5)\n");
     printf("-nod      Number of days\n");
     printf("-reffr    Reference frame (default value: 1)\n\n");

     printf("--help    This help\n"); 		
     exit (0);
     
  }

    int option_index = 0;
    c = getopt_long_only (argc, argv, "a:n:b:c:d:p:g:s:r:y:", long_options, &option_index);

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
    case 'c':
      wd = (char *) malloc (1+strlen(optarg));
      strcpy (wd, optarg);
      break;      
    case 'd':
      strcpy (dtaprefix, optarg);
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
  // The usual definition: 
      sett.fpo = 10. + 0.96875*band*(0.5/sett.dt);

  // Search settings
  search_settings(&sett);

  fprintf(stderr, "Band number is %04d\n", band);
  fprintf(stderr, "The reference frequency fpo is %f\n", sett.fpo);
  fprintf(stderr, "The data sampling time dt is %f\n", sett.dt);
  fprintf(stderr, "The reference time frame is %d\n", reffr);

  // If the -d switch is on giving a non-zero-length path 
  // somewhere, the signal will be generated using the 
  // information about lines there - will avoid known lines, 
  // as well as pole regions 

  //evade_flag = strlen(dtaprefix) ; 

  // Reading of lines data (lines related to a current band) 
  if(evade_flag) {

        fprintf(stderr, "Will evade lines\n") ;

  // Veto files
  // ----------
  // Known lines 

  sprintf (filename, "%s/veto/%03d.kl", dtaprefix, band);
  if ((data=fopen (filename, "r")) != NULL) {

	numl=0 ;   
	while(!feof(data)) {  

		fscanf(data, "%le %le", &fr, &lw) ;

		if(fr >= sett.fpo && fr <= sett.fpo + 1) { 
				
			freql[numl] = fr ;
			// Empirical 1.2 multiplication factor for the width
			linew[numl] = 1.2*lw ; 
			numl++ ; 

		}       	
 	} 

	// Info on stderr 
	fprintf (stderr, "Reading known-lines file %s for vetoing (%d found)\n", filename, numl-1) ;

  } else {
	perror (filename);
	return 1;
  }      

  fclose(data) ;

  // Forbidden declination values 
  // (signal cannot be to close to the poles)

  sprintf (filename, "%s/veto/pole_veto.d", dtaprefix);
  if ((data=fopen (filename, "r")) != NULL) {

	// Declination veto range for a given band 
	// (file pole_veto.d contains these values 
	// sorted by band)
	i=0; 
	while(i<band) {
		fscanf(data, "%le", &dvr) ;
		i++; 
	}	

	// Info on stderr 
	fprintf (stderr, "Near-pole angle veto file %s (%le rad)\n", filename, dvr) ;

	// Allowed range of angles (-dvr on the 2nd hemisphere)
	dvr = M_PI/2. - 3.*dvr ; 

  } else {
	perror (filename);
	return 1;
  }	

  fclose(data) ;  

  } // end of reading line data in case of -d 

  //#mb Changing the limits near the bands border 
  //rrn[0] = (double)(sett.nmin)/sett.nfft ; 
  //rrn[1] = (double)(sett.nmax)/sett.nfft ;  

  //#mb For O2 UL simulation run 
  rrn[0] = M_PI/20.; 
  rrn[1] = M_PI - rrn[0]; 

  // Random signal parameters 
  //------------------------- 

  // Frequency derivative 
  sgnlo[1] = sett.Smin - (sett.Smin + sett.Smax)*get_rand();
  
  if(evade_flag) { 
  // Frequency must not fall into one of know lines 
  do { 

	freq_line_check = 1 ;   

	sgnlo[0] =  get_rand()*(rrn[1] - rrn[0]) + rrn[0] ;   

	// frequency shift range in VSR1 data
	// these values are used to avoid lines
	f1 = sett.fpo + sgnlo[0] ;

  //#mb !!! reference frame !!! Here reffr (VSR1 was 67) 
	f2 = sett.fpo + sgnlo[0] - 2.*sgnlo[1]*(sett.N)*reffr ;

	for(i=0; i<numl; i++) {

		freqlp = freql[i] + linew[i] ; 
		freqlm = freql[i] - linew[i] ;

		// checks if the line is between f1 and f2
		// if so, starts again
		if(!((freqlm > f2) || (freqlp < f1))) {
			freq_line_check = 0 ;
			break; 
		} 

	}	


  } while (!freq_line_check) ; 

  // Frequency in (0, \pi) range  
  //#mb Check if this is consistent with the definition 
  // of rrn[0] and rrn[1]  
  sgnlo[0] *= 2*M_PI ;

  // Lines will not be avoided

  } else { 

	  //sgnlo[0] =  2*M_PI*get_rand()*(rrn[1] - rrn[0]) + rrn[0];	
	  //
	  //#mb this is compatible with the definition of 
	  // rrn[0] = M_PI/20., rrn[1] = M_PI - rrn[0]

	  sgnlo[0] =  get_rand()*(rrn[1] - rrn[0]) + rrn[0];	

  } 

  // Hemisphere 
  pm = round(get_rand()) + 1 ;
 
  // epsma - average inclination of Earth's axis to the ecliptic
  // (see definitions of constants in settings.h) 
  sepsm = sin(C_EPSMA);
  cepsm = cos(C_EPSMA);

  // Random sky position such that the signal is on the grid 
  // (inner loop), but is far enough from the poles (outer loop) 

  // Uniform sphere sampling algorithm 
  double x1, x2, X, Y, Z; 
  do {

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

	// breaks out of the outer loop
	// if there is no -d switch with info about known lines   
	if(!evade_flag) break; 

  } while((sgnlo[2] > dvr) || (sgnlo[2] < -dvr)) ; 

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
 		 
//  downtime = clock() / (CLOCKS_PER_SEC / 1000 ) ; 
//  time_in_seconds = (double)(downtime - uptime) / 1000 ;     
//  printf("\nTime [s]: %.3lf\n", time_in_seconds) ; 
  
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
	exit (1);
	}

	fread (&seed, sizeof (seed), 1, urandom);

	srand (seed);
	//# value for testing 
  // return 0.24; 
	return ((double)rand()/(double)(RAND_MAX));


}	
