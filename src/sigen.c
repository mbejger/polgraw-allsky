#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>

//#include <time.h>

#include "auxi.h"
#include "lvcvirgo.h"

static int help_flag=0; 

/* Default output and data directories */
#ifndef PREFIX
#define PREFIX .
#endif

#ifndef DTAPREFIX
#define DTAPREFIX .
#endif

double get_rand() ;

int main(int argc, char *argv[]) {
	
  int i, numl=0, freq_line_check, c, pm, gsize=2, band=0 ; 
  char filename[64], prefix[64], dtaprefix[64], *wd=NULL ; 
  double freql[32768], linew[32768], sgnlo[8], rrn[2], h0=2.e-3, 
	dvr, fpo, be1, be2, fr, lw, freqlp, freqlm, f1, f2,  
	sepsm, cepsm, sinalt, cosalt, sindelt, cosdelt, 
	nmin, nmax, iota, ph_o, psik, hop, hoc ; 

  FILE *data ;
  		  
//  clock_t uptime, downtime ;
//  double time_in_seconds ; 	
//  uptime = clock() / (CLOCKS_PER_SEC / 1000 ) ; 
  
  strcpy (prefix, TOSTR(PREFIX));
  strcpy (dtaprefix, TOSTR(DTAPREFIX));

  while (1) {
    static struct option long_options[] = {
      {"help", no_argument, &help_flag, 1},			
      // GW amplitude 
      {"h0", required_argument, 0, 'a'},
      // frequency band number 
      {"band", required_argument, 0, 'b'},
      // change directory parameter
      {"cwd", required_argument, 0, 'c'},
      // input data directory 
      {"data", required_argument, 0, 'd'},
      // grid search range 
      {"gsize", required_argument, 0, 's'},
      {0, 0, 0, 0}
    };

  if(help_flag) {
	  
     printf("*** Software injection parameters - cont. GW signal ***\n"); 
     printf("Usage: ./sigen -[switch1] <value1> -[switch2] <value2> ...\n") ;
     printf("Switches are:\n\n"); 
     printf("-a		GW amplitude (default value: 2.e-3)\n"); 
     printf("-b		Band number\n"); 
     printf("-c		Change to directory <dir>\n");     
     printf("-d		Data directory (default is .)\n"); 
     printf("-s		Grid search range (default value: 2)\n");
     printf("--help		This help\n"); 		
     exit (0);
     
  }

    int option_index = 0;
    c = getopt_long (argc, argv, "a:b:c:d:s:", long_options,	\
		     &option_index);
if (c == -1)
      break;

    switch (c) {
    case 'a':
      h0  = atof (optarg);
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
    case 's':
      gsize = atoi (optarg);
      break;
    case '?':
      break;
    default: break ; 
    } /* switch c */
  } /* while 1 */

  if (wd) {
    printf ("Changing working directory to %s\n", wd);
    if (chdir(wd)) {
      perror (wd);
      abort ();
    }
  }
  	
  // Starting band frequency
  fpo = 100. + 0.96875 * band;
  // Detector, ephemerides etc. constants 
  lvcvirgo (fpo);
  
  fprintf(stderr, "Starting band frequency fpo: %lf\n", fpo) ; 

  // Veto files
  // ----------
  // Known lines 

  sprintf (filename, "%s/veto/%03d.kl", dtaprefix, band);
  if ((data=fopen (filename, "r")) != NULL) {

	numl=0 ;   
	while(!feof(data)) {  

		fscanf(data, "%le %le", &fr, &lw) ;

		if(fr >= fpo && fr <= fpo + 1) { 
				
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

  // Because of frequency-domain filters, we search 
  // F-statistic in range (nmin+1, nmax) of data points 
  nmin = 2.*NAV ;
  nmax = nfft - 2.*NAV;
  rrn[0] = nmin/nfft ; 
  rrn[1] = nmax/nfft ;  



  // Random signal parameters 
  //------------------------- 

  // Frequency derivative 
  sgnlo[1] = -Smax*get_rand() ;

  // Frequency must not fall into one of know lines 
  do { 

	freq_line_check = 1 ;   

	sgnlo[0] =  get_rand()*(rrn[1] - rrn[0]) + rrn[0] ;   

	// frequency shift range in VSR1 data
	// these values are used to avoid lines
	f1 = fpo + sgnlo[0] ;
	f2 = fpo + sgnlo[0] - 2.*sgnlo[1]*Nv*67 ;

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

  // Frequency 
  sgnlo[0] *= M_PI ;

  // Hemisphere 
  pm = round(get_rand()) + 1 ;
 
  // epsma - average inclination of Earth's axis to the ecliptic
  sepsm = sin(epsma);
  cepsm = cos(epsma);

  // Random sky position such that the signal is on the grid 
  // (inner loop), but is far enough from the poles (outer loop) 
  do { 
  	do { 
		be1 = 2.*get_rand() - 1 ; 
		be2 = 2.*get_rand() - 1 ;

  	} while(sqr(be1)+sqr(be2) > 1) ; 

  	lin2ast (be1, be2, pm, sepsm, cepsm, &sinalt, &cosalt, &sindelt, &cosdelt);
  
  	// Sky position: declination
  	sgnlo[2] = asin (sindelt);

	// Right ascension
  	sgnlo[3] = fmod (atan2 (sinalt, cosalt) + 2.*M_PI, 2.*M_PI);
  
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
    
  // Output
  printf("%le\n%d\n%d\n", h0, gsize, pm) ;   		 
  printf("%.16le\n%.16le\n%.16le\n%.16le\n%.16le\n%.16le\n%.16le\n%.16le\n", 
			sgnlo[0], sgnlo[1], sgnlo[2], sgnlo[3], 
			sgnlo[4], sgnlo[5], sgnlo[6], sgnlo[7]) ;
  printf("%.16le\n%.16le\n", be1, be2) ; 			 
   
  // Testing printouts	
/*
  printf("Random number between 0 and 1: %lf\n", rand1) ; 
  printf("pm, mm, nn, spnd: %d %d %d %d\n", pm, mm, nn, spnd) ;

  printf("%le\n%d\n%d\n", h0, gsize, pm) ;
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

} /* sigen() */


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
	//return 0.476559 ; 
	return ((double)rand()/(double)(RAND_MAX));


}	
