#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>

#include "settings.h"
#include "auxi.h"
#include "struct.h"


/* Search settings: 
 * FFT lenghts & other details, bandwidth and Earth parameters
 */

void search_settings(Search_settings* sett) {

  double dt, B, oms, omr, Smin, Smax;
  int N, nfft, s, nd, interpftpad;


  dt = sett->dt;                    // data sampling time:  
                                    // set in handle_opts() from the command line
                                    // (the default value is dt=0.5)

  B = 0.5/dt;                       // Bandwidth
  oms = 2.*M_PI*(sett->fpo)*dt;     // Dimensionless angular frequency

  omr = C_OMEGA_R*dt;

  N = round (sett->nod*C_SIDDAY/dt);      // No. of data points

  nfft = 1 << (int)ceil(log(N)/log(2.));    // length of FFT
  s = 1;                                    // No. of spindowns

  /*
  Smin = 1000.*C_YEARSEC;                   // Minimum spindown time 
                                            // [sec.]

  // Maximum spindown (1000 years) [angular, dimensionless]
  Smax = 2.*M_PI*(sett->fpo + B)*dt*dt/(2.*Smin);   
  */

  //#mb ranges of spindown (RDC O1)
  double fdotmin, fdotmax; 
  fdotmin = 0.5e-8; 
  fdotmax = 0.5e-9; 

  Smax = 2.*M_PI*fdotmin*dt*dt; 
  Smin = 2.*M_PI*fdotmax*dt*dt;

  nd = 2;     // Degree of freedom, 
              // (2*nd = deg. no ofrees of freedom for chi^2)

  interpftpad = 2;

  sett->B=B;          	// bandwidth
  sett->oms=oms;      	// dimensionless angular frequency
  sett->omr=omr;      	// C_OMEGA_R * dt
  sett->N=N;          	// number of data points
  sett->nfft=nfft;    	// length of fft
  sett->s=s;          	// number of spindowns
  sett->Smin=Smin;    	// minimum spindown
  sett->Smax=Smax;    	// maximum spindown
  sett->nd=nd;        	// degrees of freedom
  sett->interpftpad=interpftpad;

  // Because of frequency-domain filters, we search
  // F-statistic in range (nmin+1, nmax) of data points
  // 
  // The value of sett->fftpad (zero padding - original grids: 2, new grids: 1) 
  // is read from the grid.bin file in read_grid() (see init.c)
  
#if 1 
  sett->nmin = sett->fftpad*NAV*sett->B;
  sett->nmax = (sett->nfft/2 - NAV*sett->B)*sett->fftpad;
#else
  sett->nmin = pow(2,-6)*sett->nfft/2*sett->fftpad;
  sett->nmax = (1 - pow(2,-6))*sett->nfft/2*sett->fftpad;
#endif
  printf("------------------------ Settings -------------------------\n");
  printf(" B         N          nfft         Fstat_nmin    Fstat_nmax\n");
  printf("%6.3f  %9d  %9d  %9d       %9d\n", sett->B, sett->N, sett->nfft, sett->nmin, sett->nmax);
  printf("-----------------------------------------------------------\n");
    
  // initial value of number of known instrumental lines in band 
  sett->numlines_band=0; 

} // search settings  



/* Network of detectors' discovery: 
 * finds subdirectories in the main input directory, 
 * which by convention should be named like V1, L1, H1 
 * and which contain input data and ephemerids; 
 * writes appropriate detector-related data into structs. 
 */ 

void detectors_settings(Search_settings* sett, 
			Command_line_opts *opts) {

  int i=0; 

  char dirname[516], x[1332];
  // Main input directory name 
  sprintf (dirname, "%s/%03d", opts->dtaprefix, opts->ident); 

  DIR *dp;
  struct dirent *ep;

  char **detnames  = malloc(MAX_DETECTORS*sizeof(char*));   
  char **xnames = malloc(MAX_DETECTORS*sizeof(char*));

  dp = opendir (dirname);
  if (dp != NULL) {
    while ((ep = readdir (dp))) { 

      // Subdirectory names checkup:
      // check if it's a dir
      // name is 2 char long
      // not a directory name of the type "./" or ".."
      // if opts->usedet is not set (length equal 0), 
      // or is set and dir name is substring of it 
      if((ep->d_type == DT_DIR) && 
        (strlen(ep->d_name)==DETNAME_LENGTH) && 
        (strncmp(&ep->d_name[0],".",1)) && 
        (!strlen(opts->usedet) || (strlen(opts->usedet) && (strstr(opts->usedet, ep->d_name))))) { 

          FILE *data;

          // Input time-domain data handling
          // 
          // We assume that in each subdirectory corresponding 
          // to the detector the input data will look as following: 
/*          sprintf(x, "%s/%03d/%s/xdatc_%03d%s.bin",
          opts->dtaprefix, opts->ident, ep->d_name,
          opts->ident, opts->label);
*/

          sprintf(x, "%s/%03d/%s/xdatsc_%03d_%04d%s.bin",
          opts->dtaprefix, opts->ident, ep->d_name,
          opts->ident, opts->band, opts->label);

          if((data = fopen(x, "r")) != NULL) {

            xnames[i]   = calloc(strlen(x)+1, sizeof(char));
            detnames[i] = calloc(DETNAME_LENGTH+1, sizeof(char));

            strncpy(xnames[i], x, strlen(x));
            strncpy(detnames[i], ep->d_name, DETNAME_LENGTH);
            i++;

          } else { 
            printf("Directory %s exists, but no input file found:\n%s missing...\n", 
              ep->d_name, x);  
            //perror (x);
          }

          fclose(data); 
          memset(x, 0, sizeof(x));
      }
    } 
      
    (void) closedir(dp);

  } else perror ("Couldn't open the input directory...");

  sett->nifo=i;      // number of detectors  
  if(sett->nifo) { 
    printf("Number of detectors: %d\n", sett->nifo); 

  } else { 
    printf("No subdirectories with detector data found. Exiting...\n"); 
    exit(EXIT_FAILURE);
  }  

  for(i=0; i<sett->nifo; i++) { 

    // Virgo detector
    if(!strcmp("V1", detnames[i])) {

      strncpy(ifo[i].xdatname, xnames[i], strlen(xnames[i]));
      strncpy(ifo[i].name, detnames[i], DETNAME_LENGTH);

      // Geographical latitude phi in radians
      ifo[i].ephi = (43.+37./60.+53.0880/3600.)/RAD_TO_DEG;
      // Geographical longitude in radians
      ifo[i].elam = (10.+30./60.+16.1885/3600.)/RAD_TO_DEG;
      // Height h above the Earth ellipsoid in meters
      ifo[i].eheight = 51.884;
      // Orientation of the detector gamma
      ifo[i].egam = (135. - (19.0+25./60.0+57.96/3600.))/RAD_TO_DEG;

      printf("Using %s IFO as detector #%d... %s as input time series data\n", 
        ifo[i].name, i, ifo[i].xdatname);

    // Hanford H1 detector
    } else if(!strcmp("H1", detnames[i])) {

      strncpy(ifo[i].xdatname, xnames[i], strlen(xnames[i]));
      strncpy(ifo[i].name, detnames[i], DETNAME_LENGTH);

      // Geographical latitude phi in radians
      ifo[i].ephi = (46+(27+18.528/60.)/60.)/RAD_TO_DEG;
      // Geographical longitude in radians
      ifo[i].elam = -(119+(24+27.5657/60.)/60.)/RAD_TO_DEG;
      // Height h above the Earth ellipsoid in meters
      ifo[i].eheight = 142.554;
      // Orientation of the detector gamma
      ifo[i].egam = 170.9994/RAD_TO_DEG;

      printf("Using %s IFO as detector #%d... %s as input time series data\n", 
        ifo[i].name, i, ifo[i].xdatname);

    // Livingston L1 detector
    } else if(!strcmp("L1", detnames[i])) {

      strncpy(ifo[i].xdatname, xnames[i], strlen(xnames[i]));
      strncpy(ifo[i].name, detnames[i], DETNAME_LENGTH);

      // Geographical latitude phi in radians
      ifo[i].ephi = (30+(33+46.4196/60.)/60.)/RAD_TO_DEG;
      // Geographical longitude in radians
      ifo[i].elam = -(90+(46+27.2654/60.)/60.)/RAD_TO_DEG;
      // Height h above the Earth ellipsoid in meters
      ifo[i].eheight = -6.574;
      // Orientation of the detector gamma
      ifo[i].egam = 242.7165/RAD_TO_DEG;

      printf("Using %s IFO as detector #%d... %s as input time series data\n", 
        ifo[i].name, i, ifo[i].xdatname);

    } else {

      printf("Meh, unknown detector %s (see settings.c) Exiting...\n", 
              detnames[i]);
      exit(EXIT_FAILURE);
    }

  } 

  // memory free for detnames and xdatnames
  for(i=0; i<sett->nifo; i++) { 
    free(detnames[i]);
    free(xnames[i]); 
  } 
   

  free(detnames); 
  free(xnames); 

} // detectors settings


  /* Coefficients of the amplitude modulation functions
   * of the Virgo detector
   */ 

void rogcvir(Detector_settings *ifo) {

  /* In the notation of Phys. Rev. D 58, 063001 (1998):
   * ephi = lambda (geographical latitude phi in radians)
   * egam = gamma (orientation of the detector)
   * 
   * (see modvir function in jobcore.c for Eqs. 12 and 13)
   */ 

  //printf("Calculating the amplitude modulation functions for %s...\n", ifo->name); 

  ifo->amod.c1 = .25*sin(2.*ifo->egam)*(1+sqr(sin(ifo->ephi)));
  ifo->amod.c2 = -.5*cos(2.*ifo->egam)*sin(ifo->ephi);
  ifo->amod.c3 = .5*sin(2.*ifo->egam)*sin(2.*ifo->ephi);
  ifo->amod.c4 = -cos(2.*ifo->egam)*cos(ifo->ephi);
  ifo->amod.c5 = .75*sin(2.*ifo->egam)*sqr(cos(ifo->ephi));
  ifo->amod.c6 = cos(2.*ifo->egam)*sin(ifo->ephi);
  ifo->amod.c7 = .5*sin(2.*ifo->egam)*(1.+sqr(sin(ifo->ephi)));
  ifo->amod.c8 = cos(2.*ifo->egam)*cos(ifo->ephi);
  ifo->amod.c9 = .5*sin(2.*ifo->egam)*sin(2.*ifo->ephi);


} // rogcvir


  /* Amplitude modulation of the signal
   */ 

void modvir(double sinal, double cosal, double sindel, double cosdel,
	    int Np, Detector_settings *ifo, Aux_arrays *aux) {

  int t;
  double cosalfr, sinalfr, c2d, c2sd, c, s, c2s, cs;

  double c1 = ifo->amod.c1,
         c2 = ifo->amod.c2,
         c3 = ifo->amod.c3,
         c4 = ifo->amod.c4,
         c5 = ifo->amod.c5,
         c6 = ifo->amod.c6,
         c7 = ifo->amod.c7,
         c8 = ifo->amod.c8,
         c9 = ifo->amod.c9;

  cosalfr = cosal*(ifo->sig.cphir) + sinal*(ifo->sig.sphir);
  sinalfr = sinal*(ifo->sig.cphir) - cosal*(ifo->sig.sphir);
  c2d = sqr(cosdel);
  c2sd = sindel*cosdel;

  // Modulation factor for every data point 
  for (t=0; t<Np; t++) { 

    c = cosalfr*aux->cosmodf[t] + sinalfr*aux->sinmodf[t];
    s = sinalfr*aux->cosmodf[t] - cosalfr*aux->sinmodf[t];
    c2s = 2.*sqr(c);
    cs = c*s;

    // modulation factors aa and bb  
    ifo->sig.aa[t] = c1*(2.-c2d)*c2s + c2*(2.-c2d)*2.*cs +
           c3*c2sd*c + c4*c2sd*s - c1*(2.-c2d) + c5*c2d;

    ifo->sig.bb[t] = c6*sindel*c2s + c7*sindel*2.*cs + 
           c8*cosdel*c + c9*cosdel*s - c6*sindel;

  } 

} // modvir

int read_lines(Search_settings *sett, 
	       Command_line_opts *opts, 
	       Detector_settings *ifo) { 
  
  int i=0, lnum, alllines, j; 
  double l[MAXL][7]; 
  char linefile[542], line[128] = {0}; 
  FILE *data; 
  
  sprintf(linefile, "%s/O2LinesToBeCleaned_%s_v1.txt", 
	  opts->dtaprefix, ifo->name);
  
  // Reading line data from the input file (data)   
  // Columns are: 
  // 1 - frequency spacing (Hz) of comb (or frequency of single line)
  // 2 - comb type (0 - singlet, 1 - comb with fixed width, 2 - comb with scaling width)
  // 3 - frequency offset of 1st visible harmonic (Hz)
  // 4 - index of first visible harmonic
  // 5 - index of last visible harmonic
  // 6 - width of left band (Hz)
  // 7 - width of right band (Hz)
  
  if ((data = fopen(linefile, "r")) != NULL) {
    
    while (fgets(line, 128, data) != NULL) {
      
      // Skip comment lines beginning with '%'
      if (*line == '%') continue;
      
      j=0; 
      char *t;  
      for(t = strtok(line, "\t "); t != NULL; t = strtok(NULL, "\t ")) {
	l[i][j] = atof(t);
	j++; 
      }

      if (++i > MAXL-1) {
	printf("Can't read %s, increase MAXL!\n", linefile); 
	exit(EXIT_FAILURE);
      }
    }
    
  } else {  
    printf("No known lines file %s for %s IFO found. Continuing without.\n", 
	   linefile, ifo->name);
    return 0; 
  } 

  // The number of data for lines
  lnum = i; 
  
  fclose(data); 
  
  // Produce line widths 
  //--------------------
  
  // Fixed line widths for each detector 
  // based on the broadening modulation 
  // due to the detector movement 
  
  double *dE; 
  dE = (double *)calloc(3*sett->N, sizeof(double));
  
  // First derivative DetSSB (velocity)  
  for(i=0; i<sett->N-1; i++) { 
    
    for(j=0; j<3; j++)  
      dE[i*3+j] = fabs(ifo->sig.DetSSB[(i+1)*3+j] - ifo->sig.DetSSB[i*3+j]); 
    
  } 

  double dEmax[3] = {0}; 

  // Find maximum absolute values 
  for(i=0; i<sett->N-1; i++) {

    for(j=0; j<3; j++) 
      if(dE[i*3+j] > dEmax[j]) 
	dEmax[j] = dE[i*3+j]; 
    
  } 
  
  //#mb test printout 
  //    printf("dEmax %f %f %f\n", dEmax[0], dEmax[1], dEmax[2]); 
 
  double *dtE;
  dtE = (double *)calloc(3*sett->N, sizeof(double));
  
  // First derivative 
  for(i=0; i<sett->N-1; i++) { 
    
    for(j=0; j<3; j++)  
      dtE[i*3+j] = fabs(ifo->sig.DetSSB[(i+1)*3+j]*(i+1) - ifo->sig.DetSSB[i*3+j]*i)*sett->dt; 
    
  } 
  
  double dtEmax[3] = {0}; 
    
  // Find maximum absolute values 
  for(i=0; i<sett->N-1; i++) {
    
    for(j=0; j<3; j++) 
      if(dtE[i*3+j] > dtEmax[j]) 
	dtEmax[j] = dtE[i*3+j]; 
    
  } 

  //#mb !!! Limits put by hand for RDC_O1 !!! 
  double fdotMax = 0.5e-8, Dfmaxmax, normdtEmax=0., normdEmax=0.;

  for(j=0; j<3; j++) { 
    
    normdtEmax += pow(dtEmax[j], 2.); 
    normdEmax  += pow(dEmax[j], 2.);
    
  } 
  
  // Free auxiliary allocs 
  free(dE); 
  free(dtE); 
  
  // Apply line widths 
  //------------------
  
  j=0; 
  for(i=0; i<lnum; i++) { 
    
    int line_type = (int)(l[i][1]);
    int indf = (int)(l[i][3]);
    int indl = (int)(l[i][4]);
    int k; 
    
    switch(line_type) { 
      
    // Singlet
    case 0: 

      // Line width from the resampling broadening 
      Dfmaxmax = 2.*fdotMax*(sett->N*sett->dt + sqrt(normdtEmax)) + l[i][0]*sqrt(normdEmax);
      // scaling with nod (relative to nod=6)
      //Dfmaxmax *= 6./sett->nod;
      Dfmaxmax *= 0.25;
      
      // The original width of a line is replaced with Dfmaxmax
      ifo->lines[j][0] = l[i][0] - Dfmaxmax; 
      ifo->lines[j][1] = l[i][0] + Dfmaxmax;
      ifo->lines[j][2] = i;
      if (++j > MAXL-1) {printf("Too many lines, increase MAXL!\n"); exit(EXIT_FAILURE);}
      break; 
      

    // Comb with fixed width. Vetoing the band 
    // [offset+index*spacing-leftwidth, offset+index*spacing+rightwidth] 
    case 1:

      for(k=indf; k<=indl; k++) {
	// Line width from the resampling broadening 
	double linefreq = l[i][2] + k*l[i][0]; 
	
	Dfmaxmax = 2.*fdotMax*(sett->N*sett->dt + sqrt(normdtEmax)) 
	         + linefreq*sqrt(normdEmax);
	// scaling with nod (relative to nod=6)
	Dfmaxmax *= 6./sett->nod;
	
	// The original width of a line is replaced with Dfmaxmax
	ifo->lines[j][0] = linefreq - Dfmaxmax;
	ifo->lines[j][1] = linefreq + Dfmaxmax;
	ifo->lines[j][2] = i;
	if (++j > MAXL-1) {printf("Too many lines, increase MAXL!\n"); exit(EXIT_FAILURE);}
      } 
      break; 
      

    // Comb with scaling-width. Vetoing the band 
    // [offset+index*spacing-index*leftwidth, offset+index*spacing+index*rightwidth]       
    case 2:
      
      for(k=indf; k<=indl; k++) { 
	// Line width from the resampling broadening 
	double linefreq = l[i][2] + k*l[i][0];
	
	Dfmaxmax = 2.*fdotMax*(sett->N*sett->dt + sqrt(normdtEmax)) 
	         + linefreq*sqrt(normdEmax);
	// scaling with nod (relative to nod=6)
	Dfmaxmax *= 6./sett->nod;
	
	ifo->lines[j][0] = linefreq - k*Dfmaxmax;
	ifo->lines[j][1] = linefreq + k*Dfmaxmax;
	ifo->lines[j][2] = i;
	if (++j > MAXL-1) {printf("Too many lines, increase MAXL!\n"); exit(EXIT_FAILURE);}
      }
      break; 
    } 
    
  }
  
  ifo->numlines = j; 
  printf("Number of known lines in %s: %d\n", ifo->name, ifo->numlines); 
  
  return 0; 
  
}


void narrow_down_band(Search_settings* sett,
		      Command_line_opts *opts) {
  
  // Adding excluding ranges near the edges 
  // to the known lines list 
  sett->lines[0][0] = 0;
  sett->lines[0][1] = M_PI_2 - opts->narrowdown;
  sett->lines[1][0] = M_PI_2 + opts->narrowdown;
  sett->lines[1][1] = M_PI;
  
  // These ranges added 
  sett->numlines_band = 2;
  printf("Band is narrowed-down to [%f, %f]\n", sett->lines[0][1], sett->lines[1][0]);
  
}


void lines_in_band(Search_settings* sett, 
		   Command_line_opts *opts) {

    int i, j, k=0; 
    double bs, be;        // Band start and end  
    double ll, lr, l, r;  // Line left and right side

    // if narrowdown option 
    if(opts->narrowdown < 0.5*M_PI) k = sett->numlines_band;

    bs = sett->fpo; 
    be = sett->fpo + sett->B; 

    for(i=0; i<sett->nifo; i++) { 

      printf("Looking for known lines for %s between %f and %f...\n", 
        ifo[i].name, bs, be); 

      for(j=0; j<ifo[i].numlines; j++) {

        double l = ifo[i].lines[j][0]; 
        double r = ifo[i].lines[j][1];

        if(!(r < bs || l > be)) { 

           if(l > bs) ll = l;
           else ll = bs;     

           if(r < be) lr = r;
           else lr = be;     

           sett->lines[k][0] = (ll-bs)/(sett->B)*M_PI;  
           sett->lines[k][1] = (lr-bs)/(sett->B)*M_PI;

	   printf("%s: %f %f %f %f [line %d in veto file]\n", ifo[i].name, ll, lr, sett->lines[k][0], sett->lines[k][1], (int)ifo[i].lines[j][2] );
           k++; 

        }      
      }
    } 

    sett->numlines_band = k;

}     


void veto_fraction(Search_settings* sett) { 
  
  int i; 
  double ll=0., gap=0.;
  
  // Sorting the excluded parts of a band (1st then 2nd column) 
  qsort(sett->lines, sett->numlines_band, 2*sizeof(double), compared2c); 
  
  /*
    printf("After sorting...\n"); 
    for(i=0; i<sett->numlines_band; i++) 
    printf("%f %f\n", sett->lines[i][0], sett->lines[i][1]);
  */
  
  for(i=0; i<sett->numlines_band; i++) { 
    
    // Looking for a gap between the excluded regions 
    if(sett->lines[i][0] >= ll) {
      gap += sett->lines[i][0] - ll;
      ll = sett->lines[i][1];
    } else {
      if (ll < sett->lines[i][1]) ll = sett->lines[i][1];
    }
  }
  if ( ll < M_PI) gap += M_PI - ll;
  
  printf("Band veto fraction = %6.4f\n", (M_PI-gap)/M_PI);

  if(gap <= 0.) { 
    printf("This band is fully vetoed. My work here is done, exiting...\n"); 
    exit(EXIT_SUCCESS); 
  }
 
}
