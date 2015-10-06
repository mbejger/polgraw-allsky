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

void search_settings(
  Search_settings* sett) {

  double dt, B, oms, omr, Smin, Smax;
  int nod, N, nfft, s, nd, interpftpad;


  dt = sett->dt;                    // data sampling time:  
                                    // set in handle_opts() from the command line
                                    // (the default value is dt=0.5)

  B = 0.5/dt;                       // Bandwidth
  oms = 2.*M_PI*(sett->fpo)*dt;     // Dimensionless angular frequency

  omr = C_OMEGA_R*dt;

  nod = 2;                          // Observation time in days
  N = round (nod*C_SIDDAY/dt);      // No. of data points

  nfft = 1 << (int)ceil(log(N)/log(2.));    // length of FFT
  s = 1;                                    // No. of spindowns

  Smin = 1000.*C_YEARSEC;                   // Minimum spindown time 
                                            // [sec.]

  // Maximum spindown (1000 years) [angular, dimensionless]
  Smax = 2.*M_PI*(sett->fpo + B)*dt*dt/(2.*Smin);   

  nd = 2;     // Degree of freedom, 
              // (2*nd = deg. no ofrees of freedom for chi^2)

  interpftpad = 2;

  sett->B=B;          	// bandwidth
  sett->oms=oms;      	// dimensionless angular frequency
  sett->omr=omr;      	// C_OMEGA_R * dt
  sett->nod=nod;      	// number of days of observation
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

  sett->nmin = sett->fftpad*NAV*sett->B;
  sett->nmax = (sett->nfft/2 - NAV*sett->B)*sett->fftpad;

} // search settings  



/* Network of detectors' discovery: 
 * finds subdirectories in the main input directory, 
 * which by convention should be named like V1, L1, H1 
 * and which contain input data and ephemerids; 
 * writes appropriate detector-related data into structs. 
 */ 

void detectors_settings(
  Search_settings* sett, 
  Command_line_opts *opts) {

  int i=0; 

  char dirname[512], x[512];
  // Main input directory name 
  sprintf (dirname, "%s/%03d", opts->dtaprefix, opts->ident); 

  DIR *dp;
  struct dirent *ep;

  char **detnames  = malloc(MAX_DETECTORS*sizeof(char*));   
  char **xnames = malloc(MAX_DETECTORS*sizeof(char*));

  dp = opendir (dirname);
  if (dp != NULL) {
    while ((ep = readdir (dp))) { 

      // Subdirectory names: 2 char long
      if((ep->d_type == DT_DIR) && 
        (strlen(ep->d_name)==DETNAME_LENGTH) && 
        strncmp(&ep->d_name[0],".",1)) { 

          FILE *data;

          // Input time-domain data handling
          // 
          // We assume that in each subdirectory corresponding 
          // to the detector the input data will look as following: 
          sprintf(x, "%s/%03d/%s/xdatc_%03d_%03d%s.bin",
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

          memset(x, 0, sizeof(x));
      }
    } 
      
    (void) closedir(dp);

  } else perror ("Couldn't open the input directory...");

  sett->nifo=i;      // number of detectors  
  if(sett->nifo) { 
    printf("Settings - number of detectors: %d\n", sett->nifo); 

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
      ifo[i].eheight = 53.238;
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

void modvir(
  double sinal, 
  double cosal, 
  double sindel, 
  double cosdel,
  int Np,
  Detector_settings *ifo, 
  Aux_arrays *aux) {

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



    /* Obtaining the trigger files list from a given directory
     * and reading their content 
     */ 

void read_trigger_files(
  Search_settings *sett, 
  Command_line_opts_coinc *opts, 
  Candidate_triggers *trig) {

  int i=0, trig_size=1048576, current_frame=0, numofframes=0;  

  char dirname[512], filename[512]; 
  // Trigger files directory name 
  sprintf (dirname, "%s", opts->dtaprefix); 

  DIR *dp;
  struct dirent *ep;
  FILE *data; 

  FLOAT_TYPE c[5]; 

  FLOAT_TYPE *t; int *ti; // testing for realloc  

  trig->f   = (FLOAT_TYPE *)calloc(trig_size, sizeof(FLOAT_TYPE));
  trig->s   = (FLOAT_TYPE *)calloc(trig_size, sizeof(FLOAT_TYPE));
  trig->d   = (FLOAT_TYPE *)calloc(trig_size, sizeof(FLOAT_TYPE));
  trig->a   = (FLOAT_TYPE *)calloc(trig_size, sizeof(FLOAT_TYPE));
  trig->snr = (FLOAT_TYPE *)calloc(trig_size, sizeof(FLOAT_TYPE));
  trig->fr  = (int *)calloc(trig_size, sizeof(int));

  dp = opendir (dirname);
  if (dp != NULL) {
    while ((ep = readdir (dp))) { 

      if((ep->d_type == DT_REG) &&
        (strstr(ep->d_name, opts->trigname) != NULL)) { 

          sprintf(filename, "%s/%s", opts->dtaprefix, ep->d_name); 

          // This part looks for the first number in the trigger file name, 
          // under the assumption that this is the frame number
          char *tmp, *epdname;
          epdname = strdup(ep->d_name);  
          while((tmp = strsep(&epdname, "_"))!=NULL) {
            if(tmp[0] >= '0' && tmp[0] <= '9') {
                current_frame = atoi(tmp);
                printf("Reading %s (frame number %d)...\n", ep->d_name, current_frame); 
                break; 
            }  
          } 

          numofframes++; 

          if((data = fopen(filename, "r")) != NULL) {

            FLOAT_TYPE finband; 

            while(fread((void *)c, sizeof(FLOAT_TYPE), 5, data)==5) {  

              //#mb
              finband = c[0] + 2.*c[1]*sett->N*(opts->refr - current_frame); 

              if((finband>0) && (finband<M_PI)) { 

                trig->f[i]   = finband; trig->s[i]  = c[1]; 
                trig->d[i]   = c[2];    trig->a[i]  = c[3]; 
                trig->snr[i] = c[4];    trig->fr[i] = current_frame;  

              i++; 

              }   

              if(i==trig_size) {

                // Doubling the triggers' array size 
                trig_size *= 2;
 
                t = (FLOAT_TYPE*)realloc(trig->f, trig_size*sizeof(FLOAT_TYPE));
                if(t!=NULL) trig->f = t;
                t = (FLOAT_TYPE*)realloc(trig->s, trig_size*sizeof(FLOAT_TYPE));
                if(t!=NULL) trig->s = t;
                t = (FLOAT_TYPE*)realloc(trig->d, trig_size*sizeof(FLOAT_TYPE));
                if(t!=NULL) trig->d = t;
                t = (FLOAT_TYPE*)realloc(trig->a, trig_size*sizeof(FLOAT_TYPE));
                if(t!=NULL) trig->a = t;
                t = (FLOAT_TYPE*)realloc(trig->snr, trig_size*sizeof(FLOAT_TYPE));
                if(t!=NULL) trig->snr = t;
                ti = (int*)realloc(trig->fr, trig_size*sizeof(int));
                if(t!=NULL) trig->fr = ti;

              }

            }

          } else { 
            printf("Problem with %s...\n", filename);  
            perror (filename);
          }

          memset(filename, 0, sizeof(filename));
          fclose(data); 

      } 
   
    } 

  } 

  (void) closedir(dp);

  trig->num_of_trig = i; 
  trig->num_of_frames = numofframes; 

}
