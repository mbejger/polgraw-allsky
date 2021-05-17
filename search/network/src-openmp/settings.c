#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <float.h>

#include "settings.h"
#include "auxi.h"
#include "struct.h"
#include <glob.h>


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

  // spindown range of NS
  // we assume minimum NS age 1000 yr
  double fdotmin, fdotmax;
  if (sett->fpo < 200.) {
      fdotmin = (sett->fpo+B)/(2.*1000.*C_YEARSEC);
      fdotmax = 0.;
  } else {
      fdotmin = 1e-10;
      fdotmax = 1e-11;
  }

  // dimensionless spindown range
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
  
  sett->nmin = sett->fftpad*NAV*sett->B;
  sett->nmax = (sett->nfft/2 - NAV*sett->B)*sett->fftpad;

  printf("------------------------ Settings --------------------------\n");
  printf(" B         N            nfft         Fstat_nmin   Fstat_nmax\n");
  printf(" %-9.3f %-12d %-12d %-12d %-12d\n", sett->B, sett->N, sett->nfft, sett->nmin, sett->nmax);
  printf(" fpo       -fdotmin     fdotmax      Smin         -Smax\n");
  printf(" %-9.3f %-12.4e %-12.4e %-12.4e %-12.4e\n", sett->fpo, fdotmin, fdotmax, sett->Smin, sett->Smax);
  printf("------------------------------------------------------------\n");
    
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

  int i=0, j=0; 
  char dirname[1024], x[1332];
  DIR *dp;
  FILE *data;
  const char *dets[] = {"H1", "L1", "V1"};
  char det[DETNAME_LENGTH];

  // Test frame input directory 
  sprintf (dirname, "%s/%03d", opts->dtaprefix, opts->ident); 
  dp = opendir(dirname);
  if (dp) {
       closedir(dp);
  } else {
       printf("Can't open the input directory: %s", dirname);
       exit(EXIT_FAILURE);
  }

  // test availability of data for detectors
  for (i=0; i<3; i++) {
       if ( !strlen(opts->usedet) || (strlen(opts->usedet) && (strstr(opts->usedet, dets[i]))) ) {
	    // detector directory
	    memset(dirname, 0, sizeof(dirname));
	    sprintf (dirname, "%s/%03d/%s", opts->dtaprefix, opts->ident, dets[i]);
	    dp = opendir(dirname);
	    if (dp) {
		 closedir(dp);
		 sprintf(x, "%s/xdatsc_%03d_%04d%s.bin", dirname, opts->ident,
			                                 opts->band, opts->label);
		 data = fopen(x, "r");
		 if (data) {
		      fclose(data);
		      strncpy(ifo[j].xdatname, x, strlen(x));
		      strncpy(ifo[j].name, dets[i], DETNAME_LENGTH);
		      memset(x, 0, sizeof(x));
		      j++;
		 } else {
		      printf("Directory %s exists, but no input file found:\n%s missing...\n", dirname, x);
		      exit(EXIT_FAILURE);
		 }
		 
	    } else {
		 if ( strlen(opts->usedet) && (strstr(opts->usedet, dets[i])) ) {
		      printf("Can't open the input directory requied by -usedet: %s", dirname);
		      exit(EXIT_FAILURE);
		 }
	    } // if dp
       } // if
  } // for i
  
  sett->nifo = j;

  for(i=0; i<sett->nifo; i++) {

       printf("Using %s IFO as detector #%d... %s as input time series data\n", 
	      ifo[i].name, i, ifo[i].xdatname);

       // Virgo detector
       if(!strcmp("V1", ifo[i].name)) {

	    // Geographical latitude phi in radians
	    ifo[i].ephi = (43.+37./60.+53.0880/3600.)/RAD_TO_DEG;
	    // Geographical longitude in radians
	    ifo[i].elam = (10.+30./60.+16.1885/3600.)/RAD_TO_DEG;
	    // Height h above the Earth ellipsoid in meters
	    ifo[i].eheight = 51.884;
	    // Orientation of the detector gamma
	    ifo[i].egam = (135. - (19.0+25./60.0+57.96/3600.))/RAD_TO_DEG;

       // Hanford H1 detector
       } else if(!strcmp("H1", ifo[i].name )) {

	    // Geographical latitude phi in radians
	    ifo[i].ephi = (46+(27+18.528/60.)/60.)/RAD_TO_DEG;
	    // Geographical longitude in radians
	    ifo[i].elam = -(119+(24+27.5657/60.)/60.)/RAD_TO_DEG;
	    // Height h above the Earth ellipsoid in meters
	    ifo[i].eheight = 142.554;
	    // Orientation of the detector gamma
	    ifo[i].egam = 170.9994/RAD_TO_DEG;

       // Livingston L1 detector
       } else if(!strcmp("L1", ifo[i].name )) {

	    // Geographical latitude phi in radians
	    ifo[i].ephi = (30+(33+46.4196/60.)/60.)/RAD_TO_DEG;
	    // Geographical longitude in radians
	    ifo[i].elam = -(90+(46+27.2654/60.)/60.)/RAD_TO_DEG;
	    // Height h above the Earth ellipsoid in meters
	    ifo[i].eheight = -6.574;
	    // Orientation of the detector gamma
	    ifo[i].egam = 242.7165/RAD_TO_DEG;
       }
  }  // for i
  
  // todo: check if there are -usedet detectors without directory match
  
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

  // Modulation factors aa, bb for every NON-ZERO data point
  for (t=0; t<Np; t++) {
      if ( fabs(ifo->sig.xDat[t]) > DBL_MIN ) {
	  c = cosalfr*aux->cosmodf[t] + sinalfr*aux->sinmodf[t];
	  s = sinalfr*aux->cosmodf[t] - cosalfr*aux->sinmodf[t];
	  c2s = 2.*sqr(c);
	  cs = c*s;
	  
	  ifo->sig.aa[t] = c1*(2.-c2d)*c2s + c2*(2.-c2d)*2.*cs +
	      c3*c2sd*c + c4*c2sd*s - c1*(2.-c2d) + c5*c2d;

	  ifo->sig.bb[t] = c6*sindel*c2s + c7*sindel*2.*cs + 
	      c8*cosdel*c + c9*cosdel*s - c6*sindel;
      } else {
	  ifo->sig.aa[t] = 0.;
	  ifo->sig.bb[t] = 0.;
      }
  } 

} // modvir


int read_lines( Search_settings *sett,
		Command_line_opts *opts ){

  int i=0, lnum, j;
  char linefile[1200], line[512], *lfile;
  FILE *data;
  struct vlines {
       double f;
       int type;
       double offset;
       int iharm1, iharm2;
       double lwidth, rwidth;
       int det;
       int nline;
       char vfile[512];
  } vline[MAXVFILEL];
  char line_aux[MAXVFILEL][512];
  double fl, fr;
  
  glob_t globbuf;
  globbuf.gl_offs = 0;
  
  double fdotMax, Dfmaxmax;
  double normdtEmax[MAX_DETECTORS], normdEmax[MAX_DETECTORS];
  
  // calculate fdotMax from sett.Smin
  fdotMax = sett->Smin/(2.*M_PI*sett->dt*sett->dt);

  //
  // for each detector
  // calculate max line broadening due to demodulation
  //

  double *dE, *dtE;
  dE = (double *)calloc(3*sett->N, sizeof(double));
  dtE = (double *)calloc(3*sett->N, sizeof(double));

  for(int det=0; det < sett->nifo; det++) {
       
       // First derivative DetSSB (velocity)  
       for(i=0; i<sett->N-1; i++) { 
	    for(j=0; j<3; j++)  
		 dE[i*3+j] = fabs(ifo[det].sig.DetSSB[(i+1)*3+j] -
				  ifo[det].sig.DetSSB[i*3+j]); 
       } 

       double dEmax[3] = {0};

       // Find maximum absolute values 
       for(i=0; i<sett->N-1; i++) {
	    for(j=0; j<3; j++) 
		 if(dE[i*3+j] > dEmax[j]) dEmax[j] = dE[i*3+j];
       }
  
       // First derivative 
       for(i=0; i<sett->N-1; i++) {
	    for(j=0; j<3; j++)
		 dtE[i*3+j] = fabs(ifo[det].sig.DetSSB[(i+1)*3+j]*(i+1) -
				   ifo[det].sig.DetSSB[i*3+j]*i)*sett->dt; 
       }
  
       double dtEmax[3] = {0};
    
       // Find maximum absolute values 
       for(i=0; i<sett->N-1; i++) {
	    for(j=0; j<3; j++) 
		 if(dtE[i*3+j] > dtEmax[j]) dtEmax[j] = dtE[i*3+j];
       }

       normdtEmax[det]=0.;
       normdEmax[det]=0.;
       for(j=0; j<3; j++) {
	    normdtEmax[det] += pow(dtEmax[j], 2.); 
	    normdEmax[det]  += pow(dEmax[j], 2.);
       } 
  
  }

  // Free auxiliary allocs 
  free(dE);
  free(dtE);

  //
  // read veto files and find lines in band
  //
  int iold = 0;
  i = 0; // veto line number (global for all veto files)

  for(int det=0; det < sett->nifo; det++) {
       // search for all veto files matching pattern <data>/lines/<det_name>lines*.csv
       sprintf(linefile, "%s/lines/%slines*.csv", opts->dtaprefix, ifo[det].name);
       
       printf("[%s] Looking for %s ... ", ifo[det].name, linefile);
       glob(linefile, GLOB_DOOFFS, NULL, &globbuf);
       printf("%ld files match\n", globbuf.gl_pathc);
    
       for (size_t ifile = 0; ifile != globbuf.gl_pathc; ++ifile){
	    lfile = globbuf.gl_pathv[ifile];
	    printf("   [%s] %s ", ifo[det].name, lfile);

	    // Reading line data from the input file (data)
	    // Columns are: 
	    // 1 - frequency spacing (Hz) of comb (or frequency of single line)
	    // 2 - comb type (0 - singlet, 1 - comb with fixed width, 2 - comb with scaling width)
	    // 3 - frequency offset of 1st visible harmonic (Hz)
	    // 4 - index of first visible harmonic
	    // 5 - index of last visible harmonic
	    // 6 - width of left band (Hz)
	    // 7 - width of right band (Hz)
	    // 8 - comments
  
	    if ((data = fopen(lfile, "r")) == NULL) {
		 printf("Can't open %s \n Aborting!\n", lfile);
		 exit(EXIT_FAILURE);
	    }
       
	    int nline = 0; // line number in file
	    while (fgets(line, 512, data) != NULL) {
		 nline++;
		 // Skip comment lines beginning with '%'
		 if (*line == '%') continue;
		 vline[i].f       = atof(strtok(line,","));
		 vline[i].type    = atoi(strtok(NULL,","));
		 vline[i].offset  = atof(strtok(NULL,","));
		 vline[i].iharm1  = atoi(strtok(NULL,","));
		 vline[i].iharm2  = atoi(strtok(NULL,","));
		 vline[i].lwidth  = atof(strtok(NULL,","));
		 vline[i].rwidth  = atof(strtok(NULL,","));
		 vline[i].det     = det;
		 vline[i].nline = nline;
		 strcpy(vline[i].vfile, lfile);
		 /*printf("%f  %d  %f  %d   %d   %f   %f   %d   %s\n",  vline[i].f, vline[i].type,
		   vline[i].offset, vline[i].iharm1, vline[i].iharm2, vline[i].lwidth,
		   vline[i].rwidth, vline[i].nline, vline[i].vfile); */
		 if (++i > MAXL-1) {
		      printf("Too many lines in file %s, increase MAXL!\n", lfile); 
		      exit(EXIT_FAILURE);
		 }
	    }

	    printf("(%d data lines) \n", i-iold);
	    iold = i;
	    fclose(data);

       }
  
       globfree(&globbuf);
       
  } // det
  

  lnum = i;
  //  printf("Total number of data lines: %d\n", lnum);

  j=0; // index of line in band
  if(opts->narrowdown < 0.5*M_PI) j = sett->numlines_band;

  // Apply line widths 
  //------------------

  for(i=0; i<lnum; i++) {
    
       int k;
       
       switch(vline[i].type) {
	    
       // Singlet
       case 0:
	    
	    // Line width from the resampling broadening 
	    Dfmaxmax = 2.*fdotMax*(sett->N*sett->dt +
				   sqrt(normdtEmax[vline[i].det])) +
		 vline[i].f*sqrt(normdEmax[vline[i].det]);
	    
	    fl = vline[i].f - vline[i].lwidth - Dfmaxmax;
	    fr = vline[i].f + vline[i].rwidth + Dfmaxmax;

	    if (line_in_band(&fl, &fr, sett)) {
		 sett->lines[j][0] = fl;
		 sett->lines[j][1] = fr;
		 sprintf(line_aux[j], "singlet          [l.%d]%s", vline[i].nline, vline[i].vfile);
		 if (++j > MAXL-1) {printf("Too many lines, increase MAXL!\n"); exit(EXIT_FAILURE);}
	    }
	    break; 
	    

	    // Comb with fixed width. Vetoing the band 
	    // [offset+index*spacing-leftwidth, offset+index*spacing+rightwidth] 
       case 1:

	    for(k=vline[i].iharm1; k<=vline[i].iharm2; k++) {
		 
		 double linefreq = vline[i].offset + k*vline[i].f; 
		 // Line width from the resampling broadening 
		 Dfmaxmax = 2.*fdotMax*(sett->N*sett->dt +
					sqrt(normdtEmax[vline[i].det])) 
		      + linefreq*sqrt(normdEmax[vline[i].det]);

		 fl = linefreq - vline[i].lwidth - Dfmaxmax;
		 fr = linefreq + vline[i].rwidth + Dfmaxmax;

		 if (line_in_band(&fl, &fr,  sett)) {
		      sett->lines[j][0] = fl;
		      sett->lines[j][1] = fr;
		      sprintf(line_aux[j], "fcomb rank %4d  [l.%d]%s", k, vline[i].nline, vline[i].vfile);
		      if (++j > MAXL-1) {printf("Too many lines, increase MAXL!\n"); exit(EXIT_FAILURE);}
		 }
		 
	    }
	    break; 
      

	    // Comb with scaling-width. Vetoing the band 
	    // [offset+index*spacing-index*leftwidth, offset+index*spacing+index*rightwidth]       
       case 2:
      
	    for(k=vline[i].iharm1; k<=vline[i].iharm2; k++) { 

		 double linefreq = vline[i].offset + k*vline[i].f; 
		 // Line width from the resampling broadening 
		 Dfmaxmax = 2.*fdotMax*(sett->N*sett->dt +
					sqrt(normdtEmax[vline[i].det])) 
		      + linefreq*sqrt(normdEmax[vline[i].det]);

		 fl = linefreq - k*vline[i].lwidth - Dfmaxmax;
		 fr = linefreq + k*vline[i].rwidth + Dfmaxmax;
		 
		 if (line_in_band(&fl, &fr, sett)) {
		      sett->lines[j][0] = fl;
		      sett->lines[j][1] = fr;
		      sprintf(line_aux[j], "scomb rank %4d  [l.%d]%s", k, vline[i].nline, vline[i].vfile);
		      if (++j > MAXL-1) {printf("Too many lines, increase MAXL!\n"); exit(EXIT_FAILURE);}
		 }
	    } //k
	    break;
	    
       } // switch
  } // i

  printf("%d veto lines in band [Hz, radians, line info]:\n", j-sett->numlines_band);


  // save veto lines to a file
  sprintf(linefile, "%s/triggers_%03d_%04d%s.vlines", 
	  opts->prefix, opts->ident, opts->band, opts->label);
  if ( !(data = fopen(linefile, "w")) ) {
       printf("Can't open %s for writing!\n", linefile);
       exit(EXIT_FAILURE);
  }

  // scale veto lines to radians (narrowdown lines are already scaled)
  for(i=sett->numlines_band; i<j; i++) {
       fl = sett->lines[i][0];
       fr = sett->lines[i][1];
       sett->lines[i][0] = (sett->lines[i][0] - sett->fpo)/(sett->B)*M_PI;
       sett->lines[i][1] = (sett->lines[i][1] - sett->fpo)/(sett->B)*M_PI;
       
       printf("   %f  %f  %f  %f  %s\n",
	      fl, fr, sett->lines[i][0], sett->lines[i][1], line_aux[i]);
       fprintf(data, "   %f  %f  %f  %f  %s\n",
	      fl, fr, sett->lines[i][0], sett->lines[i][1], line_aux[i]);
  }
  
  fclose(data);
  printf("Wrote veto lines in band to: %s\n", linefile);

  lines_veto_fraction(sett, sett->numlines_band, j, opts->veto_flag);
  
  // set number of veto lines only if veto option is given
  if (opts->veto_flag) {
       printf("Veto lines will be applied!\n");
       sett->numlines_band = j;
  } else {
       printf("Veto lines WILL NOT be applied!\n");
  }
  // printf("Number of known lines in band: %d\n", sett->numlines_band);

  printf("Excluded frequencies in band (incl. narrowdown, in radians):\n"); 
  for(i=0; i<sett->numlines_band; i++) 
       printf("   %f %f\n", sett->lines[i][0], sett->lines[i][1]);
  
  
  return 0; 
  
}


int line_in_band(double* fl, double* fr, Search_settings* sett ) {

    double bs, be;        // Band start and end  

    bs = sett->fpo; 
    be = sett->fpo + sett->B; 
    
    if (!(*fr < bs || *fl > be)) { 
	 if (*fl < bs) *fl = bs;
	 if (*fr > be) *fr = be;
	 //printf("[line in band] %f %f  %f  %f\n", bs, be, *fl, *fr);
	 return(1);
    }
    return(0);
}  



void narrow_down_band(Search_settings* sett, Command_line_opts *opts) {
  
  // Adding excluding ranges near the edges to the known lines list 
  sett->lines[0][0] = 0;
  sett->lines[0][1] = M_PI_2 - opts->narrowdown;
  sett->lines[1][0] = M_PI_2 + opts->narrowdown;
  sett->lines[1][1] = M_PI;
  
  sett->numlines_band = 2;
  printf("Band is narrowed-down to [%f, %f] (narrowdown=%f)\n", sett->lines[0][1], sett->lines[1][0],
	 opts->narrowdown/M_PI);
  
}



void lines_veto_fraction(Search_settings* sett, int lf, int le, int vflag) {
     
  // lf - index of first line, le - index of last line
  int i; 
  double ll=0., gap=0.;
  
  // Sorting veto lines in band (1st then 2nd column) 
  qsort(&sett->lines[lf], le-lf, 2*sizeof(double), compared2c); 

  for(i=lf; i<le; i++) {
    
       // Looking for a gap between lines
       if(sett->lines[i][0] >= ll) {
	    gap += sett->lines[i][0] - ll;
	    ll = sett->lines[i][1];
       } else {
	    if (ll < sett->lines[i][1]) ll = sett->lines[i][1];
       }
  }
  if ( ll < M_PI) gap += M_PI - ll;
  
  printf("Band veto fraction = %6.4f\n", (M_PI-gap)/M_PI);
  
  if( (gap <= 1.e-10) && vflag) {
       printf("This band is fully vetoed. My work here is done, exiting...\n"); 
       exit(EXIT_SUCCESS); 
  }
 
}

