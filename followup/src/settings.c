#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>

#include "settings.h"
#include "auxi.h"
#include "struct.h"


/*************************************************************** 
Search settings: number of data points, spindown range, bandwidth etc.
***************************************************************/

void search_settings(Search_settings* sett) {
	double dt, B, oms, omr, Smin, Smax;
	int N, nd;
	dt = sett->dt;				// data sampling time:  
                                    		// set in handle_opts() from the command line
                                    		// (the default value is dt=0.5)
	B = 0.5/dt;                       	// Bandwidth
	oms = 2.*M_PI*(sett->fpo)*dt;     	// Dimensionless angular frequency
	omr = C_OMEGA_R*dt;
	N = round (sett->nod*C_SIDDAY/dt);      // No. of data points
//Put by hand: #mb ranges of spindown (RDC O1) 
	double fdotmin, fdotmax; 
	fdotmin = 0.5e-8; 
	fdotmax = 0.5e-9; 
	Smax = 2.*M_PI*fdotmin*dt*dt; 
	Smin = 2.*M_PI*fdotmax*dt*dt;
	nd = 2;     				// Degree of freedom, 
	      					// (2*nd = deg. no ofrees of freedom for chi^2)

	sett->B=B;          			// bandwidth
	sett->oms=oms;      			// dimensionless angular frequency
	sett->omr=omr;      			// C_OMEGA_R * dt
	sett->N=N;          			// number of data points
	sett->Smin=Smin;    			// minimum spindown
	sett->Smax=Smax;    			// maximum spindown
	sett->nd=nd;        			// degrees of freedom
} // search settings()  



/*************************************************************** 
Detectors settings:
Network of detectors' discovery: finds subdirectories in the main input directory,
which by convention should be named like V1, L1, H1 and which contain
input data and ephemerids; writes appropriate detector-related data into structs. 
***************************************************************/ 

void detectors_settings(Search_settings* sett, Command_line_opts *opts) {

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
// Subdirectory names checkup: 
// check if it's a dir
// name is 2 char long
// not a directory name of the type "./" or ".."
// if usedef is not set (length equal 0), or is set and dir name is substring of it 
			if((ep->d_type == DT_DIR) && 
			(strlen(ep->d_name)==DETNAME_LENGTH) && 
			(strncmp(&ep->d_name[0],".",1)) && 
			(!strlen(opts->usedet) || (strlen(opts->usedet) && (strstr(opts->usedet, ep->d_name))))) { 
	 			FILE *data;
// Input time-domain data handling
// We assume that in each subdirectory corresponding to the detector the input data will look as following: 
				sprintf(x, "%s/%03d/%s/xdatc_%03d_%04d%s.bin", opts->dtaprefix, opts->ident, ep->d_name, opts->ident, opts->band, opts->label);
				if((data = fopen(x, "r")) != NULL) {
					xnames[i]   = calloc(strlen(x)+1, sizeof(char));
					detnames[i] = calloc(DETNAME_LENGTH+1, sizeof(char));
					strncpy(xnames[i], x, strlen(x));
					strncpy(detnames[i], ep->d_name, DETNAME_LENGTH);
					i++;
				} else { 
				printf("Directory %s exists, but no input file found:\n%s missing...\n", ep->d_name, x);  
				} 
				fclose(data); 
				memset(x, 0, sizeof(x));
			} //if directory is right
		}// while reading directories() 
		(void) closedir(dp);
	} 
	else{ 
		perror ("Couldn't open the input directory...");
	}
	sett->nifo=i;      // number of detectors  
	if(sett->nifo) { 
		printf("Settings - number of detectors: %d\n", sett->nifo); 
	} 
	else { 
		printf("No subdirectories with detector data found. Exiting...\n"); 
		exit(EXIT_FAILURE);
	} 
//Loop over detectors 
	for(i=0; i<sett->nifo; i++) { 
		if(!strcmp("V1", detnames[i])) {
// Virgo detector
			strncpy(ifo[i].xdatname, xnames[i], strlen(xnames[i]));
			strncpy(ifo[i].name, detnames[i], DETNAME_LENGTH);
// Geographical latitude phi in radians
			ifo[i].ephi = (43.+37./60.+53.0880/3600.)/RAD_TO_DEG;
// Geographical longitude in radians
			ifo[i].elam = (10.+30./60.+16.1885/3600.)/RAD_TO_DEG;
// Height h above the Earth ellipsoid in meters
			ifo[i].eheight = 51.238;
// Orientation of the detector gamma
			ifo[i].egam = (135. - (19.0+25./60.0+57.96/3600.))/RAD_TO_DEG;
			if (opts->gauss_flag){  
				printf("Using %s IFO as detector #%d... Gaussian noise as input time series data\n", ifo[i].name, i);
			}
			else {
				printf("Using %s IFO as detector #%d... %s as input time series data\n", 
				ifo[i].name, i, ifo[i].xdatname);
			}
	    	} 
		else if(!strcmp("H1", detnames[i])) {
// Hanford H1 detector
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
			if (opts->gauss_flag){  
				printf("Using %s IFO as detector #%d... Gaussian noise as input time series data\n", ifo[i].name, i);
			}
			else {
				printf("Using %s IFO as detector #%d... %s as input time series data\n", 
				ifo[i].name, i, ifo[i].xdatname);
			}
		} 
		else if(!strcmp("L1", detnames[i])) {
// Livingston L1 detector
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
			if (opts->gauss_flag){  
				printf("Using %s IFO as detector #%d... Gaussian noise as input time series data\n", ifo[i].name, i);
			}
			else {
				printf("Using %s IFO as detector #%d... %s as input time series data\n", 
				ifo[i].name, i, ifo[i].xdatname);
			}
		} 
		else {
			printf("Meh, unknown detector %s (see settings.c) Exiting...\n", detnames[i]);
			exit(EXIT_FAILURE);
		}
	}//loop over detectors 
// Memory free for detnames and xdatnames
	for(i=0; i<sett->nifo; i++) { 
		free(detnames[i]);
		free(xnames[i]); 
	} 
	free(detnames); 
	free(xnames); 
} // detectors settings()

/*************************************************************** 
Coefficients of the amplitude modulation functions of the Virgo detector
***************************************************************/ 
void rogcvir(Detector_settings *ifo) {
/* In the notation of Phys. Rev. D 58, 063001 (1998):
ephi = lambda (geographical latitude phi in radians)
egam = gamma (orientation of the detector)
(see modvir function in jobcore.c for Eqs. 12 and 13)
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
} // rogcvir()

/*************************************************************** 
Amplitude modulation of the signal
***************************************************************/ 
void modvir(double sinal, double cosal, double sindel, double cosdel, int Np, Detector_settings *ifo, Aux_arrays *aux, double *sigaa, double *sigbb) {
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
// Modulation factors aa and bb  
		sigaa[t] = c1*(2.-c2d)*c2s + c2*(2.-c2d)*2.*cs + c3*c2sd*c + c4*c2sd*s - c1*(2.-c2d) + c5*c2d;
		sigbb[t] = c6*sindel*c2s + c7*sindel*2.*cs + c8*cosdel*c + c9*cosdel*s - c6*sindel;
  	}//loop over data points 
} // modvir()
