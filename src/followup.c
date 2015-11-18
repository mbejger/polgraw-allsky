#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
//#include <fcntl.h>
#include <getopt.h>
#include <gsl/gsl_linalg.h>
#include <time.h>
#include <dirent.h>

#include "auxi.h"
#include "settings.h"
#include "struct.h"
#include "init.h"
//#include "timer.h"

#include <assert.h>
#if defined(SLEEF)
//#include "sleef-2.80/purec/sleef.h"
#include <sleefsimd.h>
#elif defined(YEPPP)
#include <yepLibrary.h>
#include <yepCore.h>
#include <yepMath.h>
#endif

// Default output and data directories

#ifndef PREFIX
#define PREFIX ./FSTAT_OUT
#endif

#ifndef DTAPREFIX
#define DTAPREFIX .
#endif

// Fstat function declaration
double *Fstatnet(Search_settings *sett, Command_line_opts *opts, Aux_arrays *aux, double *F, double *sgnlo){

	double nSource[3];
	double xa_real = 0., xa_imag = 0., xb_real = 0., xb_imag = 0., xasum_real = 0., xasum_imag = 0., xbsum_real = 0., xbsum_imag = 0.;
	double shft1, cosPH, sinPH, phase[sett->N];
  	double sinalt, cosalt, sindelt, cosdelt;
	int i = 0, n = 0; 
	double *fstat_out = malloc(6*sizeof(int)); //output
	double aa = 0., bb = 0., aaa = 0., bbb = 0.;
#ifdef YEPPP
//#define VLEN 1 //2048
    int VLEN = sett->N;
    yepLibrary_Init();

    Yep64f _sph[VLEN];
    Yep64f _cph[VLEN];
    enum YepStatus status;

#endif

//	printf("sgnlo: %Le %Le %Le %Le\n", sgnlo[0], sgnlo[1], sgnlo[2], sgnlo[3]);
//	printf("sett->N = %d\n", sett->N);

	sinalt = sin(sgnlo[3]);
	cosalt = cos(sgnlo[3]);
	sindelt = sin(sgnlo[2]);
	cosdelt = cos(sgnlo[2]);

	nSource[0] = cosalt*cosdelt;
	nSource[1] = sinalt*cosdelt;
	nSource[2] = sindelt;


//	printf("sinalt = %f cosalt = %f sindelt = %f cosdelt = %f\n", sinalt, cosalt, sindelt, cosdelt);

//From jobcore.c, line 237 
//Loop for each detector 
  	for(n=0; n < sett->nifo; ++n) { 
		modvir(sinalt, cosalt, sindelt, cosdelt, 
	   		sett->N, &ifo[n], aux);  

// Calculate detector positions with respect to baricenter
// Copied from jobcore.c, line 248



		xa_real = 0.;	
		xa_imag = 0.;
		xb_real = 0.;	
		xb_imag = 0.;
		aa = 0.;
		bb = 0.;
//Inside loop
    		for(i=0; i<sett->N; ++i) {

      			ifo[n].sig.shft[i] = nSource[0]*ifo[n].sig.DetSSB[i*3]
		         	+ nSource[1]*ifo[n].sig.DetSSB[i*3+1]
		         	+ nSource[2]*ifo[n].sig.DetSSB[i*3+2];
    
// Phase modulation function
// Copied from jobcore.c, line 265

			phase[i] = sgnlo[0]*(i + ifo[n].sig.shft[i]) 
				+ sgnlo[1]*i*i + (sett->oms 
				+ 2*sgnlo[1]*i)*ifo[n].sig.shft[i];
		}			

		status = yepMath_Cos_V64f_V64f(phase, _cph, VLEN);
		assert(status == YepStatusOk);
		status = yepMath_Sin_V64f_V64f(phase, _sph, VLEN);
		assert(status == YepStatusOk);

//      			exph = cosPH - I*sinPH;
// Matched filter 
// Copied from jobcore.c, line 276 and 337


		for (i = 0; i<sett->N; ++i){
			xa_real = xa_real + ifo[n].sig.xDat[i]*ifo[n].sig.aa[i]*_cph[i]/ifo[n].sig.sig2;
			xa_imag = xa_imag - ifo[n].sig.xDat[i]*ifo[n].sig.aa[i]*_sph[i]/ifo[n].sig.sig2;
			xb_real = xb_real + ifo[n].sig.xDat[i]*ifo[n].sig.bb[i]*_cph[i]/ifo[n].sig.sig2;
			xb_imag = xb_imag - ifo[n].sig.xDat[i]*ifo[n].sig.bb[i]*_sph[i]/ifo[n].sig.sig2;		
			
     			aa += sqr(ifo[n].sig.aa[i]);
      			bb += sqr(ifo[n].sig.bb[i]);
	 
		}	// End of inside loop


    		aaa += aa/ifo[n].sig.sig2; 
    		bbb += bb/ifo[n].sig.sig2;

		xasum_real += xa_real;	
		xasum_imag += xa_imag;
		xbsum_real += xb_real;
		xbsum_imag += xb_imag;

	}	// End of detector loop

// F - statistic

	fstat_out[5] = - ((( sqr(xasum_real) + sqr(xasum_imag))/aaa)
			+ ((sqr(xbsum_real) +sqr(xbsum_imag))/bbb));

// Amplitude estimates
	fstat_out[0] = 2*xasum_real/aaa;
	fstat_out[1] = 2*xbsum_real/bbb;
	fstat_out[2] = -2*xasum_imag/aaa;
	fstat_out[3] = -2*xbsum_imag/bbb;

// Signal-to-noise ratio
	fstat_out[4] = sqrt(2*(-fstat_out[5]-2));		

	printf("%Le %Le %Le %Le %Le %Le\n", fstat_out[0], fstat_out[1], fstat_out[2],
		fstat_out[3], fstat_out[4], fstat_out[5]);
	return fstat_out;


}

int main (int argc, char* argv[]) {

	Search_settings sett;	
	Command_line_opts opts;
  	Search_range s_range; 
  	Aux_arrays aux_arr;
  	double *F; 			  // F-statistic array
  	int i, j; 
	double *results;		  // Vector with results from Fstatnet function
	double s1, s2, s3, s4;
	double sgnlo[4];
//	double time_elapsed;
	double tt, tdif;
	clock_t t1, t2;

// Input file
	FILE *sig;
	sig=fopen("/home/magdalena/polgraw-allsky/allM/sigFine", "r");
	if (sig !=NULL) {
		fscanf (sig, "%lf %lf %lf %lf\n", &sgnlo[0], &sgnlo[1], &sgnlo[2], &sgnlo[3]);
	}
	else (puts("Cannot open data file"));
		fclose(sig);

// Command line options 
	handle_opts(&sett, &opts, argc, argv);  
// Output data handling
  	struct stat buffer;

  	if (stat(opts.prefix, &buffer) == -1) {
    		if (errno == ENOENT) {
// Output directory apparently does not exist, try to create one
      			if(mkdir(opts.prefix, S_IRWXU | S_IRGRP | S_IXGRP 
          			| S_IROTH	| S_IXOTH) == -1) {
	      			perror (opts.prefix);
	      			return 1;
      			}
    		} 
		else { // can't access output directory
      			perror (opts.prefix);
      			return 1;
    		}
  	}
// Search settings
  	search_settings(&sett); 
// Detector network settings
  	detectors_settings(&sett, &opts); 
// Array initialization
  	init_arrays(&sett, &opts, &aux_arr, &F);
// Amplitude modulation functions for each detector  
	for(i=0; i<sett.nifo; i++)   
		rogcvir(&ifo[i]); 
// F - statistic

	results = Fstatnet(&sett, &opts, &aux_arr, F, sgnlo);

// Cleanup & memory free 
//  	cleanup(&sett, &opts, &aux_arr, F);

	return 0;

}


