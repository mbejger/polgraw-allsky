#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <getopt.h>
#include <gsl/gsl_linalg.h>
#include <time.h>
#include <dirent.h>

#include "auxi.h"
#include "settings.h"
#include "struct.h"
#include "jobcore.h"
#include "init.h"

// Default output and data directories

#ifndef PREFIX
#define PREFIX ./candidates
#endif

#ifndef DTAPREFIX
#define DTAPREFIX .
#endif


int main (int argc, char* argv[]) {

  Command_line_opts opts;
  Search_settings sett;
  Aux_arrays aux_arr;
  double *F; 			// F-statistic array
  Signals sig; 			// signals
  Search_range s_range; 
  int i; 

  // Command line options 
  handle_opts(argc, argv, &opts, &sett);  
	
  // Output data handling
  struct stat buffer;

  if (stat(opts.prefix, &buffer) == -1) {
    if (errno == ENOENT) {
      // Output directory apparently does not exist, try to create one
      if (mkdir(opts.prefix, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH	\
		| S_IXOTH) == -1) {
	perror (opts.prefix);
	return 1;
      }
    } else { // can't access output directory
      perror (opts.prefix);
      return 1;
    }
  }
 
  // Grid data 
  read_grid(&sett, &opts);	
	
  // Search and detector network settings
  settings(&sett, &opts, &aux_arr); 

  // Amplitude modulation functions for each detector  
  for(i=0; i<sett.nifo; i++)   
    rogcvir(&ifo[i]); 

  // Array initialization
  init_arrays(&sett, &opts, &aux_arr, &F);

  // Set search range 
  set_search_range(&sett, &opts, &s_range);

  // Cleanup & memory free 
  cleanup(&sett, &opts, &s_range, 
//	&fftw_arr, 
//	&fftw_plans, 
	&aux_arr, F);

  return 0; 

/*

  // Array initialization
  init_arrays(&sig, &aux_arr, &F, &opts, &sett);

  // FFT plans 
  FFTW_plans fftw_plans;
  FFTW_arrays fftw_arr;
  plan_fftw(&fftw_plans, &fftw_arr, &sig, &aux_arr, &sett, &opts);

  // Checkpointing
  int Fnum;			// candidate signal number
  read_checkpoints(&s_range, &Fnum, &opts);
	
  // main search job
  search(
	&sett,
	&opts,
	&s_range,
	&fftw_arr,
	&sig,
	&fftw_plans,
	&aux_arr,
	&amod,
	&Fnum,
	F);

  // state file zeroed at the end of the run
  FILE *state;
  if(opts.checkp_flag) {
    state = fopen (opts.qname, "w");
    fclose (state);
  }
	
  // Cleanup & memory free 
  cleanup(
	&sett,
	&opts, 
	&s_range, 
	&fftw_arr, 
	&sig, 
	&fftw_plans, 
	&aux_arr, 
	&amod,
	F);
*/ 

  return 0; 
	
}
