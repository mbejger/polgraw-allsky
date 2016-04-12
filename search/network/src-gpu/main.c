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

#include "floats.h"
#include "auxi.h"
#include "struct.h"
#include "settings.h"
#include "jobcore.h"
#include "init.h"

// Default output and data directories

#ifndef PREFIX
#define PREFIX ./candidates
#endif

#ifndef DTAPREFIX
#define DTAPREFIX .
#endif

Detector_settings ifo[MAX_DETECTORS];

int main (int argc, char* argv[]) {

  Command_line_opts opts;
  Search_settings sett;
  Search_range s_range; 
  Aux_arrays aux_arr;
  FLOAT_TYPE *F_d; 			  // F-statistic array
  int i; 

  /* init CUDA device,
     enable mapped memory mode;
     CUDA_DEV is set in Makefile */
  if (cuinit(CUDA_DEV) == -1) {
    printf("\nGPU device initialization error!\n");
    exit(EXIT_FAILURE);
  }

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
    } else { // can't access output directory
      perror (opts.prefix);
      return 1;
    }
  }
 
  // Grid data 
  read_grid(&sett, &opts);	
	
  // Search settings
  search_settings(&sett); 

  // Detector network settings
  detectors_settings(&sett, &opts); 

  printf("detectors set!\n");

  // Array initialization
  init_arrays(&sett, &opts, &aux_arr, &F_d);


  printf("arrays initialized!\n");

  // Amplitude modulation functions for each detector  
  for(i=0; i<sett.nifo; i++)   
    rogcvir(&ifo[i]); 

  printf("after rogcvir\n");
#if 0
  // Grid search range
  if(strlen(opts.addsig))
    // If addsig switch used, add signal from file, 
    // search around this position (+- gsize)
    add_signal(&sett, &opts, &aux_arr, &s_range); 
  else 
#endif

    // Set search range from range file  
    set_search_range(&sett, &opts, &s_range);

  // FFT plans 
  FFT_plans fft_plans;
  FFT_arrays fft_arr;
  //  plan_fft( &sett, &opts, &fft_plans, &fft_arr, &aux_arr );
  plan_fft( &sett, &fft_plans, &fft_arr );
  printf("fft plan ready\n");


  // Checkpointing
  int Fnum=0;			        // candidate signal number
  read_checkpoints(&opts, &s_range, &Fnum);


  // main search job
  search(&sett, &opts, &s_range, 
	 &fft_plans, &fft_arr, &aux_arr,
	 &Fnum, F_d);

  // state file zeroed at the end of the run
  FILE *state;
  if(opts.checkp_flag) {
    state = fopen (opts.qname, "w");
    fclose (state);
  }
	
  // Cleanup & memory free 
  cleanup( &sett, &opts, &s_range, &fft_plans, &fft_arr, &aux_arr, F_d );

  return 0; 
	
}
