// ISO C behavioral defines
#define __STDC_WANT_LIB_EXT1__ 1

// Standard C includes
#include <stdio.h>      // fopen_s
#include <stdlib.h>
//#include <unistd.h>
#include <math.h>
#include <complex.h>
//#include <fftw3.h>
#include <string.h>
#include <errno.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#if defined _MSC_VER
#include <direct.h>
#endif

#include <fcntl.h>
#include <XGetopt.h>
//#include <gsl/gsl_linalg.h>
#include <time.h>

#include "floats.h"
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

Detector_settings ifo[MAX_DETECTORS];

int main (int argc, char* argv[]) {

  Command_line_opts opts;
  OpenCL_handles cl_handles;
  Search_settings sett;
  Search_range s_range; 
  Aux_arrays aux_arr;
  FLOAT_TYPE *F_d; 			  // F-statistic array
  int i; 

  // Command line options 
  handle_opts(&sett, &opts, argc, argv);

  // Initialize OpenCL
  init_opencl(&cl_handles, &opts);
	
  // Output data handling
  struct stat buffer;
  
  if (stat(opts.prefix, &buffer) == -1)
  {
      if (errno == ENOENT)
    {
      // Output directory apparently does not exist, try to create one
#ifdef WIN32
        if (_mkdir(opts.prefix) == -1)
#else
        if (mkdir(opts.prefix, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) == -1)
#endif
      {
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

  // Array initialization
  init_arrays(&sett, &opts, &aux_arr, &F_d);

  // Amplitude modulation functions for each detector  
  for(i=0; i<sett.nifo; i++)   
    rogcvir(&ifo[i]); 

  // Set search range from range file  
  set_search_range(&sett, &opts, &s_range);

  // FFT plans 
  FFT_plans fft_plans;
  FFT_arrays fft_arr;
  plan_fft(&sett, &fft_plans, &fft_arr);

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
    errno_t err = fopen_s(&state, opts.qname, "w");
    if (err)
        perror("Error zeroing out state file.");
    fclose(state);
  }
	
  // Cleanup & memory free 
  cleanup(&sett, &opts, &s_range, 
          &fft_plans, &fft_arr, &aux_arr, F_d);

  return 0; 
	
}
