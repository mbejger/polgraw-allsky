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

#ifndef CODEVER
#define CODEVER unknown
#endif

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
  Search_range s_range; 
  Aux_arrays aux_arr;
  double *F; 			  // F-statistic array
  int i; 

#define QUOTE(x) #x
#define STR(macro) QUOTE(macro)
#define CVSTR STR(CODEVER)

  printf("Code version : " CVSTR "\n");

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

  // Array initialization and reading the ephemerids 
  init_arrays(&sett, &opts, &aux_arr, &F);

  // Narrowing-down the band (excluding the edges 
  // according to the opts.narrowdown parameter) 
  if(opts.narrowdown < 0.5*M_PI)
    narrow_down_band(&sett, &opts);

  // Reading known lines data from external files 
  if(opts.veto_flag) { 
    for(i=0; i<sett.nifo; i++) {
      printf("Reading known lines data for %s from %s\n", ifo[i].name, opts.dtaprefix);
      read_lines(&sett, &opts, &ifo[i]);
    }

    // Vetoing known lines in band 
    lines_in_band(&sett, &opts); 
  } 

  // If excluded parts of band, list them
  // and check if the band isn't fully vetoed 
  if(sett.numlines_band) {     

    int k; 
    printf("list of excluded frequencies in band (in radians):\n"); 
    for(k=0; k<sett.numlines_band; k++) 
      printf("%f %f\n", sett.lines[k][0], sett.lines[k][1]);

    fraction_of_band_vetoed(&sett, &opts); 

  } else 
    printf("Veto fraction for band %04d: %f\n", opts.band, 0.);  

  // Amplitude modulation functions for each detector  
  for(i=0; i<sett.nifo; i++)   
    rogcvir(&ifo[i]); 

  // Grid search range
  if(strlen(opts.addsig)) { 
    // If addsig switch used, add signal from file, 
    // search around this position (+- gsize)
    add_signal(&sett, &opts, &aux_arr, &s_range); 
  } else 
    // Set search range from range file  
    set_search_range(&sett, &opts, &s_range);

  // FFT plans 
  FFTW_plans fftw_plans;
  FFTW_arrays fftw_arr;
  plan_fftw(&sett, &opts, &fftw_plans, &fftw_arr, &aux_arr);
  if (strlen(opts.getrange)) exit(EXIT_SUCCESS);

  // Checkpointing
  int Fnum=0;			        // candidate signal number
  read_checkpoints(&opts, &s_range, &Fnum);

  // main search job
  search(&sett, &opts, &s_range, 
        &fftw_plans, &fftw_arr, &aux_arr,
	      &Fnum, F);

  // state file zeroed at the end of the run
  FILE *state;
  if(opts.checkp_flag) {
    state = fopen (opts.qname, "w");
    fclose (state);
  }
	
  // Cleanup & memory free 
  cleanup(&sett, &opts, &s_range, 
          &fftw_plans, &fftw_arr, &aux_arr, F);

  return 0; 
	
}
