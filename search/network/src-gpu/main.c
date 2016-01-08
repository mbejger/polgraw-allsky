#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <getopt.h>
#include <time.h>

#include "auxi.h"
#include "struct.h"
#include "settings.h"
#include "jobcore.h"
#include "init.h"

/* Default output and data directories */

#ifndef PREFIX
#define PREFIX ./candidates
#endif

#ifndef DTAPREFIX
#define DTAPREFIX .
#endif


/*
 * Main
 */
int main (int argc, char* argv[]) {

  Command_line_opts opts;
  Detector_settings sett;
  FLOAT_TYPE *cu_F;
  Search_range s_range;
  Arrays arr;

  /* init */

  /*
    ############ Command line options ################
  */

  handle_opts(argc, argv, &opts, &sett); 
	
  /*
    ############ Output data handling ################
  */
  struct stat buffer;

  if (stat(opts.prefix, &buffer) == -1) {
    if (errno == ENOENT) {
      /* Output directory apparently does not exist, try to create one */
      if (mkdir(opts.prefix, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH	\
		| S_IXOTH) == -1) {
	perror (opts.prefix);
	return 1;
      }
    } else { /* can't access output directory */
      perror (opts.prefix);
      return 1;
    }
  }

  /*
    ############ Read grid data ################
  */
  read_grid(&sett, &opts);

  /*
    ############ Detector settings ################
  */
  settings (&sett, &opts, &arr); //detector settings
	
  /*
    ############ Amplitude modulation functions ################
  */

  Ampl_mod_coeff amod;
  rogcvir(&amod, &sett); //Virgo amplitude modulation function coefficients

  /*
    ############ Initialize arrays ################
  */

  init_arrays(&arr, &cu_F, &opts, &sett);

  /*
    ############ Set searching ranges ################
  */
  set_search_range(&s_range, &opts, &sett);

  /*
    ############ FFT planning ################
  */
  FFT_plans fft_plans;
  plan_fft(&fft_plans, &arr, &sett, &opts);

  /*
    ############ Read checkpoints ################
  */
  int Fnum; //candidate signal number
  read_checkpoints(&s_range, &Fnum, &opts);
	
  /*
    ############ Main job ################
  */

  search(
	 &sett,
	 &opts,
	 &s_range,
	 &arr,
	 &fft_plans,
	 &amod,
	 &Fnum,
	 cu_F
	 );
	
  FILE *state;
  // state file zeroed at the end of calculations
  if(opts.checkp_flag) {
    state = fopen (opts.qname, "w");
    fclose (state);
  }
	
  cleanup(&sett, &opts, &s_range, &arr, &fft_plans, &amod, cu_F);
	
  return 0;
}

