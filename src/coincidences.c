#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <getopt.h>
#include <gsl/gsl_linalg.h>
//#include <time.h>
#include <dirent.h>

#include "auxi.h"
#include "init.h" 
#include "settings.h"
#include "struct.h"

// Default output and data directories
#ifndef PREFIX
#define PREFIX ./coinc-results
#endif

#ifndef DTAPREFIX
#define DTAPREFIX ./candidates
#endif


int main (int argc, char* argv[]) {

  Search_settings sett;
  Command_line_opts_coinc opts;
  Candidate_triggers trig; 
  int i, j; 

  // Command line options 
  handle_opts_coinc(&sett, &opts, argc, argv);  

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
 
  // Manage grid matrix  
  manage_grid_matrix(&sett, &opts);	

  // Search settings 
  search_settings(&sett); 

  printf("Settings - dt: %f, oms: %f\n", sett.dt, sett.oms); 

  read_trigger_files(&sett, &opts, &trig); 

  // Free arrays at the end 
  free(sett.M); 

  return 0; 
	
}
