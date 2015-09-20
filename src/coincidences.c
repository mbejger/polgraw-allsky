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
#include  "init.h" 
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

  printf("These are the eigenvectors:\n"); 
  for(i=0; i<4; i++) {  
    for(j=0; j<4; j++) 
        printf("%g ", sett.eigvec[i][j]); 
    printf("\n");     
  } 

  // Search settings
  search_settings(&sett); 

  printf("dt, oms from settings.h %f %f\n", sett.dt, sett.oms); 

  

  int num_of_trig; 
  num_of_trig = read_trigger_files(&sett, &opts, &trig); 

  printf("The triggers read, in total %d\n", num_of_trig); 
  for(i=0; i<num_of_trig; i++) 
    printf("%f %f %f %f %f\n", trig.f[i], trig.s[i], trig.a[i], trig.d[i], trig.snr[i]); 

  return 0; 
	
}
