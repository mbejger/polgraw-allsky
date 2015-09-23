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
#include <gsl/gsl_linalg.h>
#include <time.h>
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

  for(i=0; i<16; i++) {
    printf("%d %le\n", i, sett.M[i]);
  }

  printf("These are the eigenvectors:\n"); 
  for(i=0; i<4; i++) {  
    for(j=0; j<4; j++) 
        printf("%.15le ", sett.eigvec[j][i]); 
    printf("\n");
  } 

  printf("These are the eigenvalues:\n");
  for(i=0; i<4; i++) printf("%.15le\n", sett.eigval[i]);

  // Search settings 
  search_settings(&sett); 

  printf("dt, oms from settings: %f %f\n", sett.dt, sett.oms); 

  read_trigger_files(&sett, &opts, &trig); 

  printf("Triggers read, in total %d\n", trig.num_of_trig); 

/* for(i=0; i<trig.num_of_trig; i++) 
    printf("%le %le %le %le %f %d\n", 
    trig.f[i], trig.s[i], trig.a[i], trig.d[i], trig.snr[i], trig.fr[i]); 
*/ 

  convert_to_linear(&sett, &opts, &trig); 

  // Cleanup: free arrays at the end
   
  free(trig.f);
  free(trig.s);
  free(trig.a);
  free(trig.d);
  free(trig.snr);
  free(trig.fr);

  free(trig.fi);
  free(trig.si);
  free(trig.ai);
  free(trig.di);


  return 0; 
	
}
