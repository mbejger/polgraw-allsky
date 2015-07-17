#ifndef __JOBCORE_H__
#define __JOBCORE_H__

#include "struct.h"

void search(
	Search_settings *sett,
	Command_line_opts *opts,
	Search_range *s_range,
  FFTW_plans *plans,
	FFTW_arrays *fftw_arr,
	Aux_arrays *aux,
	int *Fnum,
	double *F);


double* job_core(
  int pm,                     // hemisphere
  int spinpos,                // no. of sky position in the spotlight range file
  Search_settings *sett,      // search settings
  Command_line_opts *opts,    // cmd opts
  Search_range *s_range,      // range for searching
  FFTW_plans *plans,          // plans for fftw
  FFTW_arrays *fftw_arr,      // arrays for fftw
  Aux_arrays *aux,            // auxiliary arrays
  double *F,                  // F-statistics array
  int *sgnlc,                 // reference to array with the parameters
                              // of the candidate signal
                              // (used below to write to the file)
  int *FNum);                 // Candidate signal number
      

#endif
