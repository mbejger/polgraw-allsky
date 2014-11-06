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
			int pm,			// hemisphere
			int mm,			// grid 'sky position'
			int nn,			// other grid 'sky position'
			Search_settings *sett, // search settings
			Command_line_opts *opts, // cmd opts
			Search_range *s_range,	 // range for searching
      FFTW_plans *plans,       // plans for fftw
			FFTW_arrays *fftw_arr,   // arrays for fftw
			Aux_arrays *aux, 			 // auxiliary arrays
			double *F,		// F-statistics array
			int *sgnlc,		// reference to array with the parameters
								// of the candidate signal
								// (used below to write to the file)
			int *FNum		// Candidate signal number
       );


//#mb 
//void modvir (double sinal, double cosal, double sindel, double cosdel, double sphir, double cphir, double *a, double *b, int Np, Ampl_mod_coeff *amod, Aux_arrays *aux);

void modvir (double sinal, double cosal, double sindel, double cosdel,	
        int Np, Detector_settings *ifo, Aux_arrays *aux);


#endif
