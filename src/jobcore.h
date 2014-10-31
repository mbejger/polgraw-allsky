#ifndef __JOBCORE_H__
#define __JOBCORE_H__

#include "struct.h"


void search(
				Search_settings *sett,
				Command_line_opts *opts,
				Search_range *s_range,
				FFTW_arrays *fftw_arr,
				Signals *sig,
				FFTW_plans *plans,
				Aux_arrays *aux,
				Ampl_mod_coeff *amod,
				int *Fnum,
				double *F
				);


double* job_core(
			int pm,			// hemisphere
			int mm,			// grid 'sky position'
			int nn,			// other grid 'sky position'
			Search_settings *sett, // search settings
			Command_line_opts *opts, // cmd opts
			Search_range *s_range,	 // range for searching
			Signals *sig,
			FFTW_arrays *fftw_arr,   // arrays for fftw
			FFTW_plans *plans,       // plans for fftw
			Aux_arrays *aux, 			 // auxiliary arrays
			Ampl_mod_coeff *amod,    //amplitude modulation functions coefficients 
			double *F,		// F-statistics array
			int *sgnlc,		// reference to array with the parameters
								// of the candidate signal
								// (used below to write to the file)
			int *FNum		// Candidate signal number
       );


void
modvir (double sinal, double cosal, double sindel, double cosdel,	\
        double sphir, double cphir, double *a, double *b, int Np, Ampl_mod_coeff *amod,
        Aux_arrays *aux);

#endif
