#ifndef __INIT_H__
#define __INIT_H__

#include "struct.h"

void handle_opts(
	Detector_settings *sett,
	Command_line_opts *opts, 
	int argc, 
	char* argv[]); 

void init_arrays(
	Detector_settings *sett, 
	Command_line_opts *opts, 
	Signals *sig, 
	Aux_arrays *aux_arr, 
	double** F); 

void set_search_range(
	Detector_settings *sett, 
  Command_line_opts *opts, 
	Search_range *s_range); 

void read_grid(
	Detector_settings *sett, 
	Command_line_opts *opts);

void plan_fftw(
	Detector_settings *sett, 
	Command_line_opts *opts, 
	FFTW_plans *plans, 
	FFTW_arrays *fftw_arr, 
	Signals *sig, 
	Aux_arrays *aux_arr); 
 
void read_checkpoints(
	Command_line_opts *opts, 
	Search_range *s_range, 
	int *Fnum); 

void cleanup(Detector_settings *sett,
				Command_line_opts *opts,
				Search_range *s_range,
  			Signals *sig,
				FFTW_plans *plans,
				FFTW_arrays *fftw_arr,
				Aux_arrays *aux,
				Ampl_mod_coeff *amod,
				double *F);

#endif
