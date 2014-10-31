#ifndef __INIT_H__
#define __INIT_H__

#include "struct.h"

/*
char **detector_network(
	Command_line_opts *opts);
*/ 

void handle_opts(
	int argc, 
	char* argv[], 
	Command_line_opts *opts, 
	Search_settings *sett);

void init_arrays(
	Signals *sig, 
	Aux_arrays *aux_arr, 
	double** F, 
	Command_line_opts *opts, 
	Search_settings *sett);

void set_search_range(
	Search_range *s_range, 
	Command_line_opts *opts, 
	Search_settings *sett);

void read_grid(
	Search_settings *sett, 
	Command_line_opts *opts);

void plan_fftw(
	FFTW_plans *plans, 
	FFTW_arrays *fftw_arr, 
	Signals *sig, 
	Aux_arrays *aux_arr, 
	Search_settings *sett, 
	Command_line_opts *opts);

void read_checkpoints(
	Search_range *s_range, 
	int *Fnum, 
	Command_line_opts *opts);

void cleanup(
	Search_settings *sett,
	Command_line_opts *opts,
	Search_range *s_range,
	FFTW_arrays *fftw_arr,
	Signals *sig,
	FFTW_plans *plans,
	Aux_arrays *aux,
	Ampl_mod_coeff *amod,
	double *F);

#endif
