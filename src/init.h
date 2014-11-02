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

/*
//#mb older version
  void init_arrays(
	Signals *sig, 
	Aux_arrays *aux_arr, 
	double** F, 
	Command_line_opts *opts, 
	Search_settings *sett);
*/ 

void init_arrays(
  Search_settings *sett,
  Command_line_opts *opts, 
  Aux_arrays *aux_arr, 
  double** F);


void set_search_range(
	Search_settings *sett, 
	Command_line_opts *opts, 
	Search_range *s_range);  

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
//#mb 
//	FFTW_arrays *fftw_arr,
//	FFTW_plans *plans,
	Aux_arrays *aux,
	double *F);

#endif
