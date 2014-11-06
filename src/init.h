#ifndef __INIT_H__
#define __INIT_H__

#include "struct.h"

void handle_opts(
	Search_settings *sett,
	Command_line_opts *opts,
  int argc,  
	char* argv[]);  

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
	Search_settings *sett, 
	Command_line_opts *opts, 
  FFTW_plans *plans, 
	FFTW_arrays *fftw_arr, 
	Aux_arrays *aux_arr);

void read_checkpoints(
	Command_line_opts *opts, 
  Search_range *s_range,
  int *Fnum);

void cleanup(
	Search_settings *sett,
	Command_line_opts *opts,
	Search_range *s_range,
	FFTW_plans *plans,
	FFTW_arrays *fftw_arr,
	Aux_arrays *aux,
	double *F);

#endif
