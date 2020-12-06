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
		 Aux_arrays *aux_arr
		 );

void add_signal(
		Search_settings *sett,
		Command_line_opts *opts,
		Aux_arrays *aux_arr,
		Search_range *s_range);

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
	     Aux_arrays *aux
	     );

// Coincidences specific functions 
void handle_opts_coinc(
		       Search_settings *sett,
		       Command_line_opts_coinc *opts,
		       int argc,  
		       char* argv[]);  

void manage_grid_matrix_old( Search_settings *sett,
			     Command_line_opts_coinc *opts);

void manage_grid_matrix( Search_settings *sett, char *gridfile );

void convert_to_linear(
		       Search_settings *sett,
		       Command_line_opts_coinc *opts, 
		       Candidate_triggers *trig);

#endif
