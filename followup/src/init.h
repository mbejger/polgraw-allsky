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
  Aux_arrays *aux_arr);

void read_grid(
  Search_settings *sett, 
  Command_line_opts *opts);

unsigned long int random_seed();

void gauss_xdat(
  Search_settings *sett,
  double amplitude,
  double sigma,
  int i);

void add_signal(
  Search_settings *sett,
  Command_line_opts *opts,
  Aux_arrays *aux_arr);

void read_grid(
	Search_settings *sett, 
	Command_line_opts *opts);

void cleanup_followup(
	Search_settings *sett,
	Command_line_opts *opts,
	Aux_arrays *aux);

#endif
