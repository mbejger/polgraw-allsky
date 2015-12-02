
#ifndef __INIT_H__
#define __INIT_H__

#include "struct.h"
#include "floats.h"


//command line options handling
void handle_opts(int argc, char* argv[], Command_line_opts *opts, Detector_settings *sett);

//array creation
void init_arrays(Arrays *arr, FLOAT_TYPE** cu_F,  Command_line_opts *opts, Detector_settings *sett);

//determination of search range 
void set_search_range(Search_range *s_range, Command_line_opts *opts, Detector_settings *sett);

//grid reading
void read_grid(Detector_settings *sett, Command_line_opts *opts);

//planning FFTs
void plan_fft(FFT_plans *plans, Arrays *arr, Detector_settings *sett, Command_line_opts *opts);

//reading checkpoint status
void read_checkpoints(Search_range *s_range, int *Fnum, Command_line_opts *opts);

//array deallocation
void cleanup(Detector_settings *sett,
				Command_line_opts *opts,
				Search_range *s_range,
				Arrays *arr,
				FFT_plans *plans,
				Ampl_mod_coeff *amod,
				FLOAT_TYPE *cu_F);


#endif
