#ifndef __INIT_H__
#define __INIT_H__

// Polgraw includes
#include <struct.h>


/// <summary>Command line options handling: search</summary>
///
void handle_opts(Search_settings* sett,
                 OpenCL_settings* cl_sett,
		         Command_line_opts* opts,
		         int argc,  
		         char* argv[]);

/// <summary>Initialize OpenCL devices based on user preference.</summary>
/// <remarks>Currently, only a sinle platform can be selected.</remarks>
///
void init_opencl(OpenCL_handles* cl_handles,
                 OpenCL_settings* cl_sett);

/// <summary>Tries selecting the platform with the specified index.</summary>
///
cl_platform_id select_platform(cl_uint plat_id);

/// <summary>Selects all devices of the specified type.</summary>
///
cl_device_id* select_devices(cl_platform_id platform,
                             cl_device_type dev_type,
                             cl_uint* count);

/// <summary>Create a context that holds all specified devices.</summary>
///
cl_context create_standard_context(cl_device_id* devices,
                                   cl_uint count);

/// <summary>Create a set of command queues to all the devices in the context.</summary>
///
cl_command_queue* create_command_queue_set(cl_context context);

/// <summary>Load kernel file from disk.</summary>
///
char* load_program_file(const char* filename);

/// <summary>Build program file</summary>
///
cl_program build_program_source(cl_context context,
                                const char* source);

/// <summary>Create a kernel for all entry points in the program.</summary>
///
cl_kernel* create_kernels(cl_program program);

/// <summary>Obtain the name of the kernel of a given index.</summary>
///
const char* obtain_kernel_name(cl_uint i);

/// <summary>Obtain kernel with the specified index.</summary>
///
cl_kernel obtain_kernel(cl_program program, cl_uint i);

/// <summary>Generate grid from the M matrix.</summary>
/// <remarks>Processes the file 'grid.bin'</remarks>
///
void read_grid(Search_settings *sett,
               Command_line_opts *opts);

/// <summary>Initialize auxiliary and F-statistics arrays.</summary>
///
void init_arrays(Search_settings* sett,
                 OpenCL_handles* cl_handles,
		         Command_line_opts* opts, 
		         Aux_arrays* aux_arr, 
		         cl_mem* F_d);

void add_signal(
		Search_settings *sett,
		Command_line_opts *opts,
		Aux_arrays *aux_arr,
		Search_range *s_range);

/// <summary>Set search ranges based on user preference.</summary>
///
void set_search_range(Search_settings *sett, 
                      Command_line_opts *opts, 
                      Search_range *s_range);

/// <summary>Sets up BLAS internal states.</summary>
///
void init_blas(Search_settings* sett,
               OpenCL_handles* cl_handles,
               BLAS_handles* blas_handles);

/// <summary>Sets up FFT plans.</summary>
///
void plan_fft(Search_settings* sett,
              OpenCL_handles* cl_handles,
	          FFT_plans* plans, 
	          FFT_arrays* fft_arr);

void read_checkpoints(
		      Command_line_opts *opts, 
		      Search_range *s_range,
		      int *Fnum);

/// <summary>Frees all resources for termination.</summary>
///
void cleanup(Search_settings *sett,
	         Command_line_opts *opts,
	         Search_range *s_range,
             OpenCL_handles* cl_handles,
             BLAS_handles* blas_handles,
	         FFT_plans *plans,
	         FFT_arrays *fft_arr,
	         Aux_arrays *aux,
	         cl_mem F_d);

// Coincidences specific functions 
void handle_opts_coinc(
		       Search_settings *sett,
		       Command_line_opts_coinc *opts,
		       int argc,  
		       char* argv[]);  

void manage_grid_matrix(
			Search_settings *sett,
			Command_line_opts_coinc *opts);

void convert_to_linear(
		       Search_settings *sett,
		       Command_line_opts_coinc *opts, 
		       Candidate_triggers *trig);

int cuinit(int cdev);

#endif // __INIT_H__
