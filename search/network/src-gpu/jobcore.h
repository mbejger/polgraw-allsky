#ifndef __JOBCORE_H__
#define __JOBCORE_H__

// Polgraw includes
#include <struct.h>     // Search_settings, Command_line_opts, Search_range, FFT_plans, FFT_arrays, Aux_arrays
#include <floats.h>     // FLOAT_TYPE, HOST_COMPLEX_TYPE

// clBLAS includes
#include <clBLAS.h>


#define BLOCK_SIZE 256
#define BLOCK_SIZE_RED 128
#define BLOCK_DIM(n, b) ((n)/b + ((n)%b==0 ? 0 : 1))
#define NAV_THREADS 16

/// <summary>Main searching function.</summary>
/// <remarks>This function loops over hemispheres, sky positions and spin-downs.</remarks>
///
void search(Search_settings *sett,
            Command_line_opts *opts,
            Search_range *s_range,
            FFT_plans *plans,
            FFT_arrays *fft_arr,
            Aux_arrays *aux,
            int *Fnum,
            double *F);

/// <summary>Main job function.</summary>
/// <remarks>The output is stored in single or double precision. (<c>FLOAT_TYPE</c> defined in struct.h)</remarks>
///
FLOAT_TYPE *job_core(int pm,                   // hemisphere
                     int mm,                   // grid 'sky position'
                     int nn,                   // other grid 'sky position'
                     Search_settings *sett,    // search settings
                     Command_line_opts *opts,  // cmd opts
                     Search_range *s_range,    // range for searching
                     FFT_plans *plans,         // plans for fftw
                     FFT_arrays *fft_arr,      // arrays for fftw
                     Aux_arrays *aux,          // auxiliary arrays
                     double *F,                // F-statistics array
                     int *sgnlc,               // reference to array with the parameters of the candidate signal
                                               // (used below to write to the file)
                     int *FNum);               // candidate signal number
                     //cublasHandle_t scale    // handle for scaling

/// <summary>Saves the designated array into a file with the specified name.</summary>
///
void save_array(HOST_COMPLEX_TYPE *arr, int N, const char* file);

void FStat_gpu_simple(FLOAT_TYPE *F_d, int nfft, int nav);
double FStat (double *, int, int, int);
void FStat_gpu(FLOAT_TYPE *F_d, int N, int nav, FLOAT_TYPE *mu_d, FLOAT_TYPE *mu_t_d);

#endif
