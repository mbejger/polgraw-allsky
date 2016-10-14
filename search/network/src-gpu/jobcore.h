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
void search(Search_settings* sett,
            Command_line_opts* opts,
            Search_range* s_range,
            OpenCL_handles* cl_handles,
            BLAS_handles* blas_handles,
            FFT_plans* plans,
            FFT_arrays* fft_arr,
            Aux_arrays* aux,
            int* Fnum,
            cl_mem F_d);

/// <summary>Main job function.</summary>
/// <remarks>The output is stored in single or double precision. (<c>real_t</c> defined in struct.h)</remarks>
///
real_t* job_core(int pm,                        // hemisphere
                 int mm,                        // grid 'sky position'
                 int nn,                        // other grid 'sky position'
                 Search_settings *sett,         // search settings
                 Command_line_opts *opts,       // cmd opts
                 Search_range *s_range,         // range for searching
                 FFT_plans *plans,              // plans for fftw
                 FFT_arrays *fft_arr,           // arrays for fftw
                 Aux_arrays *aux,               // auxiliary arrays
                 cl_mem F,                      // F-statistics array
                 int *sgnlc,                    // reference to array with the parameters of the candidate signal
                                                // (used below to write to the file)
                 int *FNum,                     // candidate signal number
                 OpenCL_handles* cl_handles,    // handles to OpenCL resources
                 BLAS_handles* blas_handles);   // handle for scaling

/// <summary>Copies amplitude modulation coefficients to constant memory.</summary>
///
void copy_amod_coeff(cl_int nifo,
                     OpenCL_handles* cl_handles,
                     Aux_arrays* aux);

/// <summary>The purpose of this function was undocumented.</summary>
///
void modvir_gpu(real_t sinal,
                real_t cosal,
                real_t sindel,
                real_t cosdel,
                cl_int Np,
                Detector_settings* ifoi,
                OpenCL_handles* cl_handles,
                Aux_arrays* aux,
                cl_int idet);

/// <summary>The purpose of this function was undocumented.</summary>
///
void tshift_pmod_gpu(real_t shft1,
                     real_t het0,
                     real_t ns0,
                     real_t ns1,
                     real_t ns2,
                     real_t* xDat_d,
                     complex_t* xa_d,
                     complex_t* xb_d,
                     real_t* shft_d,
                     real_t* shftf_d,
                     real_t* tshift_d,
                     real_t* aa_d,
                     real_t* bb_d,
                     real_t* DetSSB_d,
                     real_t oms,
                     cl_int N,
                     cl_int nfft,
                     cl_int interpftpad);

/// <summary>Shifts frequencies and remove those over Nyquist.</summary>
///
void resample_postfft_gpu(cl_mem xa_d,
                          cl_mem xb_d,
                          cl_int nfft,
                          cl_int Ninterp,
                          cl_int nyqst,
                          OpenCL_handles* cl_handles);

/// <summary>Scales vectors with a constant.</summary>
///
void blas_scale(cl_mem xa_d,
                cl_mem xa_b,
                cl_uint n,
                real_t a,
                OpenCL_handles* cl_handles,
                BLAS_handles* blas_handles);

/// <summary>Calculates the inner product of both <c>x</c> and <c>y</c>.</summary>
/// <remarks>The function allocates an array of 2 and gives ownership to the caller.</remarks>
/// <remarks>Consider making the temporaries persistent, either providing them via function params or give static storage duration.</remarks>
///
real_t* blas_dot(cl_mem x,
                 cl_mem y,
                 cl_uint n,
                 OpenCL_handles* cl_handles,
                 BLAS_handles* blas_handles);

/// <summary>The purpose of this function was undocumented.</summary>
///
void phase_mod_1_gpu(cl_mem xa,
                     cl_mem xb,
                     cl_mem xar,
                     cl_mem xbr,
                     real_t het1,
                     real_t sgnlt1,
                     cl_mem shft,
                     cl_int N,
                     OpenCL_handles* cl_handles);

/// <summary>The purpose of this function was undocumented.</summary>
///
void phase_mod_2_gpu(cl_mem xa,
                     cl_mem xb,
                     cl_mem xar,
                     cl_mem xbr,
                     real_t het1,
                     real_t sgnlt1,
                     cl_mem shft,
                     cl_int N,
                     OpenCL_handles* cl_handles);

/// <summary>Compute F-statistics.</summary>
/// 
void compute_Fstat_gpu(cl_mem xa,
                       cl_mem xb,
                       cl_mem F,
                       cl_mem maa_d,
                       cl_mem mbb_d,
                       cl_int nmin,
                       cl_int nmax,
                       OpenCL_handles* cl_handles);

/// <summary>Compute F-statistics.</summary>
///
void FStat_gpu_simple(cl_mem F_d,
                      cl_int nfft,
                      cl_int nav,
                      OpenCL_handles* cl_handles);

/// <summary>Saves the designated array into a file with the specified name.</summary>
///
void save_array(HOST_COMPLEX_TYPE *arr, int N, const char* file);


double FStat (double *, int, int, int);
void FStat_gpu(FLOAT_TYPE *F_d, int N, int nav, FLOAT_TYPE *mu_d, FLOAT_TYPE *mu_t_d);

#endif
