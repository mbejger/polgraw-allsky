// C behavioral defines
//
// ISO: request safe versions of functions
#define __STDC_WANT_LIB_EXT1__ 1

// Polgraw includes
#include <floats.h>
#include <auxi.h>
#include <settings.h>
#include <struct.h>
#include <jobcore.h>
#include <init.h>

// Posix includes
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#ifdef WIN32
#include <direct.h>
#include <posix/dirent.h>
#include <posix/getopt.h>
#else
#include <unistd.h>
#include <dirent.h>
#include <getopt.h>
#endif // WIN32

// Standard C includes
#include <stdio.h>      // fopen_s
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <errno.h>
#include <time.h>

// Default output and data directories
#ifndef PREFIX
#define PREFIX ./candidates
#endif

#ifndef DTAPREFIX
#define DTAPREFIX .
#endif

Detector_settings ifo[MAX_DETECTORS];

int main (int argc, char* argv[])
{
    Command_line_opts opts;     // User command line options
    OpenCL_settings cl_sett;    // User settings of OpenCL resource usage
    OpenCL_handles cl_handles;  // OpenCL related resource handles
    Search_settings sett;       // User settings of the search algorithm
    Search_range s_range;       // User provided range information
    Aux_arrays aux_arr;         // Auxiliary arrays used by the application
    BLAS_handles blas_handles;  // Optimized BLAS plans
    FFT_plans fft_plans;        // Optimized FFT plans
    FFT_arrays fft_arr;         // FFT arrays
    cl_mem F_d;                 // F-statistic array
    int i;

    struct stat buffer;         // Output data handling

    // Command line options 
    handle_opts(&sett, &cl_sett, &opts, argc, argv);

    // Initialize OpenCL
    init_opencl(&cl_handles, &cl_sett);

    // Setup output buffer
    setup_output(&buffer, &opts);

    // Grid data 
    read_grid(&sett, &opts);

    // Search settings
    search_settings(&sett);

    // Detector network settings
    detectors_settings(&sett, &opts);

    // Array initialization
    init_arrays(&sett, &cl_handles, &opts, &aux_arr, &F_d);

    // Amplitude modulation functions for each detector
    for (i = 0; i<sett.nifo; i++) rogcvir(&ifo[i]);

    // Set search range from range file  
    set_search_range(&sett, &opts, &s_range);

    // BLAS init
    init_blas(&sett, &cl_handles, &blas_handles);

    // FFT init 
    plan_fft(&sett, &cl_handles, &fft_plans, &fft_arr);

    // Checkpointing
    int Fnum = 0;			        // candidate signal number
    read_checkpoints(&opts, &s_range, &Fnum);

    // main search job
    search(&sett, &opts, &s_range, &cl_handles, &blas_handles, &fft_plans, &fft_arr, &aux_arr, &Fnum, F_d);

    // state file zeroed at the end of the run
    FILE *state;
    if (opts.checkp_flag) {
        errno_t err = fopen_s(&state, opts.qname, "w");
        if (err)
            perror("Error zeroing out state file.");
        fclose(state);
    }

    // Cleanup & memory free 
    cleanup(&sett, &opts, &s_range, &cl_handles, &blas_handles,
        &fft_plans, &fft_arr, &aux_arr, F_d);

    return EXIT_SUCCESS;
}
