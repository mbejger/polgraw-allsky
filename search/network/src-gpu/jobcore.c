// C behavioral defines
//
// MSVC: macro to include constants, such as M_PI (include before math.h)
#define _USE_MATH_DEFINES
// ISO: request safe versions of functions
#define __STDC_WANT_LIB_EXT1__ 1
// GCC: hope this macro is not actually needed
//#define _GNU_SOURCE

// Polgraw includes
#include <CL/util.h>        // checkErr
#include <struct.h>
#include <jobcore.h>
#include <auxi.h>
#include <settings.h>
#include <timer.h>
#include <spline_z.h>
#include <floats.h>

// clFFT includes
#include <clFFT.h>

// OpenCL includes
#include <CL/cl.h>

// Posix includes
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#ifdef WIN32
#include <io.h>             // _chsize_s
#include <direct.h>
#include <posix/dirent.h>
#include <posix/getopt.h>
#else
#include <unistd.h>         // ftruncate
#include <dirent.h>
#include <getopt.h>
#endif // WIN32

// Standard C includes
#include <math.h>
#include <stdio.h>          // fopen/fclose, fprintf
#include <malloc.h>
#include <complex.h>
#include <string.h>         // memcpy_s
#include <errno.h>          // errno_t
#include <stdlib.h>         // EXIT_FAILURE


void save_array(HOST_COMPLEX_TYPE *arr, int N, const char* file)
{
    int i;
    FILE *fc = fopen(file, "w");
    for (i=0; i<N; i++)
        fprintf(fc, "%d %e + i %e\n", i, creal(arr[i]), cimag(arr[i]));
    fclose(fc);
}

void save_array_double(double *arr, int N, const char* file) {
  int i;
  FILE *fc = fopen(file, "w");
  for (i=0; i<N; i++) {
    fprintf(fc, "%d %e\n", i, arr[i]);
  }
  fclose(fc);
}


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
            cl_mem F_d)
{
    cl_int CL_err = CL_SUCCESS;
    int pm, mm, nn;    // hemisphere, sky positions 
    int sgnlc;         // number of candidates
    FLOAT_TYPE *sgnlv;    // array with candidates data

    char outname[512];
#ifdef WIN32
    int low_state;
#endif // WIN32
    FILE *state;

#if TIMERS>0
    struct timespec tstart = get_current_time(), tend;
#endif

    // Copy amod coefficients to device
    copy_amod_coeff(sett->nifo, cl_handles, aux);

    int cand_buffer_count = 0;

    //allocate vector for FStat_gpu
    int nav_blocks = (sett->nmax - sett->nmin) / NAV;     //number of nav-blocks
    int blocks = (sett->nmax - sett->nmin) / BLOCK_SIZE;  //number of blocks for Fstat-smoothing

    aux->mu_t_d = clCreateBuffer(cl_handles->ctx, CL_MEM_READ_WRITE, blocks * sizeof(real_t), NULL, &CL_err);
    aux->mu_d = clCreateBuffer(cl_handles->ctx, CL_MEM_READ_WRITE, nav_blocks * sizeof(real_t), NULL, &CL_err);

    state = NULL;
    if (opts->checkp_flag)
#ifdef WIN32
    {
        _sopen_s(&low_state, opts->qname,
            _O_RDWR | _O_CREAT,   // Allowed operations
            _SH_DENYNO,           // Allowed sharing
            _S_IREAD | _S_IWRITE);// Permission settings

        state = _fdopen(low_state, "w");
    }
#else
        state = fopen(opts->qname, "w");
#endif // WIN32

    // Loop over hemispheres //
    for (pm = s_range->pst; pm <= s_range->pmr[1]; ++pm)
    {
        sprintf(outname, "%s/triggers_%03d_%03d%s_%d.bin",
                opts->prefix,
                opts->ident,
                opts->band,
                opts->label,
                pm);

        // Two main loops over sky positions //
        for (mm = s_range->mst; mm <= s_range->mr[1]; ++mm)
        {
            for (nn = s_range->nst; nn <= s_range->nr[1]; ++nn)
            {
                if (opts->checkp_flag)
                {
#ifdef WIN32
                    if (_chsize(low_state, 0))
                    {
                        printf("Failed to resize file");
                        exit(EXIT_FAILURE);
                    }
#else
                    if (ftruncate(fileno(state), 0))
                    {
                        printf("Failed to resize file");
                        exit(EXIT_FAILURE);
                    }
#endif // WIN32
                    fprintf(state, "%d %d %d %d %d\n", pm, mm, nn, s_range->sst, *Fnum);
                    fseek(state, 0, SEEK_SET);
                }

                // Loop over spindowns is inside job_core() //
                sgnlv = job_core(pm,            // hemisphere
                                 mm,            // grid 'sky position'
                                 nn,            // other grid 'sky position'
                                 sett,          // search settings
                                 opts,          // cmd opts
                                 s_range,       // range for searching
                                 plans,         // fftw plans 
                                 fft_arr,       // arrays for fftw
                                 aux,           // auxiliary arrays
                                 F_d,           // F-statistics array
                                 &sgnlc,        // reference to array with the parameters
                                                // of the candidate signal
                                                // (used below to write to the file)
                                 Fnum,          // Candidate signal number
                                 cl_handles,    // handles to OpenCL resources
                                 blas_handles   // handle for scaling
                );

                // Get back to regular spin-down range
                s_range->sst = s_range->spndr[0];

                // Add trigger parameters to a file //

                // if any signals found (Fstat>Fc)
                if (sgnlc)
                {
                    FILE* fc = fopen(outname, "w");
                    if (fc == NULL) perror("Failed to open output file.");

                    size_t count = fwrite((void *)(sgnlv), sizeof(FLOAT_TYPE), sgnlc*NPAR, fc);
                    if (count < sgnlc*NPAR) perror("Failed to write output file.");

                    int close = fclose(fc);
                    if (close == EOF) perror("Failed to close output file.");

                } // if sgnlc
                free(sgnlv);
            } // for nn
            s_range->nst = s_range->nr[0];
        } // for mm
        s_range->mst = s_range->mr[0];
    } // for pm

    if (opts->checkp_flag)
        fclose(state);

    // cublasDestroy(scale);

#if TIMERS>0
    tend = get_current_time();
    // printf("tstart = %d . %d\ntend = %d . %d\n", tstart.tv_sec, tstart.tv_usec, tend.tv_sec, tend.tv_usec);
    double time_elapsed = get_time_difference(tstart, tend);
    printf("Time elapsed: %e s\n", time_elapsed);
#endif

} // end of search

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
                 cl_mem F_d,                    // F-statistics array
                 int *sgnlc,                    // reference to array with the parameters of the candidate signal
                                                // (used below to write to the file)
                 int *FNum,                     // candidate signal number
                 OpenCL_handles* cl_handles,    // handles to OpenCL resources
                 BLAS_handles* blas_handles)    // handle for scaling
{
    int smin = s_range->sst, smax = s_range->spndr[1];
    real_t al1, al2, sinalt, cosalt, sindelt, cosdelt, sgnlt[NPAR], nSource[3], het0, sgnl0, ft;

    // VLA version
    //double _tmp1[sett->nifo][sett->N];

    // Non-VLA version
    real_t** _tmp1;
    _tmp1 = (real_t**)malloc(sett->nifo * sizeof(real_t*));
    for (int x = 0; x < sett->nifo; ++x)
        _tmp1[x] = (real_t*)malloc(sett->N * sizeof(real_t));

    real_t* sgnlv;

    // Stateful function (local variable with static storage duration)
    static real_t *F;
    if (F == NULL) F = (real_t*)malloc(2 * sett->nfft * sizeof(real_t));

    /* Matrix	M(.,.) (defined on page 22 of PolGrawCWAllSkyReview1.pdf file)
    defines the transformation form integers (bin, ss, nn, mm) determining
    a grid point to linear coordinates omega, omegadot, alpha_1, alpha_2),
    where bin is the frequency bin number and alpha_1 and alpha_2 are
    defined on p. 22 of PolGrawCWAllSkyReview1.pdf file.

    [omega]                          [bin]
    [omegadot]       = M(.,.) \times [ss]
    [alpha_1/omega]                  [nn]
    [alpha_2/omega]                  [mm]

    Array M[.] is related to matrix M(.,.) in the following way;

    [ M[0] M[4] M[8]  M[12] ]
    M(.,.) =   [ M[1] M[5] M[9]  M[13] ]
    [ M[2] M[6] M[10] M[14] ]
    [ M[3] M[7] M[11] M[15] ]

    and

    M[1] = M[2] = M[3] = M[6] = M[7] = 0
    */

    // Grid positions
    al1 = nn*sett->M[10] + mm*sett->M[14];
    al2 = nn*sett->M[11] + mm*sett->M[15];

    sgnlv = NULL;
    *sgnlc = 0;

    // check if the search is in an appropriate region of the grid
    // if not, returns NULL
    if ((sqr(al1) + sqr(al2)) / sqr(sett->oms) > 1.) return NULL;

    int ss;
    real_t shft1/*, phase, cp, sp*/; // Unused variables
    // complex_t exph; // Unused variable

    // Change linear (grid) coordinates to real coordinates
    lin2ast(al1 / sett->oms, al2 / sett->oms,
            pm, sett->sepsm, sett->cepsm,
            &sinalt, &cosalt, &sindelt, &cosdelt);

    // calculate declination and right ascention
    // written in file as candidate signal sky positions
    sgnlt[2] = asin(sindelt);
    sgnlt[3] = fmod(atan2(sinalt, cosalt) + 2.*M_PI, 2.*M_PI);

    het0 = fmod(nn*sett->M[8] + mm*sett->M[12], sett->M[0]);

    // Nyquist frequency 
    int nyqst = (sett->nfft) / 2 + 1;

    // Loop for each detector 
    for (int n = 0; n<sett->nifo; ++n) {

        /* Amplitude modulation functions aa and bb
        * for each detector (in signal sub-struct
        * of _detector, ifo[n].sig.aa, ifo[n].sig.bb)
        */
        modvir_gpu(sinalt, cosalt, sindelt, cosdelt, sett->N, &ifo[n], cl_handles, aux, n);

        // Calculate detector positions with respect to baricenter
        nSource[0] = cosalt*cosdelt;
        nSource[1] = sinalt*cosdelt;
        nSource[2] = sindelt;

        shft1 = nSource[0] * ifo[n].sig.DetSSB[0] +
                nSource[1] * ifo[n].sig.DetSSB[1] +
                nSource[2] * ifo[n].sig.DetSSB[2];

        tshift_pmod_gpu(shft1, het0, nSource[0], nSource[1], nSource[2],
                        ifo[n].sig.xDat_d, fft_arr->xa_d, fft_arr->xb_d,
                        ifo[n].sig.shft_d, ifo[n].sig.shftf_d,
                        aux->tshift_d,
                        ifo[n].sig.aa_d, ifo[n].sig.bb_d,
                        ifo[n].sig.DetSSB_d,
                        sett->oms, sett->N, sett->nfft, sett->interpftpad, cl_handles);

        clfftStatus CLFFT_status = CLFFT_SUCCESS;
        cl_event fft_exec[2];
        CLFFT_status = clfftEnqueueTransform(plans->plan /*pl_int*/, CLFFT_FORWARD, 1, cl_handles->exec_queues, 0, NULL, &fft_exec[0], &fft_arr->xa_d, NULL, NULL /*May be slow, consider using tmp_buffer*/);
        checkErrFFT(CLFFT_status, "clfftEnqueueTransform(CLFFT_FORWARD)");
        CLFFT_status = clfftEnqueueTransform(plans->plan /*pl_int*/, CLFFT_FORWARD, 1, cl_handles->exec_queues, 0, NULL, &fft_exec[1], &fft_arr->xb_d, NULL, NULL /*May be slow, consider using tmp_buffer*/);
        checkErrFFT(CLFFT_status, "clfftEnqueueTransform(CLFFT_FORWARD)");

        clWaitForEvents(2, fft_exec);

        resample_postfft_gpu(fft_arr->xa_d,
                             fft_arr->xb_d,
                             sett->nfft,
                             sett->Ninterp,
                             nyqst,
                             cl_handles);

        // Backward fft (len Ninterp = nfft*interpftpad)
        clfftEnqueueTransform(plans->pl_inv, CLFFT_BACKWARD, 1, cl_handles->exec_queues, 0, NULL, &fft_exec[0], &fft_arr->xa_d, NULL, NULL /*May be slow, consider using tmp_buffer*/);
        checkErrFFT(CLFFT_status, "clfftEnqueueTransform(CLFFT_BACKWARD)");
        clfftEnqueueTransform(plans->pl_inv, CLFFT_BACKWARD, 1, cl_handles->exec_queues, 0, NULL, &fft_exec[1], &fft_arr->xb_d, NULL, NULL /*May be slow, consider using tmp_buffer*/);
        checkErrFFT(CLFFT_status, "clfftEnqueueTransform(CLFFT_BACKWARD)");

        clWaitForEvents(2, fft_exec);
        for (size_t i = 0; i < 2; ++i) clReleaseEvent(fft_exec[i]);

        //scale fft with cublas
        ft = (double)sett->interpftpad / sett->Ninterp;
        blas_scale(fft_arr->xa_d,
                   fft_arr->xb_d,
                   sett->Ninterp,
                   ft,
                   cl_handles,
                   blas_handles);

        // Spline interpolation to xDatma, xDatmb arrays
        gpu_interp(fft_arr->xa_d,       //input data
                   sett->Ninterp,       //input data length
                   aux->tshift_d,       //output time domain
                   ifo[n].sig.xDatma_d, //output values
                   sett->N,             //output data length
                   aux->diag_d,         //diagonal
                   aux->ldiag_d,        //lower diagonal
                   aux->udiag_d,        //upper diagonal
                   aux->B_d,            //coefficient matrix
                   cl_handles);

        gpu_interp(fft_arr->xb_d,       //input data
                   sett->Ninterp,       //input data length
                   aux->tshift_d,       //output time domain
                   ifo[n].sig.xDatmb_d, //output values
                   sett->N,             //output data length
                   aux->diag_d,         //diagonal
                   aux->ldiag_d,        //lower diagonal
                   aux->udiag_d,        //upper diagonal
                   aux->B_d,            //coefficient matrix
                   cl_handles);

        ft = 1. / ifo[n].sig.sig2;

        blas_scale(ifo[n].sig.xDatma_d,
                   ifo[n].sig.xDatmb_d,
                   sett->N,
                   ft,
                   cl_handles,
                   blas_handles);
    } // end of detector loop 

    real_t _maa = 0;
    real_t _mbb = 0;

    for (int n = 0; n<sett->nifo; ++n)
    {
        real_t* temp = blas_dot(ifo[n].sig.aa_d, ifo[n].sig.bb_d, sett->N, cl_handles, blas_handles);

        _maa += temp[0] / ifo[n].sig.sig2;
        _mbb += temp[1] / ifo[n].sig.sig2;

        free(temp);
    }

    // Copy sums to constant memory
    {
        cl_event write_event[2];
        clEnqueueWriteBuffer(cl_handles->write_queues[0], aux->maa_d, CL_FALSE, 0, sizeof(real_t), &_maa, 0, NULL, &write_event[0]);
        clEnqueueWriteBuffer(cl_handles->write_queues[0], aux->mbb_d, CL_FALSE, 0, sizeof(real_t), &_mbb, 0, NULL, &write_event[1]);

        clWaitForEvents(2, write_event);
        for (size_t i = 0; i < 2; ++i) clReleaseEvent(write_event[i]);
    }

    // Spindown loop //

#if TIMERS>2
    struct timespec tstart, tend;
    double spindown_timer = 0;
    int spindown_counter = 0;
#endif

    // Check if the signal is added to the data 
    // or the range file is given:  
    // if not, proceed with the wide range of spindowns 
    // if yes, use smin = s_range->sst, smax = s_range->spndr[1]  
    if (!strcmp(opts->addsig, "") && !strcmp(opts->range, "")) {

        // Spindown range defined using Smin and Smax (settings.c)  
        smin = (int)trunc((sett->Smin - nn*sett->M[9] - mm*sett->M[13]) / sett->M[5]);  // Cast is intentional and safe (silences warning).
        smax = (int)trunc(-(nn*sett->M[9] + mm*sett->M[13] + sett->Smax) / sett->M[5]); // Cast is intentional and safe (silences warning).
    }

    printf("\n>>%d\t%d\t%d\t[%d..%d]\n", *FNum, mm, nn, smin, smax);

    // No-spindown calculations
    if (opts->s0_flag) smin = smax;

    // if spindown parameter is taken into account, smin != smax
    for (ss = smin; ss <= smax; ++ss)
    {

#if TIMERS>2
        tstart = get_current_time();
#endif 

        // Spindown parameter
        sgnlt[1] = ss*sett->M[5] + nn*sett->M[9] + mm*sett->M[13];

        //    // Spindown range
        //    if(sgnlt[1] >= -sett->Smax && sgnlt[1] <= sett->Smax) { 

        int ii;
        real_t Fc, het1;

#ifdef VERBOSE
        //print a 'dot' every new spindown
        printf("."); fflush(stdout);
#endif 

        het1 = fmod(ss*sett->M[4], sett->M[0]);
        if (het1<0) het1 += sett->M[0];

        sgnl0 = het0 + het1;
        // printf("%d  %d\n", BLOCK_SIZE, (sett->N + BLOCK_SIZE - 1)/BLOCK_SIZE );

        phase_mod_1_gpu(fft_arr->xa_d,
                        fft_arr->xb_d,
                        ifo[0].sig.xDatma_d,
                        ifo[0].sig.xDatmb_d,
                        het1,
                        sgnlt[1],
                        ifo[0].sig.shft_d,
                        sett->N,
                        cl_handles);

        for (int n = 1; n<sett->nifo; ++n)
        {
            phase_mod_2_gpu(fft_arr->xa_d,
                            fft_arr->xb_d,
                            ifo[n].sig.xDatma_d,
                            ifo[n].sig.xDatmb_d,
                            het1,
                            sgnlt[1],
                            ifo[n].sig.shft_d,
                            sett->N,
                            cl_handles);
        }

        // initialize arrays to 0. with integer 0
        // assuming double , remember to change when switching to float
        {
            cl_int CL_err = CL_SUCCESS;
            complex_t pattern = {0, 0};
            cl_event fill_event[2];

            // Zero pad from offset until the end
            CL_err = clEnqueueFillBuffer(cl_handles->write_queues[0], fft_arr->xa_d, &pattern, sizeof(complex_t), sett->N * sizeof(complex_t), (sett->nfftf - sett->N) * 2 * sizeof(complex_t), 0, NULL, &fill_event[0]);
            checkErr(CL_err, "clEnqueueFillBuffer");
            CL_err = clEnqueueFillBuffer(cl_handles->write_queues[0], fft_arr->xb_d, &pattern, sizeof(complex_t), sett->N * sizeof(complex_t), (sett->nfftf - sett->N) * 2 * sizeof(complex_t), 0, NULL, &fill_event[1]);
            checkErr(CL_err, "clEnqueueFillBuffer");

            clWaitForEvents(2, fill_event);

            clReleaseEvent(fill_event[0]);
            clReleaseEvent(fill_event[1]);
        }

        // fft length fftpad*nfft
        {
            clfftStatus CLFFT_status = CLFFT_SUCCESS;
            cl_event fft_exec[2];
            clfftEnqueueTransform(plans->plan, CLFFT_FORWARD, 1, cl_handles->exec_queues, 0, NULL, &fft_exec[0], &fft_arr->xa_d, NULL, NULL /*May be slow, consider using tmp_buffer*/);
            checkErrFFT(CLFFT_status, "clfftEnqueueTransform(CLFFT_FORWARD)");
            clfftEnqueueTransform(plans->plan, CLFFT_FORWARD, 1, cl_handles->exec_queues, 0, NULL, &fft_exec[1], &fft_arr->xb_d, NULL, NULL /*May be slow, consider using tmp_buffer*/);
            checkErrFFT(CLFFT_status, "clfftEnqueueTransform(CLFFT_FORWARD)");

            clWaitForEvents(2, fft_exec);
        }

        (*FNum)++;

        compute_Fstat_gpu(fft_arr->xa_d,
                          fft_arr->xb_d,
                          F_d,
                          aux->maa_d,
                          aux->mbb_d,
                          sett->nmin,
                          sett->nmax,
                          cl_handles);

#ifdef GPUFSTAT
        if (!(opts->white_flag))  // if the noise is not white noise
            FStat_gpu(F_d + sett->nmin, sett->nmax - sett->nmin, NAV, aux->mu_d, aux->mu_t_d);

#else
        if (!(opts->white_flag))  // if the noise is not white noise
            FStat_gpu_simple(F_d, sett->nfft, NAVFSTAT, cl_handles);
#endif

        cl_int CL_err = CL_SUCCESS;

        CL_err = clEnqueueReadBuffer(cl_handles->read_queues[0], F_d, CL_TRUE, 0, 2 * sett->nfft * sizeof(real_t), F, 0, NULL, NULL);

        /*
        FILE *f1 = fopen("fstat-gpu.dat", "w");
        for(i=sett->nmin; i<sett->nmax; i++)
        fprintf(f1, "%d   %lf\n", i, F[i]);

        fclose(f1);
        printf("wrote fstat-gpu.dat | ss=%d  \n", ss);
        //exit(EXIT_SUCCESS);
        */

        for (int i = sett->nmin; i<sett->nmax; i++) {
            if ((Fc = F[i]) > opts->trl) { // if F-stat exceeds trl (critical value)
                                           // Find local maximum for neighboring signals 
                ii = i;

                while (++i < sett->nmax && F[i] > opts->trl) {
                    if (F[i] >= Fc) {
                        ii = i;
                        Fc = F[i];
                    } // if F[i] 
                } // while i 

                  // Candidate signal frequency
                sgnlt[0] = 2.*M_PI*ii / ((double)sett->fftpad*sett->nfft) + sgnl0;
                // Signal-to-noise ratio
                sgnlt[4] = sqrt(2.*(Fc - sett->nd));

                (*sgnlc)++; // increase found number

                            // Add new parameters to output array 
                sgnlv = (FLOAT_TYPE *)realloc(sgnlv, NPAR*(*sgnlc) * sizeof(FLOAT_TYPE));

                for (int j = 0; j<NPAR; ++j) // save new parameters
                    sgnlv[NPAR*(*sgnlc - 1) + j] = (FLOAT_TYPE)sgnlt[j];

#ifdef VERBOSE
                printf("\nSignal %d: %d %d %d %d %d \tsnr=%.2f\n",
                    *sgnlc, pm, mm, nn, ss, ii, sgnlt[4]);
#endif 

            } // if Fc > trl 

        } // for i


#if TIMERS>2
        tend = get_current_time();
        spindown_timer += get_time_difference(tstart, tend);
        spindown_counter++;
#endif


        //    } // if sgnlt[1] 

    } // for ss 


#ifndef VERBOSE
    printf("Number of signals found: %d\n", *sgnlc);
#endif 

#if TIMERS>2
    printf("\nTotal spindown loop time: %e s, mean spindown time: %e s (%d runs)\n",
        spindown_timer, spindown_timer / spindown_counter, spindown_counter);
#endif

    // Non-VLA free _tmp1
    for (int x = 0; x < sett->nifo; ++x)
        free(_tmp1[x]);
    free(_tmp1);

    return sgnlv;

} // jobcore

/// <summary>Copies amplitude modulation coefficients to constant memory.</summary>
///
void copy_amod_coeff(int nifo,
                     OpenCL_handles* cl_handles,
                     Aux_arrays* aux)
{
    cl_int CL_err = CL_SUCCESS;

    Ampl_mod_coeff* tmp = clEnqueueMapBuffer(cl_handles->exec_queues[0],
                                             aux->ifo_amod_d,
                                             CL_TRUE,
                                             CL_MAP_WRITE_INVALIDATE_REGION,
                                             0,
                                             nifo * sizeof(Ampl_mod_coeff),
                                             0,
                                             NULL,
                                             NULL,
                                             &CL_err);

    for (size_t i = 0; i < nifo; ++i) tmp[i] = ifo[i].amod;

    cl_event unmap_event;
    clEnqueueUnmapMemObject(cl_handles->exec_queues[0], aux->ifo_amod_d, tmp, 0, NULL, &unmap_event);

    clWaitForEvents(1, &unmap_event);

    clReleaseEvent(unmap_event);
}

/// <summary>The purpose of this function was undocumented.</summary>
///
void modvir_gpu(real_t sinal,
                real_t cosal,
                real_t sindel,
                real_t cosdel,
                int Np,
                Detector_settings* ifoi,
                OpenCL_handles* cl_handles,
                Aux_arrays* aux,
                int idet)
{
    cl_int CL_err = CL_SUCCESS;
    real_t cosalfr, sinalfr, c2d, c2sd/*, c, s, c2s, cs*/; // Unused variables

    cosalfr = cosal*(ifoi->sig.cphir) + sinal*(ifoi->sig.sphir);
    sinalfr = sinal*(ifoi->sig.cphir) - cosal*(ifoi->sig.sphir);
    c2d = sqr(cosdel);
    c2sd = sindel*cosdel;

    clSetKernelArg(cl_handles->kernels[Modvir], 0, sizeof(cl_mem), &ifoi->sig.aa_d);
    clSetKernelArg(cl_handles->kernels[Modvir], 1, sizeof(cl_mem), &ifoi->sig.bb_d);
    clSetKernelArg(cl_handles->kernels[Modvir], 2, sizeof(real_t), &cosalfr);
    clSetKernelArg(cl_handles->kernels[Modvir], 3, sizeof(real_t), &sinalfr);
    clSetKernelArg(cl_handles->kernels[Modvir], 4, sizeof(real_t), &c2d);
    clSetKernelArg(cl_handles->kernels[Modvir], 5, sizeof(real_t), &c2sd);
    clSetKernelArg(cl_handles->kernels[Modvir], 6, sizeof(cl_mem), &aux->sinmodf_d);
    clSetKernelArg(cl_handles->kernels[Modvir], 7, sizeof(cl_mem), &aux->cosmodf_d);
    clSetKernelArg(cl_handles->kernels[Modvir], 8, sizeof(real_t), &sindel);
    clSetKernelArg(cl_handles->kernels[Modvir], 9, sizeof(real_t), &cosdel);
    clSetKernelArg(cl_handles->kernels[Modvir], 10, sizeof(cl_int), &Np);
    clSetKernelArg(cl_handles->kernels[Modvir], 11, sizeof(cl_int), &idet);
    clSetKernelArg(cl_handles->kernels[Modvir], 12, sizeof(cl_mem), &aux->ifo_amod_d);

    cl_event exec;
    size_t size_Np = (size_t)Np; // Helper variable to make pointer types match. Cast to silence warning

    CL_err = clEnqueueNDRangeKernel(cl_handles->exec_queues[0], cl_handles->kernels[Modvir], 1, NULL, &size_Np, NULL, 0, NULL, &exec);

    clWaitForEvents(1, &exec);

    clReleaseEvent(exec);
}

/// <summary>The purpose of this function was undocumented.</summary>
///
void tshift_pmod_gpu(real_t shft1,
                     real_t het0,
                     real_t ns0,
                     real_t ns1,
                     real_t ns2,
                     cl_mem xDat_d,
                     cl_mem xa_d,
                     cl_mem xb_d,
                     cl_mem shft_d,
                     cl_mem shftf_d,
                     cl_mem tshift_d,
                     cl_mem aa_d,
                     cl_mem bb_d,
                     cl_mem DetSSB_d,
                     real_t oms,
                     cl_int N,
                     cl_int nfft,
                     cl_int interpftpad,
                     OpenCL_handles* cl_handles)
{
    cl_int CL_err = CL_SUCCESS;

    clSetKernelArg(cl_handles->kernels[TShiftPMod], 0, sizeof(real_t), &shft1);
    clSetKernelArg(cl_handles->kernels[TShiftPMod], 1, sizeof(real_t), &het0);
    clSetKernelArg(cl_handles->kernels[TShiftPMod], 2, sizeof(real_t), &ns0);
    clSetKernelArg(cl_handles->kernels[TShiftPMod], 3, sizeof(real_t), &ns1);
    clSetKernelArg(cl_handles->kernels[TShiftPMod], 4, sizeof(real_t), &ns2);
    clSetKernelArg(cl_handles->kernels[TShiftPMod], 5, sizeof(cl_mem), &xDat_d);
    clSetKernelArg(cl_handles->kernels[TShiftPMod], 6, sizeof(cl_mem), &xa_d);
    clSetKernelArg(cl_handles->kernels[TShiftPMod], 7, sizeof(cl_mem), &xb_d);
    clSetKernelArg(cl_handles->kernels[TShiftPMod], 8, sizeof(cl_mem), &shft_d);
    clSetKernelArg(cl_handles->kernels[TShiftPMod], 9, sizeof(cl_mem), &shftf_d);
    clSetKernelArg(cl_handles->kernels[TShiftPMod], 10, sizeof(cl_mem), &tshift_d);
    clSetKernelArg(cl_handles->kernels[TShiftPMod], 11, sizeof(cl_mem), &aa_d);
    clSetKernelArg(cl_handles->kernels[TShiftPMod], 12, sizeof(cl_mem), &bb_d);
    clSetKernelArg(cl_handles->kernels[TShiftPMod], 13, sizeof(cl_mem), &DetSSB_d);
    clSetKernelArg(cl_handles->kernels[TShiftPMod], 14, sizeof(real_t), &oms);
    clSetKernelArg(cl_handles->kernels[TShiftPMod], 15, sizeof(cl_int), &N);
    clSetKernelArg(cl_handles->kernels[TShiftPMod], 16, sizeof(cl_int), &nfft);
    clSetKernelArg(cl_handles->kernels[TShiftPMod], 17, sizeof(cl_int), &interpftpad);

    cl_event exec;
    size_t size_nfft = (size_t)nfft; // Helper variable to make pointer types match. Cast to silence warning

    CL_err = clEnqueueNDRangeKernel(cl_handles->exec_queues[0], cl_handles->kernels[TShiftPMod], 1, NULL, &size_nfft, NULL, 0, NULL, &exec);

    clWaitForEvents(1, &exec);

    clReleaseEvent(exec);
}

/// <summary>Shifts frequencies and remove those over Nyquist.</summary>
///
void resample_postfft_gpu(cl_mem xa_d,
                          cl_mem xb_d,
                          cl_int nfft,
                          cl_int Ninterp,
                          cl_int nyqst,
                          OpenCL_handles* cl_handles)
{
    cl_int CL_err = CL_SUCCESS;

    clSetKernelArg(cl_handles->kernels[ResamplePostFFT], 0, sizeof(cl_mem), &xa_d);
    clSetKernelArg(cl_handles->kernels[ResamplePostFFT], 1, sizeof(cl_mem), &xb_d);
    clSetKernelArg(cl_handles->kernels[ResamplePostFFT], 2, sizeof(real_t), &nfft);
    clSetKernelArg(cl_handles->kernels[ResamplePostFFT], 3, sizeof(real_t), &Ninterp);
    clSetKernelArg(cl_handles->kernels[ResamplePostFFT], 4, sizeof(real_t), &nyqst);

    cl_event exec;
    size_t size_Ninterp = (size_t)Ninterp; // Helper variable to make pointer types match. Cast to silence warning

    CL_err = clEnqueueNDRangeKernel(cl_handles->exec_queues[0], cl_handles->kernels[ResamplePostFFT], 1, NULL, &size_Ninterp, NULL, 0, NULL, &exec);

    clWaitForEvents(1, &exec);

    clReleaseEvent(exec);
}

/// <summary>Scales vectors with a constant.</summary>
///
void blas_scale(cl_mem xa_d,
                cl_mem xb_d,
                cl_uint n,
                real_t a,
                OpenCL_handles* cl_handles,
                BLAS_handles* blas_handles)
{
    clblasStatus status[2];
    cl_event blas_exec[2];
#ifdef COMP_FLOAT
    status[0] = clblasSscal(n, a, xa_d, 0, 1, 1, cl_handles->exec_queues, 0, NULL, &blas_exec[0]);
    status[1] = clblasSscal(n, a, xb_d, 0, 1, 1, cl_handles->exec_queues, 0, NULL, &blas_exec[1]);
#else
    status[0] = clblasDscal(n, a, xa_d, 0, 1, 1, cl_handles->exec_queues, 0, NULL, &blas_exec[0]);
    status[1] = clblasDscal(n, a, xb_d, 0, 1, 1, cl_handles->exec_queues, 0, NULL, &blas_exec[1]);
#endif // COMP_FLOAT

    clWaitForEvents(2, blas_exec);

    for (size_t i = 0; i < 2; ++i) clReleaseEvent(blas_exec[i]);
}

/// <summary>Calculates the inner product of both <c>x</c> and <c>y</c>.</summary>
/// <remarks>The function allocates an array of 2 and gives ownership to the caller.</remarks>
/// <remarks>Consider making the temporaries persistent, either providing them via function params or give static storage duration.</remarks>
///
real_t* blas_dot(cl_mem x,
                 cl_mem y,
                 cl_uint n,
                 OpenCL_handles* cl_handles,
                 BLAS_handles* blas_handles)
{
    cl_int CL_err = CL_SUCCESS;
    clblasStatus status[2];
    cl_event blas_exec[2], unmap;
    cl_mem result_buf, scratch_buf[2];

    real_t* result = (real_t*)malloc(2 * sizeof(real_t));

    result_buf = clCreateBuffer(cl_handles->ctx, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, 2 * sizeof(real_t), NULL, &CL_err);
    scratch_buf[0] = clCreateBuffer(cl_handles->ctx, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, n * sizeof(real_t), NULL, &CL_err);
    scratch_buf[1] = clCreateBuffer(cl_handles->ctx, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, n * sizeof(real_t), NULL, &CL_err);

#ifdef COMP_FLOAT
    status[0] = clblasSdot(n, result_buf, 0, x, 0, 1, x, 0, 1, scratch_buf[0], 1, cl_handles->exec_queues, 0, NULL, &blas_exec[0]);
    status[0] = clblasSdot(n, result_buf, 1, x, 0, 1, x, 0, 1, scratch_buf[1], 1, cl_handles->exec_queues, 0, NULL, &blas_exec[1]);
#else
    status[0] = clblasDdot(n, result_buf, 0, x, 0, 1, x, 0, 1, scratch_buf[0], 1, cl_handles->exec_queues, 0, NULL, &blas_exec[0]);
    status[0] = clblasDdot(n, result_buf, 1, x, 0, 1, x, 0, 1, scratch_buf[1], 1, cl_handles->exec_queues, 0, NULL, &blas_exec[1]);
#endif // COMP_FLOAT

    void* res = clEnqueueMapBuffer(cl_handles->read_queues[0], result_buf, CL_TRUE, CL_MAP_READ, 0, 2 * sizeof(real_t), 2, blas_exec, NULL, &CL_err);
    checkErr(CL_err, "clEnqueueMapMemObject(result_buf)");

    errno_t CRT_err = memcpy_s(result, 2 * sizeof(real_t), res, 2 * sizeof(real_t));
    if (CRT_err != 0)
        exit(EXIT_FAILURE);

    CL_err = clEnqueueUnmapMemObject(cl_handles->read_queues[0], result_buf, res, 0, NULL, &unmap);
    checkErr(CL_err, "clEnqueueUnmapMemObject(result_buf)");

    clWaitForEvents(1, &unmap);

    // cleanup
    clReleaseEvent(blas_exec[0]);
    clReleaseEvent(blas_exec[1]);
    clReleaseEvent(unmap);
    clReleaseMemObject(result_buf);
    clReleaseMemObject(scratch_buf[0]);
    clReleaseMemObject(scratch_buf[1]);

    return result;
}

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
                     OpenCL_handles* cl_handles)
{
    cl_int CL_err = CL_SUCCESS;

    clSetKernelArg(cl_handles->kernels[PhaseMod1], 0, sizeof(cl_mem), &xa);
    clSetKernelArg(cl_handles->kernels[PhaseMod1], 1, sizeof(cl_mem), &xb);
    clSetKernelArg(cl_handles->kernels[PhaseMod1], 2, sizeof(cl_mem), &xar);
    clSetKernelArg(cl_handles->kernels[PhaseMod1], 3, sizeof(cl_mem), &xbr);
    clSetKernelArg(cl_handles->kernels[PhaseMod1], 4, sizeof(real_t), &het1);
    clSetKernelArg(cl_handles->kernels[PhaseMod1], 5, sizeof(real_t), &sgnlt1);
    clSetKernelArg(cl_handles->kernels[PhaseMod1], 6, sizeof(cl_mem), &shft);
    clSetKernelArg(cl_handles->kernels[PhaseMod1], 7, sizeof(cl_int), &N);

    cl_event exec;
    size_t size_N = (size_t)N; // Helper variable to make pointer types match. Cast to silence warning

    CL_err = clEnqueueNDRangeKernel(cl_handles->exec_queues[0], cl_handles->kernels[PhaseMod1], 1, NULL, &size_N, NULL, 0, NULL, &exec);
    checkErr(CL_err, "clEnqueueNDRangeKernel(PhaseMod1)");

    clWaitForEvents(1, &exec);

    clReleaseEvent(exec);
}

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
                     OpenCL_handles* cl_handles)
{
    cl_int CL_err = CL_SUCCESS;

    clSetKernelArg(cl_handles->kernels[PhaseMod2], 0, sizeof(cl_mem), &xa);
    clSetKernelArg(cl_handles->kernels[PhaseMod2], 1, sizeof(cl_mem), &xb);
    clSetKernelArg(cl_handles->kernels[PhaseMod2], 2, sizeof(cl_mem), &xar);
    clSetKernelArg(cl_handles->kernels[PhaseMod2], 3, sizeof(cl_mem), &xbr);
    clSetKernelArg(cl_handles->kernels[PhaseMod2], 4, sizeof(real_t), &het1);
    clSetKernelArg(cl_handles->kernels[PhaseMod2], 5, sizeof(real_t), &sgnlt1);
    clSetKernelArg(cl_handles->kernels[PhaseMod2], 6, sizeof(cl_mem), &shft);
    clSetKernelArg(cl_handles->kernels[PhaseMod2], 7, sizeof(cl_int), &N);

    cl_event exec;
    size_t size_N = (size_t)N; // Helper variable to make pointer types match. Cast to silence warning

    CL_err = clEnqueueNDRangeKernel(cl_handles->exec_queues[0], cl_handles->kernels[PhaseMod2], 1, 0, &size_N, NULL, 0, NULL, &exec);
    checkErr(CL_err, "clEnqueueNDRangeKernel(PhaseMod2)");

    clWaitForEvents(1, &exec);

    clReleaseEvent(exec);
}

/// <summary>Compute F-statistics.</summary>
/// 
void compute_Fstat_gpu(cl_mem xa,
                       cl_mem xb,
                       cl_mem F,
                       cl_mem maa_d,
                       cl_mem mbb_d,
                       cl_int nmin,
                       cl_int nmax,
                       OpenCL_handles* cl_handles)
{
    cl_int CL_err = CL_SUCCESS;
    cl_int N = nmax - nmin;

    clSetKernelArg(cl_handles->kernels[ComputeFStat], 0, sizeof(cl_mem), &xa);
    clSetKernelArg(cl_handles->kernels[ComputeFStat], 1, sizeof(cl_mem), &xb);
    clSetKernelArg(cl_handles->kernels[ComputeFStat], 2, sizeof(cl_mem), &F);
    clSetKernelArg(cl_handles->kernels[ComputeFStat], 3, sizeof(cl_mem), &maa_d);
    clSetKernelArg(cl_handles->kernels[ComputeFStat], 4, sizeof(cl_mem), &mbb_d);
    clSetKernelArg(cl_handles->kernels[ComputeFStat], 5, sizeof(cl_int), &N);

    cl_event exec;
    size_t size_N = (size_t)N,
           size_nmin = (size_t)nmin; // Helper variable to make pointer types match. Cast to silence warning

    CL_err = clEnqueueNDRangeKernel(cl_handles->exec_queues[0], cl_handles->kernels[ComputeFStat], 1, &size_nmin, &size_N, NULL, 0, NULL, &exec);
    checkErr(CL_err, "clEnqueueNDRangeKernel(ComputeFStat)");

    clWaitForEvents(1, &exec);

    clReleaseEvent(exec);
}

/// <summary>Compute F-statistics.</summary>
///
void FStat_gpu_simple(cl_mem F_d,
    cl_uint nfft,
    cl_uint nav,
    OpenCL_handles* cl_handles)
{
    cl_int CL_err = CL_SUCCESS;
    cl_uint N = nfft / nav;
    size_t max_wgs;             // maximum supported wgs on the device (limited by register count)
    cl_ulong local_size;        // local memory size in bytes
    cl_uint ssi;                // shared size in num gentypes

    CL_err = clGetKernelWorkGroupInfo(cl_handles->kernels[FStatSimple], cl_handles->devs[0], CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &max_wgs, NULL);
    checkErr(CL_err, "clGetKernelWorkGroupInfo(FStatSimple, CL_KERNEL_WORK_GROUP_SIZE)");

    CL_err = clGetDeviceInfo(cl_handles->devs[0], CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &local_size, NULL);
    checkErr(CL_err, "clGetDeviceInfo(FStatSimple, CL_DEVICE_LOCAL_MEM_SIZE)");

    // How long is the array of local memory
    ssi = (cl_uint)(local_size / sizeof(real_t)); // Assume integer is enough to store gentype count (well, it better)

    clSetKernelArg(cl_handles->kernels[FStatSimple], 0, sizeof(cl_mem), &F_d);
    clSetKernelArg(cl_handles->kernels[FStatSimple], 1, ssi * sizeof(real_t), NULL);
    clSetKernelArg(cl_handles->kernels[FStatSimple], 2, sizeof(cl_uint), &ssi);
    clSetKernelArg(cl_handles->kernels[FStatSimple], 3, sizeof(cl_uint), &nav);

    size_t gsi = N;

    cl_event exec;
    CL_err = clEnqueueNDRangeKernel(cl_handles->exec_queues[0], cl_handles->kernels[FStatSimple], 1, NULL, &gsi, &max_wgs, 0, NULL, &exec);
    checkErr(CL_err, "clEnqueueNDRangeKernel(FStatSimple)");

    clWaitForEvents(1, &exec);

    clReleaseEvent(exec);
}

#ifdef GPUFSTAT
/* WARNING
   This won't necessarily work for other values than:
   NAV = 4096
   N = nmax - nmin = 507904
   For those given above it works fine.
   Reason is the "reduction depth", i.e. number of required \
   reduction kernel calls.
*/
void FStat_gpu(FLOAT_TYPE *F_d, int N, int nav, FLOAT_TYPE *mu_d, FLOAT_TYPE *mu_t_d) {

  int nav_blocks = N/nav;           //number of blocks
  int nav_threads = nav/BLOCK_SIZE; //number of blocks computing one nav-block
  int blocks = N/BLOCK_SIZE;

  //    CudaSafeCall ( cudaMalloc((void**)&cu_mu_t, sizeof(float)*blocks) );
  //    CudaSafeCall ( cudaMalloc((void**)&cu_mu, sizeof(float)*nav_blocks) );

  //sum fstat in blocks
  reduction_sum<BLOCK_SIZE_RED><<<blocks, BLOCK_SIZE_RED, BLOCK_SIZE_RED*sizeof(FLOAT_TYPE)>>>(F_d, mu_t_d, N);
  CudaCheckError();

  //sum blocks computed above and return 1/mu (number of divisions: blocks), then fstat_norm doesn't divide (potential number of divisions: N)
  reduction_sum<<<nav_blocks, nav_threads, nav_threads*sizeof(FLOAT_TYPE)>>>(mu_t_d, mu_d, blocks);
  CudaCheckError();
  
  //divide by mu/(2*NAV)
  fstat_norm<<<blocks, BLOCK_SIZE>>>(F_d, mu_d, N, nav);
  CudaCheckError();
  
}
#endif
