// MSVC macro to include constants, such as M_PI (include before math.h)
#define _USE_MATH_DEFINES

// Standard C includes
//#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
//#include <unistd.h>
#include <malloc.h>
#include <complex.h>

/* JobCore file */
#include "struct.h"
#include "jobcore.h"
#include "auxi.h"
#include "settings.h"
#include "timer.h"
//#include "kernels.h"
#include "spline_z.h"
//#include "cuda_error.h"
#include <assert.h>
#include <floats.h>

// __constant__ double maa_d, mbb_d;  // mean Amplitude modulation functions

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
    int fd;
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

    aux->mu_t_d = clCreateBuffer(cl_handles->ctx, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, blocks * sizeof(real_t), NULL, &CL_err);
    aux->mu_d = clCreateBuffer(cl_handles->ctx, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, nav_blocks * sizeof(real_t), NULL, &CL_err);

    state = NULL;
    if (opts->checkp_flag)
        state = fopen(opts->qname, "w");

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
                    ftruncate(fileno(state), 0);
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
                 cl_mem F,                      // F-statistics array
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
    real_t shft1, phase, cp, sp;
    complex_t exph;

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
        modvir_gpu(sinalt, cosalt, sindelt, cosdelt,
                   sett->N, &ifo[n], aux, n);

        // Calculate detector positions with respect to baricenter
        nSource[0] = cosalt*cosdelt;
        nSource[1] = sinalt*cosdelt;
        nSource[2] = sindelt;

        shft1 = nSource[0] * ifo[n].sig.DetSSB[0]
            + nSource[1] * ifo[n].sig.DetSSB[1]
            + nSource[2] * ifo[n].sig.DetSSB[2];

        tshift_pmod_gpu(shft1, het0, nSource[0], nSource[1], nSource[2],
                        ifo[n].sig.xDat_d, fft_arr->xa_d, fft_arr->xb_d,
                        ifo[n].sig.shft_d, ifo[n].sig.shftf_d,
                        aux->tshift_d,
                        ifo[n].sig.aa_d, ifo[n].sig.bb_d,
                        ifo[n].sig.DetSSB_d,
                        sett->oms, sett->N, sett->nfft, sett->interpftpad, cl_handles);

        cl_event fft_exec[2];
        clfftEnqueueTransform(plans->pl_int, CLFFT_FORWARD, 1, cl_handles->exec_queues, 0, NULL, &fft_exec[0], fft_arr->xa_d, NULL, NULL /*May be slow, consider using tmp_buffer*/);
        clfftEnqueueTransform(plans->pl_int, CLFFT_FORWARD, 1, cl_handles->exec_queues, 0, NULL, &fft_exec[1], fft_arr->xb_d, NULL, NULL /*May be slow, consider using tmp_buffer*/);

        clWaitForEvents(2, fft_exec);

        resample_postfft_gpu(fft_arr->xa_d,
                             fft_arr->xb_d,
                             sett->nfft,
                             sett->Ninterp,
                             nyqst,
                             cl_handles);

        // Backward fft (len Ninterp = nfft*interpftpad)
        clfftEnqueueTransform(plans->pl_inv, CLFFT_BACKWARD, 1, cl_handles->exec_queues, 0, NULL, &fft_exec[0], fft_arr->xa_d, NULL, NULL /*May be slow, consider using tmp_buffer*/);
        clfftEnqueueTransform(plans->pl_inv, CLFFT_BACKWARD, 1, cl_handles->exec_queues, 0, NULL, &fft_exec[1], fft_arr->xb_d, NULL, NULL /*May be slow, consider using tmp_buffer*/);

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
                   aux->B_d);           //coefficient matrix

        gpu_interp(fft_arr->xb_d,       //input data
                   sett->Ninterp,       //input data length
                   aux->tshift_d,       //output time domain
                   ifo[n].sig.xDatmb_d, //output values
                   sett->N,             //output data length
                   aux->diag_d,         //diagonal
                   aux->ldiag_d,        //lower diagonal
                   aux->udiag_d,        //upper diagonal
                   aux->B_d);           //coefficient matrix

        ft = 1. / ifo[n].sig.sig2;
        //    cublasZdscal( scale, sett->N, &ft, ifo[n].sig.xDatma_d, 1);
        //    cublasZdscal( scale, sett->N, &ft, ifo[n].sig.xDatmb_d, 1);
        //    CudaCheckError();

    } // end of detector loop 

    double _maa = 0.;
    double _mbb = 0.;

    for (int n = 0; n<sett->nifo; ++n)
    {
        //    double aatemp = 0., bbtemp = 0.;
        //    // square sums of modulation factors
        //    cublasDdot (scale, sett->N , 
        //      (const double *)ifo[n].sig.aa_d, 1,
        //      (const double *)ifo[n].sig.aa_d, 1,
        //      &aatemp);
        //    cublasDdot (scale, sett->N , 
        //      (const double *)ifo[n].sig.bb_d, 1,
        //      (const double *)ifo[n].sig.bb_d, 1,
        //      &bbtemp);
        //    
        //    /* or use sqr( cublasSnrm2()) */
        //    _maa += aatemp/ifo[n].sig.sig2;
        //    _mbb += bbtemp/ifo[n].sig.sig2;
    }
    //  CudaCheckError();

    //  printf("maa_d=%f", _maa);
    //  cudaMemcpyToSymbol(maa_d, &_maa, sizeof(double), 0, cudaMemcpyHostToDevice);
    //  cudaMemcpyToSymbol(mbb_d, &_mbb, sizeof(double), 0, cudaMemcpyHostToDevice);
    //  CudaCheckError();  


    /* Spindown loop */

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
        smin = trunc((sett->Smin - nn*sett->M[9] - mm*sett->M[13]) / sett->M[5]);
        smax = trunc(-(nn*sett->M[9] + mm*sett->M[13] + sett->Smax) / sett->M[5]);
    }

    printf("\n>>%d\t%d\t%d\t[%d..%d]\n", *FNum, mm, nn, smin, smax);

    // No-spindown calculations
    if (opts->s0_flag) smin = smax;

    // if spindown parameter is taken into account, smin != smax
    for (ss = smin; ss <= smax; ++ss) {

#if TIMERS>2
        tstart = get_current_time();
#endif 

        // Spindown parameter
        sgnlt[1] = ss*sett->M[5] + nn*sett->M[9] + mm*sett->M[13];

        //    // Spindown range
        //    if(sgnlt[1] >= -sett->Smax && sgnlt[1] <= sett->Smax) { 

        int ii;
        double Fc, het1;

#ifdef VERBOSE
        //print a 'dot' every new spindown
        printf("."); fflush(stdout);
#endif 

        het1 = fmod(ss*sett->M[4], sett->M[0]);
        if (het1<0) het1 += sett->M[0];

        sgnl0 = het0 + het1;
        //      printf("%d  %d\n", BLOCK_SIZE, (sett->N + BLOCK_SIZE - 1)/BLOCK_SIZE );


        //      phase_mod_1<<<(sett->N + BLOCK_SIZE - 1)/BLOCK_SIZE, BLOCK_SIZE>>>
        //        ( fft_arr->xa_d, fft_arr->xb_d,
        //          ifo[0].sig.xDatma_d, ifo[0].sig.xDatmb_d,
        //          het1, sgnlt[1], ifo[0].sig.shft_d,
        //          sett->N );
        //      
        //      cudaDeviceSynchronize();

        for (int n = 1; n<sett->nifo; ++n) {
            //	phase_mod_2<<<(sett->N + BLOCK_SIZE - 1)/BLOCK_SIZE, BLOCK_SIZE>>>
            //	  ( fft_arr->xa_d, fft_arr->xb_d,
            //	    ifo[n].sig.xDatma_d, ifo[n].sig.xDatmb_d,
            //	    het1, sgnlt[1], ifo[n].sig.shft_d,
            //	    sett->N );
        }

        // initialize arrays to 0. with integer 0
        // assuming double , remember to change when switching to float
        // cuMemsetD32Async((CUdeviceptr) (fft_arr->xa_d + sett->N), 0,
        //       (sett->nfftf - sett->N)*2*(sizeof(double)/4), NULL);
        // cuMemsetD32Async((CUdeviceptr) (fft_arr->xb_d + sett->N), 0,
        //       (sett->nfftf - sett->N)*2*(sizeof(double)/4), NULL);
        // CudaCheckError();

        // fft length fftpad*nfft
        //cufftExecZ2Z(plans->plan, fft_arr->xa_d, fft_arr->xa_d, CUFFT_FORWARD);
        //cufftExecZ2Z(plans->plan, fft_arr->xb_d, fft_arr->xb_d, CUFFT_FORWARD);

        (*FNum)++;

        // compute_Fstat<<<(sett->nmax-sett->nmin + BLOCK_SIZE - 1)/BLOCK_SIZE, BLOCK_SIZE>>>
        //   ( fft_arr->xa_d + sett->nmin,
        //     fft_arr->xb_d + sett->nmin,
        //     F_d + sett->nmin,
        //     sett->nmax - sett->nmin );
        // CudaCheckError();

#define GPUFSTAT_NO
#ifdef GPUFSTAT
        if (!(opts->white_flag))  // if the noise is not white noise
            FStat_gpu(F_d + sett->nmin, sett->nmax - sett->nmin, NAV, aux->mu_d, aux->mu_t_d);

#else
        if (!(opts->white_flag))  // if the noise is not white noise
            FStat_gpu_simple(F_d + sett->nmin, sett->nmax - sett->nmin, NAVFSTAT);
#endif

        // CudaSafeCall ( cudaMemcpy(F, F_d, 2*sett->nfft*sizeof(FLOAT_TYPE), cudaMemcpyDeviceToHost));
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
    real_t cosalfr, sinalfr, c2d, c2sd, c, s, c2s, cs;

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
    CL_err = clEnqueueNDRangeKernel(cl_handles->exec_queues[0], cl_handles->kernels[Modvir], 1, NULL, &Np, NULL, 0, NULL, &exec);

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
                     int N,
                     int nfft,
                     int interpftpad,
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
    CL_err = clEnqueueNDRangeKernel(cl_handles->exec_queues[0], cl_handles->kernels[TShiftPMod], 1, NULL, &nfft, NULL, 0, NULL, &exec);

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
    CL_err = clEnqueueNDRangeKernel(cl_handles->exec_queues[0], cl_handles->kernels[ResamplePostFFT], 1, NULL, &Ninterp, NULL, 0, NULL, &exec);

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

/* just simple patch - to be replaced */
void FStat_gpu_simple(FLOAT_TYPE *F_d, int nfft, int nav) {
  
  // int blocks = nfft/nav;
  // fstat_norm_simple<<<blocks, 1>>>(F_d, nav);
  // CudaCheckError();

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
