// Polgraw includes
#include <CL/util.h>        // checkErr
#include <spline_z.h>

// Standard C includes
#include <stdlib.h>          // calloc, free


/// <summary>Initialize the spline matrices.</summary>
/// <remarks>PCI Should replace it with kernels that initialize on the device.</remarks>
///
void init_spline_matrices(OpenCL_handles* cl_handles,
                          cl_mem* cu_d,  // buffer of complex_devt
                          cl_mem* cu_dl, // buffer of complex_devt
                          cl_mem* cu_du, // buffer of complex_devt
                          cl_mem* cu_B,  // buffer of complex_devt
                          cl_int N)
{
    cl_int CL_err = CL_SUCCESS;
    N -= 1; // N is number of intervals here

    complex_devt *d, *du, *dl;

    d = (complex_devt*)calloc(N - 1, sizeof(complex_devt));
    du = (complex_devt*)calloc(N - 1, sizeof(complex_devt));
    dl = (complex_devt*)calloc(N - 1, sizeof(complex_devt));

    for (int i = 0; i<N - 2; i++)
    {
        dl[i + 1].s[0] = 1;
        du[i].s[0] = 1;
        d[i].s[0] = 4;

        dl[i].s[1] = 0;
        du[i].s[1] = 0;
        d[i].s[1] = 0;
    }

    // dl[0] is 0 and du[N-2]=0
    dl[0].s[0] = 0;
    du[N - 2].s[0] = 0;
    d[N - 2].s[0] = 4;

    dl[N - 2].s[1] = 0;
    du[N - 2].s[1] = 0;
    d[N - 2].s[1] = 0;

    // copy to gpu
    *cu_d = clCreateBuffer(cl_handles->ctx,
                           CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                           (N - 1) * sizeof(complex_devt),
                           d,
                           &CL_err);
    checkErr(CL_err, "clCreateBuffer(cu_d)");

    *cu_dl = clCreateBuffer(cl_handles->ctx,
                            CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                            (N - 1) * sizeof(complex_devt),
                            d,
                            &CL_err);
    checkErr(CL_err, "clCreateBuffer(cu_dl)");

    *cu_du = clCreateBuffer(cl_handles->ctx,
                            CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                            (N - 1) * sizeof(complex_devt),
                            d,
                            &CL_err);
    checkErr(CL_err, "clCreateBuffer(cu_du)");

    // allocate B (or z) vector
    *cu_B = clCreateBuffer(cl_handles->ctx,
                           CL_MEM_READ_WRITE,
                           (N + 1) * sizeof(complex_devt),
                           NULL,
                           &CL_err);
    checkErr(CL_err, "clCreateBuffer(cu_B)");

    //clean up
    free(d);
    free(du);
    free(dl);
}

/// <summary>Spline interpolation to xDatma, xDatmb arrays.</summary>
///
void gpu_interp(cl_mem cu_y,                // buffer of complex_t
                cl_int N,
                cl_mem cu_new_x,            // buffer of real_t
                cl_mem cu_new_y,            // buffer of complex_t
                cl_int new_N,
                cl_mem cu_d,                // buffer of complex_t
                cl_mem cu_dl,               // buffer of complex_t
                cl_mem cu_du,               // buffer of complex_t
                cl_mem cu_B,                // buffer of complex_t
                OpenCL_handles* cl_handles) // handles to OpenCL resources
{
    N-=1; // N is number of intervals here
  
    // allocate and compute vector B=z (replaced on gtsv)
    // z has size N+1 (i=0..N), but we solve only for (i=1..N-1)
    // (z[0] = z[N] = 0) because of `natural conditions` of spline
#ifdef WIN32
    complex_t pattern = {(real_t)0, (real_t)0};
#else
    complex_t pattern = 0;
#endif
    cl_event fill_event;

    clEnqueueFillBuffer(cl_handles->write_queues[0], cu_B, &pattern, sizeof(complex_t), 0, (N + 1) * sizeof(complex_t), 0, NULL, &fill_event);

    clWaitForEvents(1, &fill_event);
    clReleaseEvent(fill_event);

    computeB_gpu(cu_y, cu_B, N, cl_handles); // TODO almost certainly wrong indexing

    tridiagMul_gpu(cu_dl,
                   cu_d,
                   cu_du,
                   cu_B,
                   N + 1,
                   cl_handles); // TODO almost certainly wrong indexing of cu_B + use persistent tmp buffer


    // here, vector cu_z=cu_B is computed
    // time to interpolate


  interpolate_gpu(cu_new_x,
                  cu_new_y,
                  cu_B,
                  cu_y,
                  N,
                  new_N,
                  cl_handles);
}

/// <summary>The purpose of this function was undocumented.</summary>
///
void computeB_gpu(cl_mem y,
                  cl_mem B,
                  cl_int N,
                  OpenCL_handles* cl_handles)
{
    cl_int CL_err = CL_SUCCESS;

    clSetKernelArg(cl_handles->kernels[ComputeB], 0, sizeof(cl_mem), &y);
    clSetKernelArg(cl_handles->kernels[ComputeB], 1, sizeof(cl_mem), &B);
    clSetKernelArg(cl_handles->kernels[ComputeB], 2, sizeof(cl_int), &N);

    cl_event exec;
    size_t size_N = (size_t)N; // Helper variable to make pointer types match. Cast to silence warning

    CL_err = clEnqueueNDRangeKernel(cl_handles->exec_queues[0], cl_handles->kernels[ComputeB], 1, NULL, &size_N, NULL, 0, NULL, &exec);

    clWaitForEvents(1, &exec);

    clReleaseEvent(exec);
}

/// <summary>Multiplies the tridiagonal matrix specified by <c>{dl, d, du}</c> with dense vector <c>x</c>.</summary>
///
void tridiagMul_gpu(cl_mem dl,
                    cl_mem d,
                    cl_mem du,
                    cl_mem x,
                    cl_int length,
                    OpenCL_handles* cl_handles)
{
    cl_int CL_err = CL_SUCCESS;

    cl_mem y = clCreateBuffer(cl_handles->ctx, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, length * sizeof(complex_t), NULL, &CL_err);

    clSetKernelArg(cl_handles->kernels[TriDiagMul], 0, sizeof(cl_mem), &dl);
    clSetKernelArg(cl_handles->kernels[TriDiagMul], 1, sizeof(cl_mem), &d);
    clSetKernelArg(cl_handles->kernels[TriDiagMul], 2, sizeof(cl_mem), &du);
    clSetKernelArg(cl_handles->kernels[TriDiagMul], 3, sizeof(cl_mem), &x);
    clSetKernelArg(cl_handles->kernels[TriDiagMul], 4, sizeof(cl_mem), &y);

    cl_event events[2];
    size_t size_length = (size_t)length; // Helper variable to make pointer types match. Cast to silence warning

    CL_err = clEnqueueNDRangeKernel(cl_handles->exec_queues[0], cl_handles->kernels[TriDiagMul], 1, NULL, &size_length, NULL, 0, NULL, &events[0]);
    CL_err = clEnqueueCopyBuffer(cl_handles->write_queues[0], y, x, 0, 0, length * sizeof(complex_t), 1, &events[0], &events[1]);

    clWaitForEvents(2, events);

    for (size_t i = 0 ; i < 2; ++i) clReleaseEvent(events[i]);
    clReleaseMemObject(y);
}

/// <summary>The purpose of this function was undocumented.</summary>
///
void interpolate_gpu(cl_mem new_x,
                     cl_mem new_y,
                     cl_mem z,
                     cl_mem y,
                     cl_int N,
                     cl_int new_N,
                     OpenCL_handles* cl_handles)
{
    cl_int CL_err = CL_SUCCESS;

    clSetKernelArg(cl_handles->kernels[Interpolate], 0, sizeof(cl_mem), &new_x);
    clSetKernelArg(cl_handles->kernels[Interpolate], 1, sizeof(cl_mem), &new_y);
    clSetKernelArg(cl_handles->kernels[Interpolate], 2, sizeof(cl_mem), &z);
    clSetKernelArg(cl_handles->kernels[Interpolate], 3, sizeof(cl_mem), &y);
    clSetKernelArg(cl_handles->kernels[Interpolate], 4, sizeof(cl_int), &N);
    clSetKernelArg(cl_handles->kernels[Interpolate], 5, sizeof(cl_int), &new_N);

    cl_event exec;
    size_t size_new_N = (size_t)new_N; // Helper variable to make pointer types match. Cast to silence warning
    CL_err = clEnqueueNDRangeKernel(cl_handles->exec_queues[0], cl_handles->kernels[Interpolate], 1, NULL, &size_new_N, NULL, 0, NULL, &exec);

    clWaitForEvents(1, &exec);

    clReleaseEvent(exec);
}
