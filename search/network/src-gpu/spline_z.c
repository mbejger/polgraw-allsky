// Polgraw includes
#include <spline_z.h>

// clSPARSE includes
#include <clSPARSE.h>

// Standard C includes
#include <stdlib.h>          // calloc, free


/// <summary>Initialize the spline matrices.</summary>
/// <remarks>PCI Should replace it with kernels that initialize on the device.</remarks>
///
void init_spline_matrices(OpenCL_handles* cl_handles,
                          cl_mem cu_d,  // buffer of complex_devt
                          cl_mem cu_dl, // buffer of complex_devt
                          cl_mem cu_du, // buffer of complex_devt
                          cl_mem cu_B,  // buffer of complex_devt
                          int N)
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
    cu_d = clCreateBuffer(cl_handles->ctx,
                          CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                          (N - 1) * sizeof(complex_devt),
                          d,
                          &CL_err);
    checkErr(CL_err, "clCreateBuffer(cu_d)");

    cu_dl = clCreateBuffer(cl_handles->ctx,
                           CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                           (N - 1) * sizeof(complex_devt),
                           d,
                           &CL_err);
    checkErr(CL_err, "clCreateBuffer(cu_dl)");

    cu_du = clCreateBuffer(cl_handles->ctx,
                           CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                           (N - 1) * sizeof(complex_devt),
                           d,
                           &CL_err);
    checkErr(CL_err, "clCreateBuffer(cu_du)");

    // allocate B (or z) vector
    cu_B = clCreateBuffer(cl_handles->ctx,
                          CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR,
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

    complex_t pattern = {(real_t)0, (real_t)0};
    cl_event fill_event;

    clEnqueueFillBuffer(cl_handles->write_queues[0], cu_B, &pattern, sizeof(complex_t), 0, (N + 1) * sizeof(complex_t), 0, NULL, &fill_event);

    clWaitForEvents(1, &fill_event);
    clReleaseEvent(fill_event);

    computeB_gpu(cu_y, cu_B, N, cl_handles); // TODO almost certainly wrong indexing



  //create handle for sparse
  //cusparseHandle_t handle=0;
  //cusparseCreate(&handle);
  
  //compute vector Z
//  cusparseStatus_t status = cusparseZgtsv(
//					  handle,
//					  N-1,			//size of linear system (N+1 points - 2 = N-1)
//					  1,				//number of right-hand sides
//					  cu_dl,			//lower diagonal
//					  cu_d,			//diagonal
//					  cu_du,			//upper diagonal
//					  cu_B+1,		//right-hand-side vector
//					  N-1);			//leading dimension
//  if (status != CUSPARSE_STATUS_SUCCESS) {
//    printf("sparse status: %d\n", status);
//    printf("Blad cusparse!\n");
//  }
//  CudaCheckError();

  //here, vector cu_z=cu_B is computed
  //time to interpolate
//  interpolate<<< new_N/SPLINE_BLOCK_SIZE + 1 , SPLINE_BLOCK_SIZE >>>(
//								     cu_new_x,
//								     cu_new_y,
//								     cu_B,
//								     cu_y,
//								     N,
//								     new_N);
//  CudaCheckError();

}

//__global__ void computeB(cufftDoubleComplex *cu_y, cufftDoubleComplex *cu_B, int N);

//__global__ void test_vector(cufftDoubleComplex *a, int len, char id);

//__global__ void interpolate(double *new_x, cufftDoubleComplex *new_y, cufftDoubleComplex *z, cufftDoubleComplex *y,
//									 int N, int new_N);

//__global__ void test_vector(cufftDoubleComplex *a, int len, char id) {
//  int idx = threadIdx.x + blockIdx.x * blockDim.x;
//  if (idx < len) {
//    printf("%c %d(%d.%d): %lf\n", id, idx,  blockIdx.x, threadIdx.x, a[idx]);
//  }
//}
//

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
    CL_err = clEnqueueNDRangeKernel(cl_handles->exec_queues[0], cl_handles->kernels[ComputeB], 1, NULL, &N, NULL, 0, NULL, &exec);

    clWaitForEvents(1, &exec);

    clReleaseEvent(exec);
}


//__global__ void interpolate(double *new_x, cufftDoubleComplex *new_y, 
//			    cufftDoubleComplex *z, cufftDoubleComplex *y,
//			    int N, int new_N) {
//  double alpha = 1./6.;
//  int idx = threadIdx.x + blockIdx.x * blockDim.x;
//  if (idx < new_N) {
//    double x = new_x[idx];
//    
//    //get index of interval
//    int i = floor(x);
//    //compute value:
//    // S[i](x) = z[i+1]/6 * (x-x[i])**3 + z[i]/6 *	(x[i+1]-x)**3 + C[i]*(x-x[i]) + D[i]*(x[i+1]-x)
//    // C[i] = y[i+1] - z[i+1]/6
//    // D[i] = y[i] - z[i]/6
//    // x[i] = i
//    // x = new_x
//    double dist1 = x-i;
//    double dist2 = i+1-x;
//    
//    new_y[idx].x = dist1*( z[i+1].x*alpha*(dist1*dist1 - 1) + y[i+1].x ) +
//      dist2*( z[i].x*alpha*(dist2*dist2 - 1) + y[i].x );
//    new_y[idx].y = dist1*( z[i+1].y*alpha*(dist1*dist1 - 1) + y[i+1].y ) +
//      dist2*( z[i].y*alpha*(dist2*dist2 - 1) + y[i].y );
//    
//    // that change makes kernel ~2.3x faster
//    /*
//      new_y[idx].x = dist1*( z[i+1].x/6*dist1*dist1 +  (y[i+1].x-z[i+1].x/6) ) +
//      dist2*( z[i].x/6*dist2*dist2 + (y[i].x-z[i].x/6) );
//      new_y[idx].y = dist1*( z[i+1].y/6*dist1*dist1 +  (y[i+1].y-z[i+1].y/6) ) +
//      dist2*( z[i].y/6*dist2*dist2 + (y[i].y-z[i].y/6) );
//    */
//  }
//}

