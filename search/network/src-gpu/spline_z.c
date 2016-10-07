// Polgraw includes
#include <spline_z.h>

// Standard C includes
#include <stdio.h>
//#include <cuda.h>
//#include <cuda_runtime_api.h>
//#include "cusparse_v2.h"
//#include <cufft.h>

//__global__ void computeB(cufftDoubleComplex *cu_y, cufftDoubleComplex *cu_B, int N);

//__global__ void test_vector(cufftDoubleComplex *a, int len, char id);

//__global__ void interpolate(double *new_x, cufftDoubleComplex *new_y, cufftDoubleComplex *z, cufftDoubleComplex *y,
//									 int N, int new_N);


void gpu_interp(COMPLEX_TYPE* cu_y,
                int Np,
                double *cu_new_x,
                COMPLEX_TYPE *cu_new_y,
                int new_N,
                COMPLEX_TYPE *cu_d,
                COMPLEX_TYPE *cu_dl,
                COMPLEX_TYPE *cu_du,
                COMPLEX_TYPE *cu_B) {


//  N-=1; //N is number of intervals here
  
	//allocate and compute vector B=z (replaced on gtsv)
	// z has size N+1 (i=0..N), but we solve only for (i=1..N-1)
	// (z[0] = z[N] = 0) because of `natural conditions` of spline
  //CudaSafeCall( cudaMemset(cu_B, 0, sizeof(cufftDoubleComplex)*(N+1))); //set values to zero
  //computeB<<< (N-1)/SPLINE_BLOCK_SIZE + 1 , SPLINE_BLOCK_SIZE >>>(cu_y, cu_B+1, N);
  //CudaCheckError();

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


//__global__ void test_vector(cufftDoubleComplex *a, int len, char id) {
//  int idx = threadIdx.x + blockIdx.x * blockDim.x;
//  if (idx < len) {
//    printf("%c %d(%d.%d): %lf\n", id, idx,  blockIdx.x, threadIdx.x, a[idx]);
//  }
//}
//
//
//__global__ void computeB(cufftDoubleComplex *y, cufftDoubleComplex *B, int N) {
//  int idx = threadIdx.x + blockIdx.x * blockDim.x;
//  if (idx < N-1) {
//    B[idx].x = 6*(y[idx+2].x - 2*y[idx+1].x + y[idx].x);
//    B[idx].y = 6*(y[idx+2].y - 2*y[idx+1].y + y[idx].y);
//  }
//}


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

    d  = (complex_devt*)calloc(N - 1, sizeof(complex_devt));
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
