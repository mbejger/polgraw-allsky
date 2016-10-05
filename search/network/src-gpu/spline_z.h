// OpenCL includes
#include <CL/cl.h>

// Polgraw includes
#include <floats.h>     // COMPLEX_TYPE

#define SPLINE_BLOCK_SIZE 256

void init_spline_matrices(cl_mem cu_d,  // buffer of COMPLEX_TYPE
                          cl_mem cu_dl, // buffer of COMPLEX_TYPE
                          cl_mem cu_du, // buffer of COMPLEX_TYPE
                          cl_mem cu_B,  // buffer of COMPLEX_TYPE
                          int Np);

void gpu_interp(COMPLEX_TYPE* cu_y,
                int Np,
                double *cu_new_x,
                COMPLEX_TYPE *cu_new_y,
                int new_N,
                COMPLEX_TYPE *cu_d,
                COMPLEX_TYPE *cu_dl,
                COMPLEX_TYPE *cu_du,
                COMPLEX_TYPE *cu_B);
