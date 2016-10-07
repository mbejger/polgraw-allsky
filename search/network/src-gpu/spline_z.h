// Polgraw includes
#include <floats.h>     // COMPLEX_TYPE
#include <struct.h>     // OpenCL_handles

// OpenCL includes
#include <CL/cl.h>      // cl_mem

#define SPLINE_BLOCK_SIZE 256

/// <summary>Initialize the spline matrices.</summary>
/// <remarks>PCI Should replace it with kernels that initialize on the device.</remarks>
///
void init_spline_matrices(OpenCL_handles* cl_handles, 
                          cl_mem cu_d,  // buffer of complex_devt
                          cl_mem cu_dl, // buffer of complex_devt
                          cl_mem cu_du, // buffer of complex_devt
                          cl_mem cu_B,  // buffer of complex_devt
                          int N);

void gpu_interp(COMPLEX_TYPE* cu_y,
                int Np,
                double *cu_new_x,
                COMPLEX_TYPE *cu_new_y,
                int new_N,
                COMPLEX_TYPE *cu_d,
                COMPLEX_TYPE *cu_dl,
                COMPLEX_TYPE *cu_du,
                COMPLEX_TYPE *cu_B);
