#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cufft.h>

#include "cuda_error.h"

#define SPLINE_BLOCK_SIZE 256

void init_spline_matrices(cufftDoubleComplex **cu_d, cufftDoubleComplex **cu_dl,
								  cufftDoubleComplex **cu_du, cufftDoubleComplex **cu_B,
								  int Np);

void gpu_interp(cufftDoubleComplex *cu_y, int Np, double *cu_new_x, cufftDoubleComplex *cu_new_y,
				int new_N, cufftDoubleComplex *cu_d, cufftDoubleComplex *cu_dl, cufftDoubleComplex *cu_du,
				cufftDoubleComplex *cu_B);

