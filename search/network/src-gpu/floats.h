#ifndef __FLOATS_H__
#define __FLOATS_H__

// OpenCL includes
#include <CL/cl.h>      // cl_float, cl_double

// Standard C includes
#include <complex.h>    // _Dcomplex


#undef COMP_FLOAT

//changing computations in spindown loop to single-precision arithmetic
#ifdef COMP_FLOAT //if single-precision
#define CLFFT_TRANSFORM_PRECISION CLFFT_SINGLE
#define CLFFT_TRANSFORM_LAYOUT CLFFT_REAL
#define COMPLEX_TYPE cufftComplex
#define FLOAT_TYPE float
#define HOST_COMPLEX_TYPE complex float
#else //if double-precision
#define CLFFT_TRANSFORM_PRECISION CLFFT_DOUBLE
#define CLFFT_TRANSFORM_LAYOUT CLFFT_COMPLEX_INTERLEAVED
typedef cl_double real_t;
typedef cl_double2 complex_devt;
#define FLOAT_TYPE cl_double
#define COMPLEX_TYPE cl_double2
#ifdef _WIN32
typedef _Dcomplex complex_t;
#define HOST_COMPLEX_TYPE _Dcomplex
#else
typedef complex double complex_t;
#define HOST_COMPLEX_TYPE complex double
#endif // WIN32
#endif


#endif
