#ifndef __UTIL_H__
#define __UTIL_H__

// OpenCL behavioral defines
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS 1

// clFFT includes
#include <clFFT.h>

// OpenCL includes
#include <CL/cl.h>


/// <summary>OpenCL error handling function.</summary>
///
void checkErr(cl_int err, const char * name);

/// <summary>clFFT error handling function.</summary>
///
void checkErrFFT(clfftStatus stat, const char * name);

#endif // __UTIL_H__
