#ifndef __UTIL_H__
#define __UTIL_H__

// OpenCL includes
#include <CL/cl.h>


/// <summary>OpenCL error handling function.</summary>
///
void checkErr(cl_int err, const char * name);

#endif // __UTIL_H__
