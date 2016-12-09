// Polgraw includes
#include <CL/util.h>    // function declarations

// Standard C includes
#include <stdio.h>      // printf_s
#include <stdlib.h>     // exit

/// <summary>OpenCL error handling function.</summary>
/// 
void checkErr(cl_int err, const char * name)
{
    if (err != CL_SUCCESS)
    {
#ifdef _WIN32
        int count = printf_s("ERROR: %s (%i)\n", name, err);
#else
        printf("ERROR: %s (%i)\n", name, err);
#endif
        exit(err);
    }
}

/// <summary>clFFT error handling function.</summary>
///
void checkErrFFT(clfftStatus stat, const char * name)
{
    if (stat != CLFFT_SUCCESS)
    {
        printf("ERROR: %s (%i)\n", name, stat);

        exit(stat);
    }
}
