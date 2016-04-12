#ifndef __FLOATS_H__
#define __FLOATS_H__

#undef COMP_FLOAT

//changing computations in spindown loop to single-precision arithmetic
#ifdef COMP_FLOAT //if single-precision
        #define CUFFT_TRANSFORM_TYPE CUFFT_C2C
        #define CUFFT_EXEC_FFT cufftExecC2C
        #define COMPLEX_TYPE cufftComplex
        #define FLOAT_TYPE float
#else //if double-precision
        #define CUFFT_TRANSFORM_TYPE CUFFT_Z2Z
        #define CUFFT_EXEC_FFT cufftExecZ2Z
        #define COMPLEX_TYPE cufftDoubleComplex
        #define FLOAT_TYPE double
#endif


#endif
