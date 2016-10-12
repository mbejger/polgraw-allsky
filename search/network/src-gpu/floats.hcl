#ifndef __FLOATS_HCL__
#define __FLOATS_HCL__

#ifdef COMP_FLOAT
typedef float real_t;
typedef float2 complex_t;
#else
typedef double real_t;
typedef double2 complex_t;
#endif // COMP_FLOAT


#endif // __FLOATS_HCL__
