#define NANO_INV 1000000000L

// Custom include
#include <timer.h>

// C90 includes
#include <stdio.h>          // perror

// C99 includes
#include <stdbool.h>        // _Bool


struct timespec get_current_time()
{
  struct timespec t;

  _Bool success = timespec_get(&t, TIME_UTC) == TIME_UTC;

  if (!success)
      perror("Failed to invoke CRT timer function.");

  return t;
}


double get_time_difference(struct timespec t0, struct timespec t1)
{
  return (t1.tv_sec-t0.tv_sec) + (double)(t1.tv_nsec-t0.tv_nsec) / (double)NANO_INV;
}
