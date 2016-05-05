#define NANO_INV 1000000000L
#include "timer.h"

struct timespec get_current_time() {
  struct timespec t;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t);
  return t;
}

double get_time_difference(struct timespec t0, struct timespec t1) { // in s
  return (t1.tv_sec-t0.tv_sec)+(double)(t1.tv_nsec-t0.tv_nsec)/(double)NANO_INV;
}
