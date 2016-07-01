#define NANO_INV 1000000000L
#include "timer.h"

struct timespec get_current_time(clockid_t cid) {
  struct timespec t;
  // CLOCK_PROCESS_CPUTIME_ID or CLOCK_REALTIME
  clock_gettime(cid, &t);
  return t;
}

double get_time_difference(struct timespec t0, struct timespec t1) { // in s
  return (t1.tv_sec-t0.tv_sec)+(double)(t1.tv_nsec-t0.tv_nsec)/(double)NANO_INV;
}
