#ifndef __TIMER_H__
#define __TIMER_H__

#include <time.h>

struct timespec get_current_time(clockid_t cid);

double get_time_difference(struct timespec t0, struct timespec t1);

#endif
