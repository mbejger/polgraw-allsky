#ifndef __TIMER_H__
#define __TIMER_H__

#include <sys/time.h>

struct timespec get_current_time();

double get_time_difference(struct timespec t0, struct timespec t1);

#endif
