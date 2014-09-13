#ifndef __TIMER_H__
#define __TIMER_H__

#include <sys/time.h>


struct timeval get_current_time();

double get_time_difference(struct timeval tstart, struct timeval tend);

#endif
