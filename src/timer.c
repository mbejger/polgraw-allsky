

#include "timer.h"

struct timeval get_current_time() {
	struct timeval t;
	gettimeofday(&t, 0);
	return t;
}

double get_time_difference(struct timeval tstart, struct timeval tend) { // in s
	return (tend.tv_sec - tstart.tv_sec) + 1e-6 * (tend.tv_usec - tstart.tv_usec);
}
