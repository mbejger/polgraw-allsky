/* 
Ver. 1.0 
Function takes candidate parameters and 
number of bins (as arguments) and creates 
grid around it.
MS
*/

#include <stdio.h>
#include <stdlib.h>
#include "neigh.h"

float** neigh(float s1, float s2, int bins){ 
	float **arr;
	int rows, cols = 2;
	rows = bins*bins;
	int k;
// Allocation of memory for martix
  	arr = (float **)malloc(rows*sizeof(float *));
  
  	for (k=0; k<rows; k++) arr[k] = (float *)calloc(cols, sizeof(float));


	float perc = 0.1;	// % define neighbourhood around each parameter
	float beg1, beg2;
	float width1, width2;
	float m1, m2;
	int i1, i2, i;
	beg1 = s1 - (s1*perc);
	width1 = 2*perc*s1/bins;
	width2 = 2*perc*s2/bins;
	i = 0;
	for(i1 = 0; i1 < bins; i1++){
		beg2 = s2 - (s2*perc);
		for(i2 = 0; i2 < bins; i2++, i++){
			arr[i][0] = beg1;
			arr[i][1] = beg2;
			beg2 = beg2 + width2;
		}
		beg1 = beg1 + width1;
	}
	return arr;

}
