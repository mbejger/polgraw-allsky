/* Ver. 1.0
Main followup programm. Required files
(all should be in directory 'data' - given 
as argument) :
'candidates.coi' - File with candidates from coincidences
'list.txt'
Ver. 2.0
Functions glue and neigh added!
Ver. 3.0
Mesh adaptive direct search (MADS) added (to find real
maximum)
Ver. 4.0
Simplex added
Ver. 5.0
Function neigh moved to followup.c
All vectors/arrays declarated using yeppp! library 
 = FASTER!

MS
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
//#include <fcntl.h>
#include <getopt.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <time.h>
#include <dirent.h>
#include <omp.h>

#include "auxi.h"
#include "settings.h"
#include "struct.h"
#include "init.h"
//#include "timer.h"

//#include "glue.h"
//#include "neigh.h"

#include <assert.h>
#if defined(SLEEF)
//#include "sleef-2.80/purec/sleef.h"
#include <sleefsimd.h>
#elif defined(YEPPP)
#include <yepLibrary.h>
#include <yepCore.h>
#include <yepMath.h>
#endif

// Default output and data directories

#ifndef PREFIX
#define PREFIX ./FSTAT_OUT
#endif

#ifndef DTAPREFIX
#define DTAPREFIX .
#endif

#define ZEPS 1e-10

//Function neigh takes candidate parameters and number of bins (as arguments) and creates grid around it.
void neigh(double *m, double *perc, int b, double **arr){ 
//	double **array;
	int rows, cols = 4;
	rows = pow((b+1),4);
	int k;
// Allocation of memory for martix
/*#ifdef YEPPP
    	yepLibrary_Init();
	Yep64f **array = (Yep64f**)malloc(rows*sizeof(Yep64f));
	for (k=0; k < rows; k++) array[k] = (Yep64f*)calloc(rows,sizeof(Yep64f));
    	enum YepStatus status;
#else
  	array = (double **)malloc(rows*sizeof(double *));
  	for (k=0; k < rows; k++) array[k] = (double *)calloc(cols, sizeof(double));
#endif */
	double beg[4];
	double width[4];
	int i1, i2, i3, i4, j, i;
	for (j = 0; j < 4; j++) {
		width[j] = 2*perc[j]*m[j]/b;
	}
	i = 0;
	beg[0] = m[0]*(1 - perc[0]);
	for (i1 = 0; i1 < (b + 1); i1++){
		beg[1] = m[1]*(1 - perc[1]);
		for (i2 = 0; i2 < (b + 1); i2++){
			beg[2] = m[2]*(1 - perc[2]);
			for (i3 = 0; i3 < (b + 1); i3++){
				beg[3] = m[3]*(1 - perc[3]);
				for (i4 = 0; i4 < (b + 1); i4++){
					for (j = 0; j < 4; j++) {
						arr[i][j] = beg[j];
					}
					beg[3] = beg[3] + width[3];
					i++;
				}
				beg[2] = beg[2] + width[2];
			}
			beg[1] = beg[1] + width[1];
		}
		beg[0] = beg[0] + width[0];
	}

//	return array;

}

//Function takes calculated range around point (from grid.bin file) and creates grid around it.

void neigh_from_range(double **range, int b, double **arr){

	int rows, cols = 4;
	rows = pow((b+1),4);
	int k;
// Allocation of memory for martix
/*#ifdef YEPPP
    	yepLibrary_Init();
	Yep64f **array = (Yep64f**)malloc(rows*sizeof(Yep64f));
	for (k=0; k < rows; k++) array[k] = (Yep64f*)calloc(rows,sizeof(Yep64f));
    	enum YepStatus status;
#else 
  	array = (double **)malloc(rows*sizeof(double *));
  	for (k=0; k < rows; k++) array[k] = (double *)calloc(cols, sizeof(double));
#endif */
	double beg[4];
	double width[4];
	int i1, i2, i3, i4, j, i;
	for (j = 0; j < 4; j++) {
		width[j] = (range[j][1] - range[j][0])/b;
	}
	i = 0;
	beg[0] = range[0][0];
	for (i1 = 0; i1 < (b + 1); i1++){
		beg[1] = range[1][0];
		for (i2 = 0; i2 < (b + 1); i2++){
			beg[2] = range[2][0];
			for (i3 = 0; i3 < (b + 1); i3++){
				beg[3] = range[3][0];
				for (i4 = 0; i4 < (b + 1); i4++){
					for (j = 0; j < 4; j++) {
						arr[i][j] = beg[j];
					}
					beg[3] = beg[3] + width[3];
					i++;
				}
				beg[2] = beg[2] + width[2];
			}
			beg[1] = beg[1] + width[1];
		}
		beg[0] = beg[0] + width[0];
	}

} 

// Allocation of memory for double martix with given number of rows and columns
double** matrix(int rows, int cols) {

  	int k;
	double **m;
  	m = (double **)malloc(rows*sizeof(double *));
  
  	for (k=0; k < rows; k++)
    		m[k] = (double *)calloc(cols, sizeof(double));
  
  	return m;
}
// Allocation of memory for complex double martix with given number of rows and columns
complex double** matrix_complex(int rows, int cols) {

  	int k;
	complex double **m;
  	m = (complex double **)malloc(rows*sizeof(complex double *));
  
  	for (k=0; k < rows; k++)
    		m[k] = (complex double *)calloc(cols, sizeof(complex double));
  
  	return m;
}

// Allocation of memory for vector and martix with given number of rows and columns
double * alloc_vector(int cols){
	return (double *) malloc(sizeof(double) * cols);
}

int * alloc_vector_int(int cols){
	return (int*) malloc(sizeof(int) * cols);
}

void free_vector(double * vector, int cols){
	free(vector);
}

void free_vector_int(int * vector, int cols){
	free(vector);
}

double ** alloc_matrix(int rows, int cols){
	int i;
	double ** matrix = (double **) malloc(sizeof(double *) * rows);
	for (i = 0; i < rows; i++) matrix[i] = alloc_vector(cols);
	return matrix;
}
void free_matrix(double **matrix, int rows, int cols){
	int i;
	for (i = 0; i < rows; i++) free(matrix[i]);
	free(matrix);
}
void free_matrix_complex(complex double **matrix, int rows, int cols){
	int i;
	for (i = 0; i < rows; i++) free(matrix[i]);
	free(matrix);
}

// Function calculates F-statistics in given point
double* Fstatnet(Search_settings *sett, double *sgnlo, double *nSource, double **sigaa, double **sigbb){

	int i = 0, n = 0; 
	double aatemp, bbtemp, aa = 0., bb = 0.;
	complex double exph, xasum, xbsum;
	double **shft; 			// old ifo[i].sig.shft[n]
	complex double **xDatma, **xDatmb;	// old ifo[n].sig.xDatma[i], ifo[n].sig.xDatmb[i]
        shft = matrix(sett->nifo, sett->N);
	xDatma = matrix_complex(sett->nifo, sett->N);
	xDatmb = matrix_complex(sett->nifo, sett->N);		


#ifdef YEPPP
	int VLEN = sett->N;
	yepLibrary_Init();
//Yep64f *x = (Yep64f*)calloc(ARRAY_SIZE, sizeof(Yep64f));
	Yep64f *_sph = (Yep64f*)malloc(sizeof(Yep64f)*VLEN);
	Yep64f *_cph = (Yep64f*)malloc(sizeof(Yep64f)*VLEN);
	Yep64f *phase = (Yep64f*)malloc(sizeof(Yep64f)*VLEN); 
	Yep64f *fstat_out = (Yep64f*)malloc(sizeof(Yep64f)*11); 
	enum YepStatus status;
#endif

//From jobcore.c, line 237 
//Loop for each detector 1

	xasum = 0 - I * 0;
	xbsum = 0 - I * 0;
  	for(n=0; n < sett->nifo; ++n) { 

// Calculate detector positions with respect to baricenter
// Copied from jobcore.c, line 248

//shft & phase loop
    		for(i=0; i < sett->N; ++i) {

      			shft[n][i] = nSource[0]*ifo[n].sig.DetSSB[i*3]
		         	+ nSource[1]*ifo[n].sig.DetSSB[i*3+1]
		         	+ nSource[2]*ifo[n].sig.DetSSB[i*3+2];
    
// Phase modulation function
// Copied from jobcore.c, line 265

			phase[i] = sgnlo[0]*(double)(i + shft[n][i]) 
				+ (sgnlo[1]*i*i) + ((sett->oms 
				+ 2.*sgnlo[1]*i)*shft[n][i]);
		} //shft & phase			

//Sin & Cos calculations using Yeppp!

		status = yepMath_Cos_V64f_V64f(phase, _cph, VLEN);
		assert(status == YepStatusOk);
		status = yepMath_Sin_V64f_V64f(phase, _sph, VLEN);
		assert(status == YepStatusOk);

// Matched filter 
// Copied from jobcore.c, line 276 and 337

		for (i = 0; i < sett->N; ++i){

	      		exph = _cph[i] - I * _sph[i];

	      		xDatma[n][i] = ifo[n].sig.xDat[i]*sigaa[n][i]*exph;
	      		xDatmb[n][i] = ifo[n].sig.xDat[i]*sigbb[n][i]*exph;

		}
	} //End of detector loop 1

//Loop for each detector 2

  	for(n=0; n < sett->nifo; ++n) { 
		aatemp = 0.;
		bbtemp = 0.;
		for(i = 0; i < sett->N; i++){

			aatemp += sqr(sigaa[n][i]);
			bbtemp += sqr(sigbb[n][i]);
		}
		for(i=0; i < sett->N; ++i) {
			xDatma[n][i] /= ifo[n].sig.sig2;
			xDatmb[n][i] /= ifo[n].sig.sig2;

			xasum += xDatma[n][i];
			xbsum += xDatmb[n][i];
		}
		aa += aatemp/ifo[n].sig.sig2; 
		bb += bbtemp/ifo[n].sig.sig2;   

	}// End of detector loop 2

// F - statistic
	fstat_out[5] = - ((( sqr(creal(xasum)) + sqr(cimag(xasum)))/aa)
			+ ((sqr(creal(xbsum)) + sqr(cimag(xbsum)))/bb));

// Amplitude estimates
	fstat_out[0] = 2*creal(xasum)/aa;
	fstat_out[1] = 2*creal(xbsum)/bb;
	fstat_out[2] = -2*cimag(xasum)/aa;
	fstat_out[3] = -2*cimag(xbsum)/bb;

// Signal-to-noise ratio
	fstat_out[4] = sqrt(2*(-fstat_out[5]-2));

	fstat_out[6] = sgnlo[0];
	fstat_out[7] = sgnlo[1];
	fstat_out[8] = sgnlo[2];
	fstat_out[9] = sgnlo[3];

// Signal-to-noise ratio from estimated amplitudes (for h0 = 1)
	fstat_out[10] = sqrt(sqr(2*creal(xasum)) + sqr(2*creal(xbsum)) + sqr(2*cimag(xasum)) + sqr(2*cimag(xbsum))); 	

	free(_sph);
	free(_cph);
	free(phase);
	free_matrix(shft, sett->nifo, sett->N);
	free_matrix_complex(xDatma, sett->nifo, sett->N);
	free_matrix_complex(xDatmb, sett->nifo, sett->N);

	return fstat_out;
}



//mesh adaptive direct search (MADS) maximum search declaration
Yep64f* MADS(Search_settings *sett, Aux_arrays *aux, double* in, double *start, double delta, double *pc, int bins){

	int i, j, k, l, m, n, o, r, a = 0;
  	double sinalt, cosalt, sindelt, cosdelt;
	double param[4]; 			//initial size of mesh
	double smallest = 100;
        double **sigaa, **sigbb;   // aa[nifo][N]
        sigaa = matrix(sett->nifo, sett->N);
	sigbb = matrix(sett->nifo, sett->N);
	for(r = 0; r < 4;r++){
		param[r] = pc[r]/bins;
	}
	for(r = 0; r < 4; r++){
		if(param[r] < smallest) smallest = param[r];
	}

#ifdef YEPPP
    yepLibrary_Init();
    Yep64f *p = (Yep64f*)malloc(sizeof(Yep64f)*4);
    Yep64f *out = (Yep64f*)malloc(sizeof(Yep64f)*11);
    Yep64f *nSource = (Yep64f*)malloc(sizeof(Yep64f)*3); 
    Yep64f *extr = (Yep64f*)malloc(sizeof(Yep64f)*11); 
    Yep64f *res = (Yep64f*)malloc(sizeof(Yep64f)*11);
    enum YepStatus status;

#endif
//	puts("MADS");

	for(i = 0; i < 11; i++) extr[i] = in[i];
	while(smallest >= delta){  	//when to end computations
		k = 0;
		for(j = 0; j < 8; j++){
			for(i = 0; i < 4; i++) p[i] = in[6+i];
			if(j < 4){
				p[k] = in[6+k] - start[k]*(param[k] - delta);
				k++;
			}
			else { 
				k--;
				p[k] = in[6+k] + start[k]*(param[k] - delta);
			}
			sinalt = sin(p[3]);
			cosalt = cos(p[3]);
			sindelt = sin(p[2]);
			cosdelt = cos(p[2]);

			nSource[0] = cosalt*cosdelt;
			nSource[1] = sinalt*cosdelt;
			nSource[2] = sindelt;

			for (o = 0; o < sett->nifo; ++o){
				modvir(sinalt, cosalt, sindelt, cosdelt, 
			   		sett->N, &ifo[o], aux, sigaa[o], sigbb[o]);  
			}
			res = Fstatnet(sett, p, nSource, sigaa, sigbb); //Fstat function for mesh points

			if (res[5] < extr[5]){
				for (l = 0; l < 11; l++){ 
					extr[l] = res[l];
				}
			}

		} //j
		if (extr[5] == in[5]){
			smallest = smallest - delta;
			for(r = 0; r < 4; r++){
				param[r] = param[r] - delta;
			}
		}
		else{
			for(m = 0; m < 4; m++) p[m] = extr[6+m];
			for(m = 0; m < 11; m++) in[m] = extr[m];
		}
	} // while
	for (n= 0; n < 11; n++){
		out[n] = extr[n];
	}

	free(p);
	free(nSource);
	free(extr);
	free(res);
	free_matrix(sigaa, sett->nifo, sett->N);
	free_matrix(sigbb, sett->nifo, sett->N);
	return out;
}

// Few functions for Nelder-Mead (simplex) algorithm  

double ** make_simplex(double * point, int dim, double *pc2){
	int i, j;
	double ** simplex = alloc_matrix(dim + 1, dim);
	for (i = 0; i < dim + 1; i++){
		for (j = 0; j < dim; j++){
			simplex[i][j] = point[j];
		}
	}
	for (i = 0; i < dim; i++){
		simplex[i][i] = (1 + pc2[i])*simplex[i][i];
	}
	return simplex;
}

void evaluate_simplex(double ** simplex, int dim, double ** fx, Search_settings *sett, Aux_arrays *aux, double *nS, double **sigaa_max, double **sigbb_max){
	double sinalt, cosalt, sindelt, cosdelt;
	double *out = alloc_vector(11);
	int i, o, j;
	for (i = 0; i < dim + 1; i++){
			sinalt = sin(simplex[i][3]);
			cosalt = cos(simplex[i][3]);
			sindelt = sin(simplex[i][2]);
			cosdelt = cos(simplex[i][2]);

			nS[0] = cosalt*cosdelt;
			nS[1] = sinalt*cosdelt;
			nS[2] = sindelt;

			for (o = 0; o < sett->nifo; ++o){
				modvir(sinalt, cosalt, sindelt, cosdelt, 
			   		sett->N, &ifo[o], aux, sigaa_max[o], sigbb_max[o]);  
			}
			out = Fstatnet(sett, simplex[i], nS, sigaa_max, sigbb_max);
			for (j = 0; j < 11; j++) fx[i][j] = out[j];
	}
	free_vector(out, 11);
}

int* simplex_extremes(double ** fx, int dim){
	int i;
	int ihi, ilo, inhi;
	static int ihe[3];

	if (fx[0][5] > fx[1][5]){ 
		ihi = 0; ilo = inhi = 1; 
	}
	else { 
		ihi = 1; 
		ilo = inhi = 0; 
	}
	for (i = 2; i < dim + 1; i++){
		if (fx[i][5] <= fx[ilo][5]) {
			ilo = i;
		}		
		else if (fx[i][5] > fx[ihi][5]){ 
			inhi = ihi; 
			ihi = i; 
		}
		else if (fx[i][5] > fx[inhi][5]){
			inhi = i;
		}
	}
	ihe[0] = ihi;
	ihe[1] = ilo;
	ihe[2] = inhi;

	return ihe;
}

void simplex_bearings(double ** simplex, int dim, double * midpoint, double * line, int ihi){
	int i, j;
	for (j = 0; j < dim; j++) midpoint[j] = 0.0;
	for (i = 0; i < dim + 1; i++){
		if (i != ihi){
			for (j = 0; j < dim; j++) midpoint[j] += simplex[i][j];
		}
	}
	for (j = 0; j < dim; j++){
		midpoint[j] /= dim;
		line[j] = simplex[ihi][j] - midpoint[j];
	}
}
int update_simplex(double ** simplex, int dim, double  fmax, double ** fx, int ihi, double * midpoint, double * line, double scale, Search_settings *sett, Aux_arrays *aux, double *nS, double **sigaa_max, double **sigbb_max){
	int i, o, j, update = 0; 
	double * next = alloc_vector(dim);
	double * fx2 = alloc_vector(11)	;
	double sinalt, cosalt, sindelt, cosdelt;

	for (i = 0; i < dim; i++) next[i] = midpoint[i] + scale * line[i];

	sinalt = sin(next[3]);
	cosalt = cos(next[3]);
	sindelt = sin(next[2]);
	cosdelt = cos(next[2]);

	nS[0] = cosalt*cosdelt;
	nS[1] = sinalt*cosdelt;
	nS[2] = sindelt;

	for (o = 0; o < sett->nifo; ++o){
		modvir(sinalt, cosalt, sindelt, cosdelt, 
	   		sett->N, &ifo[o], aux, sigaa_max[o], sigbb_max[o]);  
	}
	fx2 = Fstatnet(sett, next, nS, sigaa_max, sigbb_max);
	if (fx2[5] < fmax){
		for (i = 0; i < dim; i++) simplex[ihi][i] = next[i];
		for (j = 0; j < 11; j++) fx[ihi][j] = fx2[j];
		update = 1;
	}

	free_vector(next, dim);
	free_vector(fx2, 11);
	return update;
}

void contract_simplex(double ** simplex, int dim, double ** fx, int ilo, int ihi, Search_settings *sett, Aux_arrays *aux, double *nS, double **sigaa_max, double **sigbb_max){
  	double sinalt, cosalt, sindelt, cosdelt;
	double *fx3= alloc_vector(11);
	int i, j, k, o;

	for (i = 0; i < dim + 1; i++){
		if (i != ilo){
			for (j = 0; j < dim; j++) simplex[i][j] = (simplex[ilo][j]+simplex[i][j])*0.5;
			sinalt = sin(simplex[i][3]);
			cosalt = cos(simplex[i][3]);
			sindelt = sin(simplex[i][2]);
			cosdelt = cos(simplex[i][2]);

			nS[0] = cosalt*cosdelt;
			nS[1] = sinalt*cosdelt;
			nS[2] = sindelt;

			for (o = 0; o < sett->nifo; ++o){
				modvir(sinalt, cosalt, sindelt, cosdelt, 
			   		sett->N, &ifo[o], aux, sigaa_max[o], sigbb_max[o]);  
			}

			fx3 = Fstatnet(sett, simplex[i], nS, sigaa_max, sigbb_max);
			for (k = 0; k < 11; k++) fx[i][k] = fx3[k];
		}
	}
	free_vector(fx3, 11);

}

int check_tol(double fmax, double fmin, double ftol){
	double delta = fabs(fmax - fmin);
	double accuracy = (fabs(fmax) + fabs(fmin)) * ftol;
	return (delta < (accuracy + ZEPS));
}

double * amoeba(Search_settings *sett, Aux_arrays *aux, double *point, double *nS, double *res_max, int dim, double tol, double *pc2, double **sigaa_max, double **sigbb_max){
	int ihi, ilo, inhi;
// ihi = ih[0], ilo = ih[1], inhi = ih[2];
	int *ih;
	int j, i;
	static double NM_out[11];
	double ** fx = alloc_matrix(dim + 1, 11);
	double * midpoint = alloc_vector(dim);
	double * line = alloc_vector(dim);
	double ** simplex = make_simplex(point, dim, pc2);

	evaluate_simplex(simplex, dim, fx, sett, aux, nS, sigaa_max, sigbb_max);

	while (true)
	{
		ih = simplex_extremes(fx, dim);
		ihi = ih[0];
		ilo = ih[1]; 
		inhi = ih[2];
		simplex_bearings(simplex, dim, midpoint, line, ihi);
		if(check_tol(fx[ihi][5], fx[ilo][5], tol)) break;
		update_simplex(simplex, dim, fx[ihi][5], fx, ihi, midpoint, line, -1.0, sett, aux, nS, sigaa_max, sigbb_max);
		if (fx[ihi][5] < fx[ilo][5]){
			update_simplex(simplex, dim, fx[ihi][5], fx, ihi, midpoint, line, -2.0, sett, aux, nS, sigaa_max, sigbb_max);
		}
		else if (fx[ihi][5] > fx[inhi][5]){
			if (!update_simplex(simplex, dim, fx[ihi][5], fx, ihi, midpoint, line, 0.5, sett, aux, nS, sigaa_max, sigbb_max)){
				contract_simplex(simplex, dim, fx, ilo, ihi, sett, aux, nS, sigaa_max, sigbb_max);
			}
		}
	}

	for (j = 0; j < dim; j++) point[j] = simplex[ilo][j];
	for (j = 0; j < 11; j++) NM_out[j] = fx[ilo][j];
	free_matrix(fx, dim + 1, 10);
	free_vector(midpoint, dim);
	free_vector(line, dim);
//	free_vector_int(ih, 3);
	free_matrix(simplex, dim + 1, dim);


/*	free(fx);
	free(midpoint);
	free(line);
	free(simplex);
*/
	return NM_out;
}

// Main programm
int main (int argc, char *argv[]) {

	Search_settings sett;	
	Command_line_opts opts;
  	Aux_arrays aux_arr;
//  	double *F; 			// F-statistic array
  	int i, j, r, c, a, b, g, flag; 	
	int d, o, m, k, s;
	int bins = 2, ROW, dim = 4;	// neighbourhood of point will be divide into defined number of bins
	int gsize = 2;			// grid size where followup will be searching maximum
	int spndr[2], nr[2], mr[2];	// range in linear unities
	int hemi; 			// hemisphere
	double pc[4];			// % define neighbourhood around each parameter for initial grid
	double pc2[4];			// % define neighbourhood around each parameter for direct maximum search (MADS & Simplex)
	double tol = 1e-10;
	double cof, al1, al2;
//	double delta = 1e-5;		// initial step in MADS function
//	double *results;		// Vector with results from Fstatnet function
//	double *maximum;		// True maximum of Fstat
//	double results_max[11];	
	double s1, s2, s3, s4;
	double sgnlo[4];
	double **sgnlo_range;  
	double sgnlol[4]; 
	double be[2];	
  	double *MM ; 	 		// Auxiliary array for grid points	 
	double **arr;			// arr[ROW][COL] - array of points where Fstat will be calculated
	double nSource[3];
  	double sinalt, cosalt, sindelt, cosdelt;
	double F_min;
	double au;
//	char path[512];
	double x, y;
	ROW = pow((bins+1),4);

	arr = matrix(ROW, 4);

#ifdef YEPPP
    	yepLibrary_Init();
   	Yep64f *results_max = (Yep64f*)malloc(sizeof(Yep64f)*11); 
    	Yep64f *results_first = (Yep64f*)malloc(sizeof(Yep64f)*11);
    	Yep64f *results = (Yep64f*)malloc(sizeof(Yep64f)*11);
	Yep64f *maximum = (Yep64f*)malloc(sizeof(Yep64f)*11);
//    Yep64f *sgnlo = (Yep64f*)malloc(sizeof(Yep64f)*4);  
//    Yep64f *nSource = (Yep64f*)malloc(sizeof(Yep64f)*3); 
	Yep64f *mean = (Yep64f*)malloc(sizeof(Yep64f)*4); 

    enum YepStatus status;

#endif

	double sgnlo_max[4];
	double nSource_max[3];
	double **sigaa_max, **sigbb_max;   // aa[nifo][N]

// Time tests
	double tdiff;
	clock_t tstart, tend;

// Command line options 
	handle_opts(&sett, &opts, argc, argv); 
// Search settings
  	search_settings(&sett); 
// Detector network settings
  	detectors_settings(&sett, &opts); 

// Array initialization
  	init_arrays(&sett, &opts, &aux_arr);

	if(opts.simplex_flag){
		sigaa_max = matrix(sett.nifo, sett.N);
		sigbb_max = matrix(sett.nifo, sett.N);
	}
// Grid data
	if(!opts.neigh_flag){
		sgnlo_range = matrix(4, 2);
		read_grid(&sett, &opts);
	}
// Amplitude modulation functions for each detector  
				for(i=0; i<sett.nifo; i++) rogcvir(&ifo[i]); 
// Adding signal from file
  				if(strlen(opts.addsig)) { 
    					add_signal(&sett, &opts, &aux_arr);
  				}

// Output data handling
/*  struct stat buffer;

  if (stat(opts.prefix, &buffer) == -1) {
    if (errno == ENOENT) {
      // Output directory apparently does not exist, try to create one
      if(mkdir(opts.prefix, S_IRWXU | S_IRGRP | S_IXGRP 
          | S_IROTH	| S_IXOTH) == -1) {
	      perror (opts.prefix);
	      return 1;
      }
    } else { // can't access output directory
      perror (opts.prefix);
      return 1;
    }
  }
*/
//	sprintf(path, "%s/candidates_%03d.coi", opts.dtaprefix, opts.ident);
//	sprintf(path, "candidates.coi", opts.ident);

//Glue function
/*	if(strlen(opts.glue)) {
		glue(&opts);

		sprintf(opts.dtaprefix, "./data_total");
		sprintf(opts.dtaprefix, "%s/followup_total_data", opts.prefix); 
		opts.ident = 000;
	}
*/	
	FILE *coi;
	int z;
	if ((coi = fopen(opts.candidates, "r")) != NULL) {
//		while(!feof(coi)) {

/*			if(!fread(&w, sizeof(unsigned short int), 1, coi)) { break; } 
		  	fread(&mean, sizeof(float), 5, coi);
		  	fread(&fra, sizeof(unsigned short int), w, coi); 
		  	fread(&ops, sizeof(int), w, coi);

			if((fread(&mean, sizeof(float), 4, coi)) == 4){
*/
			while(fscanf(coi, "%le %le %le %le", &mean[0], &mean[1], &mean[2], &mean[3]) == 4){
//Time test
//			tstart = clock();

//Frequency shifted from the reference frame to the current frame
				if (opts.refr > 0){
					mean[0] += -2.*mean[1]*(sett.N)*(opts.refr - opts.ident); 
				}

//Function neighbourhood - generating grid around point
//Area around starting point is calculating as a percent from initial values
				if(opts.neigh_flag){
					pc[0] = 0.0003;
					pc[1] = 0.01;
					pc[2] = 0.01;
					pc[3] = 0.01;

					for (i = 0; i < 4; i++){
						pc2[i] = 2*pc[i]/bins;
					}

					neigh(mean, pc, bins, arr);

				} //if neigh
				else{

printf("%d %le %le %le %le\n", sett.N, mean[0], mean[1], mean[2], mean[3]);

					cof = sett.oms + mean[0]; 
					for(i=0; i<2; i++) sgnlol[i] = mean[i]; 
					hemi = ast2lin(mean[3], mean[2], C_EPSMA, be);
 
					sgnlol[2] = be[0]*cof; 
					sgnlol[3] = be[1]*cof;

// solving a linear system in order to translate 
// sky position, frequency and spindown (sgnlo parameters) 
// into the position in the grid


					MM = (double *) calloc (16, sizeof (double));
					for(i=0; i<16; i++) MM[i] = sett.M[i];
					gsl_vector *x = gsl_vector_alloc (4);
					gsl_matrix_view m = gsl_matrix_view_array (MM, 4, 4);
					gsl_matrix_transpose (&m.matrix) ; 
					gsl_vector_view b = gsl_vector_view_array (sgnlol, 4);
					gsl_permutation *p = gsl_permutation_alloc (4);
					 
					gsl_linalg_LU_decomp (&m.matrix, p, &s);
					gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
					  
					spndr[0] = round(gsl_vector_get(x,1)); 
					nr[0] 	= round(gsl_vector_get(x,2));
					mr[0] 	= round(gsl_vector_get(x,3));
					  
					gsl_permutation_free (p);
					gsl_vector_free (x);
					free (MM);

// Check if calculated point is in an appropriate region

					al1 = nr[0]*sett.M[10] + mr[0]*sett.M[14];
					al2 = nr[0]*sett.M[11] + mr[0]*sett.M[15];
  					if ((sqr(al1)+sqr(al2))/sqr(sett.oms) > 1.){ 
						puts("Current candidate is in an inappropriate region! -> going to the next candidate");
						continue;
					}
					else{

// Define the grid range in which the signal will be looked for

						spndr[1] = spndr[0] + gsize; 
						spndr[0] -= gsize;
						nr[1] = nr[0] + gsize; 
						nr[0] -= gsize;
						mr[1] = mr[0] + gsize; 
						mr[0] -= gsize;

// Go back to astrophysical units - calculate range for calculations

						flag = 0;

// Loop over min/max values
						for(i = 0; i < 2; i++){
//printf("nr = %d %d; mr = %d %d; i = %d\n", nr[0], nr[1], mr[0], mr[1], i);

							al1 = nr[i]*sett.M[10] + mr[i]*sett.M[14];
							al2 = nr[i]*sett.M[11] + mr[i]*sett.M[15];
lin2ast(al1/sett.oms, al2/sett.oms, hemi, sett.sepsm, sett.cepsm, &sinalt, &cosalt, &sindelt, &cosdelt);
printf("%le %le %le\n", spndr[i]*sett.M[5] + nr[i]*sett.M[9] + mr[i]*sett.M[13], asin(sindelt), fmod(atan2(sinalt, cosalt) + 2.*M_PI, 2.*M_PI));

// check if the range is in an appropriate region of the grid
// if not - go to the next candidate 

	  						if ((sqr(al1)+sqr(al2))/sqr(sett.oms) > 1.){ 
								puts("Inappropriate region of the range! -> changing range");
								if (i == 0){
									nr[i] = nr[i] + 1;
									mr[i] = mr[i] + 1;
									i = -1;
								}
								if (i == 1){
									nr[i] = nr[i] - 1;
									mr[i] = mr[i] - 1;
									i = 0;
								}
								if ((nr[1] == nr[0])&&(mr[1] == mr[0])){
									flag = 1;
									continue;
								}
							}
							else{

								lin2ast(al1/sett.oms, al2/sett.oms, hemi, sett.sepsm, 
									sett.cepsm, &sinalt, &cosalt, &sindelt, &cosdelt);

								sgnlo_range[2][i] = asin(sindelt);
								sgnlo_range[3][i] = fmod(atan2(sinalt, cosalt) + 2.*M_PI, 2.*M_PI);
								sgnlo_range[1][i] = spndr[i]*sett.M[5] + nr[i]*sett.M[9] + mr[i]*sett.M[13];
							}
						
						} //for min/max values

						if (flag == 0){
//							sgnlo_range[0][0] = 0.9997*mean[0];
//							sgnlo_range[0][1] = 1.0003*mean[0];
							sgnlo_range[0][0] = mean[0] - 0.03;
							sgnlo_range[0][1] = mean[0] + 0.03;
							for (i = 0; i < 4; i++){ 
								if (sgnlo_range[i][0] > sgnlo_range[i][1]){
									au = sgnlo_range[i][0];
									sgnlo_range[i][0] = sgnlo_range[i][1];
									sgnlo_range[i][1] = au;
								}

								printf("range[%d][0] = %le, range[%d][1] = %le\n", i, sgnlo_range[i][0], i, sgnlo_range[i][1]);
							}
							neigh_from_range(sgnlo_range, bins, arr);

//							pc2[0] = 0.0006/bins;
							for (j = 0; j < 4; j++) pc2[j] = (sgnlo_range[j][1] - sgnlo_range[j][0])/bins;
						}

					} //if candidate is in good region
				} //if not neigh

				if(flag == 1) continue;

// Output data handling
/*  				struct stat buffer;

  				if (stat(opts.prefix, &buffer) == -1) {
    					if (errno == ENOENT) {
// Output directory apparently does not exist, try to create one
      						if(mkdir(opts.prefix, S_IRWXU | S_IRGRP | S_IXGRP 
          						| S_IROTH	| S_IXOTH) == -1) {
	      						perror (opts.prefix);
	      						return 1;
      						}
    					} 
					else { // can't access output directory
			      			perror (opts.prefix);
			      			return 1;
			    		}
  				}
*/
		
// Setting number of using threads (not required)
//omp_set_num_threads(1);


				results_max[5] = 0.;

// Main loop - over all parameters + parallelisation
#pragma omp parallel default(shared) private(d, i, sgnlo, sinalt, cosalt, sindelt, cosdelt, nSource, results, maximum)
{

                       		double **sigaa, **sigbb;   // old aa[nifo][N], bb[nifo][N]
              			sigaa = matrix(sett.nifo, sett.N);
				sigbb = matrix(sett.nifo, sett.N);

#pragma omp for  
				for (d = 0; d < ROW; ++d){

					for (i = 0; i < 4; i++){
						sgnlo[i] = arr[d][i];
					}
 
					sinalt = sin(sgnlo[3]);
					cosalt = cos(sgnlo[3]);
					sindelt = sin(sgnlo[2]);
					cosdelt = cos(sgnlo[2]);

					nSource[0] = cosalt*cosdelt;
					nSource[1] = sinalt*cosdelt;
					nSource[2] = sindelt;

					for (i = 0; i < sett.nifo; ++i){
						modvir(sinalt, cosalt, sindelt, cosdelt, 
					   		sett.N, &ifo[i], &aux_arr, sigaa[i], sigbb[i]);  
					}
// F-statistic in given point

					results = Fstatnet(&sett, sgnlo, nSource, sigaa, sigbb);
//printf("%le %le %le %le %le %le\n", results[6], results[7], results[8], results[9], results[5], results[4]);

// Check is it the biggest found value of F-statistic
#pragma omp critical
					if(results[5] < results_max[5]){
						for (i = 0; i < 11; i++){
							results_max[i] = results[i];
						}
						if(opts.simplex_flag){
							for (i = 0; i < 4; i ++){
								sgnlo_max[i] = sgnlo[i];
							}
							for (i = 0; i < 3; i++){
								nSource_max[i] = nSource[i];
							}
							for (g = 0; g < sett.nifo; g++){
								for (i = 0; i < sett.N; i++){
									sigaa_max[g][i] = sigaa[g][i];
									sigbb_max[g][i] = sigbb[g][i];
								}
							}
						} //if --simplex

					} // if results < results_max


				} // d - main outside loop
				free_matrix(sigaa, sett.nifo, sett.N);
				free_matrix(sigbb, sett.nifo, sett.N);

} //pragma

				for(g = 0; g < 11; g++) results_first[g] = results_max[g];

// Maximum search using MADS algorithm
  				if(opts.mads_flag) {
//					puts("MADS");
					maximum = MADS(&sett, &aux_arr, results_max, mean, tol, pc2, bins);
				}

// Maximum search using simplex algorithm
				if(opts.simplex_flag){
					puts("Simplex");
					maximum = amoeba(&sett, &aux_arr, sgnlo_max, nSource_max, results_max, dim, tol, pc2, sigaa_max, sigbb_max);
//printf("Amoeba: %le %le %le %le %le %le\n", maximum[6], maximum[7], maximum[8], maximum[9], maximum[5], maximum[4]);

// Maximum value in points searching
					if(maximum[5] < results_max[5]){
						for (i = 0; i < 11; i++){
							results_max[i] = maximum[i];
						}
					}

				} //simplex

//Time test
//				tend = clock();
//				tdiff = (tend - tstart)/(double)CLOCKS_PER_SEC;
				puts("Maximum:");
				printf("%le %le %le %le %le %le\n", results_max[6], results_max[7], results_max[8], results_max[9], results_max[5], results_max[4]);



			} // while fread coi
//		}
	} //if coi
	else {
		
		printf ("Problem with %s\n", opts.candidates);
		return 1;
	}

// Output information
/*	puts("**********************************************************************");
	printf("***	Maximum value of F-statistic for grid is : (-)%.8le	***\n", -results_first[5]);
	printf("Sgnlo: %.8le %.8le %.8le %.8le\n", results_first[6], results_first[7], results_first[8], results_first[9]);
	printf("Amplitudes: %.8le %.8le %.8le %.8le\n", results_first[0], results_first[1], results_first[2], results_first[3]);
	printf("Signal-to-noise ratio: %.8le\n", results_first[4]); 
	printf("Signal-to-noise ratio from estimated amplitudes (for h0 = 1): %.8le\n", results_first[10]);
	puts("**********************************************************************");
if((opts.mads_flag)||(opts.simplex_flag)){
	printf("***	True maximum is : (-)%.8le				***\n", -maximum[5]);
	printf("Sgnlo for true maximum: %.8le %.8le %.8le %.8le\n", maximum[6], maximum[7], maximum[8], maximum[9]);
	printf("Amplitudes for true maximum: %.8le %.8le %.8le %.8le\n", maximum[0], maximum[1], maximum[2], maximum[3]);
	printf("Signal-to-noise ratio for true maximum: %.8le\n", maximum[4]); 
	printf("Signal-to-noise ratio from estimated amplitudes (for h0 = 1) for true maximum: %.8le\n", maximum[10]);
	puts("**********************************************************************");
}*/
// Cleanup & memory free 
	free(results_max);
	free(results_first);
	free(results);
//	free(maximum);
	free(mean);
	free_matrix(sgnlo_range, 4, 2);
	free_matrix(arr, ROW, 4);
	if(opts.simplex_flag){
		free_matrix(sigaa_max, sett.nifo, sett.N);
		free_matrix(sigbb_max, sett.nifo, sett.N);
	}
  	cleanup_followup(&sett, &opts, &aux_arr);
	return 0;
}


