/* 

***************************************************************
Main followup programm. Code starts from the given point
and investigate space around it (calculating F-statistics
for a given frequency, spindown and sky position), looking
for maximum of F-statistics.
***************************************************************
One (best so far) option is calculating F-statistics on the
very precise, very dense grid (Fisher matrix is reading from
file), choose point with the higher value of F-statistics and
from this point start invertedMADS algorithm (other options 
of maximum search are 'normal' MADS and Simplex).
***************************************************************
Big loops (on grid points and inside invertedMADS algorithm) 
use OpenMP parallelisation.
***************************************************************

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
#include <getopt.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>
#include <time.h>
#include <dirent.h>
#include <omp.h>

#include "auxi.h"
#include "settings.h"
#include "struct.h"
#include "init.h"

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

//Accuracy in Simplex algorithm
#define ZEPS 1e-14


/***************************************************************
Allocation of memory for vector with given number of columns
***************************************************************/
double * alloc_vector(int cols){
	return (double *) malloc(sizeof(double) * cols);
} //end alloc_vector()

/***************************************************************
Allocation of memory for double martix with given number of rows and columns
***************************************************************/
double** matrix(int rows, int cols) {
  	int k;
	double **m;
  	m = (double **)malloc(rows*sizeof(double *));
  
  	for (k=0; k < rows; k++){
    		m[k] = (double *)calloc(cols, sizeof(double));  
	}
  	return m;
} //end matrix()

/***************************************************************
Allocation of memory for complex double martix with given number of rows and columns
***************************************************************/
complex double** matrix_complex(int rows, int cols) {
  	int k;
	complex double **m;
  	m = (complex double **)malloc(rows*sizeof(complex double *));
  
  	for (k=0; k < rows; k++){
    		m[k] = (complex double *)calloc(cols, sizeof(complex double));
	}  
  	return m;
} //end matrix_complex()

/***************************************************************
Free memory for allocated vectors 
***************************************************************/
void free_vector(double * vector, int cols){
	free(vector);
} //end free_vector()

/***************************************************************
Free memory for allocated matricies 
***************************************************************/
void free_matrix(double **matrix, int rows, int cols){
	int i;
	for (i = 0; i < rows; i++){ 
		free(matrix[i]);
	}
	free(matrix);
} //end free_matrix

/***************************************************************
Free memory for allocated complex matricies 
***************************************************************/
void free_matrix_complex(complex double **matrix, int rows, int cols){
	int i;
	for (i = 0; i < rows; i++){ 
		free(matrix[i]);
	}
	free(matrix);
}//end free_matrix_complex()

/***************************************************************
Generate matrix which will be used for 4D optimal grid
***************************************************************/
void A4opt(double minimalm, int no, double *gam, double **Mopt){
	int i, j, k, l;
	double r; // covering radius
	r = sqrt(1.-(minimalm*minimalm));

	double V[4][4], D[4];
  	double Mg[4][4];
  	double Mg_t[4][4];
	double Moptn[4][4];
	double temp[16];

	double datam[4][4] = { 	{r*sqrt(5.),		0.,			0.,			0.},
        			{r*sqrt(5.)/2.,		r*sqrt(15.)/2.,		0.,			0.},
        			{r*sqrt(5.)/2.,		r*sqrt(5./3.)/2.,	r*sqrt(10./3.),		0.},
       				{-r*sqrt(5.)/2.,	-r*sqrt(5./3.)/2.,	-r*sqrt(5./6.)/2.,	-r/(2.*sqrt(2.))}};

//Fill with zeros
	for (i = 0; i < 4; i++){ 
		D[i] = 0.;
		for (j = 0; j < 4; j++){ 
			V[i][j] = 0.;
			Mg[i][j] = 0.;
			Mg_t[i][j] = 0.;
			Moptn[i][j] = 0.;
		}
	}

//few GSL function to find eigenvalues and eigenvectors  	
  	gsl_matrix_view gamgsl = gsl_matrix_view_array (gam, 4, 4);
  	gsl_vector *eval = gsl_vector_alloc (4);
  	gsl_matrix *evec = gsl_matrix_alloc (4, 4);
  	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (4);
	gsl_eigen_symmv (&gamgsl.matrix, eval, evec, w);
	gsl_eigen_symmv_free (w);

	for (i = 0; i < 4; i++){
// D is diagonal matrix with eigenvalues of the Fisher matrix 
		D[i] = sqrt(fabs(1./gsl_vector_get(eval, i)));
		gsl_vector_view evec_i = gsl_matrix_column(evec, i);
		for (k = 0; k < 4; k++){
// V is matrix in which columns are eigenvectors of the Fisher matrix
			V[k][i] = gsl_vector_get(&evec_i.vector, k);
		}

	}

//Free memory
	gsl_vector_free (eval);
	gsl_matrix_free (evec);

// V * sqrt(abs(D^{-1})) -> Mg

  	for(i=0; i<4; i++) for(j=0; j<4; j++){ 
		Mg[i][j] = V[i][j]*D[j]; 
	}
	k=0;
	for(i=0; i<4; i++){ 
		for(j=0; j<4; j++){
			temp[k] = Mg[i][j];
			k++;
		}
	}
// Transpose of the Mg matrix with GSL -> Mg^T
  	gsl_matrix_view Mggsl = gsl_matrix_view_array (temp, 4, 4);
	gsl_matrix_transpose(&Mggsl.matrix);
	for(i=0; i<4; i++){ 
		for(j=0; j<4; j++){ 
			Mg_t[i][j] = gsl_matrix_get(&Mggsl.matrix, i, j); 
		}
	}
  	for(i=0; i<4; i++){ 
    		for(j=0; j<4; j++) { 
      			for(k=0; k<4; k++){
//Final matrix: datam * Mg^T -> Mopt
				Mopt[i][j] += datam[i][k]*Mg_t[k][j];
			}
//Normalisation
     			Mopt[i][j] /= no; 
    		}
	}
// Special normalisation for spindown 
  	for(i=0; i<4; i++) Mopt[i][1] /= no; 

} //end A4opt()

/***************************************************************
Use matrix to generate points on the optimal 4D grid
***************************************************************/
void InjLoc(int *gri1, double *sl, double ** mx){
	double vec[16];
	double inva[16];
	double nmso[4];
	double Mp_t[4][4], Mp_i[4][4];
	int s, i, j, k=0;
	for (i = 0; i < 4; i++){
		for (j = 0; j < 4; j++){ 
			vec[k] = mx[i][j];
			k++; 
		}
	}
//Mp = Mopt
//Transposition Mp^T=Mopt^T with GSL
	gsl_matrix_view Mpgsl = gsl_matrix_view_array (vec, 4, 4);
	gsl_matrix_transpose(&Mpgsl.matrix);
	for(i=0; i<4; i++){ 
		for(j=0; j<4; j++){ 
			Mp_t[i][j] = gsl_matrix_get(&Mpgsl.matrix, i, j);
		} 
	}
	gsl_matrix_view inv = gsl_matrix_view_array(inva,4,4);
	gsl_permutation * p = gsl_permutation_alloc (4);
//Inversion Mp^{-1}
	gsl_linalg_LU_decomp (&Mpgsl.matrix, p, &s);
	gsl_linalg_LU_invert (&Mpgsl.matrix, p, &inv.matrix);
	for(i=0; i<4; i++){ 
		for(j=0; j<4; j++){ 
			Mp_i[i][j] = gsl_matrix_get(&inv.matrix,i,j); 
		}
	}
	for (i = 0; i < 4; i++){ 
		nmso[i] = 0.;
	}
//sl=sgnlol - this is starting point, but the sky position are rescaled to the
//linear coordinates
//nmso = Mp^{-1}*sgnlol
	for (i = 0; i < 4; i++){
		for (j = 0; j < 4; j++){
			nmso[i] += inva[4*i+j]*sl[j];
		}
	}
//gri1 is the closest point on the optimal grid to the starting point sgnlol
	for (i = 0; i < 4; i++){ 
		gri1[i] = round(nmso[i]);
	}
	gsl_permutation_free (p);
} //end InjLoc()

/***************************************************************
//Function neigh takes candidate parameters and number of bins (as arguments) 
//and creates uniform grid around it.
//Range is defined as a % value from the candidate's parameters
***************************************************************/
void neigh(double *m, double *perc, int b, double **arr){ 
	int rows, cols = 4;
	rows = pow((b+1),4);
	int k;
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
} //end neigh()

/***************************************************************
Similar as neigh() function, but takes +/- values instead of %
***************************************************************/
void neigh2(double *m, double *perc, int b, double **arr){ 
	int rows, cols = 4;
	rows = pow((b+1),4);
	int k;
	double beg[4];
	double width[4];
	int i1, i2, i3, i4, j, i;
	for (j = 0; j < 4; j++) {
		width[j] = 2*perc[j]/b;
	}
	i = 0;
	beg[0] = m[0] - perc[0];
	for (i1 = 0; i1 < (b + 1); i1++){
		beg[1] = m[1] - perc[1];
		for (i2 = 0; i2 < (b + 1); i2++){
			beg[2] = m[2] - perc[2];
			for (i3 = 0; i3 < (b + 1); i3++){
				beg[3] = m[3] - perc[3];
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
} //end neigh2()

/***************************************************************
Function takes calculated range around point (from grid.bin file) and creates grid around it.
This grid is not such precise as in the case of functions A4opt and InjLoc due to the
linear approximation used to generate files grid.bin
***************************************************************/
void neigh_from_range(double **range, int b, double **arr){

	int rows, cols = 4;
	rows = pow((b+1),4);
	int k;
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
} // end neigh_from_range()

/***************************************************************
Function calculates F-statistics for a given point (frequency, spindown, sky position)
***************************************************************/
double* Fstatnet(Search_settings *sett, double *sgnlo, double *nSource, double **sigaa, double **sigbb){
	int i = 0, n = 0; 
	double aatemp, bbtemp, aa = 0., bb = 0.;
	complex double exph, xasum, xbsum;
	double **shft; 						// Shift of the phase
        shft = matrix(sett->nifo, sett->N);
	complex double **xDatma, **xDatmb;			// Matched filter
	xDatma = matrix_complex(sett->nifo, sett->N);
	xDatmb = matrix_complex(sett->nifo, sett->N);		


#ifdef YEPPP
	int VLEN = sett->N;					// Number of data points
	yepLibrary_Init();
	Yep64f *_sph = (Yep64f*)malloc(sizeof(Yep64f)*VLEN);	// Sines
	Yep64f *_cph = (Yep64f*)malloc(sizeof(Yep64f)*VLEN);	// Cosines
	Yep64f *phase = (Yep64f*)malloc(sizeof(Yep64f)*VLEN); 	// Phase
	Yep64f *fstat_out = (Yep64f*)malloc(sizeof(Yep64f)*11); // Result vector
	enum YepStatus status;
#endif

	xasum = 0 - I * 0;
	xbsum = 0 - I * 0;
//Loop for each detector 1
  	for(n=0; n < sett->nifo; ++n) { 

// Calculate detector positions with respect to baricenter
    		for(i=0; i < sett->N; ++i) {
// Shift function
      			shft[n][i] = nSource[0]*ifo[n].sig.DetSSB[i*3]
		         	+ nSource[1]*ifo[n].sig.DetSSB[i*3+1]
		         	+ nSource[2]*ifo[n].sig.DetSSB[i*3+2];
    
// Phase modulation function
			phase[i] = sgnlo[0]*(double)(i + shft[n][i]) 
				+ (sgnlo[1]*i*i) + ((sett->oms 
				+ 2.*sgnlo[1]*i)*shft[n][i]);
		} //shft & phase loop			

//Sin & Cos calculations using Yeppp!

		status = yepMath_Cos_V64f_V64f(phase, _cph, VLEN);
		assert(status == YepStatusOk);
		status = yepMath_Sin_V64f_V64f(phase, _sph, VLEN);
		assert(status == YepStatusOk);

// Matched filter 
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

// F - statistic value
	fstat_out[5] = - ((( sqr(creal(xasum)) + sqr(cimag(xasum)))/aa)
			+ ((sqr(creal(xbsum)) + sqr(cimag(xbsum)))/bb));
// Amplitude estimates
	fstat_out[0] = 2*creal(xasum)/aa;
	fstat_out[1] = 2*creal(xbsum)/bb;
	fstat_out[2] = -2*cimag(xasum)/aa;
	fstat_out[3] = -2*cimag(xbsum)/bb;
// Signal-to-noise ratio
	fstat_out[4] = sqrt(2*(-fstat_out[5]-2));
// Initial signal
	fstat_out[6] = sgnlo[0];
	fstat_out[7] = sgnlo[1];
	fstat_out[8] = sgnlo[2];
	fstat_out[9] = sgnlo[3];
// Signal-to-noise ratio from estimated amplitudes (for h0 = 1)
	fstat_out[10] = sqrt(sqr(2*creal(xasum)) + sqr(2*creal(xbsum)) + sqr(2*cimag(xasum)) + sqr(2*cimag(xbsum))); 	
//Free memory
	free(_sph);
	free(_cph);
	free(phase);
	free_matrix(shft, sett->nifo, sett->N);
	free_matrix_complex(xDatma, sett->nifo, sett->N);
	free_matrix_complex(xDatmb, sett->nifo, sett->N);
	return fstat_out;
	free(fstat_out);
}//end Fstatnet()

/***************************************************************
Inverted MADS (mesh adaptive direct search) maximum search algorithm: 
starts from initial point and create tiny 4D grid around it (points 
of the grid are in vertices, edges, faces and cells of the 4D hypercube). 
Meshes of the grid increase in every step and shrink if new local maximum 
is found. This method gives better results than 'normal' MADS and it is faster. 
***************************************************************/
Yep64f* invertedMADS(Search_settings *sett, Aux_arrays *aux, double* sgnl, double* rslts, double *p2){
	int i, j, k, l, m, n, o, r, a = 0;
	int count=0, limit = 4000;		//maximal number of steps to avoid infinite loops
  	double sinalt, cosalt, sindelt, cosdelt;
	double fp, fm, sp, sm, dp, dm, ap, am;
	double paramfin[4]; 			//final size of mesh
	double paramstart[4];			//initial size of mesh
	double param[4];
	double array[81][4];			//points on the hypercube
	double p[4];
	double nSource[3];
	double incr[4];				//how much meshes of the grid increase size
	incr[0] = 5e-5;				//in every step, in every direction
	incr[1] = 5e-11;			//values were put by hand; based on tests
	incr[2] = incr[3] = 5e-5;
	for(r = 0; r < 4;r++){
		paramfin[r] = p2[r];
		paramstart[r] = param[r] = incr[r];
	}
#ifdef YEPPP
	yepLibrary_Init();
	Yep64f *out = (Yep64f*)malloc(sizeof(Yep64f)*11);
	Yep64f *extr = (Yep64f*)malloc(sizeof(Yep64f)*11); 
	Yep64f *res = (Yep64f*)malloc(sizeof(Yep64f)*11);
	Yep64f *maxres = (Yep64f*)malloc(sizeof(Yep64f)*11);
	enum YepStatus status;

#endif

	double **saa, **sbb;   		// amplitude modulation factors
	saa = matrix(sett->nifo, sett->N);
	sbb = matrix(sett->nifo, sett->N);

//Main invertedMADS loop: when to end computations
	while((param[0] <= paramfin[0]) && (param[1] <= paramfin[1]) && (param[2] <= paramfin[2]) && (param[3] <= paramfin[3])){  	
		count++;
		for (i = 0; i < 11; i++){ 
			maxres[i] = extr[i] = rslts[i];
		}
//Seed of the step
		for (k = 0; k < 4; k++){
			array[0][k] = extr[6+k];
		}
		fp = extr[6] + param[0];
		fm = extr[6] - param[0];
		sp = extr[7] + param[1];
		sm = extr[7] - param[1];
		dp = extr[8] + param[2];
		dm = extr[8] - param[2];
		ap = extr[9] + param[3];
		am = extr[9] - param[3];
//Construction of the points on the 4D hypercube
//Vertices 
		for (j = 1; j < 9; j++) array[j][0] = fp;
		for (j = 9; j < 17; j++) array[j][0] = fm;

		for (j = 1; j < 5; j++) array[j][1] = array[j+8][1] = sp;
		for (j = 5; j < 9; j++) array[j][1] = array[j+8][1] = sm;

		for (j = 1; j < 3; j++) array[j][2] = array[j+4][2] = array[j+8][2] = array[j+12][2] = dp;
		for (j = 3; j < 5; j++) array[j][2] = array[j+4][2] = array[j+8][2] = array[j+12][2] = dm;

		for (j = 1; j < 17; j = j + 2) array[j][3] = ap;
		for (j = 2; j < 17; j = j + 2) array[j][3] = am;

//Edges
		for (j = 17; j < 25; j++) array[j][0] = extr[6];
		for (j = 25; j < 33; j++) array[j][1] = extr[7];
		for (j = 33; j < 41; j++) array[j][2] = extr[8];
		for (j = 41; j < 49; j++) array[j][3] = extr[9];

		for (j = 25; j < 29; j++) array[j][0] = array[j+8][0] = array[j+16][0] = fp;
		for (j = 29; j < 33; j++) array[j][0] = array[j+8][0] = array[j+16][0] = fm;

		for (j = 17; j < 21; j++) array[j][1] = sp;
		for (j = 21; j < 25; j++) array[j][1] = sm;
		for (j = 33; j < 35; j++) array[j][1] = array[j+4][1] = array[j+8][1] = array[j+12][1] = sp;
		for (j = 35; j < 37; j++) array[j][1] = array[j+4][1] = array[j+8][1] = array[j+12][1] = sm;

		for (j = 17; j < 19; j++) array[j][2] = array[j+4][2] = array[j+8][2] = array[j+12][2] = dp;
		for (j = 19; j < 21; j++) array[j][2] = array[j+4][2] = array[j+8][2] = array[j+12][2] = dm;
		for (j = 41; j < 49; j = j + 2) array[j][2] = dp;
		for (j = 42; j < 49; j = j + 2) array[j][2] = dm;

		for (j = 17; j < 41; j = j + 2) array[j][3] = ap;
		for (j = 18; j < 41; j = j + 2) array[j][3] = am;
		
//Faces
		for (j = 49; j < 61; j++) array[j][0] = extr[6];
		for (j = 49; j < 53; j++) array[j][1] = array[j+12][1] = array[j+18][1] = extr[7];
		for (j = 53; j < 55; j++) array[j][2] = array[j+4][2] = array[j+8][2] = array[j+12][2] = array[j+14][2] = array[j+18][2] = extr[8];
		for (j = 55; j < 57; j++) array[j][3] = array[j+4][3] = array[j+8][3] = array[j+10][3] = array[j+14][3] = array[j+16][3] = extr[9];

		for (j = 61; j < 67; j++) array[j][0] = fp;
		for (j = 67; j < 73; j++) array[j][0] = fm;
		
		for (j = 53; j < 57; j++) array[j][1] = sp;
		for (j = 57; j < 61; j++) array[j][1] = sm;
		array[65][1] = array[71][1] = sp;
		array[66][1] = array[72][1] = sm;

		for (j = 49; j < 51; j++) array[j][2] = dp;
		for (j = 51; j < 53; j++) array[j][2] = dm;
		for (j = 55; j < 65; j = j + 4) array[j][2] = dp;
		for (j = 56; j < 65; j = j + 4) array[j][2] = dm;
		array[69][2] = dp;
		array[70][2] = dm;

		for (j = 49; j < 55; j = j + 2) array[j][3] = ap;
		for (j = 50; j < 55; j = j + 2) array[j][3] = am;
		array[57][3] = array[61][3] = array[67][3] = ap; 
		array[58][3] = array[62][3] = array[68][3] = am; 
//Cells
		for (j = 73; j < 79; j++) array[j][0] = extr[6];
		for (j = 75; j < 81; j++) array[j][1] = extr[7];
		for (j = 73; j < 77; j++) array[j][2] = extr[8];
		for (j = 79; j < 81; j++) array[j][2] = extr[8];
		for (j = 73; j < 75; j++) array[j][3] = extr[9];
		for (j = 77; j < 81; j++) array[j][3] = extr[9];

		array[79][0] = fp;
		array[80][0] = fm;

		array[73][1] = sp;
		array[74][1] = sm;

		array[77][2] = dp;
		array[78][2] = dm;

		array[75][3] = ap;
		array[76][3] = am;

//Calculate F-statistics for the seed of the step
		for(i = 0; i < 4; i++){ 
			p[i] = array[0][i];
		}
		sinalt = sin(p[3]);
		cosalt = cos(p[3]);
		sindelt = sin(p[2]);
		cosdelt = cos(p[2]);

		nSource[0] = cosalt*cosdelt;
		nSource[1] = sinalt*cosdelt;
		nSource[2] = sindelt;
		for (o = 0; o < sett->nifo; ++o){
			modvir(sinalt, cosalt, sindelt, cosdelt, sett->N, &ifo[o], aux, saa[o], sbb[o]);  
		}
		extr = Fstatnet(sett, p, nSource, saa, sbb);
		for (i = 0; i < 11; i++) maxres[i] = extr[i];

//Calculate F-statistics for the rest points - parallel
		#pragma omp parallel default(shared) private(j, o, p, i, sinalt, cosalt, sindelt, cosdelt, nSource, res)
		{
			double **sigaa, **sigbb;   
			sigaa = matrix(sett->nifo, sett->N);
			sigbb = matrix(sett->nifo, sett->N);
			#pragma omp for
			for(j = 1; j < 81; j++){
				for(i = 0; i < 4; i++){ 
					p[i] = array[j][i];
				}

				sinalt = sin(p[3]);
				cosalt = cos(p[3]);
				sindelt = sin(p[2]);
				cosdelt = cos(p[2]);

				nSource[0] = cosalt*cosdelt;
				nSource[1] = sinalt*cosdelt;
				nSource[2] = sindelt;
				for (o = 0; o < sett->nifo; ++o){
					modvir(sinalt, cosalt, sindelt, cosdelt, sett->N, &ifo[o], aux, sigaa[o], sigbb[o]);  
				}
				res = Fstatnet(sett, p, nSource, sigaa, sigbb); 
//Compare points and find one with maximal value of F-statistics
				#pragma omp critical
				if (res[5] < maxres[5]){
					for (i = 0; i < 11; i++){ 
						maxres[i] = res[i];
					}
				}
			} //loop over points
//Free memory
			free_matrix(sigaa, sett->nifo, sett->N);
			free_matrix(sigbb, sett->nifo, sett->N);
		} //pragma
//Compare with the result from the previous step
		if(extr[5] <= maxres[5]){
//If old result is better - increase size of the meshes
			for(r = 0; r < 4; r++){ 
				param[r] = param[r] + incr[r];
			}
			for(j = 0; j < 11; j++){ 
				rslts[j] = extr[j];
			}
		}
		else{
//If new result is better - take it as initial point to the next step
//and shrink meshes to the original size
			for(j = 0; j < 11; j++){ 
				rslts[j] = maxres[j];
			}
			for(r = 0; r < 4; r++){ 
				param[r] = paramstart[r];
			}
		}
//Check if number of steps is not exceeded - if yes, stop while loop
		if(count >= limit){
			for(r = 0; r < 4; r++) param[r] = 100;
		}

	} // while loop

	for (n= 0; n < 11; n++){
		out[n] = rslts[n];
	}

//Free memory
	free(extr);
	free(res);
	free_matrix(saa, sett->nifo, sett->N);
	free_matrix(sbb, sett->nifo, sett->N);
//Print info about number of steps
	printf("invertedMADS ends in %d steps.\n", count);
	return out;
}// end invertedMADS()

/***************************************************************
Directed (inverted) MADS (mesh adaptive direct search) maximum search algorithm: 
starts from initial point and create tiny 2D grid around it (only in frequency
and spindown parameters - keep sky position same as in initial point). 
Meshes of the grid increase in every step and shrink if new local maximum 
is found. 
***************************************************************/

Yep64f* skyMADS(Search_settings *sett, Aux_arrays *aux, double* sgnl, double* rslts, double *p0){
	int i, j, k, l, m, n, o, r, a = 0;
	int count=0, limit = 10000;		//maximal number of steps to avoid infinite loops
  	double sinalt, cosalt, sindelt, cosdelt;
	double fp, fm, sp, sm, fp2, fm2, sp2, sm2;
	double paramfin[2]; 			//final size of mesh
	double paramstart[2];			//initial size of mesh
	double param[2];
	double array[17][2];			//points on the square
	double p[4];
	double nSource[3];
	double incr[2];				//how much meshes of the grid increase size
	incr[0] = 1e-5;				//in every step, in every direction
	incr[1] = 1e-12;			//values were put by hand; based on tests
	for(r = 0; r < 2;r++){
		paramfin[r] = p0[r];
		paramstart[r] = param[r] = incr[r];
	}
#ifdef YEPPP
	yepLibrary_Init();
	Yep64f *out = (Yep64f*)malloc(sizeof(Yep64f)*11);
	Yep64f *extr = (Yep64f*)malloc(sizeof(Yep64f)*11); 
	Yep64f *res = (Yep64f*)malloc(sizeof(Yep64f)*11);
	Yep64f *maxres = (Yep64f*)malloc(sizeof(Yep64f)*11);
	enum YepStatus status;

#endif
	double **saa, **sbb;   		// amplitude modulation factors
	saa = matrix(sett->nifo, sett->N);
	sbb = matrix(sett->nifo, sett->N);

//Main invertedMADS loop: when to end computations
	while((param[0] <= paramfin[0]) && (param[1] <= paramfin[1])){  	
		count++;
		for (i = 0; i < 11; i++){ 
			maxres[i] = extr[i] = rslts[i];
		}
//Seed of the step
		for (k = 0; k < 2; k++){
			array[0][k] = extr[6+k];
		}
		fp = extr[6] + param[0];
		fm = extr[6] - param[0];
		sp = extr[7] + param[1];
		sm = extr[7] - param[1];
		fp2 = extr[6] + 0.5*param[0];
		fm2 = extr[6] - 0.5*param[0];
		sp2 = extr[7] + 0.5*param[1];
		sm2 = extr[7] - 0.5*param[1];

//Construction of the points on the square
		for (j = 1; j < 6; j++) array[j][0] = fp;
		for (j = 6; j < 11; j++) array[j][0] = fm;
		for (j = 11; j < 13; j++) array[j][0] = fp2;
		for (j = 13; j < 15; j++) array[j][0] = fm2;
		for (j = 15; j < 17; j++) array[j][0] = extr[6];

		for (j = 1; j < 6; j++) array[j][1] = sp;
		for (j = 6; j < 11; j++) array[j][1] = sm;
		for (j = 11; j < 13; j++) array[j][1] = sp2;
		for (j = 13; j < 15; j++) array[j][1] = sm2;
		for (j = 15; j < 17; j++) array[j][1] = extr[7];

//Calculate F-statistics for the seed of the step
		for(i = 0; i < 2; i++){ 
			p[i] = array[0][i];
		}
		p[2] = rslts[8];
		p[3] = rslts[9];		

		sinalt = sin(p[3]);
		cosalt = cos(p[3]);
		sindelt = sin(p[2]);
		cosdelt = cos(p[2]);

		nSource[0] = cosalt*cosdelt;
		nSource[1] = sinalt*cosdelt;
		nSource[2] = sindelt;
		for (o = 0; o < sett->nifo; ++o){
			modvir(sinalt, cosalt, sindelt, cosdelt, sett->N, &ifo[o], aux, saa[o], sbb[o]);  
		}
		extr = Fstatnet(sett, p, nSource, saa, sbb);
		for (i = 0; i < 11; i++) maxres[i] = extr[i];

//Calculate F-statistics for the rest points - parallel
		#pragma omp parallel default(shared) private(j, o, i, p, res)
		{
/*			double **sigaa, **sigbb;   
			sigaa = matrix(sett->nifo, sett->N);
			sigbb = matrix(sett->nifo, sett->N);

			sinalt = sin(p[3]);
			cosalt = cos(p[3]);
			sindelt = sin(p[2]);
			cosdelt = cos(p[2]);

			nSource[0] = cosalt*cosdelt;
			nSource[1] = sinalt*cosdelt;
			nSource[2] = sindelt;
			for (o = 0; o < sett->nifo; ++o){
				modvir(sinalt, cosalt, sindelt, cosdelt, sett->N, &ifo[o], aux, sigaa[o], sigbb[o]);  
			}*/
			#pragma omp for
			for(j = 1; j < 17; j++){
				for(i = 0; i < 2; i++){ 
					p[i] = array[j][i];
				}
				p[2] = rslts[8];
				p[3] = rslts[9];
//				res = Fstatnet(sett, p, nSource, sigaa, sigbb); 
				res = Fstatnet(sett, p, nSource, saa, sbb);
//Compare points and find one with maximal value of F-statistics
				#pragma omp critical
				if (res[5] < maxres[5]){
					for (i = 0; i < 11; i++){ 
						maxres[i] = res[i];
					}
				}
			} //loop over points
//			free_matrix(sigaa, sett->nifo, sett->N);
//			free_matrix(sigbb, sett->nifo, sett->N);
		} //pragma
//Compare with the result from the previous step
		if(extr[5] <= maxres[5]){
//If old result is better - increase size of the meshes
			for(r = 0; r < 2; r++){ 
				param[r] = param[r] + incr[r];
			}
			for(j = 0; j < 11; j++){ 
				rslts[j] = extr[j];
			}
		}
		else{
//If new result is better - take it as initial point to the next step
//and shrink meshes to the original size
			for(j = 0; j < 11; j++){ 
				rslts[j] = maxres[j];
			}
			for(r = 0; r < 2; r++){ 
				param[r] = paramstart[r];
			}
		}
//Check if number of steps is not exceeded - if yes, stop while loop
		if(count >= limit){
			for(r = 0; r < 2; r++) param[r] = 100;
		}

	} // while loop
	for (n= 0; n < 11; n++){
		out[n] = rslts[n];
	}

//Free memory
	free(extr);
	free(res);
	free_matrix(saa, sett->nifo, sett->N);
	free_matrix(sbb, sett->nifo, sett->N);
//Print info about number of steps
	printf("skyMADS ends in %d steps.\n", count);
	return out;
}// end skyMADS()


/***************************************************************
Mesh adaptive direct search (MADS) maximum search declaration - classical
version in which function starts from relatively loose and extended grid,
which shrinks in the direction of maximum. Here not all points (in 
comparision with invertedMADS) on the 4D hypercube are used.
Idea and names of variables are very similar as introduced in
invertedMADS function.
***************************************************************/
Yep64f* MADS(Search_settings *sett, Aux_arrays *aux, double* sgnl, double* rslts, double tolmads, double *p2){
	int i, j, k, l, m, n, o, r, a = 0;
	int count=0, limit = 600;		//Maximal number of steps
  	double sinalt, cosalt, sindelt, cosdelt;
	double param[4]; 			//Initial size of mesh
	double smallest = 25.0;
	double array[25][4];
	double scale = 0.997;
	double p[4];
	double nSource[3];
	for(r = 0; r < 4;r++){
		param[r] = p2[r];
	}
#ifdef YEPPP
    yepLibrary_Init();
    Yep64f *out = (Yep64f*)malloc(sizeof(Yep64f)*11);
    Yep64f *extr = (Yep64f*)malloc(sizeof(Yep64f)*11); 
    Yep64f *res = (Yep64f*)malloc(sizeof(Yep64f)*11);
    Yep64f *maxres = (Yep64f*)malloc(sizeof(Yep64f)*11);
    enum YepStatus status;

#endif
	while(smallest >= tolmads){  	//when to end computations
		count ++;
		k = 0;
		for(i = 0; i < 11; i++){ 
			maxres[i] = extr[i] = rslts[i];
		}
		for(j = 0; j < 8; j++){
			for(k = 0; k < 4; k++){
				array[j][k] = extr[6+k];
			}
		}
		for (k = 0; k < 4; k++){
			array[k][k] = array[k][k] + param[k];
		}
		for (k = 0; k < 4; k++){
			j = 4 + k;
			array[j][k] = array[j][k] - param[k];
		}
		for(k = 8; k < 16; k++){
			array[k][0] = extr[6] + param[0];
			for (j = 8; j < 12; j++){
				array[j][1] = extr[7] + param[1];
				for (l = 8; l < 10; l++){
					array[l][2] = extr[8] + param[2];
				}
				array[8][3] = extr[9] + param[3];
				array[9][3] = extr[9] - param[3];
				for (l = 10; l < 12; l++){
					array[l][2] = extr[8] - param[2];
				}
				array[10][3] = extr[9] + param[3];
				array[11][3] = extr[9] - param[3];
			}
			for (j = 12; j < 16; j++){
				array[j][1] = extr[7] - param[1];
				for (l = 12; l < 14; l++){
					array[l][2] = extr[8] + param[2];
				}
				array[12][3] = extr[9] + param[3];
				array[13][3] = extr[9] - param[3];
				for (l = 14; l < 16; l++){
					array[l][2] = extr[8] - param[2];
				}
				array[14][3] = extr[9] + param[3];
				array[15][3] = extr[9] - param[3];
			}
		}
		for(k = 16; k < 24; k++){
			array[k][0] = extr[6] - param[0];
			for (j = 16; j < 20; j++){
				array[j][1] = extr[7] + param[1];
				for (l = 16; l < 18; l++){
					array[l][2] = extr[8] + param[2];
				}
				array[16][3] = extr[9] + param[3];
				array[17][3] = extr[9] - param[3];
				for (l = 18; l < 20; l++){
					array[l][2] = extr[8] - param[2];
				}
				array[18][3] = extr[9] + param[3];
				array[19][3] = extr[9] - param[3];
			}
			for (j = 20; j < 24; j++){
				array[j][1] = extr[7] - param[1];
				for (l = 20; l < 22; l++){
					array[l][2] = extr[8] + param[2];
				}
				array[20][3] = extr[9] + param[3];
				array[21][3] = extr[9] - param[3];
				for (l = 22; l < 24; l++){
					array[l][2] = extr[8] - param[2];
				}
				array[22][3] = extr[9] + param[3];
				array[23][3] = extr[9] - param[3];
			}
		}
		for (k = 0; k < 4; k++){
			array[24][k] = extr[6+k];
		}
		double **saa, **sbb;   
		saa = matrix(sett->nifo, sett->N);
		sbb = matrix(sett->nifo, sett->N);
		for(i = 0; i < 4; i++){ 
			p[i] = array[24][i];
		}
		sinalt = sin(p[3]);
		cosalt = cos(p[3]);
		sindelt = sin(p[2]);
		cosdelt = cos(p[2]);

		nSource[0] = cosalt*cosdelt;
		nSource[1] = sinalt*cosdelt;
		nSource[2] = sindelt;
		for (o = 0; o < sett->nifo; ++o){
			modvir(sinalt, cosalt, sindelt, cosdelt, sett->N, &ifo[o], aux, saa[o], sbb[o]);  
		}
		extr = Fstatnet(sett, p, nSource, saa, sbb);
		#pragma omp parallel default(shared) private(j, o, p, i, sinalt, cosalt, sindelt, cosdelt, nSource, res)
		{
			double **sigaa, **sigbb;   
			sigaa = matrix(sett->nifo, sett->N);
			sigbb = matrix(sett->nifo, sett->N);
			#pragma omp for
			for(j = 0; j < 24; j++){
				for(i = 0; i < 4; i++){ 
					p[i] = array[j][i];
				}

				sinalt = sin(p[3]);
				cosalt = cos(p[3]);
				sindelt = sin(p[2]);
				cosdelt = cos(p[2]);

				nSource[0] = cosalt*cosdelt;
				nSource[1] = sinalt*cosdelt;
				nSource[2] = sindelt;
				for (o = 0; o < sett->nifo; ++o){
					modvir(sinalt, cosalt, sindelt, cosdelt, sett->N, &ifo[o], aux, sigaa[o], sigbb[o]);  
				}
//Fstat function for mesh points
				res = Fstatnet(sett, p, nSource, sigaa, sigbb); 

				#pragma omp critical
				if (res[5] < extr[5]){
					for (i = 0; i < 11; i++){ 
						maxres[i] = res[i];
					}
				}

			} //j
			free_matrix(sigaa, sett->nifo, sett->N);
			free_matrix(sigbb, sett->nifo, sett->N);
		} //pragma

		if(extr[5] <= maxres[5]){
			smallest = smallest*scale;
			for(r = 0; r < 4; r++){ 
				param[r] = scale*param[r];
			}
			for(j = 0; j < 11; j++){ 
				rslts[j] = extr[j];
			}
		}
		else{
			for(j = 0; j < 11; j++) rslts[j] = maxres[j];
		}
	if(count >= limit) smallest = 0.0;

	} // while loop

	for (n= 0; n < 11; n++){
		out[n] = rslts[n];
	}
//Free memory
	free(extr);
	free(res);
	return out;
} //end MADS()

/***************************************************************
Few functions for Nelder-Mead (simplex) algorithm - 'standard'
maximum search algorithm 
***************************************************************/
double ** make_simplex(double * point, int dim, double *pc2){
	int i, j;
	double ** simplexm = matrix(dim + 1, dim);
	for (i = 0; i < dim + 1; i++){
		for (j = 0; j < dim; j++){
			simplexm[i][j] = point[j];
		}
	}
	for (i = 0; i < dim; i++){
		simplexm[i][i] = simplexm[i][i] + pc2[i];
	}
	return simplexm;
	free_matrix(simplexm, dim + 1, dim);
}//end make_simplex()

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
				modvir(sinalt, cosalt, sindelt, cosdelt, sett->N, &ifo[o], aux, sigaa_max[o], sigbb_max[o]);  
			}
			out = Fstatnet(sett, simplex[i], nS, sigaa_max, sigbb_max);
			for (j = 0; j < 11; j++) fx[i][j] = out[j];
	}
	free_vector(out, 11);
} //end evaluate_simplex()

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
}//end simplex_extremes()

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
}// end simplex_bearings()

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
		modvir(sinalt, cosalt, sindelt, cosdelt, sett->N, &ifo[o], aux, sigaa_max[o], sigbb_max[o]);  
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
}//end update_simplex()

void contract_simplex(double ** simplex, int dim, double ** fx, int ilo, int ihi, Search_settings *sett, Aux_arrays *aux, double *nS, double **sigaa_max, double **sigbb_max, double scal){
  	double sinalt, cosalt, sindelt, cosdelt;
	double *fx3= alloc_vector(11);
	int i, j, k, o;

	for (i = 0; i < dim + 1; i++){
		if (i != ilo){
			for (j = 0; j < dim; j++) simplex[i][j] = (simplex[ilo][j]+simplex[i][j])*scal;
			sinalt = sin(simplex[i][3]);
			cosalt = cos(simplex[i][3]);
			sindelt = sin(simplex[i][2]);
			cosdelt = cos(simplex[i][2]);

			nS[0] = cosalt*cosdelt;
			nS[1] = sinalt*cosdelt;
			nS[2] = sindelt;

			for (o = 0; o < sett->nifo; ++o){
				modvir(sinalt, cosalt, sindelt, cosdelt, sett->N, &ifo[o], aux, sigaa_max[o], sigbb_max[o]);  
			}
			fx3 = Fstatnet(sett, simplex[i], nS, sigaa_max, sigbb_max);
			for (k = 0; k < 11; k++) fx[i][k] = fx3[k];
		}
	}
	free_vector(fx3, 11);

}//end contract_simplex()

int check_tol(double fmax, double fmin, double ftol){
	double delta = fabs(fmax - fmin);
	double accuracy = (fabs(fmax) + fabs(fmin)) * ftol;
	return (delta < (accuracy + ZEPS));
}//end check_tol

/***************************************************************
Main function of the simplex (Nelder-Mead) algorithm
***************************************************************/
double * amoeba(Search_settings *sett, Aux_arrays *aux, double *point, double *nS, double *res_max, int dim, double tol, double *pc2, double **sigaa_max, double **sigbb_max){
	int ihi, ilo, inhi;
	int *ih;
	int j, i;
	static double NM_out[11];
	double ** fx = matrix(dim + 1, 11);
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
				contract_simplex(simplex, dim, fx, ilo, ihi, sett, aux, nS, sigaa_max, sigbb_max, 0.5);
			}
		}
	}

	for (j = 0; j < dim; j++) point[j] = simplex[ilo][j];
	for (j = 0; j < 11; j++) NM_out[j] = fx[ilo][j];
	free_matrix(fx, dim + 1, 10);
	free_vector(midpoint, dim);
	free_vector(line, dim);
	free_matrix(simplex, dim + 1, dim);

	return NM_out;
}//end amoeba()

/***************************************************************
Main part of the code
***************************************************************/
int main (int argc, char *argv[]) {

//variables declaration
	Search_settings sett;	
	Command_line_opts opts;
  	Aux_arrays aux_arr;
  	int i, j, r, c, a, b, g, flag=0; 	
	int k1, k2, k3, k4;
	int d, o, m, k, s;
	int bins = 17;				// neighbourhood of point will be divide into defined number of bins for neigh and naive flags
	int ROW, dim = 4;
	int gsize = 2;				// grid size where followup will be searching maximum for naive flag
	int hemi; 				// hemisphere
	int spndr[2], nr[2], mr[2], fo[2];	// range in linear unities
	int gri1[4];
	int nof[4];				// Number of points in optimal grid where Fstat will be calculated (in every dim)
	double minm;				// minimal match used in optimal 4d grid generation
	double tol = 1e-7;			// accuracy of Simplex
	double tolmads = 1e-2;			// accuracy of MADS (but not invertedMADS)
	double cof, al1, al2;
	double s1, s2, s3, s4;
  	double sinalt, cosalt, sindelt, cosdelt;
	double F_min;
	double au;
	double temp1, temp2;
	double x, y;
	double sc;				// scaling factor
	double pc[4];				// % define neighbourhood around each parameter for initial grid for neigh flag
	double pc2[4];				// % define neighbourhood around each parameter for direct maximum search (invertedMADS, MADS & Simplex)
	double pc0[2];				// define neighbourhood around frequency and spindown parameter (skyMADS)
	double sgnlo[4];
	double nearest_point[4];
	double sgnlol[4]; 
	double be[2];		 
	double nSource[3];
	double sgnlo_max[4];
	double nSource_max[3];
	double **sigaa_max, **sigbb_max; 	// Auxiliary arrays to storage amplitude modulation factors for the maximum on the grid
  	double *MM ; 	 			// Auxiliary array for grid points in naive grid
	double **arr;				// arr[ROW][COL] - array of points where Fstat will be calculated
	double **sgnlo_range;  
	sgnlo_range = matrix(4, 2); 

#ifdef YEPPP
    	yepLibrary_Init();
   	Yep64f *results_max = (Yep64f*)malloc(sizeof(Yep64f)*11); 
    	Yep64f *results_first = (Yep64f*)malloc(sizeof(Yep64f)*11);
    	Yep64f *results = (Yep64f*)malloc(sizeof(Yep64f)*11);
	Yep64f *maximum = (Yep64f*)malloc(sizeof(Yep64f)*11);
	Yep64f *mean = (Yep64f*)malloc(sizeof(Yep64f)*4); 

    enum YepStatus status;

#endif

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
// Scaling factor
//	sc = 6./sett.nod;
	sc = 1.;
//Auxiliary arrays
	sigaa_max = matrix(sett.nifo, sett.N);
	sigbb_max = matrix(sett.nifo, sett.N);
// Define on how many grid points Fstat will be calculated on (optimal, naive, neigh) grid - in every dimension
	if ((!opts.naive_flag)&&(!opts.neigh_flag)){
		nof[0] = 5;
		nof[1] = 5;
		nof[2] = 5;
		nof[3] = 5;
		ROW = (2*nof[0]+1)*(2*nof[1]+1)*(2*nof[2]+1)*(2*nof[3]+1);
	}
	else{
		ROW = pow((bins+1),4);
	}
	arr = matrix(ROW, 4);
// Grid data
	if(!opts.neigh_flag){
		read_grid(&sett, &opts);
	}
// Amplitude modulation functions for each detector  
	for(i=0; i<sett.nifo; i++){ 
		rogcvir(&ifo[i]); 
	}
// Adding signal from file
  	if(strlen(opts.addsig)) { 
    		add_signal(&sett, &opts, &aux_arr);
  	}
// Reading initial point from file	
	FILE *coi;
	int z;
	if ((coi = fopen(opts.candidates, "r")) != NULL) {
		while(fscanf(coi, "%le %le %le %le", &mean[0], &mean[1], &mean[2], &mean[3]) == 4){
// Time test
//			tstart = clock();

//Frequency shifted from the reference frame to the current frame
			if (opts.refr > 0){
				mean[0] += -2.*mean[1]*(sett.N)*(opts.refr - opts.ident); 
			}
// In the case of --onepoint flag code doesn't calculate ANY grid,
//just computes F-statistics for the initial point and (in the case of
//--mads or --simplex flags) runs maximum search
			if(!opts.onepoint_flag){
//Function neighbourhood - generating grid around point
//Area around starting point is calculating as a percent from initial values
				if(opts.neigh_flag){
					pc[0] = 5e-5;
					pc[1] = 5e-10;
					pc[2] = 0.015;
					pc[3] = 0.015;
					for (i = 0; i < 4; i++){
						pc2[i] = 2*pc[i]/bins;
					}
					neigh2(mean, pc, bins, arr);
				} //if neigh
//Naive grid
				else if(opts.naive_flag){
//if naive or optimal grid -> find closest points on the grid from grid.bin file
// and uniformly divide the space between them
					cof = sett.oms + mean[0]; 
					for(i=0; i<2; i++) sgnlol[i] = mean[i]; 
					hemi = ast2lin(mean[3], mean[2], C_EPSMA, be);	 
					sgnlol[2] = be[0]*cof; 
					sgnlol[3] = be[1]*cof;
// solving a linear system in order to translate sky position, frequency 
//and spindown (sgnlo parameters) into the position in the grid
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

					lin2ast(al1/sett.oms, al2/sett.oms, hemi, sett.sepsm, sett.cepsm, &sinalt, &cosalt, &sindelt, &cosdelt);
// Find nearest point on the grid to the initial point
					nearest_point[0] = mean[0];	
					nearest_point[1] = spndr[0]*sett.M[5] + nr[0]*sett.M[9] + mr[0]*sett.M[13];
					nearest_point[2] = asin(sindelt);
					nearest_point[3] = fmod(atan2(sinalt, cosalt) + 2.*M_PI, 2.*M_PI);
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
							al1 = nr[i]*sett.M[10] + mr[i]*sett.M[14];
							al2 = nr[i]*sett.M[11] + mr[i]*sett.M[15];
// Check if the range is in an appropriate region of the grid
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
							sgnlo_range[0][0] = mean[0] - 0.005;
							sgnlo_range[0][1] = mean[0] + 0.005;
							for (i = 0; i < 4; i++){ 
								if (sgnlo_range[i][0] > sgnlo_range[i][1]){
									au = sgnlo_range[i][0];
									sgnlo_range[i][0] = sgnlo_range[i][1];
									sgnlo_range[i][1] = au;
								}
								printf("range[%d][0] = %le, range[%d][1] = %le\n", i, sgnlo_range[i][0], i, sgnlo_range[i][1]);
							}
							neigh_from_range(sgnlo_range, bins, arr);
							for (j = 0; j < 4; j++) pc2[j] = (sgnlo_range[j][1] - sgnlo_range[j][0])/(2*bins);
						}

					} //if candidate is in good region
				} //if naive
// Optimal, precise grid calculations
				else if ((!opts.naive_flag)&&(!opts.neigh_flag)){

					char filename[512];
					FILE *fish;
					double wg1,wg2;
					double mg[16];	
					double mg1[16];				
					double mg2[16];	
					double **Mopt;	
					Mopt = matrix(4, 4);
// In the case of two detectors: weighted matrices
					wg1=wg2=0.5;
// Convert coordinates						
					cof = sett.oms + mean[0]; 
					for(i=0; i<2; i++) sgnlol[i] = mean[i]; 
					hemi = ast2lin(mean[3], mean[2], C_EPSMA, be);	 
					sgnlol[2] = be[0]*cof; 
					sgnlol[3] = be[1]*cof;
// Read file with the Fisher matrix
					if(strlen(opts.usedet)==2){
						sprintf (filename, "%s/%03d/%s/%s", opts.dtaprefix, opts.ident, opts.usedet, opts.fisher);
						if ((fish=fopen (filename, "r")) != NULL) {
							for(i=0; i<16; i++) fscanf(fish, "%le",i+mg); 
    							fclose (fish);
						}
						else { 
							perror(filename);
							puts("Problem with reduced Fisher matrix file. Exiting...");
							exit(0);
						}
					}
					else { 
						sprintf (filename, "%s/%03d/H1/%s", opts.dtaprefix, opts.ident, opts.fisher);
						if ((fish=fopen (filename, "r")) != NULL) {
							for(i=0; i<16; i++) fscanf(fish, "%le",i+mg1); 
    							fclose (fish);
						}
						else { 
							perror(filename);
							puts("Problem with reduced Fisher matrix file (H1). Exiting...");
							exit(0);
						}
						sprintf (filename, "%s/%03d/L1/%s", opts.dtaprefix, opts.ident, opts.fisher);
						if ((fish=fopen (filename, "r")) != NULL) {
							for(i=0; i<16; i++) fscanf(fish, "%le",i+mg2); 
    							fclose (fish);
						}
						else { 
							perror(filename);
							puts("Problem with reduced Fisher matrix file (L1). Exiting...");
							exit(0);
						}
//Weigh matrices in the case of 2 detectors
						for(i=0; i<16; i++){ 
							mg[i] = wg1*mg1[i]+wg2*mg2[i];
						}
					}
// Values of minimal match put by hand (according to our tests)
					if (sett.nod == 6){
						minm = 0.999;
					}
					else {
						minm = 0.9999;
					}
//Calculate grid points
					A4opt(minm, sett.N, mg, Mopt);
					InjLoc(gri1, mean, Mopt);
// Define the grid range in which the signal will be looked for
					fo[0] = gri1[0];
					spndr[0] = gri1[1];
					nr[0] = gri1[2];
					mr[0] = gri1[3];
					fo[1] = fo[0] + nof[0]; 
					fo[0] -= nof[0];
					spndr[1] = spndr[0] + nof[1]; 
					spndr[0] -= nof[1];
					nr[1] = nr[0] + nof[2]; 
					nr[0] -= nof[2];
					mr[1] = mr[0] + nof[3]; 
					mr[0] -= nof[3];
					printf("Ranges on the grid: fo = %d %d; spndr = %d %d; nr = %d %d; mr = %d %d\n", fo[0], fo[1], spndr[0], spndr[1], nr[0], nr[1], mr[0], mr[1]);
// Translate to the astrophysical coordinates
					for (i = 0; i < 2; i++){
						sgnlo_range[0][i] = fo[i]*Mopt[0][0] + spndr[i]*Mopt[1][0] + nr[i]*Mopt[2][0] + mr[i]*Mopt[3][0];
						sgnlo_range[1][i] = fo[i]*Mopt[0][1] + spndr[i]*Mopt[1][1] + nr[i]*Mopt[2][1] + mr[i]*Mopt[3][1];
						sgnlo_range[2][i] = fo[i]*Mopt[0][2] + spndr[i]*Mopt[1][2] + nr[i]*Mopt[2][2] + mr[i]*Mopt[3][2];
						sgnlo_range[3][i] = fo[i]*Mopt[0][3] + spndr[i]*Mopt[1][3] + nr[i]*Mopt[2][3] + mr[i]*Mopt[3][3];
					}
					for (i = 0; i < 4; i++){ 
						if (sgnlo_range[i][0] > sgnlo_range[i][1]){
							au = sgnlo_range[i][0];
							sgnlo_range[i][0] = sgnlo_range[i][1];
							sgnlo_range[i][1] = au;
						}
					}
// Put optimal points into array
					i = 0;
					for (k1 = fo[0]; k1 <= fo[1]; k1++){
						for (k2 = spndr[0]; k2 <= spndr[1]; k2++){
							for (k3 = nr[0]; k3 <= nr[1]; k3++){
								for (k4 = mr[0]; k4 <= mr[1]; k4++){
									for (j = 0; j < 4; j++) {
										arr[i][j] = k1*Mopt[0][j] + k2*Mopt[1][j] + k3*Mopt[2][j] + k4*Mopt[3][j];
									}
									i++;
								}							
							}							
						}
					}
// Prepare area of search for invertedMADS, MADS and Simplex
					if(opts.simplex_flag||opts.mads_flag){
// Values put by hand, according to our tests
						pc2[0] = 0.05*sc;
						pc2[1] = 5e-8*sc;
						pc2[2] = 0.15*sc;
						pc2[3] = 0.10*sc;
					} //if simplex or mads flag
// Prepare area of search for skyMADS
					if(opts.skymads_flag){
						pc0[0] = 0.2;
						pc0[1] = 1e-7;
					}//if skymads flag
				} //if optimal grid
				if(flag == 1) continue;
				results_max[5] = 0.;

// Main loop - over all parameters saved in array + OpenMP parallelisation
				#pragma omp parallel default(shared) private(d, i, sgnlo, sinalt, cosalt, sindelt, cosdelt, nSource, results, maximum)
				{
// This needs to be inside pragma, because every thread needs its own copy
			       		double **sigaa, **sigbb;  
		      			sigaa = matrix(sett.nifo, sett.N);
					sigbb = matrix(sett.nifo, sett.N);
// Proper loop
					#pragma omp for  
					for (d = 0; d < ROW; ++d){
						if (arr[d][0] < 0.0){ 
							puts("Inappropriate region of the range! -> going to the next point");
							continue;
						}
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
							modvir(sinalt, cosalt, sindelt, cosdelt, sett.N, &ifo[i], &aux_arr, sigaa[i], sigbb[i]);  
						}
// Calculate F-statistic in given point
						results = Fstatnet(&sett, sgnlo, nSource, sigaa, sigbb);
// Check is it the biggest found value of F-statistic so far
						#pragma omp critical
						if(results[5] < results_max[5]){
							for (i = 0; i < 11; i++){
								results_max[i] = results[i];
							}
// Save information about this point for Simplex and/or invertedMADS
							if(opts.simplex_flag||opts.mads_flag||opts.skymads_flag){
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
// Free memory
					free_matrix(sigaa, sett.nifo, sett.N);
					free_matrix(sigbb, sett.nifo, sett.N);
				} //pragma 
				for(g = 0; g < 11; g++) results_first[g] = results_max[g];
			} // if not onepoint
			else{
//In the case of --onepoint calculate F-statistics only for initial point
				double **sigaa, **sigbb;   
	      			sigaa = matrix(sett.nifo, sett.N);
				sigbb = matrix(sett.nifo, sett.N);

// Prepare area of search for invertedMADS, MADS and Simplex
				if(opts.simplex_flag||opts.mads_flag){
// Values put by hand, according to our tests
					pc2[0] = 0.05*sc;
					pc2[1] = 5e-8*sc;
					pc2[2] = 0.15*sc;
					pc2[3] = 0.10*sc;
				} //if simplex or mads flag
// Prepare area of search for skyMADS
				if(opts.skymads_flag){
					pc0[0] = 0.2;
					pc0[1] = 1e-7;
				}//if skymads flag

				for (i = 0; i < 4; i++){
					sgnlo[i] = mean[i];
				}
				sinalt = sin(sgnlo[3]);
				cosalt = cos(sgnlo[3]);
				sindelt = sin(sgnlo[2]);
				cosdelt = cos(sgnlo[2]);

				nSource[0] = cosalt*cosdelt;
				nSource[1] = sinalt*cosdelt;
				nSource[2] = sindelt;

				for (i = 0; i < sett.nifo; ++i){
					modvir(sinalt, cosalt, sindelt, cosdelt, sett.N, &ifo[i], &aux_arr, sigaa[i], sigbb[i]);  
				}
				results = Fstatnet(&sett, sgnlo, nSource, sigaa, sigbb);
				puts("Maximum from --onepoint:");
				printf("%le %le %le %le %le %le\n", results[6], results[7], results[8], results[9], results[5], results[4]);
				for (i = 0; i < 11; i++){
					results_max[i] = results[i];
				}
				if((opts.simplex_flag)||(opts.mads_flag)||opts.skymads_flag){
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
				} //if --simplex or --mads
			} // if onepoint

// Maximum search using invertedMADS/MADS algorithm
  			if(opts.mads_flag) {
				puts("MADS");
//					maximum = MADS(&sett, &aux_arr, sgnlo_max, results_max, tolmads, pc2);
				maximum = invertedMADS(&sett, &aux_arr, sgnlo_max, results_max, pc2);
				if(maximum[5] < results_max[5]){
					for (i = 0; i < 11; i++){
						results_max[i] = maximum[i];
					}
				}
			} //mads
// Maximum search using skyMADS algorithm - only in frequency and spindown parameters
  			if(opts.skymads_flag) {
				puts("skyMADS");
//					maximum = MADS(&sett, &aux_arr, sgnlo_max, results_max, tolmads, pc2);
				maximum = skyMADS(&sett, &aux_arr, sgnlo_max, results_max, pc0);
				if(maximum[5] < results_max[5]){
					for (i = 0; i < 11; i++){
						results_max[i] = maximum[i];
					}
				}
			} //mads

// Maximum search using Simplex algorithm
			if(opts.simplex_flag){
				puts("Simplex");
				maximum = amoeba(&sett, &aux_arr, sgnlo_max, nSource_max, results_max, dim, tol, pc2, sigaa_max, sigbb_max);
				if(maximum[5] < results_max[5]){
					for (i = 0; i < 11; i++){
						results_max[i] = maximum[i];
					}
				}
			} //simplex
//Time test
//			tend = clock();
//			tdiff = (tend - tstart)/(double)CLOCKS_PER_SEC;
			puts("Maximum from grid:");
			printf("%le %le %le %le %le %le\n", results_first[6], results_first[7], results_first[8], results_first[9], results_first[5], results_first[4]);
			puts("Maximum from Simplex/invertedMADS/MADS:");
			printf("%le %le %le %le %le %le\n", results_max[6], results_max[7], results_max[8], results_max[9], results_max[5], results_max[4]);
		} // while fread coi
	} //if coi
	else {
		
		printf ("Problem with %s\n", opts.candidates);
		return 1;
	}

// Additional output information
/*	puts("**********************************************************************");
	printf("***	Maximum value of F-statistic for grid is : (-)%.8le	***\n", -results_first[5]);
	printf("Sgnlo: %.8le %.8le %.8le %.8le\n", results_first[6], results_first[7], results_first[8], results_first[9]);
	printf("Amplitudes: %.8le %.8le %.8le %.8le\n", results_first[0], results_first[1], results_first[2], results_first[3]);
	printf("Signal-to-noise ratio: %.8le\n", results_first[4]); 
	printf("Signal-to-noise ratio from estimated amplitudes (for h0 = 1): %.8le\n", results_first[10]);
	puts("**********************************************************************");
if((opts.mads_flag)||(opts.simplex_flag)||(opts.skymads_flag)){
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
	free(mean);
	if(!opts.neigh_flag){
		free_matrix(sgnlo_range, 4, 2);
	}
	free_matrix(arr, ROW, 4);
	if(opts.simplex_flag){
		free_matrix(sigaa_max, sett.nifo, sett.N);
		free_matrix(sigbb_max, sett.nifo, sett.N);
	}
  	cleanup_followup(&sett, &opts, &aux_arr);
	return 0;
}
