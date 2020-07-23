#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics.h>
#include <time.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#ifdef DEBUG
#include <stdio.h>
#endif

double zcritical (double alpha, int n) {
     /* Computes the critical z value for rejecting outliers (GRUBBS TEST) */
     double tcrit;
     
     tcrit = gsl_cdf_tdist_Pinv (alpha/(2*n), n-2);
     return (n-1)/sqrt(n)*(sqrt(tcrit*tcrit/(n-2+tcrit*tcrit)));
}

int GaussOutliersOne (double *buf, int len, double alpha, gsl_rng *r) {
     /*
       Replace outliers in the input buffer with zeros using the Grubbs
       test (at the significance level alpha). The buffer is
       overwritten. Return the number of outliers found.
     */
     int noutliers = 0, idx, n, outlier = 1;
     double meanval, maxval, sdval, tn, critval, v, *tmp;

     tmp = (double *) malloc (len*sizeof (double));
     while (outlier) {
	  n = 0;
	  for (idx=0; idx<len; idx++)
	       if ((v=*(buf+idx)))
		    *(tmp+n++) = v;
	  if (n>2) {
	       meanval = gsl_stats_mean (tmp, 1, n);
	       maxval = meanval;
	       for (idx=0; idx<n; idx++) {
		    v = *(tmp+idx);
		    if (fabs(v-meanval) > fabs(maxval-meanval))
			 maxval = v;
	       }
	       sdval = gsl_stats_sd_m (tmp, 1, n, meanval);
	       tn = fabs ((maxval-meanval)/sdval);
	       critval = zcritical (alpha, n);
	       if (tn>critval) {
		    for (idx=0; idx<len; idx++)
			 if (*(buf+idx) == maxval) {
			      *(buf+idx) = gsl_ran_gaussian_ziggurat(r, sdval);
			      //printf("g %f ",*(buf+idx)); fflush(stdout);
			      noutliers++;
			 }
	       } else 
		    outlier = 0;
	  } else
	       outlier = 0;
     }
     free (tmp);
     return noutliers;
}


int ZeroOutliersOne (double *buf, int len, double alpha) {
     /*
       Replace outliers in the input buffer with zeros using the Grubbs
       test (at the significance level alpha). The buffer is
       overwritten. Return the number of outliers found.
     */
     int noutliers = 0, idx, n, outlier = 1;
     double meanval, maxval, sdval, tn, critval, v, *tmp;

     tmp = (double *) malloc (len*sizeof (double));
     while (outlier) {
	  n = 0;
	  for (idx=0; idx<len; idx++)
	       if ((v=*(buf+idx)))
		    *(tmp+n++) = v;
	  if (n>2) {
	       meanval = gsl_stats_mean (tmp, 1, n);
	       maxval = meanval;
	       for (idx=0; idx<n; idx++) {
		    v = *(tmp+idx);
		    if (fabs(v-meanval) > fabs(maxval-meanval))
			 maxval = v;
	       }
	       sdval = gsl_stats_sd_m (tmp, 1, n, meanval);
	       tn = fabs ((maxval-meanval)/sdval);
	       critval = zcritical (alpha, n);
	       if (tn>critval) {
		    for (idx=0; idx<len; idx++)
			 if (*(buf+idx) == maxval) {
			      *(buf+idx) = 0.;
			      noutliers++;
			 }
	       } else 
		    outlier = 0;
	  } else
	       outlier = 0;
     }
     free (tmp);
     return noutliers;
}


		  
int GrubbsOutliersMany (double *inp, double *out, int len, int bufsize,	\
			double alpha, char *out_replace) {
     /*
       Split the input into sequences of the length bufsize, omitting
       zeros. Run ZeroOutliersOne on any of short sequences and
       concatenate the results. Return the number of outliers removed.
     */
     double *buffer, *xpos;
     int i, ncopied, nout, nremoved = 0, *nonzero;
     
     // initialize random numbers
     const gsl_rng_type *T;
     gsl_rng *r = NULL;
     unsigned long mySeed;
     mySeed = time(NULL);
     gsl_rng_env_setup();
     T = gsl_rng_default;
     r = gsl_rng_alloc (T);
     gsl_rng_set(r, mySeed);
     
     buffer = (double *) malloc (bufsize*sizeof (double));
     nonzero = (int *) malloc (bufsize*sizeof (int));
     memset (out, 0, len*sizeof (double));
     xpos = inp;
     while (xpos<inp+len) {
	  ncopied = 0;
	  while (ncopied<bufsize && xpos<inp+len) {
	       if ((*(buffer+ncopied) = *xpos++) != 0.)
		    *(nonzero+ncopied++) = xpos-inp-1;
	  }
#ifdef DEBUG
	  fprintf (stderr, "%d samples in buffer, ", ncopied);
#endif
	  if (! strcmp(out_replace, "gauss"))
	       nout = GaussOutliersOne (buffer, ncopied, alpha, r);
	  else
	       nout = ZeroOutliersOne (buffer, ncopied, alpha);
    
	  nremoved += nout;
#ifdef DEBUG
	  fprintf (stderr, "%d outlier(s) cleaned, %4.1f%% done\n",	\
		   nout, 100.0*(xpos-inp)/(double)len);
#endif
	  for (i=0; i<ncopied; i++)
	       *(out+nonzero[i]) = buffer[i];
     }
     free (nonzero);
     free (buffer);
     
     return nremoved;
}
