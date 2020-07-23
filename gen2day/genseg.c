#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <fftw3.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctype.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

int GrubbsOutliersMany (double *, double *, int, int, double, char *);
#include "iniparser.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

/* Default configuration file. */
/* Cann be overridden by the first command-line arg */
#define INI_FNAME "input_s6.txt"

/* sft filename extension */
#define SFT_SUFFIX ".out"

/* Buffer size for the outlier removal algorithm */
#define BUFSIZE 1<<15

#define EPSILON 0.40909280422232891

#define MAX_LINE 8192
#define MAX_FILES 8192

#ifdef USE_LAL
/* Compile the ephemeris generation code */
/* Requires the LALPulsar library */
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/LALBarycenter.h>
#include "EphemerisDetector.h"
#endif



int main (int argc, char *argv[]) {
     FILE *data, *td_stream;
     char inp[MAX_LINE], buf[MAX_LINE], td_fname[MAX_LINE],		\
	  tdc_fname[MAX_LINE], td_dir[MAX_LINE], ini_fname[MAX_LINE],	\
	  date_fname[MAX_LINE], sft_fname[MAX_LINE];
     double fpo, tmp1, tmp2, gpsd, gpsd1, gpsdn, gps1, *gpsdv, *gpsd1v,	\
	  gpsdRest = 0;
     int N, lfft, lfftr, lfftm, lenx, ldat, n, i, j, yndx, notempty, \
	  bufsize = BUFSIZE, dd = 1, nSeg, nfiles = 0,		\
	  nxRest = 0, nxall = 0, flidx, offset, gen_sci, nout;
     double *rdt, *rtmp, *xtime, *x0, alpha, dt, othr, *xRest=NULL,	\
	  *xall, w_taper_dt;
     int fftsize[] = {0};
     fftw_complex *Xin_array, *Xout_array;
     fftw_plan plan;
     struct stat st = {0};
     dictionary *ini;
     char *site, *plsr, *DataDir, *SftDir, *flsum, *fnptr,	\
	  *tptr, *flnames[MAX_FILES], *segments_fname, *out_replace ;
     
     double mingps=0., maxgps=0.;
     double start, end;
     int numbersegan=0;
     
#ifdef USE_LAL
     int gen_eph;
     char *EphDir, *efile, *sfile, eFname[MAX_LINE], sFname[MAX_LINE],	\
	  eph_fname[MAX_LINE];
     double *DetSSB, *rDet, *rSSB, mjd1, phir, elam = 0, position[4];
     EphemerisData *edat = NULL;
     Detectors detector = 0;
#endif
     double startgps;

     fprintf (stderr, "%s: extracting time-domain data from sft's\n", argv[0]);
     if (argc > 1) strcpy (ini_fname, argv[1]); else strcpy (ini_fname, INI_FNAME);
     fprintf (stderr, "Loading config file %s\n", ini_fname);
     if ((ini = iniparser_load (ini_fname)) == NULL) {
	  perror (ini_fname);
	  return 1;
     }

     startgps = iniparser_getdouble(ini, "general:startgps", 0.0); // write time sequences starting from this time

     site = iniparser_getstring (ini, "data:site", NULL);
     plsr = iniparser_getstring (ini, "data:plsr", NULL);
     DataDir = iniparser_getstring (ini, "data:datadir", NULL);
     SftDir = iniparser_getstring (ini, "data:sftdir", NULL);
     flsum = iniparser_getstring (ini, "data:flsum", NULL);
     N = iniparser_getint (ini, "general:n", 344656);
     lfft = iniparser_getint (ini, "general:lfft", 2048);
     *fftsize = lfft;
     lfftr = lfft/4;
     alpha = iniparser_getdouble (ini, "general:alpha", 0.1);
     dt = iniparser_getdouble (ini, "general:dt", 0.5);
     othr = iniparser_getdouble (ini, "general:othresh", 250.);   // threshold to clean large outliers
     out_replace = iniparser_getstring (ini, "data:out_replace", "gauss");  // replace outliers with zero or random gaussian value
     gen_sci = iniparser_getboolean (ini, "general:scientific", 1); // use science data only
     w_taper_dt = iniparser_getdouble (ini, "general:w_taper_dt", 1200.);   // Tukey window tapering size in seconds; if <=0 do not apply window to science regions


     const gsl_rng_type *T;
     gsl_rng *r = NULL;
     unsigned long mySeed;
     if (! strcmp(out_replace,"gauss")) {
	  // initialize random numbers
	  mySeed = time(NULL);
	  gsl_rng_env_setup();
	  T = gsl_rng_default;
	  r = gsl_rng_alloc (T);
	  gsl_rng_set(r, mySeed);
     }


     double *segar = NULL;

     if(gen_sci){
	  segments_fname = iniparser_getstring (ini, "general:segments", NULL);

	  FILE  *segment_file = fopen(segments_fname, "r");
	  if(segment_file == NULL){
	       perror("Error! Cannot open segments file \n" );
	       return 1;
	  } else {
	       printf("Reading file %s\n", segments_fname);
	       fscanf(segment_file, "%lf %lf", &start, &end);
	       mingps = start;
	       maxgps = end;
	       while(fscanf(segment_file, "%lf %lf", &start, &end) != EOF){
		    if (maxgps < end)   maxgps = end;
		    if (mingps > start) mingps = start;
	       }
	       
	       fprintf(stderr, "[sci] mingps = %f, maxgps = %f, dt = %f\n", mingps, maxgps, dt);
	       numbersegan = (int)ceil((maxgps-mingps)/dt);
	       fprintf(stderr, "[sci] numbersegan = %d \n", numbersegan);
	       fprintf(stderr, "[sci] w_taper_dt = %f [s]\n", w_taper_dt);
	       segar = (double *)malloc(numbersegan*sizeof(double));
	       // initialize filter to 0
	       for(i=0; i<numbersegan; i++)
		    segar[i] = 0.;
	       
	       rewind(segment_file);
	       // set mas to 1 in scientific regions
	       while(fscanf(segment_file, "%lf %lf", &start, &end) != EOF){
		    int starti = (int)floor((start-mingps)/dt);
		    int endi   = (int)ceil((end-mingps)/dt);
		    fprintf(stderr, "starti = %d endi = %d\n", starti, endi);
		    for(i=starti; i<=endi; i++)
			 segar[i] = 1.;
	       }
	  }
     }//end if(gen_sci)

     
#ifdef USE_LAL
     gen_eph = iniparser_getboolean (ini, "general:ephemeris", 0);
     if (gen_eph) {
	  detector = get_detector (site);
	  get_position (detector, position);
	  elam = position[1];
	  fprintf (stderr, "Detector ephemeris for %s will be created\n",	names[detector]);
	  EphDir = iniparser_getstring (ini, "general:EphDir", NULL);
	  efile = iniparser_getstring (ini, "general:efile", NULL);
	  sfile = iniparser_getstring (ini, "general:sfile", NULL);
	  sprintf (eFname, "%s/%s", EphDir, efile);
	  sprintf (sFname, "%s/%s", EphDir, sfile);
	  if ((edat = XLALInitBarycenter (eFname, sFname)) == NULL) {
	       perror (eFname);
	       return 1;
	  };
     }
     DetSSB = (double *) calloc (3*N+2, sizeof(double));
     rDet = (double *) calloc (3*N, sizeof(double));
     rSSB = (double *) calloc (3*N, sizeof(double));
#endif

     if (startgps != 0) gpsdRest = startgps; 
     printf (" startgps = %15.5f gpsdRest = %15.5f\n", startgps, gpsdRest);

     fnptr = flsum;
     while (*fnptr != '\0') {
	  flnames[nfiles] = (char *) malloc (MAX_LINE);
	  sprintf (flnames[nfiles], "%s_", plsr);
	  tptr = flnames[nfiles]+strlen (flnames[nfiles]);
	  while (isgraph((*tptr = *fnptr))) {
	       tptr++;
	       fnptr++;
	  }
	  *tptr = '\0';
	  strcat (flnames[nfiles], SFT_SUFFIX);
	  while (isblank(*fnptr)) fnptr++;
	  nfiles++; // count number of files in worklist
     }

     fprintf (stderr, "Time sequences for %s will be created. ", plsr);
     fprintf (stderr, "Found %d sft files in the worklist\n", nfiles);


     gpsdv = (double *) calloc (nfiles, sizeof (double)); // array for starting GPS time for files
     gpsd1v = (double *) calloc (nfiles, sizeof (double));

     // Read "nfile" number of *.out files
     for (flidx=0; flidx < nfiles; flidx++) {
	  fprintf (stderr, ">>> Reading data file %s ...", flnames[flidx]);
	  sprintf (sft_fname, "%s/%s", SftDir, flnames[flidx]);
	  // read header and check the file
	  if ((data=fopen (sft_fname, "r")) != NULL) {
	       if (fgets (inp, MAX_LINE, data) == NULL) {
		    perror ("Incorrect data file format");
		    return 1;
	       }
	       if (fgets (inp, MAX_LINE, data) == NULL) {
		    perror ("Incorrect data file format");
		    return 1;
	       }
	       sscanf (inp, "%s %lf", buf, &fpo);
	       if (fgets (inp, MAX_LINE, data) == NULL) {
		    perror ("Incorrect data file format");
		    return 1;
	       }
	       if (fgets (inp, MAX_LINE, data) == NULL) {
		    perror ("Incorrect data file format");
		    return 1;
	       }
	       sscanf (inp, "%s %lf %lf %lf %lf", buf, &tmp1, &tmp2, &gpsd, &gpsdn);
	       gpsd += gpsdn/1.0e-9;  // staring GPS time for a file
	       gpsdv[flidx] = gpsd;
	       ldat = 0;
	       rdt = (double *) malloc (2*sizeof(double));
	       rtmp = rdt;

	       while (fgets (inp, MAX_LINE, data)) {
		    if (*inp != '%' && sscanf (inp, "%lf %lf", rtmp, rtmp+1) == 2) { 
			 ldat++;
			 rdt = (double *) realloc (rdt, 2*sizeof (double)*(ldat+1));
			 rtmp = rdt+ldat*2;
		    }
	       }
	       fclose (data);

	  } else {
	       perror (sft_fname);
	       return 1;
	  }
	  fprintf (stderr, " Done.\n");

	  n = 2*ldat/lfft;
	  lfftm = lfft-2*lfftr;
	  lenx = n*lfftm+2*lfftr;
	  xtime = (double *) malloc (lenx*sizeof (double));

	  /* Reconstruct full FFT */

	  fprintf (stderr, "--- Reconstructing time domain sequence... \n");

	  Xin_array = (fftw_complex *) fftw_malloc (n*lfft*sizeof (fftw_complex));
	  Xout_array = (fftw_complex *) fftw_malloc (n*lfft*sizeof (fftw_complex));
	  for (i=0; i<n; i++) {
	       memcpy (Xin_array+i*lfft, rdt+i*lfft, lfft/2*sizeof (fftw_complex));
	       memset (Xin_array+(2*i+1)*lfft/2, 0, lfft/2*sizeof (fftw_complex));
	  }

	  plan = fftw_plan_many_dft (1, fftsize, n, Xin_array, NULL, 1, lfft, \
				     Xout_array, NULL, 1, lfft,		\
				     FFTW_BACKWARD, FFTW_ESTIMATE);
	  fftw_execute (plan);
	  fftw_destroy_plan (plan);

	  /* Set up the time domain sequence.		*/
	  /* Remove overlapping, add borders.		*/
	  for (j=0; j<lfftr; j++)
	       *(xtime+j) = (*(Xout_array+j))[0]/lfft;
	  for (i=0; i<n; i++)
	       for (j=0; j<lfftm; j++)
		    *(xtime+lfftr+i*lfftm+j) = (*(Xout_array+i*lfft+lfftr+j))[0]/lfft;
	  for (j=0; j<lfftr; j++)
	       *(xtime+lfftr+n*lfftm+j) = (*(Xout_array+(n-1)*lfft+lfftr+lfftm+j))[0]/lfft;
	  fftw_free (Xin_array);
	  fftw_free (Xout_array);

	  if (gpsdRest == 0) gpsdRest = gpsd; //relevant for the first file

	  offset = (gpsd-gpsdRest)/dt;
	  nxall = lenx+offset;
	  printf("nxall = %d\n", nxall);
	  xall = (double *) calloc (nxall, sizeof (double));

	  fprintf (stderr, " Done.\n");
#ifdef DEBUG
	  fprintf (stderr, "\nLength of data = %d\n", lenx);
	  fprintf (stderr, "gpsd = %f\n", gpsd);
	  fprintf (stderr, "gpsdRest = %f\n", gpsdRest);
	  fprintf (stderr, "offset = %d\n", offset);
	  fprintf (stderr, "lenx = %d\n", lenx);
	  fprintf (stderr, "nxall = %d\n", nxall);
#endif

	  if (offset <= 0) {
	       //printf("offset = %d\n", offset);
	       if (nxRest > 0) {
		    //printf("nxRest = %d\n", nxRest);

#ifdef DEBUG
		    fprintf (stderr, "%s: long xRest", flnames[flidx]);
#else 
		    fprintf (stderr, " It was long xRest");
#endif
	       }
	       fprintf (stderr, "\n");
	       memcpy (xall, xtime-offset, nxall*sizeof(double));
	  } else {
	       printf("nxRest > offset\n");
	       if (nxRest > offset) {
#ifdef DEBUG
		    fprintf (stderr, "%s: short xRest\n", flnames[flidx]);
#else
		    fprintf (stderr, " It was short xRest\n");
#endif
	       } else {
#ifdef DEBUG
		    fprintf (stderr, "%s: all xRest\n", flnames[flidx]);
#else
		    fprintf (stderr, " It was all xRest\n");
#endif
	       }
	       memcpy (xall, xRest, nxRest*sizeof(double));
	       printf("offset(xtime) = %d, lenx = %d, nxall = %d\n", offset, lenx, nxall);
	       memcpy (xall+offset, xtime, lenx*sizeof(double));
	  }
	  gpsd1 = gpsdRest;
	  gpsd1v[flidx] = gpsd1;
	  fprintf (stderr, "gpsd1 = %f\n", gpsd1);

	  /* Divide data into N*dt [s] long segments */

	  nSeg = nxall/N;
	  nxRest = nxall-N*nSeg;
	  gpsdRest = gpsd1+N*nSeg*dt;
	  xRest = (double *) realloc (xRest, nxRest*sizeof(double));
	  memcpy (xRest, xall+N*nSeg, nxRest*sizeof(double));
#ifdef DEBUG
	  fprintf (stderr, "xRest: %d samples\n", nxRest);
#endif

	  // if not enough data for a full segment
	  if (nSeg < 1) continue;
	  
	  fprintf (stderr, "--- Preparing td segments...\n");

	  if (stat (DataDir, &st) == -1) {
	       mkdir (DataDir, 0755);
	       fprintf (stderr, "Data directory created: %s/\n", DataDir);
	  }

	  double *seg_sci_mask = malloc(N*sizeof(double)); // science data mask
	  //double *seg_out_mask = malloc(N*sizeof(double)); // outliers groups mask
	  int seg, i1, i2, l;

	  // loop over all segments ready
	  for (seg=0; seg<nSeg; seg++) {
	       gps1 = gpsd1+N*seg*dt;
	       
	       fprintf (stderr, "    Segment %d (local %d) ,  gps1=%f\n", dd, seg, gps1);
	       
	       // init mask of scientific regions in segment
	       int offsetgps=(int)round((gps1 - mingps)/dt);
	       for (i=0; i < N; i++){
		    if ((offsetgps+i) < 0  || (offsetgps+i) > numbersegan) {
			 seg_sci_mask[i] = 0.;
		    }else{
			 seg_sci_mask[i] = segar[offsetgps+i];
		    }
	       }
	       
	       // create Tukey windows in scientific regions
#define PI    3.141592653589793

	       int lobe_isize = (int)ceil(w_taper_dt/dt);
	       int li;

	       i = 0;
	       while(i < N){
		    if (seg_sci_mask[i] < 0.5) {i++; continue;}
		    i1 = i; // first point in sci region
		    i2 = i1;
		    while(++i)
			 if ((seg_sci_mask[i] < 0.5) || (i == N)) break;
		    
		    i2 = i-1; // last point in sci region

#if 1
		    // exclude large outliers groups at edges of sci regions
		    int is;
		    fprintf(stderr,"Sci reg: %d %d   ", i1, i2);
		    for (is=i1; is <= i2; is++){
			 if ( fabs(xall[seg*N+is]) < othr &&
			      fabs(xall[seg*N+is+1]) < othr ) {
			      break;
			 }
			 seg_sci_mask[is] = 0.;
		    }
		    i1 = is; 
		    for (is=i2; is >= i1; is--){
			 if ( fabs(xall[seg*N+is]) < othr &&
			      fabs(xall[seg*N+is-1]) < othr ) {
			      break;
			 }
			 seg_sci_mask[is] = 0.;
		    }
		    i2 = is;
		    fprintf(stderr,"changed to: %d %d\n", i1, i2);
#endif
		    
#if 1
		    // in case lobe_isize*2 < L ; +2 accounts for even and odd cases
		    li = MIN(lobe_isize, (i2-i1+2)/2 );
		    fprintf(stderr, "sci i1=%d  i2=%d  lobe_isize=%d\n", i1, i2, li);
		    for(l=0; l<li; l++){
			 double lobe = 0.5*(1.-cos(PI*l/li));
			 seg_sci_mask[i1 + l] = lobe;
			 seg_sci_mask[i2 - l] = lobe;
		    }
#endif		    
	       }

	       
	       // apply scientific mask
	       for (i=0; i < N; i++) {
		    xall[seg*N+i] *= seg_sci_mask[i];
	       }
	       
#if 1
               // remove ramaining large outliers above othr
	       fprintf (stderr, "--- Cleaning large outliers - replace with %s \n", out_replace);
	       
	       double sdval;
	       nout = 0;
	       sdval = 0.;
	       for (i=0; i<N; i++){
		    double v = xall[seg*N+i];
		    // do not count zeros
		    if (fabs(v) < othr && fabs(v) > 1.e-30){
			 sdval += v*v;
			 nout++;
		    }
	       }
	       sdval = sqrt(sdval/(double)(nout-1));
	       fprintf(stderr, "    stdval = %f ,  noutzero = %d\n", sdval, nout);
	       
	       nout = 0;
	       if (! strcmp(out_replace, "gauss")) {
		    for (i=0; i<N; i++)
			 // if (fabs(xall[seg*N+i]) > othr){
			 if (fabs(xall[seg*N+i]) > 6.*sdval){
			      xall[seg*N+i] = gsl_ran_gaussian_ziggurat(r, sdval);
			      nout++;
			 }
	       } else {
		    for (i=0; i<N; i++)
			 // if (fabs(xall[seg*N+i]) > othr){
			 if (fabs(xall[seg*N+i]) > 6.*sdval){			      
			      xall[seg*N+i] = 0.;
			      nout++;
			 }
	       }
	       fprintf(stderr,"    Large outliers removed : %d / %d = %d%% \n", nout, N, 100*nout/N);
#endif
	       // remove remaining outliers using Grubbs test
	       x0 = (double *) malloc (N*sizeof (double));
#if 1
	       nout = GrubbsOutliersMany(xall+seg*N, x0, N, bufsize, alpha, out_replace);
	       fprintf(stderr,"    Grubbs outliers removed : %d / %d = %d%% \n", nout, N, 100*nout/N);
#else
	       memcpy (x0, xall, N*sizeof(double));
#endif
	       // write output files
	       sprintf (td_dir, "%s/%03d", DataDir, dd);
	       if (stat (td_dir, &st) == -1) {
		    mkdir (td_dir, 0755);
	       }
	       sprintf (td_dir, "%s/%03d/%s", DataDir, dd, site);
	       if (stat (td_dir, &st) == -1) {
		    mkdir (td_dir, 0755);
	       }
	       if(gen_sci){
		    sprintf (td_fname, "%s/xdats_%03d_%s.bin", td_dir, dd, plsr);
		    sprintf (tdc_fname, "%s/xdatsc_%03d_%s.bin", td_dir, dd, plsr);}
	       else{
		    sprintf (td_fname, "%s/xdat_%03d_%s.bin", td_dir, dd, plsr);
		    sprintf (tdc_fname, "%s/xdatc_%03d_%s.bin", td_dir, dd, plsr);
	       }
	       sprintf (date_fname, "%s/starting_date", td_dir);
	       fprintf (stderr, "Segment %03d\t", dd);
	       // Do "xall" contains only zeros?
	       notempty=0;
	       for (yndx=0; yndx < N; yndx++){
		    if(xall[seg*N+yndx] != 0.0){
			 notempty=1;// not empty
			 break;
		    }
	       }
	       // write xdat_*.bin file if not contains only zeros
	       if(notempty){
		    if ((td_stream=fopen (td_fname, "w")) != NULL) {
			 fwrite ((void *)(xall+seg*N), sizeof(double), N, td_stream);
			 fclose (td_stream);
		    }
	       }
	       // Do this part of "x0" contains only zeros?
	       notempty=0;
	       for (yndx=0; yndx < N; yndx++){
		    if(x0[seg*N+yndx] !=0.0){
			 notempty=1;
			 break;
		    }
	       }
	       // write xdatc_*.bin file
	       if(notempty){
		    if ((td_stream=fopen (tdc_fname, "w")) != NULL) {
			 fwrite ((void *)(x0+seg*N), sizeof(double), N, td_stream);
			 fclose (td_stream);
		    }
	       }
	       if ((td_stream=fopen (date_fname, "w")) != NULL) {
		    fprintf (td_stream, "%.10e\n", gps1);
		    fclose (td_stream);
	       }

	       /* Generate detector ephemeris */
#ifdef USE_LAL
	       if (gen_eph) {
		    get_barycenter (gps1, detector, edat, DetSSB, rDet, dt, N);
		    for (j=0; j<3*N; j++) rSSB[j] = DetSSB[j] - rDet[j];
		    mjd1 = gps2mjd (gps1);
		    phir = sid (mjd1, elam);
		    DetSSB[3*N] = phir;
		    DetSSB[3*N+1] = EPSILON;

		    sprintf (eph_fname, "%s/DetSSB.bin", td_dir);
		    if ((td_stream=fopen (eph_fname, "w")) != NULL) {
			 fwrite ((void *)DetSSB, sizeof(double), 3*N+2, td_stream);
			 fclose (td_stream);
		    }
		    sprintf (eph_fname, "%s/rDet.bin", td_dir);
		    if ((td_stream=fopen (eph_fname, "w")) != NULL) {
			 fwrite ((void *)rDet, sizeof(double), 3*N, td_stream);
			 fclose (td_stream);
		    }
		    sprintf (eph_fname, "%s/rSSB.bin", td_dir);
		    if ((td_stream=fopen (eph_fname, "w")) != NULL) {
			 fwrite ((void *)rSSB, sizeof(double), 3*N, td_stream);
			 fclose (td_stream);
		    }
	       }
#endif
	       dd++;
	  }
	  free (xall);
	  free (xtime);
	  free (rdt);
	  free (x0);
	  free (flnames[flidx]);

	  fprintf (stderr, "Done.\n");
     }

     iniparser_freedict (ini);
     free (gpsdv);
     free (gpsd1v);
#ifdef USE_LAL
     free (rSSB);
     free (rDet);
     free (DetSSB);
#endif
     return 0;  
}
