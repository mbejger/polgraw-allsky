#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <getopt.h>
#include <dirent.h>
#include <time.h>
#include <unistd.h>
#include <limits.h>

#include "auxi.h"
#include "init.h" 
#include "settings.h"


// Default output and data directories
#ifndef PREFIX
#define PREFIX ./coinc-results
#endif

#ifndef DTAPREFIX
#define DTAPREFIX ./candidates
#endif


int main (int argc, char* argv[]) {

  Search_settings sett;
  Command_line_opts_coinc opts;
  Candidate_triggers trig; 

  // Command line options 
  handle_opts_coinc(&sett, &opts, argc, argv);  

  // Output data handling
  struct stat buffer;

  if (stat(opts.prefix, &buffer) == -1) {
    if (errno == ENOENT) {
      // Output directory apparently does not exist, try to create one
      if(mkdir(opts.prefix, 
	       S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) == -1) {
	perror (opts.prefix);
	return 1;
      }
    } else { // can't access output directory
      perror (opts.prefix);
      return 1;
    }
  }
  
  // Manage grid matrix  
  manage_grid_matrix(&sett, &opts);	

  // Search settings 
  search_settings(&sett); 

  printf("Settings dt: %f, oms: %f\n", sett.dt, sett.oms); 
  printf("Reference frame number: %d\n", opts.refr); 

  read_trigger_files(&sett, &opts, &trig); 

  // Free arrays at the end 
  free(sett.M); 

  return 0; 
	
}


static int way2compare_4c(const void *a, const void *b){
  /* compare first 4 columns; array allocated like this: int (*candi)[CLEN] */

  const int* x = (const int *)a;
  const int* y = (const int *)b;

  int diff = *y - *x;
  if (diff) return diff;
  diff = *(y+1) - *(x+1);
  if (diff) return diff;
  diff = *(y+2) - *(x+2);
  if (diff) return diff;
  return *(y+3) - *(x+3);
}

static int way2compare_c1(const void *a, const void *b){
  /* compare column 1 */

  const int* x = *((const int **)a);
  const int* y = *((const int **)b);
  
  return *(y+1) - *(x+1);
}



inline int min ( int a, int b ) { return a < b ? a : b; }

// Allocation of memory for martix with given number of rows and columns
int** matrix(int rows, int cols) {

  int k, **m;
  m = (int **)malloc(rows*sizeof(int *));
  
  for (k=0; k<rows; k++)
    m[k] = (int *)calloc(cols, sizeof(int));
  
  return m;
}


/* 
 *  Obtaining the trigger files list from a given directory
 * and reading their content 
 */ 

void read_trigger_files(Search_settings *sett, 
			Command_line_opts_coinc *opts, 
			Candidate_triggers *trig) {
  
  long i, j;
  int current_frame=0, frcount=0;
  int val, shift[4], scale[4]; 
  long candsize=INICANDSIZE, goodcands=0;
  double sqrN, omsN, v[4][4], be[2];
  FLOAT_TYPE tmp[4];

  char dirname[512], filename[1024];
  // Trigger files directory name 
  sprintf (dirname, "%s", opts->dtaprefix); 

  DIR *dp;
  struct dirent *ep;
  FILE *data; 

  typedef struct _frame_params {
    int  num;              // frame number
    long ninband;          // number of candidates in band
    long ncands;           // number of good candidates
    long i1, i2;           // read candidates between indices i1 and i2
    int  fi2;              // fi value at i2
    char candi_fname[512]; // candi temporary file name
    FILE *candi_fh;        // file handle to candi file
    FILE *trig_fh;         // trigger file handle
  } Frame_params;

  
  // 256 is max no of frames used in struct.h ...
  Frame_params fpar[256];
  int id2frcount[256];

  long ic, filelen, maxfilelen=0;

  // random prefix for temp. file names
  char datestr[15];
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  sprintf(datestr, "ctmp_%09ld", t.tv_nsec);
  printf("tmp prefix : %s\n", datestr);

  //trigger record length
#define TRLEN 5
  // candi record length
#define CLEN 6
  // allcandi record length
#define ACLEN 6

  FLOAT_TYPE (*candf)[TRLEN] = NULL;
  // candi 7 columns: fi, si, di, ai, orgpos, fr, i
  int (*candi)[CLEN] = NULL;
  // allcandi 6 columns in a row: fi, si, di, ai, orgpos, fr
  int (*allcandi)[ACLEN] = NULL;


  dp = opendir (dirname);
  if (dp != NULL) {

    // Calculating the shifts from opts->shift
    val = opts->shift;
    for(i=0; i<4; i++) shift[i] = 0; // Initial value: no shift
    i=3;
    while (val > 0) {
      if(val%10) shift[i] = val%10;
      i--; val /= 10;
    }

    printf("Cell shifts  in f, s, d, a directions: ");
    for(i=0; i<4; i++) printf("%d ", shift[i]);
    printf("\n");

    // Calculating the scaling from opts->scale
    val = opts->scale;
    for(i=0; i<4; i++) scale[i] = 1; // Initial value: no scaling
    i=3;
    while (val > 0) {
      if(val%10) scale[i] = val%10;
      i--; val /= 10;
    }

    printf("Cell scaling in f, s, d, a directions: ");
    for(i=0; i<4; i++) printf("%d ", scale[i]);
    printf("\n");

    // Transformation matrix elements: division by the corresponding 
    // scale factors outside the main loop over candidates    
    for(j=0; j<4; j++)
      for(i=0; i<4; i++)
        v[i][j] = (sett->vedva[i][j])/scale[j];

    sqrN = pow(sett->N, 2);
    omsN = sett->oms*sett->N;


    while ((ep = readdir (dp))) {

      if( !((ep->d_type == DT_REG) || (ep->d_type == DT_LNK)) ||
	  (strstr(ep->d_name, opts->trigname) == NULL) )
	continue;

      sprintf(filename, "%s/%s", opts->dtaprefix, ep->d_name);

      // This part looks for the first number in the trigger file name,
      // under the assumption that this is the frame number
      char *fr, *epdname;
      epdname = strdup(ep->d_name);
      while((fr = strsep(&epdname, "_"))!=NULL) {
	if(fr[0] >= '0' && fr[0] <= '9') {
	  current_frame = atoi(fr);
	  printf("Reading %s... ", ep->d_name);
	  break;
	}
      }
      
      if((data = fopen(filename, "r")) != NULL) {
	
	// not used FLOAT_TYPE finband;
	// Original candidate number (position in the trigger file)
	// not needed int orgpos=0;
	// Counter for 'good' candidates i.e. these that are in band
	
	i = 0;
	frcount++;

	// reverse lookup table: frcount for given frame id
	id2frcount[current_frame] = frcount;	  
	fpar[frcount].trig_fh = data;

	// Each candidate is represented by 5 FLOAT_TYPE (double or float) numbers
	// c[0]=f, c[1]=s, c[2]=d, c[3]=a, c[4]=snr

	// get candsize from file size
	fseek(data, 0, SEEK_END);
	filelen = ftell(data);
	rewind(data);
	candsize = filelen/(TRLEN*sizeof(FLOAT_TYPE));
	printf("Frame[%04d]=%04d  ", frcount, current_frame);

	if (filelen > maxfilelen) {
	  maxfilelen = filelen;
	  //printf("\n[debug] realloc candi and candf!!!\n");
	  //candf = realloc(candf, sizeof(FLOAT_TYPE[candsize][TRLEN]));
	  //candi = realloc(candi, sizeof(int[candsize][CLEN]));
	  free(candf); free(candi);
	  candf = malloc(sizeof(FLOAT_TYPE[candsize][TRLEN]));
	  candi = malloc(sizeof(int[candsize][CLEN]));
	}

	if( fread((void *)candf, sizeof(FLOAT_TYPE[candsize][TRLEN]), 1, data) != 1 ) {
	  printf("Problem while reading from %s\n", filename);
	  exit(EXIT_FAILURE);
	}
	

	for(ic=0; ic<candsize; ++ic){
	  
	  if((candf[ic][0] > M_PI_2 - opts->narrowdown) && (candf[ic][0] < M_PI_2 + opts->narrowdown)) {
	    
	    // shifting c[0] (=frequency) to opts->refr reference frame
	    candf[ic][0] = candf[ic][0] + 2.*candf[ic][1]*(sett->N)*(opts->refr - current_frame);
	    
	    if(((candf[ic][0]>0) && (candf[ic][0]<M_PI)) && (candf[ic][4] > opts->snrcutoff)) {
	      
	      // Conversion to linear parameters
	      tmp[0] = candf[ic][0]*sett->N;
	      tmp[1] = candf[ic][1]*sqrN;
	      
	      // Transformation of astronomical to linear coordinates;
	      // C_EPSMA, an average value of epsm, is defined in settings.h
	      ast2lin(candf[ic][3], candf[ic][2], C_EPSMA, be);
	      
	      // tmp[2] corresponds to declination (d), tmp[3] to right ascension (a)
	      tmp[2] = omsN*be[0];
	      tmp[3] = omsN*be[1];
	      
	      for(j=0; j<4; j++) {
		// Integer values (0=fi, 1=si, 2=di, 3=ai)
		candi[i][j] = round(tmp[0]*v[0][j] + tmp[1]*v[1][j]
				  + tmp[2]*v[2][j] + tmp[3]*v[3][j] 
				  + 0.5*shift[j]);
	      } 
	      
	      // Saving the original position, frame number and current index
	      candi[i][4] = ic;
	      candi[i][5] = current_frame;
	      ++i; 
	      
	    } // if finband 
	  } // narrowdown
	} // ic - condidate index in trigger file
	
	// Frame number
	fpar[frcount].num = current_frame;
	// Number of candidates in band for a given frame
	fpar[frcount].ninband = i;

	// Looking for duplicates and selecting the one with highest SNR
	//--------------------------------------------------------------
	
	// Sorting the first 4 columns of candi
	qsort(candi, fpar[frcount].ninband, sizeof(int)*CLEN, way2compare_4c);

	// if frame==1 or trigger file is larger - allocate larger allcandi
	// (it will be reused later)
	if(frcount==1 || fpar[frcount].ninband > fpar[frcount-1].ninband ){
	  // if (allcandi) free(allcandi);
	  // allcandi = malloc(sizeof(int[ fpar[frcount].ninband ][ ACLEN ]));
	  //printf("\n[debug] realloc allcandi to %ld for frame %d !!!\n", fpar[frcount].ninband, frcount);
	  allcandi = realloc(allcandi, sizeof(int[ fpar[frcount].ninband ][ACLEN]));
	}
	
	long maxsnridx=0, frgoodcands=0;
	double candsnr=0.;
	for (i=0; i<fpar[frcount].ninband; ++i) {

	  int idx, idx1, diff=1;
	  for(j=0; j<4; ++j)
	    // using XOR: !(a^b) equals 1 for a=b
	    diff *= !((candi[i][j])^(candi[i+1][j]));

	  if(!diff) {
	    int k=i;
	    if(maxsnridx) k=maxsnridx;
	    for(j=0; j<6; ++j)
	      allcandi[frgoodcands][j] = candi[k][j];

 	    maxsnridx=0;
	    ++frgoodcands;
	    ++goodcands;

	  } else { // not unique - select candidate with highest SNR

	    idx = candi[i][4];
	    idx1 = candi[i+1][4];
	    if(!maxsnridx) {
	      maxsnridx = (candf[idx][4] > candf[idx1][4] ? i : i+1);
	    } else {
	      if(candf[idx][4] > candsnr) {
		maxsnridx = i;
	      } else if(candf[idx1][4] > candsnr) {
		maxsnridx = i+1;
	      }
	    }
	    candsnr = candf[candi[maxsnridx][4]][4];
	  }
	} // i - good candidates

	fpar[frcount].ncands = frgoodcands;
	printf("good/ninband/all = %10ld/%10ld/%10ld\n", frgoodcands, fpar[frcount].ninband, candsize);

	// save allcandi to disk
	sprintf(fpar[frcount].candi_fname, "%s_%s", datestr, ep->d_name); 
	printf("Writing candi to: %s\n", fpar[frcount].candi_fname);
	if((data = fopen(fpar[frcount].candi_fname, "w")) != NULL) {
	  fwrite(allcandi, sizeof(int[frgoodcands][ACLEN]), 1, data);
	  fclose(data);
	} else {
	  printf("Problem while opening file %s\n", fpar[frcount].candi_fname);
	  exit(EXIT_FAILURE);
	}

      } else {
	printf("Problem while opening file %s\n", filename);
      }
      printf("---------------------\n");
      
    } // readdir

  } // dp

  
  free(candf);
  free(candi);
  closedir(dp);
  
  printf("Total number of candidates from all frames: %ld\n", goodcands);

  // allocate allcandi 
  // available memory = max(n*GB, candf_size + candi_size)
  // candf_size+candi_size = candf_size*(1+(6*4)/(5*8)) = 1.6*candf_size
  // where candf_size is for the biggest trigger file.
  // This is max memory we had to allocate up to this point. 
  // If it's less then 4GB let's allocate 4GB because it's still reasonably small :)

#define GB 1073741824L
  long il, allcandi_size, chunk_size, cand_left, offset=0;
  int fr;

  allcandi_size = ( 0.1*GB > 1.6*maxfilelen) ? 4*GB : 1.7*maxfilelen;
  allcandi_size /= ACLEN*sizeof(int);
  chunk_size = allcandi_size/frcount;
  printf("[debug] maxfilelen=%ld   allcandi_size=%ld   chunk_size=%ld\n", maxfilelen, allcandi_size, chunk_size);
  allcandi = realloc(allcandi, sizeof(int[allcandi_size][ACLEN]));
  printf("\n[debug] realloc allcandi to %ld = %f GB!!!\n", allcandi_size, sizeof(int[allcandi_size][ACLEN])/(float)GB);

  // open files and reset i1 to 0
  for(fr=1; fr<=frcount; ++fr) {
    fpar[fr].i1 = 0;
    if( (fpar[fr].candi_fh = fopen(fpar[fr].candi_fname, "r")) != NULL) {
      printf("Opened %s\n", fpar[fr].candi_fname);
    } else {
      printf("Error opening %s\n", fpar[fr].candi_fname);
    }
  }

  // determine min fi for all frames
  int fimax_, chunk=0;
  long acoffset=0;
  //read candidates in chunks of given fi range
  while (++chunk) {
    printf("\n#########################\nchunk=%d ", chunk);
    cand_left=0;
    for(fr=1; fr<=frcount; ++fr) {
      fpar[fr].i2 = fpar[fr].ncands - 1; // the last record
      cand_left += fpar[fr].i2 - fpar[fr].i1 + 1;
    }
    if (cand_left > allcandi_size) {
      // find fi_min, set indices
      fimax_ = INT_MIN;
      for(fr=1; fr<=frcount; ++fr){
	fpar[fr].i2 = min(fpar[fr].i1 + chunk_size - 1, fpar[fr].ncands -1);
	offset = sizeof(int[fpar[fr].i2][ACLEN]);
	fseek(fpar[fr].candi_fh, offset, SEEK_SET);
	fread(&fpar[fr].fi2, sizeof(int), 1, fpar[fr].candi_fh);
	if (fpar[fr].fi2 > fimax_) fimax_ = fpar[fr].fi2;
      }
      printf("\nrewinding to max fi2 for all frames = %d\n", fimax_);
      // rewind to first record with fi > fimin_ , set i2
      for(fr=1; fr<=frcount; ++fr){
	printf(" frame[%d] i1=%ld  i2=%ld fi2=%d -> ", fr, fpar[fr].i1, fpar[fr].i2, fpar[fr].fi2);
	for (il=fpar[fr].i2; il>=fpar[fr].i1; --il){
	  offset = sizeof(int[il][ACLEN]);
	  fseek(fpar[fr].candi_fh, offset, SEEK_SET);
	  fread(&fpar[fr].fi2, sizeof(int), 1, fpar[fr].candi_fh);
	  //printf("\nfi2[%ld]=%d\n", il, fpar[fr].fi2); getchar();
	  if (fpar[fr].fi2 > fimax_){
	    fpar[fr].i2 = il;
	    break;
	  }
	}
	printf("i2=%ld fi2=%d\n", fpar[fr].i2, fpar[fr].fi2);
      }
	
    } else { // cand_left > chunk_size
      printf("cand_left = %ld < allcandi_size\n", cand_left);
      if (cand_left==0) break;
    }

    printf("frame[range] :");
    for(fr=1; fr<=frcount; ++fr)
      printf("%d[%ld,%ld] ", fpar[fr].num, fpar[fr].i1, fpar[fr].i2);

    // read candi in into allcandi
    acoffset = 0;
    for(fr=1; fr<=frcount; ++fr){
      fseek(fpar[fr].candi_fh, sizeof(int[fpar[fr].i1][ACLEN]), SEEK_SET);
      int chunk_ = fpar[fr].i2-fpar[fr].i1+1;
      //printf("\n[debug] frame %d: i1=%ld i2=%ld  size=%ld" , fr, fpar[fr].i1,  fpar[fr].i2, sizeof(int[chunk_][ACLEN]) );
      fread((void *)(&(allcandi[acoffset][0])), sizeof(int[chunk_][ACLEN]), 1, fpar[fr].candi_fh);
      printf("\ni = %10ld - %10ld  fi = %10d - %10d", acoffset,  acoffset+chunk_-1, allcandi[acoffset][0], allcandi[acoffset+chunk_-1][0]);
      acoffset += chunk_;
    }
    
    // prepare for the next chunk
    for(fr=1; fr<=frcount; ++fr) fpar[fr].i1 = fpar[fr].i2 + 1;
    
  } // chunk


  // remove temp files
  printf("\n[debug] removing temporary files\n");
  for(fr=1; fr<=frcount; ++fr){
    fclose(fpar[fr].candi_fh);
    unlink(fpar[fr].candi_fname);
  }


  // Looking for coincidences (the same integer values) between frames
  //------------------------------------------------------------------

  // Sorting the first 4 columns of allcandi
  qsort(allcandi, acoffset, sizeof(int[ACLEN]), way2compare_4c);

  int **imtr;
  int coindx=0, numc;
  unsigned short int weight=1;
  
  // Maximal possible amount of coincidences, given a threshold for 
  // a minimum number of interesting coincidences (opts->mincoin)  
  numc = acoffset/(opts->mincoin); 
  
  imtr = matrix(numc, 2); 
  
  // Coincidences: counting rows in a sorted table 
  for (i=0; i<(acoffset-1); i++) {
    
    int diff=1; 
    for(j=0; j<4; j++) 
      // using XOR: !(a^b) equals 1 for a=b
      // if((allcandi[i][j])^(allcandi[i+1][j])) { diff=0; break; } 
      diff *= !((allcandi[i][j])^(allcandi[i+1][j]));
    
    if(diff) { 
      
      weight++; 
      if(weight==opts->mincoin) // threshold value  
        coindx++; 
      
    } else {
      
      if(weight>=opts->mincoin) { 
        imtr[coindx][0] = i;
        imtr[coindx][1] = weight;
      }  

      weight=1; 
    }     
   
  } 

  // Sorting the coincidences table
  //-------------------------------
  qsort(imtr, numc, sizeof(int *), way2compare_c1);

  char outname[1040]; 
  sprintf(outname, "%s/%04d_%s.coi", opts->prefix, opts->shift, opts->trigname);
  data = fopen(outname, "w"); 
  FLOAT_TYPE trigf[TRLEN];

  FILE *coifh = fopen("coi2.txt", "w"); 

  int q; 
  for(q=0; q<coindx; q++) {  
 
    j = imtr[q][0];
    long int ops[256];
    unsigned short int l, w=imtr[q][1], fra[256];  
    double mean[5]; 
    float meanf[5]; 

    for(l=0; l<5; l++) mean[l]=0; 

    for(i=0; i<w; i++) {   
      int l, k = j-i; 
      fr = id2frcount[allcandi[k][5]];
      // ops[i]: position in trigger file #fra[i]
      ops[i] = allcandi[k][4];
      // frame id (not frcount!)
      fra[i] = (unsigned short int)allcandi[k][5]; 

      fseek(fpar[fr].trig_fh, sizeof(FLOAT_TYPE[ops[i]][TRLEN]), SEEK_SET);
      fread(trigf, sizeof(trigf), 1, fpar[fr].trig_fh);

      trigf[0] = trigf[0] + 2.*trigf[1]*(sett->N)*(opts->refr - allcandi[k][5]);

      for(l=0; l<4; l++)
        mean[l] += trigf[l]; 

//#mb definition for alpha (mean[3]) in MDC Stage4 
//      for(l=0; l<3; l++)  
//        mean[l] += allcandf[f][l];      
//
//     if(allcandf[f][3]>M_PI) 
//       mean[3] += 2*M_PI - arccos(cos(allcandf[f][3])); 
//         mean[3] += 2*M_PI - allcandf[f][3]; 
//       else 
//         mean[3] += allcandf[f][3];     

      mean[4] += trigf[4]*trigf[4];

    }
 
    for(l=0; l<4; l++) 
      mean[l] /= w; 

    mean[4] = sqrt(mean[4]); // SNR mean: sqrt of sum of squares  

    for(l=0; l<5; l++) 
      meanf[l] = (float)mean[l]; 

    // writing to binary file 

    fwrite(&w, sizeof(unsigned short int), 1, data); 
    fwrite(&meanf, sizeof(float), 5, data);          
    fwrite(&fra, sizeof(unsigned short int), w, data); 
    fwrite(&ops, sizeof(int), w, data); 

    fprintf(coifh, "%d %15.8e %15.8e %15.8e %15.8e %15.8e", 
	   w, mean[0], mean[1], mean[2], mean[3], mean[4]);
    for(i=0; i<w; ++i)
      fprintf(coifh, " %d", fra[i]); 
    for(i=0; i<w; ++i)
      fprintf(coifh, " %ld", ops[i]); 
    fprintf(coifh, "\n"); 


    // Maximal coincidence (first row of imtr[][])
    if(!q) {
      fprintf(stderr, "%s %04d %5f %5hu %5d %15.8le %5.8le %5.8le %5.8le %5le ", 
	      opts->trigname, opts->shift, sett->fpo, frcount, w,   
	      mean[0], mean[1], mean[2], mean[3], mean[4]);

      int ii, jj;
      // Number of candidates from frames that participated in the coincidence 
      for(ii=0; ii<=frcount; ii++)
        for(jj=0; jj<w; jj++)
          if(fpar[ii+1].num == fra[jj]) {
	    fprintf(stderr, "%d %ld %ld ",
		    fra[jj], fpar[ii+1].ninband, fpar[ii+1].ncands);
	    break;
          }
      fprintf(stderr, "\n");  
    }
  } // q

  // close files
  fclose(data);
  fclose(coifh);
  for(fr=1; fr<=frcount; ++fr) fclose(fpar[fr].trig_fh);

  // free memory
  free(allcandi);
  for(i=0; i<numc; i++) 
    free(imtr[i]); 
  free(imtr); 

  exit(EXIT_SUCCESS);

}
