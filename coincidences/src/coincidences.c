#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <getopt.h>
#include <dirent.h>

#include "auxi.h"
#include "init.h" 
#include "settings.h"


// Default output and data directories
#ifndef PREFIX
#define PREFIX ./coinc-results
#endif

#define REALLOC_FACTOR 1.2

int main (int argc, char* argv[]) {

  Search_settings sett;
  Command_line_opts_coinc opts;
  Candidate_triggers trig; 
  int i, j; 

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
  manage_grid_matrix(&sett, opts.refgrid);
  
  // Search settings 
  search_settings(&sett); 

  printf("Settings dt: %f, oms: %f\n", sett.dt, sett.oms); 
  printf("Reference frame number: %d\n", opts.refr); 

  read_trigger_files(&sett, &opts, &trig); 

  // Free arrays at the end 
  free(sett.M); 
  printf("\nEND\n");

  return 0; 
	
}


static int way2compare_4c(const void *a, const void *b){
  /* compare 4 columns (0,1,2,3) */

  const int* x = *((const int **)a);
  const int* y = *((const int **)b);
  
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
  
  int i, j, candsize=INICANDSIZE, allcandsize=INICANDSIZE,
       goodcands=0, current_frame=0, frcount=0;  
  int val, shift[4], scale[4]; 
  double sqrN, omsN, v[4][4], be[2];
  FLOAT_TYPE tmp[4], c[5];
  char dirname[512], filename[1024], outname[1200];
  FILE *data, *infile, *vfile;

  int **candi = malloc(candsize*sizeof(int *)); 
  int **ti; 
  // 7 columns in a row: fi, si, di, ai, orgpos, fr, i
  for(i=0; i<candsize; i++) 
       candi[i] = (int*)calloc(7, sizeof(int));

  int **allcandi = malloc(allcandsize*sizeof(int *)); 
  // 7 columns in a row: fi, si, di, ai, orgpos, fr, i
  for(i=0; i<allcandsize; i++) 
       allcandi[i] = malloc(7*sizeof(int));

  FLOAT_TYPE **candf = malloc(candsize*sizeof(FLOAT_TYPE *)); 
  FLOAT_TYPE **tf; 
  // 5 columns in a row: f, s, d, a, SNR
  for(i=0; i<candsize; i++) 
       candf[i] = malloc(5*sizeof(FLOAT_TYPE));

  FLOAT_TYPE **allcandf = malloc(allcandsize*sizeof(FLOAT_TYPE *));
  // 5 columns in a row: f, s, d, a, SNR
  for(i=0; i<allcandsize; i++)
       allcandf[i] = malloc(5*sizeof(FLOAT_TYPE));


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

  // Scaling values from opts->scale[f, s, d, a]
  scale[0] = opts->scalef; 
  scale[1] = opts->scales; 
  scale[2] = opts->scaled; 
  scale[3] = opts->scalea; 

  sqrN = pow(sett->N, 2);
  omsN = sett->oms*sett->N; 

  printf("Cell scaling in f, s, d, a directions: ");
  printf("%d %d %d %d\n", scale[0], scale[1], scale[2], scale[3]);


  // Transformation matrix elements: division by the corresponding 
  // scale factors outside the main loop over candidates    
  for(j=0; j<4; j++)  
       for(i=0; i<4; i++)
	    v[i][j] = (sett->vedva[i][j])/scale[j];

  printf("Reading trigger filenames from %s\n", opts->infile);
  
  if((infile = fopen(opts->infile, "r")) == NULL) {
       printf("Error opening %s", opts->infile);
       exit(EXIT_FAILURE);
  }
  char line[1024], vfilename[1040];
  char bandstr[5], frstr[4], bndstr[5], hstr[2], *basename_pos;
  int curerent_frame, current_band, current_hemi, nhemi;
  
  while(fgets(line, sizeof(line), infile)) {
       if( line[0] == '%') continue;
       strcpy(filename, strtok(line," \r\n"));
       strcpy(vfilename, strtok(NULL," \r\n"));

       printf("%s", filename);
       
       // determine current band, frame and hemisphere from filename
       // current_band will be used in future to handle out of band candidates
       // all files must have the same hemisphere
       basename_pos = strrchr(filename, '/');
       strncpy(frstr, basename_pos+10, 3);
       current_frame = atoi(frstr);
       
       strncpy(bndstr, basename_pos+14, 4);
       current_band = atoi(bndstr);

       char trigfilename[2][1024];
       strncpy(hstr, basename_pos+19, 1);
       if (hstr[0]=='*') {
	    nhemi = 2;
	    opts->hemi = 0;
	    (basename_pos+19)[0] = '1';
	    strcpy(trigfilename[0], filename);
	    (basename_pos+19)[0] = '2';
	    strcpy(trigfilename[1], filename);
       } else {
	    nhemi = 1;
	    current_hemi = atoi(hstr);
	    if (opts->hemi == -1) opts->hemi = current_hemi;
	    if (current_hemi != opts->hemi) {
		 printf("Change of hemisphere detected for %s \nExiting...\n", filename);
		 exit(EXIT_FAILURE);
	    }
	    strcpy(trigfilename[0], filename);
       }

       // Read veto lines file
       int nlines=0;
       if((vfile = fopen(vfilename, "r")) == NULL) {
	    printf(" [-veto] ");
       } else {
	    printf(" [+veto] ");
	    while(fgets(line, sizeof(line), vfile)) {
		 if( line[0] == '%') continue;
		 double f1, f2;
		 sscanf(line, "%lf %lf %lf %lf", &f1, &f2,
			&sett->lines[nlines][0], &sett->lines[nlines][1]);
		 //printf("Line: %f  %f\n", sett->lines[nlines][0], sett->lines[nlines][1]);
		 nlines++;
	    }
	    fclose(vfile);	    
       }


       // Counter for 'good' candidates i.e. these that are in band
       i=0;
       frcount++;

       for(int it=0; it<nhemi; it++){
	    
       	    current_hemi = it+1;
	    if((data = fopen(trigfilename[it], "r")) == NULL) {
		 printf("Can't open %s\n", filename);
		 exit(EXIT_FAILURE);
	    }

	    printf("\nfilename[%d]=%s\n", it, trigfilename[it]);

	    // Original candidate number (position in the trigger file)
	    int orgpos=-1;

       // Each candidate is represented by 5 FLOAT_TYPE (double or float) numbers
       // c[0]=f, c[1]=s, c[2]=d, c[3]=a, c[4]=snr
       while(fread((void *)c, sizeof(FLOAT_TYPE), 5, data)==5) {  
	    
	    orgpos++;
			
            // Narrowing-down the band around center  
            if( (c[0] <= M_PI_2 - opts->narrowdown) ||
		(c[0] >= M_PI_2 + opts->narrowdown) )  continue;

	    // test if trigger is in line
	    int isline=0;
	    for(int il=0; il<nlines; il++){
		 if( c[0] >= sett->lines[il][0] && c[0] <= sett->lines[il][1] ) {
		      isline = 1;
		      break;
		 }
	    }
	    if (isline) continue;
	    
	    // shifting c[0] (=frequency) to opts->refr reference frame 
	    c[0] = c[0] + 2.*c[1]*(sett->N)*(opts->refr - current_frame);
	    //printf("shifting: %d  %d  %g\n", opts->refr, current_frame, 2.*c[1]*(sett->N)*(opts->refr - current_frame));

	    // #mb todo: deal with the out-of-band candidates 
	    // c[4] = rho = \sqrt{2(F-2)}

	    if(((c[0]>0) && (c[0]<M_PI)) && (c[4] > opts->snrcutoff)) {

		 // Conversion to linear parameters
		 //--------------------------------
		      
		 tmp[0] = c[0]*sett->N;
		 tmp[1] = c[1]*sqrN;

		 // Transformation of astronomical to linear coordinates;
		 // C_EPSMA, an average value of epsm, is defined in settings.h
		 int hemi;
		 hemi = ast2lin(c[3], c[2], C_EPSMA, be);
		 // pci: we assume that hemisphere do not change when moving to the ref frame

		 // tmp[2] corresponds to declination (d), tmp[3] to right ascension (a)
		 tmp[2] = omsN*be[0];
		 tmp[3] = omsN*be[1];

		 // Saving candidate values
		 for(j=0; j<4; j++) {
			   
		      // Integer values (0=fi, 1=si, 2=di, 3=ai)
		      candi[i][j] = round(tmp[0]*v[0][j] + tmp[1]*v[1][j]
		      		  + tmp[2]*v[2][j] + tmp[3]*v[3][j] 
		      		  + 0.5*shift[j]);
    
		      // Astrophysical values (0=f, 1=s, 2=d, 3=a)
		      // f is shifted to opts->refr time frame
		      candf[i][j] = c[j];

		 }
		 //if (candi[i][3]==-2 && fabs(c[2])<0.2) printf("ra! %g %g  %d\n", c[3], tmp[3]*v[3][3], candi[i][3]);
		 //if (c[3]>6. || c[3]<0.2) printf("ra! %g %g  %d\n", c[3], tmp[3]*v[3][3], candi[i][3]);
		 // Saving the original position, frame number and current index
		 if (current_hemi == 1)
		      candi[i][4] = orgpos;
		 else
		      candi[i][4] = -orgpos;
		 
		 candi[i][5] = current_frame;
		 candi[i][6] = i;
		 // Saving the SNR value
		 candf[i][4] = c[4];
		 i++;

	    } // if inband
		 



            // Resizing the candidates' array, if the previous limit is reached
            // (realloc by a factor of 2)
            if( i==candsize ) {

		 candsize *= REALLOC_FACTOR;

		 ti = realloc(candi, candsize*sizeof(int *)); 
		 if(ti!=NULL) { 
		      
		      candi = ti; 
		      for(j=i; j<candsize; j++)
			   candi[j] = malloc(7*sizeof(int));
		      
		 } else { 
		      printf("Problem with memory realloc for candidates array (int)... exiting...\n");
		      exit(EXIT_FAILURE);
		 }
		 
		 tf = realloc(candf, candsize*sizeof(FLOAT_TYPE *)); 
		 if(tf!=NULL) { 
		      
		      candf = tf; 
		      
		      for(j=i; j<candsize; j++)
			   candf[j] = malloc(5*sizeof(FLOAT_TYPE));
		      
		 } else { 
		      printf("Problem with memory realloc for candidates array (astro)... exiting...\n"); 
		      exit(EXIT_FAILURE);
		 }
		 
            } // candsize realloc 

       } // while fread
       
       fclose(data);
       printf("Read %d inband candidates\n", i);
#if 0
       char crfname[36];
       sprintf(crfname, "candi_%03d_%04d_%1d.bin", current_frame, current_band, current_hemi);
       data = fopen(crfname, "w");
       double cc[5];
       for (int ii=0; ii<i; ii++) {
	    cc[0] = (double)(candi[ii][0]);
	    cc[1] = (double)(candi[ii][1]);
	    cc[2] = (double)(candi[ii][2]);
	    cc[3] = (double)(candi[ii][3]);
	    cc[4] = (double)(candf[ii][4]);	    
	    fwrite( cc, sizeof(double), 5, data);
       }
       fclose(data);
#endif           
       } // for over hemispheres

// save candi (linear coords, shifted to ref frame)
#if 0
       char crfname[36];
       sprintf(crfname, "candi_%03d_%04d_%d.bin", current_frame, current_band, opts->hemi);
       data = fopen(crfname, "w");
       double cc[5];
       for (int ii=0; ii<i; ii++) {
	    cc[0] = (double)(candi[ii][0]);
	    cc[1] = (double)(candi[ii][1]);
	    cc[2] = (double)(candi[ii][2]);
	    cc[3] = (double)(candi[ii][3]);
	    cc[4] = (double)(candf[ii][4]);
	    fwrite( cc, sizeof(double), 5, data);
       }
       fclose(data);
#endif           
       
       // Frame number  
       trig->frameinfo[frcount][0] = current_frame;
       // Number of candidates in band for a given frame 
       trig->frameinfo[frcount][1] = i;  
       
       
       // Looking for duplicates and selecting the one with highest SNR
       //--------------------------------------------------------------
       
       // Sorting the first 4 columns of candi
       qsort(candi, trig->frameinfo[frcount][1], sizeof(int *), way2compare_4c);
       
       int maxsnridx=0, frgoodcands=0;
       double candsnr=0;  
       for (i=0; i<trig->frameinfo[frcount][1]; i++) {
	    
            int idx, idx1, maxi, diff=1;  
            for(j=0; j<4; j++)
		 // using XOR: !(a^b) equals 1 for a=b
		 //if((candi[i][j])^(candi[i+1][j])) { diff=0; break; }             
		 diff *= !((candi[i][j])^(candi[i+1][j]));
	    
            idx = candi[i][6];
	    
            if(!diff) {

		 int k=i, kidx=idx;
		 if(maxsnridx) { k=maxi; kidx=maxsnridx; }
		 
		 // Writing to array containing all candidates 
		 for(j=0; j<6; j++)
		      allcandi[goodcands][j] = candi[k][j];
		 allcandi[goodcands][6] = goodcands; 
		 
		 for(j=0; j<5; j++)
		      allcandf[goodcands][j] = candf[kidx][j];

		 
		 maxsnridx=0;
		 goodcands++;
		 frgoodcands++;
		 
		 if(goodcands==allcandsize) {
		      
		      allcandsize *= REALLOC_FACTOR;
		      
		      ti = realloc(allcandi, allcandsize*sizeof(int *)); 
		      if(ti!=NULL) {
			   allcandi = ti;
			   for(j=goodcands; j<allcandsize; j++)
				allcandi[j] = malloc(7*sizeof(int));
		      } else { 
			   printf("Problem with memory realloc for ALL candidates array (int)... exiting...\n");
			   exit(EXIT_FAILURE);
		      }
		      
		      tf = realloc(allcandf, allcandsize*sizeof(FLOAT_TYPE *)); 
		      if(tf!=NULL) { 
			   allcandf = tf; 
			   for(j=goodcands; j<allcandsize; j++)
				allcandf[j] = malloc(5*sizeof(FLOAT_TYPE));
		      } else { 
			   printf("Problem with memory realloc for ALL candidates array (astro)... exiting...\n");
			   exit(EXIT_FAILURE);
		      } 
		      
		 }
		 // The candidate is not unique, selecting the one with the highest SNR
            } else {
		 
		 idx1 = candi[i+1][6];

		 if(!maxsnridx) {  

		      maxsnridx = (candf[idx][4] > candf[idx1][4] ? idx : idx1);  
		      maxi = (candf[idx][4] > candf[idx1][4] ? i : i+1);
		      candsnr = candf[maxsnridx][4];    
		      
		 } else {
		      
	 	      if(candf[idx][4] > candsnr) {
			   maxsnridx = idx; maxi = i; 
			   candsnr = candf[idx][4]; 
		      } else if(candf[idx1][4] > candsnr) {
			   maxsnridx = idx1; maxi = i+1; 
			   candsnr = candf[idx1][4];
		      }
		 }
            }
       }
       
       // Number of unique candidates in a given frame
       trig->frameinfo[frcount][2] = frgoodcands;
       printf("%d/%d  [unique/inband]\n", trig->frameinfo[frcount][2], trig->frameinfo[frcount][1]);

       memset(filename, 0, sizeof(filename)); 
       
  } //   while(fgets(line, sizeof(line), infile)) 

  trig->frcount = frcount;
  trig->goodcands = goodcands;

  printf("Total number of candidates from all frames: %d\n", trig->goodcands);

    
  // Looking for coincidences (the same integer values) among different frames
  //--------------------------------------------------------------------------

  // Sorting the first 4 columns of allcandi
  qsort(allcandi, trig->goodcands, sizeof(int *), way2compare_4c);

  int **imtr;
  int coindx=0, numc;
  unsigned short int weight=1, maxweight=0;   
  
  // Maximal possible amount of coincidences, given a threshold  
  // for a minimum number of interesting coincidences (opts->mincoin)  
  numc = (trig->goodcands)/(opts->mincoin); 
  
  imtr = matrix(numc, 2); 
  
  // Coincidences: counting rows in a sorted table 
  for (i=0; i<(trig->goodcands-1); i++) {

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

  // Coincidences above opts->mincoin threshold 
  //-------------------------------------------
  memset(outname, 0, sizeof(outname));
  sprintf(outname, "%s/%04d_%d-%d-%d-%d_%04d_%1d.coi", 
	  opts->prefix, opts->shift, opts->scalef, opts->scales, opts->scaled, opts->scalea, opts->band, opts->hemi);
  data = fopen(outname, "w"); 
  printf("coindx=%d\n", coindx);
  int q, maxcoin = imtr[0][1], maxcoinindex=0; 
  double maxsnr=0; 
  for(q=0; q<coindx; q++) {
       
    j = imtr[q][0];
    int ops[256];
    unsigned short int l, w=imtr[q][1], fra[256];  
    double mean[5]; 
    float meanf[5];
    //if ( j == 1425307 ) printf("Ha5! %g %g %g \n", allcandf[j][0], allcandf[j][1], allcandf[j][4]);     
    //if ( allcandf[j][4] > 5) printf("Ha5! %g %g %g \n", allcandf[j][0], allcandf[j][1], allcandf[j][4]);     

    for(l=0; l<5; l++) mean[l]=0; 

    for(i=0; i<w; i++) {   
      int l, k = j-i; 
      int f = allcandi[k][6]; 

//#mb 
      for(l=0; l<4; l++)  
        mean[l] += allcandf[f][l]; 

//#mb definition for alpha (mean[3]) in MDC Stage4 
//      for(l=0; l<3; l++)  
//        mean[l] += allcandf[f][l];      
//
//     if(allcandf[f][3]>M_PI) 
//       mean[3] += 2*M_PI - arccos(cos(allcandf[f][3])); 
//         mean[3] += 2*M_PI - allcandf[f][3]; 
//       else 
//         mean[3] += allcandf[f][3];     

      //if ( allcandf[f][4] > 10) printf("Ha6! %g %g \n", allcandf[f][0], allcandf[f][4]);     

      mean[4] += allcandf[f][4]*allcandf[f][4];  

      // ops[i]: position in trigger file #fra[i]
      ops[i] = allcandi[k][4]; 
      fra[i] = (unsigned short int)allcandi[k][5]; 

    }
 
    for(l=0; l<4; l++) mean[l] /= w;  

    mean[4] = sqrt(mean[4]); // SNR mean: sqrt of sum of squares  

    for(l=0; l<5; l++) meanf[l] = (float)mean[l]; 

    // writing to binary file 
    fwrite(&w, sizeof(unsigned short int), 1, data); 
    fwrite(&meanf, sizeof(float), 5, data);          
    fwrite(&fra, sizeof(unsigned short int), w, data); 
    fwrite(&ops, sizeof(int), w, data); 

    // Looking for a maximal coincidence with a maximal snr
    if(!q) 
      maxsnr = mean[4]; 

    if(w==maxcoin && mean[4] > maxsnr) { 
      maxsnr = mean[4]; 
      maxcoinindex = q; 
    } 

  }

  fclose(data); 

  if(maxcoin >= opts->mincoin) { // Writing out the maximal coincidence with maximal snr to stderr 

    int ops[256];
    unsigned short int l, fra[256];  
    double mean[5]; 

    for(l=0; l<5; l++) mean[l]=0; 

    for(i=0; i<maxcoin; i++) {   
      int l, k = imtr[maxcoinindex][0] - i; 
      int f = allcandi[k][6]; 
  
      for(l=0; l<4; l++)  
        mean[l] += allcandf[f][l]; 

      mean[4] += allcandf[f][4]*allcandf[f][4];  

      // ops[i]: position in trigger file #fra[i]
      ops[i] = allcandi[k][4]; 
      fra[i] = (unsigned short int)allcandi[k][5]; 

    }
 
    for(l=0; l<4; l++) mean[l] /= maxcoin;  

    mean[4] = sqrt(mean[4]); // SNR mean: sqrt of sum of squares  
 
    fprintf(stderr, "%04d_%1d %04d %5f %5hu %5d %15.8le %5.8le %5.8le %5.8le %5le ", 
	    opts->band, opts->hemi, opts->shift, sett->fpo, trig->frcount, maxcoin,   
	    mean[0], mean[1], mean[2], mean[3], mean[4]);

    // info about time segments: frame number, all, unique candidates 
    for(i=1; i<=trig->frcount; i++) 
	 fprintf(stderr, "%d %d %d ", trig->frameinfo[i][0],
		 trig->frameinfo[i][1], trig->frameinfo[i][2]); 

    // Frame numbers participating in the coincidence 
    for(i=0; i<maxcoin; i++)
      fprintf(stderr, "%d ", fra[i]); 

    fprintf(stderr, "\n"); 
 
  } 

  // Freeing auxiliary arrays at the end 
  for(i=0; i<candsize; i++) { 
    free(candi[i]);
    free(candf[i]); 
  } 
  free(candi);
  free(candf);
   
  for(i=0; i<allcandsize; i++) { 
    free(allcandi[i]);
    free(allcandf[i]);
  }
  free(allcandi);
  free(allcandf);

  for(i=0; i<numc; i++) 
    free(imtr[i]); 
  free(imtr); 

}
