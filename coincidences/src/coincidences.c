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

#ifndef DTAPREFIX
#define DTAPREFIX ./candidates
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
  manage_grid_matrix(&sett, &opts);	

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
  
  int i, j, candsize=INICANDSIZE, allcandsize=INICANDSIZE, goodcands=0, current_frame=0, frcount=0;  
  int val, shift[4], scale[4]; 
  int hemi;
  double sqrN, omsN, v[4][4], be[2];
  FLOAT_TYPE tmp[4], c[5];

  char dirname[512], filename[512], outname[512];  
  // Trigger files directory name 
  sprintf (dirname, "%s", opts->dtaprefix); 

  DIR *dp;
  struct dirent *ep;
  FILE *data; 

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

      if(((ep->d_type == DT_REG) || (ep->d_type == DT_LNK)) &&
        (strstr(ep->d_name, opts->trigname) != NULL)) {

        sprintf(filename, "%s/%s", opts->dtaprefix, ep->d_name);

        // This part looks for the first number in the trigger file name,
        // under the assumption that this is the frame number
        char *fr; 
        char *epdname = ep->d_name;  

        while((fr = strsep(&epdname, "_"))!=NULL) {
          if(fr[0] >= '0' && fr[0] <= '9') {
            current_frame = atoi(fr);
            printf("Reading %s... Frame %d: ", ep->d_name, current_frame);
            break; 
          }
        }

        if((data = fopen(filename, "r")) != NULL) {

          FLOAT_TYPE finband;
          // Original candidate number (position in the trigger file)
          int orgpos=0;
          // Counter for 'good' candidates i.e. these that are in band
          i=0; 
          frcount++;

          // Each candidate is represented by 5 FLOAT_TYPE (double or float) numbers
          // c[0]=f, c[1]=s, c[2]=d, c[3]=a, c[4]=snr
          while(fread((void *)c, sizeof(FLOAT_TYPE), 5, data)==5) {  

            //Narrowing-down the band around center  
            if((c[0] > M_PI_2 - opts->narrowdown) && (c[0] < M_PI_2 + opts->narrowdown)) {

             // shifting c[0] (=frequency) to opts->refr reference frame 
             c[0] = c[0] + 2.*c[1]*(sett->N)*(opts->refr - current_frame); 

              // #mb todo: deal with the out-of-band candidates 
              // if frequency is in band 
              // c[4] = 4.0620192023179804 corresponds to F-stat = 10.25
              // c[4] = 4.1231056256176606 corresponds to F-stat = 10.5
              // c[4] = 4.2426406871192848 corresponds to F-stat = 11 
              // c[4] = 4.4721359549995796 corresponds to F-stat = 12
              // c[4] = 4.5825756949558398 corresponds to F-stat = 12.5
              // c[4] = 5.0990195135927845 corresponds to F-stat = 15
              // because
              // c[4] = rho = \sqrt{2(F-2)}
              if(((c[0]>0) && (c[0]<M_PI)) && (c[4] > opts->snrcutoff)) {

                // Conversion to linear parameters
                //--------------------------------

                tmp[0] = c[0]*sett->N;
                tmp[1] = c[1]*sqrN;

                // Transformation of astronomical to linear coordinates;
                // C_EPSMA, an average value of epsm, is defined in settings.h
                hemi = ast2lin(c[3], c[2], C_EPSMA, be);

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

                // Saving the original position, frame number and current index
                candi[i][4] = orgpos;
                candi[i][5] = current_frame;
                candi[i][6] = i;
                // Saving the SNR value
                candf[i][4] = c[4];
                i++;

             } // if finband

            } // if narrowdown

            orgpos++;

            // Resizing the candidates' array, if the previous limit is reached
            // (realloc by a factor of 2)
            if(i==candsize) {

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
          printf("%d/%d\n", trig->frameinfo[frcount][2], trig->frameinfo[frcount][1]);

        } else { 
          printf("Problem with %s...\n", filename);  
          perror (filename);
        }

      memset(filename, 0, sizeof(filename)); 
      fclose(data); 

      } // if(((ep->d_type == DT_REG) ...

    } // while ((ep = readdir (dp))) 

  } // if (dp != NULL)  

  (void) closedir(dp);

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
  sprintf(outname, "%s/%04d_%04d_%s.coi", 
    opts->prefix, opts->shift, opts->scale, opts->trigname);
  data = fopen(outname, "w"); 

  int q, maxcoin = imtr[0][1], maxcoinindex=0; 
  double maxsnr=0; 
  for(q=0; q<coindx; q++) {  
 
    j = imtr[q][0];
    int ops[256];
    unsigned short int l, w=imtr[q][1], fra[256];  
    double mean[5]; 
    float meanf[5]; 

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
 
    fprintf(stderr, "%s %04d %5f %5hu %5d %15.8le %5.8le %5.8le %5.8le %5le ", 
      opts->trigname, opts->shift, sett->fpo, trig->frcount, maxcoin,   
      mean[0], mean[1], mean[2], mean[3], mean[4]);

    // info about time segments: frame number, all, unique candidates 
    for(i=1; i<=trig->frcount; i++) 
      fprintf(stderr, "%d %d %d ", trig->frameinfo[i][0], trig->frameinfo[i][1], trig->frameinfo[i][2]); 

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
