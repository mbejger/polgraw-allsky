// gcc -o read_coi_and_triggers read_coi_and_triggers.c -lm 
// $ ./read_coi_and_triggers coi_file path_to_triggers band reffr N

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv) { 
	
	int i, j, N, ops[256]; 
  unsigned short int w, whigh=0, refr, fra[256];  
  char x[512]; 
	FILE *cdata, *tdata;  
	float mean[5]; 
  double cand[5], meanf[5]; 
    	
  refr = atoi(argv[4]); 
  N    = atoi(argv[5]); 

  // reading from the input datafile   
	if ((cdata = fopen(argv[1], "rb")) != NULL) {

    while(!feof(cdata)) {

      // end if no more data 
      if(!fread(&w, sizeof(unsigned short int), 1, cdata)) { break; }
    
      // first pass, set up the value of the highest coincidence 
      if(!whigh) { whigh = w; }  
    
      // go out of the loop if coincidences lower than whigh appear  
      if(w < whigh) { break; } 

      // Reading the highest coincidence (first one) 
      fread(&mean, sizeof(float), 5, cdata);
      fread(&fra, sizeof(unsigned short int), w, cdata); 
      fread(&ops, sizeof(int), w, cdata);

      // clearing the test array
      for(i=0; i<5; i++) meanf[i] = 0; 

      printf("# Individual triggers (frame, f, s, d, a, snr):\n"); 

      // Reading the candidates that entered the coincidence   
      for(i=0; i<w; i++) { 

        sprintf(x, "%s/triggers_%03d_%04d_2.bin", argv[2], fra[i], atoi(argv[3]));
 
        if((tdata = fopen(x, "rb")) != NULL) {

          fseek(tdata, 5*sizeof(double)*ops[i], SEEK_SET); 
          fread(&cand, sizeof(double), 5, tdata);

          meanf[0] += cand[0] + 2*cand[1]*N*(refr - fra[i]);
          for(j=1; j<4; j++) meanf[j] += cand[j]; 
          meanf[4] += cand[4]*cand[4]; 

          printf("    %03d ", fra[i]); 
          for(j=0; j<5; j++) printf("%le ", cand[j]); 
          printf("\n"); 

        } else { 

          perror (x);
          return 1;

        }  

        fclose(tdata); 
    
      }

      for(i=0; i<4; i++) meanf[i] /= w;
      meanf[4] = sqrt(meanf[4]); 

      printf("\033[32mmean    ");  
      for(i=0; i<5; i++) printf("%le ", meanf[i]);  
      printf("\033[0m\n"); 

      printf("\033[31mcoin\033[0m%3hu \033[31m", w);
      for(i=0; i<5; i++) printf("%le ", mean[i]);
      printf("\033[0m"); 
      for(i=0; i<w; i++) printf("%hu:%d ", fra[i], ops[i]);
      printf("\n\n"); 


    } 

	  } else {
		
		  perror (argv[1]);
		  return 1;
	  }

  fclose(cdata); 

	return 0; 

}
