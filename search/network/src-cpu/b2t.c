#include <stdio.h>
#include <stdlib.h>
#define FLOAT double

int main(int argc, char **argv) { 

  
  int i, j=0, data=5 ; 
  FILE *datain ;  
  FLOAT *cand ; 
  
  cand 		= (FLOAT *) calloc (data, sizeof(FLOAT));
  
  // reading from the input datafile   
  if ((datain = fopen (argv[1], "rb")) != NULL) {
    
    while (!feof(datain)) {
      
      if((fread ((void *)(cand), sizeof (FLOAT), data, datain))==data)
	
	for(i=0; i<data; i++) 
	  printf("%.8le ", cand[i]) ;
      
      printf("\n") ;
      
    }
    
  } else {
    
    perror (argv[1]);
    return 1;
  }
  
  fclose(datain) ; 
    
  return 0 ; 
}
