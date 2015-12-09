#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) { 
	
	int i, pari[4], ops[256]; 
    unsigned short int w, fra[256];  
	FILE *data;  
	float mean[5]; 
    	
    // reading from the input datafile   
	if ((data = fopen(argv[1], "rb")) != NULL) {

		 while(!feof(data)) {

          if(!fread(&w, sizeof(unsigned short int), 1, data)) { break; } 
//          fread(&pari, sizeof(int), 4, data);
          fread(&mean, sizeof(float), 5, data);
          fread(&fra, sizeof(unsigned short int), w, data); 
          fread(&ops, sizeof(int), w, data);

          printf("%hu ", w);
//          for(i=0; i<4; i++) printf("%d ", pari[i]); 
          for(i=0; i<5; i++) printf("%le ", mean[i]);
          for(i=0; i<w; i++) printf("%hu:%d ", fra[i], ops[i]);
          printf("\n"); 

 }  


	} else {
		
		perror (argv[1]);
		return 1;
	}

    fclose(data); 
    
	return 0; 
}
