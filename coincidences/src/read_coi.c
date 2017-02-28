#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) { 
	
	int i, ops[256]; 
    unsigned short int w, fra[256];  
	FILE *data;  
	float mean[5]; 
    	
  // reading from the input datafile   
	if ((data = fopen(argv[1], "rb")) != NULL) {

    printf("# num_of_coincidences    mean_val_of_pars (f, s, d, a), snr    frame_num:trigger_num_in_trigger_file\n");
    printf("# (see http://mbejger.github.io/polgraw-allsky/coincidences for details)\n");  

    while(!feof(data)) {

      if(!fread(&w, sizeof(unsigned short int), 1, data)) { break; } 
        fread(&mean, sizeof(float), 5, data);
        fread(&fra, sizeof(unsigned short int), w, data); 
        fread(&ops, sizeof(int), w, data);

        printf("%hu ", w);
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
