/*

Code simply joins together any bin files.
There should be list_file with list of 
frames that will be glued.
It takes arguments from command line:
list_path, datapath, outpath, name_of_detector name_of_file_to_glue

For example:
./glue2 list.txt /work/psk/irods/Resc/home/public/groups/polgraw/rdc_6d_xdat_0.25 /home/msieniawska/tests/newglue/data/glued/glued_data/002 H1 rSSB.bin

Where list.txt contains:
001
002
003
004

MS
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv){

	char listpath[512];
	char datapath[512];
	char outpath[512];
	char detector[512];
	char filetoglue[512];

	strcpy(listpath, argv[1]);
	strcpy(datapath, argv[2]);
	strcpy(outpath, argv[3]);
	strcpy(detector, argv[4]);
	strcpy(filetoglue, argv[5]);
	FILE *list; 
	FILE *fp;
	FILE *output;
	int line;
	int j;
	int data = 1;
	char path[512];
	char out[512];
	double *cand;
	cand = (double *) calloc (data, sizeof(double));
	sprintf(out, "%s/%s/%s", outpath, detector, filetoglue);

	if ((list = fopen (listpath, "r")) != NULL) {
		while (fscanf(list, "%d\n", &line) == 1){
			sprintf(path, "%s/%03d/%s/%s", datapath, line, detector, filetoglue);
			if ((fp = fopen (path, "rb")) != NULL) {
				if ((output = fopen (out, "ab")) != NULL) {			
					while (!feof(fp)) {
						if((fread ((void *)(cand), sizeof (double), data, fp))==data){
							for(j = 0; j < data; j++) fwrite((void *)(cand), sizeof (double), data, output);
						}
					}
				}
			}

		}
	}

}




