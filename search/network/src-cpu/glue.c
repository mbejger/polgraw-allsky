/* 
Ver. 1.0
Function joins together DetSSB files
(last 2 values from the first file are moved
at the end of final file) and xdatc* with given 
label name. Programm uses list.txt (from data dir), 
where should exist list of interesting times. Final 
files are in output_dir/followup_total_data. 
MS
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "glue.h"

int glue(char prfx[512], char dataprfx[512], char band[512], int nod){
	FILE *list; 
	FILE *fp;
	FILE *output;
	int i, j, k, l, flag = 0;
	int data = 1;
	int line;
	int base_size = 258493; //size of DetSSB.bin (in bytes) for nod = 1
	char xdat[200];
	char listpath[200];
	char path[200];
	char part[100];
	char out[200];
	char output_tot[512];
	char output_L[512];
	char output_H[512];
	char output_0[512];
	float *cand;
	float *ssb1, *ssb2;
	cand = (float *) calloc (data, sizeof(float));
	ssb1 = (float *) calloc (data, sizeof(float));
	ssb2 = (float *) calloc (data, sizeof(float));
	sprintf(listpath, "%s/list.txt", dataprfx);

	
	if ((list = fopen (listpath, "r")) != NULL) {
		l = 0;
		while (fscanf(list, "%d\n", &line) == 1){
			if (l != 0) flag = 1;
			l++;
			sprintf(xdat, "xdatc_%02d%s.bin", line, band);
			mkdir(prfx, 0777);
			sprintf(output_tot, "%s/followup_total_data/", prfx);
			mkdir(output_tot, 0777);
			sprintf(output_0, "%s/followup_total_data/000", prfx);
			mkdir(output_0, 0777);
			sprintf(output_H, "%s/H1/", output_0);
			mkdir(output_H, 0777);
			sprintf(output_L, "%s/L1/", output_0);
			mkdir(output_L, 0777);
			for(i = 0; i < 4; i++){
				if(i == 0){
					sprintf(part, "/L1/DetSSB.bin");
					sprintf(out, "%s%s", output_0, part);
				}
				if(i == 1){
					sprintf(part, "/H1/DetSSB.bin");
					sprintf(out, "%s%s", output_0, part);
				}
				if(i == 2){
					sprintf(part, "/L1/%s", xdat);
					sprintf(out, "%s/L1/xdatc_000%s.bin", output_0, band);
				}
				if(i == 3){
					sprintf(part, "/H1/%s", xdat);
					sprintf(out, "%s/H1/xdatc_000%s.bin", output_0, band);
				}
				sprintf(path, "%s/%02d%s", dataprfx, line, part);
				if ((fp = fopen (path, "rb")) != NULL) {
					k = 0;
					if ((output = fopen (out, "ab")) != NULL) {			
						if((i == 2)||(i == 3)) {
							while (!feof(fp)) {
								if((fread ((void *)(cand), sizeof (float), data, fp))==data){
									for(j = 0; j < data; j++) fwrite((void *)(cand), sizeof (float), data, output);
								}
							}
						}
	
						else {
							while (!feof(fp)) {
								if((fread ((void *)(cand), sizeof (float), data, fp))==data){
									if(k < (base_size*nod)){ //for nod = 2: k < 516986
										fwrite((void *)(cand), sizeof (float), data, output);
									}
									if((k == base_size*nod)&&(flag == 0)) ssb1[0] = cand[0]; //for nod = 2: k == 516986
									if((k == (base_size*nod)+1)&&(flag == 0)) ssb2[0] = cand[0]; //for nod = 2: k == 516987

									k++;
								}
							}
							fwrite((void *)(ssb1), sizeof (float), data, output);
							fwrite((void *)(ssb2), sizeof (float), data, output);

						}
					}
					else {		
						printf("Problem with %s file!\n", out);
					}
				}
				else {		
					printf("Problem with %s file!\n", path);
				}
				fclose(fp);
				fclose(output);

			}
		}
	}
	else {
		
		perror (listpath);
		return 1;
	}
//	puts("END OF GLUE FUNCTION");
	fclose(list);
	return 0;


}
