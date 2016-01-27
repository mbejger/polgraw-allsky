/* 
Ver. 1.0
Function join together DetSSB files
(last 2 values from the first file are moved
at the end of final file) and xdatc* with given 
band name. Programm uses list.txt, where
should be list of interesting times. Final 
files are in data_total directory. 
MS
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "glue.h"

void glue(char band[10]){
	FILE *list; 
	FILE *fp;
	FILE *output;
	int i, j, k, l, flag = 0;
	int data = 1;
	int line;
	char xdat[24];
	char path[200];
	char part[100];
	char out[200];
	float *cand;
	float *ssb1, *ssb2;
	cand = (float *) calloc (data, sizeof(float));
	ssb1 = (float *) calloc (data, sizeof(float));
	ssb2 = (float *) calloc (data, sizeof(float));
	if ((list = fopen ("list.txt", "r")) != NULL) {
		l = 0;
		while (fscanf(list, "%d\n", &line) == 1){
			if (l != 0) flag = 1;
			l++;
			sprintf(xdat, "xdatc_%03d_%s.bin", line, band);
			mkdir("data_total", 0777);
			mkdir("data_total/H1/", 0777);
			mkdir("data_total/L1/", 0777);
			for(i = 0; i < 4; i++){
				if(i == 0){
					sprintf(part, "/L1/DetSSB.bin");
					sprintf(out, "data_total%s", part);
				}
				if(i == 1){
					sprintf(part, "/H1/DetSSB.bin");
					sprintf(out, "data_total%s", part);
				}
				if(i == 2){
					sprintf(part, "/L1/%s", xdat);
					sprintf(out, "data_total/L1/xdatc_%s.bin", band);
				}
				if(i == 3){
					sprintf(part, "/H1/%s", xdat);
					sprintf(out, "data_total/H1/xdatc_%s.bin", band);
				}
				sprintf(path, "%03d%s", line, part);
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
									if(k < 516986){
										fwrite((void *)(cand), sizeof (float), data, output);
									}
									if((k == 516986)&&(flag == 0)) ssb1[0] = cand[0];
									if((k == 516987)&&(flag == 0)) ssb2[0] = cand[0];
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
	puts("END OF GLUE FUNCTION");
	fclose(list);



}
