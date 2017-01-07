/* 
Ver. 1.0
Function joins together DetSSB files
(last 2 values from the first file are moved
at the end of final file) and xdatc* with given 
label name. Programm uses list.txt (from data dir), 
where should exist list of interesting times. Final 
files are in output_dir/followup_total_data. 
Ver. 2.0
The same, but outside followup. Programm takes 
arguments from command line. Problem with size solved.
MS
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
//#include "glue.h"

//int glue(char prfx[512], char dataprfx[512], char band[512], int nod){ //Ver. 1.0
int main(int argc, char **argv){ //Ver 2.0, Arguments: output data label 

	char prfx[512];
	char dataprfx[512];
	char band[512];
	char newident[512];
	strcpy(prfx, argv[1]);
	strcpy(dataprfx, argv[2]);
	strcpy(band, argv[3]);
	strcpy(newident, argv[4]);

	FILE *list; 
	FILE *fp;
	FILE *output;
	int i, j, k, l, flag;
	int data = 1;
	int line;
	int base_size = 6203824; //size of DetSSB.bin (in bytes); for 2days: 2067952, for 6days: 6203824
	char xdat[200];
	char listpath[200];
	char path[200];
	char part[100];
	char out[200];
	char output_tot[512];
	char output_L[512];
	char output_H[512];
	char output_0[512];
	double *cand;
	double *ssb1_l, *ssb2_l, *ssb1_h, *ssb2_h;
	cand = (double *) calloc (data, sizeof(double));
	ssb1_l = (double *) calloc (data, sizeof(double));
	ssb2_l = (double *) calloc (data, sizeof(double));
	ssb1_h = (double *) calloc (data, sizeof(double));
	ssb2_h = (double *) calloc (data, sizeof(double));
	sprintf(listpath, "%s/list.txt", dataprfx); //list.txt should be in data dir

	
	if ((list = fopen (listpath, "r")) != NULL) {
		flag = 0;
		l = 0;
		while (fscanf(list, "%d\n", &line) == 1){
			if (l != 0) flag = 1;
			l++;
			sprintf(xdat, "xdatg_%03d_%s.bin", line, band);
			mkdir(prfx, 0777);
			sprintf(output_tot, "%s/glued_data", prfx);
			mkdir(output_tot, 0777);
			sprintf(output_0, "%s/glued_data/%s", prfx, newident);
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
					sprintf(out, "%s/L1/xdatg_%s_%s.bin", output_0, newident, band);
				}
				if(i == 3){
					sprintf(part, "/H1/%s", xdat);
					sprintf(out, "%s/H1/xdatg_%s_%s.bin", output_0, newident, band);
				}
				sprintf(path, "%s/%03d%s", dataprfx, line, part);
				k = 1;
				if ((fp = fopen (path, "rb")) != NULL) {
					if ((output = fopen (out, "ab")) != NULL) {			
						if((i == 2)||(i == 3)) {
							while (!feof(fp)) {
								if((fread ((void *)(cand), sizeof (double), data, fp))==data){
									for(j = 0; j < data; j++) fwrite((void *)(cand), sizeof (double), data, output);
								}
							}
						}
						else {
							while (!feof(fp)) {
								if((fread ((void *)(cand), sizeof (double), data, fp))==data){
									if(k < ((base_size-8)/8)){ 
										fwrite((void *)(cand), sizeof (double), data, output);
									}
									if((k == ((base_size-8)/8))&&(flag == 0)&&(i == 0)){ 
										ssb1_l[0] = cand[0]; 									
									}
									if((k == (base_size/8))&&(flag == 0)&&(i == 0)){ 
										ssb2_l[0] = cand[0]; 
									}	
									if((k == ((base_size-8)/8))&&(flag == 0)&&(i == 1)){ 
										ssb1_h[0] = cand[0]; 
									}
									if((k == (base_size/8))&&(flag == 0)&&(i == 1)){ 
										ssb2_h[0] = cand[0]; 
									}
									k++;
								}
							}
							
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
	sprintf(out, "%s/DetSSB.bin",output_L);
	if ((output = fopen (out, "ab")) != NULL) {
		fwrite((void *)(ssb1_l), sizeof (double), data, output);
		fwrite((void *)(ssb2_l), sizeof (double), data, output);
	}
	else {		
		printf("Problem with %s file - at the end!\n", out);
	}
	sprintf(out, "%s/DetSSB.bin",output_H);
	if ((output = fopen (out, "ab")) != NULL) {
		fwrite((void *)(ssb1_h), sizeof (double), data, output);
		fwrite((void *)(ssb2_h), sizeof (double), data, output);
	}
	else {		
		printf("Problem with %s file - at the end!\n", out);
	}
//	puts("END OF GLUE FUNCTION");
	fclose(list);
	return 0;


}
