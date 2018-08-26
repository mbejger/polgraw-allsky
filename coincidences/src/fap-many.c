#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_combination.h>
#include <getopt.h>

#include "auxi.h"
#include "settings.h"
#include "struct.h"

#define min(x, y) (((x) < (y)) ? (x) : (y))

static int help_flag=0;

int nchoosek(int, int);
int *FalseAlarm(int, int, int, double, int*, double*);
int *FalseAlarmCFast(int, int, int, double, int*, double*);  

int main (int argc, char *argv[]) {

  Command_line_opts opts;
  Search_settings sett;

  size_t i;
  short int noc;  
  int status, band=0, cellsize=4, To; 
  char filename[512], datafile[512], griddir[512];
  double gamrn[16], vetofrac=0, threshold=0.1; 
  double f_min, f_max, fdotmin, fdotmax, detgamrn, vol4, vl, vc, Nc;
  FILE *data;

  // Default initial value of the data sampling time 
  sett.dt = 2 ; 

  // Initial value of the number of days is set to 0
  sett.nod = 0; 

  // Initial value of the band number is set to 0 
  opts.band = 0; 

  opts.help_flag=0;
  static int help_flag=0;  

  // Reading arguments 

  while (1) {
    static struct option long_options[] = {
      {"help", no_argument, &help_flag, 1},
      // frequency band number
      {"band", required_argument, 0, 'b'},
      // cell size 
      {"cellsize", required_argument, 0, 'c'},
      // input data directory
      {"data", required_argument, 0, 'd'},
      // grid matrix data directory
      {"grid", required_argument, 0, 'g'},
      // FAP threshold
      {"threshold", required_argument, 0, 't'},
      // data sampling time 
      {"dt", required_argument, 0, 's'},
      // veto fraction 
      {"vetofrac", required_argument, 0, 'v'},
      // number of days in the time-domain segment 
      {"nod", required_argument, 0, 'y'},
      {0, 0, 0, 0}
    };

    if (help_flag) {

      printf("polgraw-allsky periodic GWs FAP of coincidence calculator\n");
      printf("Usage: ./fap -[switch1] <value1> -[switch2] <value2> ...\n") ;
      printf("Switches are:\n\n");
      printf("-band         Band number\n");
      printf("-cellsize     Cell size (default value: 4)\n");
      printf("-data         Coincidence summary file\n");
      printf("-grid         Grid matrix directory (default value: .)\n");
      printf("-dt           Data sampling time dt (default value: 2)\n");
      printf("-threshold    FAP threshold (default value: 0.1)\n");
      printf("-nod          Number of days\n");
      printf("-vetofrac     Vetoed fraction of the band (default value: 0)\n\n");

      printf("Also:\n\n");
      printf("--help            This help\n");

      exit(EXIT_SUCCESS);
    }

    int option_index = 0;
    int c = getopt_long_only(argc, argv, "b:c:d:g:s:t:y:v", 
			     long_options, &option_index);
    if (c == -1)
      break;

    switch (c) {
    case 'b': // band 
      opts.band = atoi(optarg);
      break;
    case 'c': // cellsize 
      cellsize = atoi(optarg);
      break;
    case 'd': // data file  
      strcpy(datafile, optarg);
      break;
    case 'g': // grid dir 
      strcpy(griddir, optarg);
      break;
    case 's': // sampling time 
      sett.dt = atof(optarg);
      break;
    case 't': // FAP threshold 
      threshold = atof(optarg);
      break;
    case 'y': // number of days
      sett.nod = atoi(optarg);
      break;
    case 'v': // veto fraction 
      vetofrac = atof(optarg);
      break;
    case '?':
      break;
    default:
      break ;
    } /* switch c */
  } /* while 1 */


  // Check if sett.nod was set up, if not, exit
  if(!(sett.nod)||!(opts.band)) { 
    printf("Number of days or band number not set... Exiting\n"); 
    exit(EXIT_FAILURE); 
  } 

  printf("Number of days in time segments: %d\n", sett.nod); 

  printf("Input data: %s\n", datafile);
  printf("Grid matrix data directory: %s\n", griddir);

  printf("Band number: %04d (veto fraction: %f)\n", opts.band, vetofrac);

  // Starting band frequency:
  sett.fpo = 10. + 0.96875*opts.band*(0.5/sett.dt);

  printf("The reference frequency fpo: %f\n", sett.fpo);
  printf("The data sampling time dt: %f\n", sett.dt); 
  printf("FAP threshold: %f\n", threshold); 

  printf("Cell size: %d\n", cellsize); 

	
  // Search settings
  //----------------

  search_settings(&sett); 

  // Observation time  
  To = sett.N*sett.dt;

  // Read grid: gamrn the normalized Fisher matrix  
  sprintf (filename, "%s/grid.bin", griddir);

  if ((data=fopen (filename, "rb")) != NULL) {

    // skipping fftpad (1 int) and the M matrix (16 doubles) 
    fseek(data, sizeof(int) + 16*sizeof(double), SEEK_SET);
    // gamrn: normalized Fisher matrix  
    status = fread ((void *)gamrn, sizeof (double), 16, data);

    fclose (data);
  } else {
    perror (filename);
    exit(EXIT_FAILURE);
  }

/*
  printf("gamrn matrix from grid.bin:\n");

	printf("%e %e %e %e\n", gamrn[0], gamrn[1], gamrn[2], gamrn[3]);
	printf("%e %e %e %e\n", gamrn[4], gamrn[5], gamrn[6], gamrn[7]);
	printf("%e %e %e %e\n", gamrn[8], gamrn[9], gamrn[10], gamrn[11]);
	printf("%e %e %e %e\n", gamrn[12], gamrn[13], gamrn[14], gamrn[15]);   
*/ 

  // Nc: Number of cells 
  //-------------------- 

  // Determinant of the Fisher matrix 
  int signum;
  gsl_matrix_view m = gsl_matrix_view_array (gamrn, 4, 4);
  gsl_permutation *p = gsl_permutation_alloc (4);
  gsl_linalg_LU_decomp (&m.matrix, p, &signum);

  detgamrn = gsl_linalg_LU_det(&m.matrix, signum);

  gsl_permutation_free (p); 

  // Hypervolume vol4 of a 4-dimensional sphere 
  // (for an n-dimensional sphere = pi^(n/2)/gamma(n/2+1) )  
  vol4 = M_PI*M_PI/gsl_sf_gamma(3);  

  vc = vol4/sqrt(detgamrn); 

  // Values of fdotmin and fdotmax recovered from settings.c 
  fdotmin = sett.Smax/(2*M_PI*sett.dt*sett.dt); 
  fdotmax = sett.Smin/(2*M_PI*sett.dt*sett.dt); 

  f_max = sett.fpo + sett.B; f_min = sett.fpo; 

  // vl =  16/3*pi^5*To^5*(fdotmin+fdotmax)*(f_max^3-f_min^3) 
  vl = 16./3.*pow(M_PI*To, 5)*(fdotmin + fdotmax)
     *(pow(f_max, 3) - pow(f_min, 3));

  Nc = round(vl/vc); 
 
  // Taking into acount the vetoed fraction of the band 
  if(vetofrac) 
    Nc = round((1.0 - vetofrac)*Nc/pow(cellsize, 4));
  else 
    Nc = round(Nc/pow(cellsize, 4));

  // Read the coincidence data 
  //--------------------------  

  if ((data=fopen (datafile, "r")) != NULL) {

    // until the end-of-file 
    do { 

      short int i, shift, nof, band, hemi;
      double fpofile;  

      status = fscanf(data, "%hu_%hu %hu %lf %hu %hu",  
                 &band, &hemi, &shift, &fpofile, &nof, &noc);   

      int Nk[nof], Nku[nof], frn[nof], frnc[noc], Nkall=0;
      double sigpar[5]; 

      // Coincidence signal parameters (f, s, d, a, snr) 
      for(i=0; i<5; i++) {  
        status = fscanf(data, "%le", &sigpar[i]); 
      } 

      // Frames information: frame number, no. of candidates, no. of unique candidates 
      for(i=0; i<nof; i++) { 
        status = fscanf(data, "%d %d %d", &frn[i], &Nk[i], &Nku[i]); 
        Nkall += Nku[i]; 
      } 

      // Numbers of frames participating in the coincidence 
      for(i=0; i<noc; i++) { 
        status = fscanf(data, "%d", &frnc[i]); 
      } 

      double FAP, PFce[2*noc-2]; 
      FalseAlarmCFast(2, noc, nof, Nc, &Nku[0], &PFce[0]); 
      FAP = PFce[2*noc-3]; 

      // Final result: output to stderr cases when FAP threshold is reached  
      if( (FAP < threshold) && (status != EOF) ) { 
   
        fprintf(stderr, "%04d %le %le %le %d %d ", 
          band, f_min, f_max, FAP, noc, Nkall); 
        
        for(i=0; i<5; i++) fprintf(stderr, "%le ", sigpar[i]);
   
        fprintf(stderr, "%d\n", hemi); 

      }  

    } while (status != EOF);   

  } else {
    perror (filename);
    exit(EXIT_FAILURE);
  }

  fclose (data);


  return 0;
 
}

// Binomial coeficient n choose k 
//-------------------------------

int nchoosek(int n, int k) {
  return gsl_sf_fact(n)/(gsl_sf_fact(k)*gsl_sf_fact(n-k));
} 


int *FalseAlarm(int Cmax, int noc, int L, double Nc, int *Nk, double *r) { 

   int i, j, k, Nmax;
   double ee[L];

   gsl_combination *cp, *cq;

   for(i=0; i<L; i++)    //#mb length(Nk)  
      ee[i] = Nk[i]/Nc;

  Nmax = min(noc, L);  

  double C[Nmax], pf=0; 

  for(i=Cmax; i<=Nmax; i++) {  
    
    cp = gsl_combination_calloc (L, i);
    
    double P[nchoosek(L, i)];

    k=0; 
    do {
      
      P[k] = 1;   
      for(j=0; j<gsl_combination_k(cp); j++)
        P[k] *= ee[gsl_combination_get(cp, j)];

      k++;     

    } while (gsl_combination_next (cp) == GSL_SUCCESS);

    gsl_combination_free (cp);

    cq = gsl_combination_calloc (L, L-i);
    
    double Q[nchoosek(L, L-i)];

    k=0;   
    do {
      
      Q[k] = 1;   
      for(j=0; j<gsl_combination_k(cq); j++)
        Q[k] *= (1. - ee[gsl_combination_get(cq, j)]); 

      k++;     

    } while (gsl_combination_next (cq) == GSL_SUCCESS);

    gsl_combination_free (cq);

    C[i-1] = 0; 
    for(k=0; k<nchoosek(L, i); k++)
      C[i-1] += P[k]*Q[nchoosek(L, L-i) - (k+1)];

    // Probability that a cell cointains Cmax or more coincidences
    pf += C[i-1];

  } 

  // False alarm probability PF 
  // Probability that there is Cmax or more coincidences in one or more cells
  r[0] = 1. - pow(1. - pf, Nc);

  // Expected number of false alarms NF
  r[1] = Nc*pf;

  r[2] = pf; 
  
  k=1; 
  for(i=Cmax-1; i<Nmax; i++) {  
    r[2+k] = C[i]; 
    k++; 
  } 

  // The return order of r is: PF, NF, pf, C[] array

  return 0; 

}


int *FalseAlarmCFast(int Cmax, int noc, int L, double Nc, int *Nk, double *r) { 

  int i, k;

  // Arrays of noc + 2 doubles: PF, NF, pf, C[] array  
  double C0[noc+2], C1[noc+2], C2[noc+2], C3[noc+2], C4[noc+2]; 

  for(i=0; i<noc+3; i++) { 
    C0[i] = 0; C1[i] = 0; C2[i] = 0; C3[i] = 0; C4[i] = 0; 
  } 

  FalseAlarm(Cmax, noc, L, Nc, &Nk[0], &C0[0]);
  FalseAlarm(Cmax, noc, L, 2*Nc, &Nk[0], &C1[0]);
  FalseAlarm(Cmax, noc, L, 2*2*Nc, &Nk[0], &C2[0]);
  FalseAlarm(Cmax, noc, L, 2*2*2*Nc, &Nk[0], &C3[0]);
  FalseAlarm(Cmax, noc, L, 2*2*2*2*Nc, &Nk[0], &C4[0]);

  double pfe0[noc], pfe1[noc], pfe2[noc], pfe3[noc], pfe4[noc];

  for(k=0; k<noc; k++) { 

    pfe0[k] = 0; pfe1[k] = 0; pfe2[k] = 0; pfe3[k] = 0; pfe4[k] = 0;

    for(i=k; i<noc; i++) { 

      pfe0[k] += C0[3+i]; 
      pfe1[k] += C1[3+i];
      pfe2[k] += C2[3+i]; 
      pfe3[k] += C3[3+i];
      pfe4[k] += C4[3+i];

    }

  } 

  double PF0[noc-1], PF04cor[noc-1]; 

  // PF0 = 1 - (1 - pfe).^Nc
  for(i=0; i<noc-1; i++) {  

    PF0[i] = 1. - pow(1. - pfe0[i], Nc);

    PF04cor[i] = pow(2,4)*pfe0[i]  
               - ( nchoosek(4,1)*pfe1[i] 
                 + nchoosek(4,2)*pfe2[i] 
                 + nchoosek(4,3)*pfe3[i] 
                 + nchoosek(4,4)*pfe4[i] )
               - ( nchoosek(4,2)*pfe2[i] 
                 + nchoosek(4,3)*pfe3[i] 
                 + nchoosek(4,4)*pfe4[i] )
               - ( nchoosek(4,3)*pfe3[i] 
                 + nchoosek(4,4)*pfe4[i] )
               - pfe4[i]; 

    PF04cor[i] = 1. - pow(1. - PF04cor[i], Nc);
  
    // return order of r is (2*noc-2) doubles: [PF0, PF04cor] 
    r[i] = PF0[i]; 
    r[i+noc-1] = PF04cor[i]; 

  } 

  return 0; 

} 
