#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <getopt.h>
#include <gsl/gsl_linalg.h>
#include <time.h>
#include <dirent.h>

#include "auxi.h"
#include "settings.h"
#include "struct.h"
#include "init.h"

#include "arb.h"
#include "arb_mat.h" 

#ifndef CODEVER
#define CODEVER unknown
#endif

#ifndef DTAPREFIX
#define DTAPREFIX .
#endif


int fisher (
  Search_settings *sett,
  Command_line_opts *opts,
  Aux_arrays *aux);  

void cleanup_fisher(
	Search_settings *sett,
	Command_line_opts *opts,
	Aux_arrays *aux,
  double *F);


void handle_opts_fisher (
  Search_settings *sett,
  Command_line_opts *opts,
  int argc,
  char* argv[]); 


int main (int argc, char* argv[]) {

  Command_line_opts opts;
  Search_settings sett;
  Search_range s_range; 
  Aux_arrays aux_arr;
  double *F; 			  // F-statistic array
  int i; 

#define QUOTE(x) #x
#define STR(macro) QUOTE(macro)
#define CVSTR STR(CODEVER)

  printf("Code version : " CVSTR "\n");

  // Command line options 
  handle_opts_fisher(&sett, &opts, argc, argv);  
	
  // Search settings
  search_settings(&sett); 

  // Detector network settings
  detectors_settings(&sett, &opts); 

  // Array initialization and reading the ephemerids 
  init_arrays(&sett, &opts, &aux_arr, &F);

  // Amplitude modulation functions for each detector  
  for(i=0; i<sett.nifo; i++)   
    rogcvir(&ifo[i]); 


  // Main job - Fisher matrix calculations   
  for(i=0; i<sett.nifo; i++) 
    fisher(&sett, &opts, &aux_arr);


  // Cleanup & memory free 
  cleanup_fisher(&sett, &opts, &aux_arr, F);

  return 0; 
	
}

int fisher(Search_settings *sett,
           Command_line_opts *opts, 
           Aux_arrays *aux) { 


  int i, j, k, l, m, n, dim=8, dim2=4; 

  // Signal parameters: f, fdot, delta, alpha, a1, a2, a3, a4
  // (see Phys. Rev. D 82, 022005 2010, Eqs. 2.13a-d) 
  double sgnlo[dim]; 

  FILE *data; 

  // Reading signal parameters 
  if ((data=fopen (opts->addsig, "r")) != NULL) {

    for(i=0; i<dim; i++)
      fscanf(data, "%le",i+sgnlo); 
    
    fclose (data);
                 
  } else {
    perror (opts->addsig);
  }

  
  double ma[dim][dim], mFl[dim][dim]; 
  double a[dim2], F[dim2][dim2][dim2], S[dim2][dim2][dim2][dim2];
  double omega0, omega1, domega, sumhsq;  
  double sindelt, cosdelt, sinalt, cosalt; 

  omega0 = sgnlo[0]; // /sett->dt; 
  omega1 = sgnlo[1]; 

  domega = sett->oms; // 2*M_PI*sett->fpo*sett->dt; 

  sindelt = sin(sgnlo[2]); 
  cosdelt = cos(sgnlo[2]); 
  sinalt  = sin(sgnlo[3]); 
  cosalt  = cos(sgnlo[3]); 
  
  // a1 = a[0], a2 = a[1], a3 = a[2], a4 = a[3] 
  // in the signal amplitude model 
  // h = (a1*a + a2*b)*cos(psi) + (a3*a + a4*b)*sin(psi)
  a[0] = sgnlo[4]; 
  a[1] = sgnlo[5]; 
  a[2] = sgnlo[6]; 
  a[3] = sgnlo[7]; 

  // Loop for each detector 
  for(n=0; n<sett->nifo; ++n) { 

    /* Amplitude modulation functions aa and bb 
     * for each detector (in signal sub-struct 
     * of _detector, ifo[n].sig.aa, ifo[n].sig.bb) 
     */

    modvir(sinalt, cosalt, sindelt, cosdelt, 
        sett->N, &ifo[n], aux);

    for(k=0; k<dim; k++) 
      for(j=0; j<dim; j++) { 
        ma[k][j] = 0; 
        mFl[k][j] = 0; 
      }

    for(l=0; l<dim2; l++)
      for(k=0; k<dim2; k++) 
        for(j=0; j<dim2; j++) 
          F[l][k][j] = 0; 

    for(m=0; m<dim2; m++)
      for(l=0; l<dim2; l++)
        for(k=0; k<dim2; k++) 
          for(j=0; j<dim2; j++) 
            S[m][l][k][j] = 0; 

    double sumhsq = 0; 

    for(i=0; i<sett->N; ++i) {

      double psi, xet, yet, zet, aa, bb, amph, dpdA[dim2], h[dim2], dh[dim2]; 

      // Detector ephemerids at time step i
      xet = ifo[n].sig.DetSSB[i*3]; 
      yet = ifo[n].sig.DetSSB[i*3+1]; 
      zet = ifo[n].sig.DetSSB[i*3+2]; 

      // phase at time step i 
      psi = omega0*i + omega1*aux->t2[i] 
          + (cosalt*cosdelt*xet + sinalt*cosalt*yet + sindelt*zet)*(omega0 + 2*omega1*i + domega); 

      /* Phase derivatives w.r.t. intrinsic parameters A = (omega0, omega1, delta, alpha) 
       * freq = f = omega0, spindown = s = omega1, delta = d, alpha = a
       * 
       * dpdA[0] = dpdf, dpdA[1] = dpds, dpdA[2] = dpdd, dpdA[3] = dpda
       */     

      dpdA[0] = i + xet*cosalt*cosdelt + yet*cosdelt*sinalt + zet*sindelt;  
      dpdA[1] = aux->t2[i] + 2*i*(cosalt*cosdelt*xet + sinalt*cosalt*yet + sindelt*zet);   
      dpdA[2] = (domega + omega0 + 2*omega1*i)*(zet*cosdelt - (xet*cosalt + yet*sinalt)*sindelt); 
      dpdA[3] = (domega + omega0 + 2*omega1*i)*cosdelt*(yet*cosalt - xet*sinalt); 

      /* amplitude h = (a1*a + a2*b)*cos(psi) + (a3*a + a4*b)*sin(psi)
       * aa = a, bb = b 
       * 
       * h1 = h[0] = a*cos(psi)
       * h2 = h[1] = b*cos(psi) 
       * h3 = h[2] = a*sin(psi) 
       * h4 = h[3] = b*sin(psi) 
       */ 

      // amplitude modulation function at time step i 
      aa = ifo[n].sig.aa[i]; 
      bb = ifo[n].sig.bb[i]; 

      h[0] = aa*cos(psi); 
      h[1] = bb*cos(psi); 
      h[2] = aa*sin(psi); 
      h[3] = bb*sin(psi); 
    
      amph = a[0]*h[0] + a[1]*h[1] + a[2]*h[2] + a[3]*h[3];     

      // sum of amplitude h squares 
      sumhsq += amph*amph; 

      dh[0] = -aa*sin(psi); 
      dh[1] = -bb*sin(psi); 
      dh[2] =  aa*cos(psi); 
      dh[3] =  bb*cos(psi); 
  
      // F(A)
      for(l=0; l<dim2; l++) 
        for(k=0; k<dim2; k++) 
          for(j=0; j<dim2; j++)  
            F[l][k][j] += h[k]*dh[j]*dpdA[l]; 

      // S(A, B) 
      for(m=0; m<dim2; m++)
        for(l=0; l<dim2; l++)
          for(k=0; k<dim2; k++) 
            for(j=0; j<dim2; j++) 
              S[m][l][k][j] += dh[k]*dpdA[m]*dh[j]*dpdA[l]; 
 
      // Fisher matrix elements 
      ma[0][0] = (aa*aa*a[0]*a[0]*dpdA[0]*dpdA[0])/2. + (aa*aa*a[2]*a[2]*dpdA[0]*dpdA[0])/2. + aa*a[0]*a[1]*bb*dpdA[0]*dpdA[0] + aa*a[2]*a[3]*bb*dpdA[0]*dpdA[0] + (a[1]*a[1]*bb*bb*dpdA[0]*dpdA[0])/2. + (a[3]*a[3]*bb*bb*dpdA[0]*dpdA[0])/2.; 
 
      ma[0][1] = (aa*aa*a[0]*a[0]*dpdA[1]*dpdA[0])/2 + (aa*aa*a[2]*a[2]*dpdA[1]*dpdA[0])/2. + aa*a[0]*a[1]*bb*dpdA[1]*dpdA[0] + aa*a[2]*a[3]*bb*dpdA[1]*dpdA[0] + (a[1]*a[1]*bb*bb*dpdA[1]*dpdA[0])/2. + (a[3]*a[3]*bb*bb*dpdA[1]*dpdA[0])/2.;
 
      ma[0][2] = (aa*aa*a[0]*a[0]*dpdA[0]*dpdA[2])/2 + (aa*aa*a[2]*a[2]*dpdA[0]*dpdA[2])/2. + aa*a[0]*a[1]*bb*dpdA[0]*dpdA[2] + aa*a[2]*a[3]*bb*dpdA[0]*dpdA[2] + (a[1]*a[1]*bb*bb*dpdA[0]*dpdA[2])/2. + (a[3]*a[3]*bb*bb*dpdA[0]*dpdA[2])/2.; 

      ma[0][3] = (aa*aa*a[0]*a[0]*dpdA[0]*dpdA[3])/2. + (aa*aa*a[2]*a[2]*dpdA[0]*dpdA[3])/2. + aa*a[0]*a[1]*bb*dpdA[0]*dpdA[3] + aa*a[2]*a[3]*bb*dpdA[0]*dpdA[3] + (a[1]*a[1]*bb*bb*dpdA[0]*dpdA[3])/2. + (a[3]*a[3]*bb*bb*dpdA[0]*dpdA[3])/2.;  

      ma[0][4] = (aa*aa*a[2]*dpdA[0])/2. + (aa*a[3]*bb*dpdA[0])/2.;  

      ma[0][5] = (aa*a[2]*bb*dpdA[0])/2. + (a[3]*bb*bb*dpdA[0])/2.;
 
      ma[0][6] = -(aa*aa*a[0]*dpdA[0])/2. - (aa*a[1]*bb*dpdA[0])/2.;
 
      ma[0][7] = -(aa*a[0]*bb*dpdA[0])/2. - (a[1]*bb*bb*dpdA[0])/2.; 

      ma[1][1] = (aa*aa*a[0]*a[0]*dpdA[1]*dpdA[1])/2. + (aa*aa*a[2]*a[2]*dpdA[1]*dpdA[1])/2. + aa*a[0]*a[1]*bb*dpdA[1]*dpdA[1] + aa*a[2]*a[3]*bb*dpdA[1]*dpdA[1] + (a[1]*a[1]*bb*bb*dpdA[1]*dpdA[1])/2. + (a[3]*a[3]*bb*bb*dpdA[1]*dpdA[1])/2.; 

      ma[1][2] = (aa*aa*a[0]*a[0]*dpdA[1]*dpdA[2])/2. + (aa*aa*a[2]*a[2]*dpdA[1]*dpdA[2])/2. + aa*a[0]*a[1]*bb*dpdA[1]*dpdA[2] + aa*a[2]*a[3]*bb*dpdA[1]*dpdA[2] + (a[1]*a[1]*bb*bb*dpdA[1]*dpdA[2])/2. + (a[3]*a[3]*bb*bb*dpdA[1]*dpdA[2])/2.;
 
      ma[1][3] = (aa*aa*a[0]*a[0]*dpdA[1]*dpdA[3])/2. + (aa*aa*a[2]*a[2]*dpdA[1]*dpdA[3])/2. + aa*a[0]*a[1]*bb*dpdA[1]*dpdA[3] + aa*a[2]*a[3]*bb*dpdA[1]*dpdA[3] + (a[1]*a[1]*bb*bb*dpdA[1]*dpdA[3])/2. + (a[3]*a[3]*bb*bb*dpdA[1]*dpdA[3])/2.;

      ma[1][4] = (aa*aa*a[2]*dpdA[1])/2. + (aa*a[3]*bb*dpdA[1])/2.;
 
      ma[1][5] = (aa*a[2]*bb*dpdA[1])/2. + (a[3]*bb*bb*dpdA[1])/2.; 
 
      ma[1][6] = -(aa*aa*a[0]*dpdA[1])/2. - (aa*a[1]*bb*dpdA[1])/2.;
 
      ma[1][7] = -(aa*a[0]*bb*dpdA[1])/2. - (a[1]*bb*bb*dpdA[1])/2.;
 
      ma[2][2] = (aa*aa*a[0]*a[0]*dpdA[2]*dpdA[2])/2. + (aa*aa*a[2]*a[2]*dpdA[2]*dpdA[2])/2. + aa*a[0]*a[1]*bb*dpdA[2]*dpdA[2] + aa*a[2]*a[3]*bb*dpdA[2]*dpdA[2] + (a[1]*a[1]*bb*bb*dpdA[2]*dpdA[2])/2. + (a[3]*a[3]*bb*bb*dpdA[2]*dpdA[2])/2.;
 
      ma[2][3] = (aa*aa*a[0]*a[0]*dpdA[2]*dpdA[3])/2. + (aa*aa*a[2]*a[2]*dpdA[2]*dpdA[3])/2. + aa*a[0]*a[1]*bb*dpdA[2]*dpdA[3] + aa*a[2]*a[3]*bb*dpdA[2]*dpdA[3] + (a[1]*a[1]*bb*bb*dpdA[2]*dpdA[3])/2. + (a[3]*a[3]*bb*bb*dpdA[2]*dpdA[3])/2.;
 
      ma[2][4] = (aa*aa*a[2]*dpdA[2])/2. + (aa*a[3]*bb*dpdA[2])/2.;
 
      ma[2][5] = (aa*a[2]*bb*dpdA[2])/2. + (a[3]*bb*bb*dpdA[2])/2.;
 
      ma[2][6] = -(aa*aa*a[0]*dpdA[2])/2. - (aa*a[1]*bb*dpdA[2])/2.;
 
      ma[2][7] = -(aa*a[0]*bb*dpdA[2])/2. - (a[1]*bb*bb*dpdA[2])/2.;
 
      ma[3][3] = (aa*aa*a[0]*a[0]*dpdA[3]*dpdA[3])/2. + (aa*aa*a[2]*a[2]*dpdA[3]*dpdA[3])/2. + aa*a[0]*a[1]*bb*dpdA[3]*dpdA[3] + aa*a[2]*a[3]*bb*dpdA[3]*dpdA[3] + (a[1]*a[1]*bb*bb*dpdA[3]*dpdA[3])/2. + (a[3]*a[3]*bb*bb*dpdA[3]*dpdA[3])/2.;
 
      ma[3][4] = (aa*aa*a[2]*dpdA[3])/2. + (aa*a[3]*bb*dpdA[3])/2.;
 
      ma[3][5] = (aa*a[2]*bb*dpdA[3])/2. + (a[3]*bb*bb*dpdA[3])/2.;
 
      ma[3][6] = -(aa*aa*a[0]*dpdA[3])/2. - (aa*a[1]*bb*dpdA[3])/2.;
 
      ma[3][7] = -(aa*a[0]*bb*dpdA[3])/2. - (a[1]*bb*bb*dpdA[3])/2.;
 
      ma[4][4] = aa*aa/2.;
 
      ma[4][5] = (aa*bb)/2.;
 
      ma[4][6] = 0;
 
      ma[4][7] = 0;
 
      ma[5][5] = bb*bb/2.;
 
      ma[5][6] = 0;
 
      ma[5][7] = 0;
 
      ma[6][6] = aa*aa/2.;
 
      ma[6][7] = (aa*bb)/2.;
 
      ma[7][7] = bb*bb/2.;


      // mFl 
      for(k=0; k<dim; k++) 
        for(j=0; j<dim; j++) 
          mFl[k][j] += ma[k][j];

    } 

    // Symmetrize mFl 
    for(k=1; k<dim; k++) 
      for(j=0; j<k; j++)
        mFl[k][j] = mFl[j][k]; 

    printf("A0=omega0, A1=omega1, A2=delta, A3=alpha\n"); 

    for(l=0; l<dim2; l++) { 
      printf("\nF(A%d):\n", l); 
      for(k=0; k<dim2; k++) { 
        for(j=0; j<dim2; j++)  
          printf("%.16e ", F[l][k][j]);
        printf("\n");
      }
    } 

    for(m=0; m<dim2; m++)
      for(l=0; l<=m; l++) { 
        printf("\nS(A%d, A%d):\n", m, l); 
        for(k=0; k<dim2; k++) {  
          for(j=0; j<dim2; j++) 
            printf("%.16e ", S[m][l][k][j]);
          printf("\n");
        } 
      }

    //#mb printf("\nsumhsq, sqrt(sumhsq), N: %f %f %d\n", sumhsq, sqrt(sumhsq), sett->N); 

    printf("\nThe Fisher matrix:\n"); 

    int r; 
    arb_mat_t A, T; 
    arb_t t; 

    arb_mat_init(A, dim, dim); 
    arb_mat_init(T, dim, dim); 

    arb_init(t); 

    for(k=0; k<dim; k++) { 
//      printf("[");
      for(j=0; j<dim; j++) { 
        printf("%.16e ", mFl[k][j]);
        arb_set_d(arb_mat_entry(A, k, j), mFl[k][j]);
      } 
      printf("\n");
    }

    printf("Inverting the Fisher matrix...\n"); 

    r = arb_mat_spd_inv(T, A, 128);

    if(r) {
      
      printf("Diagonal elements of the covariance matrix:\n");  
      for(k=0; k<dim; k++) { 
        arb_set(t, arb_mat_entry(T, k, k)); 
        printf("%le ", arf_get_d(arb_midref(t), ARF_RND_DOWN));
      }
    } else { 
    
      printf("Failed to invert the Fisher matrix.\n"); 
    } 

    printf("\n"); 

    arb_mat_clear(A);
    arb_mat_clear(T); 


  } 

  return 0; 

} 


  /* Cleanup & memory free (calculation 
   * of the Fisher matrix) 
	 */

void cleanup_fisher(
	Search_settings *sett,
	Command_line_opts *opts,
	Aux_arrays *aux, 
  double *F) {

  int i; 

  for(i=0; i<sett->nifo; i++) {
    free(ifo[i].sig.xDat);
    free(ifo[i].sig.xDatma);
    free(ifo[i].sig.xDatmb); 
    free(ifo[i].sig.DetSSB);
    free(ifo[i].sig.aa);
    free(ifo[i].sig.bb);
    free(ifo[i].sig.shftf);
    free(ifo[i].sig.shft);
  } 
	
  free(aux->sinmodf);
  free(aux->cosmodf);
  free(aux->t2);
  free(F);  


} // end of cleanup & memory free 



/*  Command line options handling: fisher 
 */ 

void handle_opts_fisher( Search_settings *sett, 
		  Command_line_opts *opts,
		  int argc, 
		  char* argv[]) {

  strcpy (opts->prefix, TOSTR(PREFIX));
  strcpy (opts->dtaprefix, TOSTR(DTAPREFIX));

  opts->label[0]    = '\0';
  opts->range[0]    = '\0';
  opts->getrange[0] = '\0';
  opts->usedet[0]   = '\0';
  opts->addsig[0]   = '\0';
	
  // Initial value of starting frequency set to a negative quantity. 
  // If this is not changed by the command line value, fpo is calculated 
  // from the band number b (fpo = fpo = fstart + 0.96875*b/(2dt))
  sett->fpo = -1;

  // Default initial value of the data sampling time 
  sett->dt = 0.5; 

  // Initial value of the number of days is set to 0
  sett->nod = 0; 

  opts->help_flag=0;

  static int help_flag=0;

  // Reading arguments 

  while (1) {
    static struct option long_options[] = {
      {"help", no_argument, &help_flag, 1},
      // frame number
      {"ident", required_argument, 0, 'i'},
      // frequency band number
      {"band", required_argument, 0, 'b'},
      // input data directory
      {"data", required_argument, 0, 'd'},
      // fpo value
      {"fpo", required_argument, 0, 'p'},
      // add signal parameters
      {"addsig", required_argument, 0, 'x'},
      // number of days in the time-domain segment 
      {"nod", required_argument, 0, 'y'},
      // which detectors to use
      {"usedet", required_argument, 0, 'u'}, 
      // data sampling time 
      {"dt", required_argument, 0, 's'},
      {0, 0, 0, 0}
    };

    if (help_flag) {

      printf("polgraw-allsky periodic GWs: Fisher matrix calculation\n");
      printf("Usage: ./fisher -[switch1] <value1> -[switch2] <value2> ...\n") ;
      printf("Switches are:\n\n");
      printf("-data         Data directory (default is .)\n");
      printf("-ident        Frame number\n");
      printf("-band         Band number\n");
      printf("-fpo          Reference band frequency fpo value\n");
      printf("-dt           Data sampling time dt (default value: 0.5)\n");
      printf("-usedet       Use only detectors from string (default is use all available)\n");
      printf("-addsig       Add signal with parameters from <file>\n");
      printf("-nod          Number of days\n\n");


      printf("Also:\n\n");
      printf("--help            This help\n");

      exit(EXIT_SUCCESS);
    }

    int option_index = 0;
    int c = getopt_long_only(argc, argv, "i:b:dp:x:y:s:u:", 
			     long_options, &option_index);
    if (c == -1)
      break;

    switch (c) {
    case 'i':
      opts->ident = atoi (optarg);
      break;
    case 'b':
      opts->band = atoi(optarg);
      break;
    case 'd':
      strcpy(opts->dtaprefix, optarg);
      break;
    case 'p':
      sett->fpo = atof(optarg);
      break;
    case 'x':
      strcpy(opts->addsig, optarg);
      break;
    case 'y':
      sett->nod = atoi(optarg);
      break;
    case 's':
      sett->dt = atof(optarg);
      break;
    case 'u':
      strcpy(opts->usedet, optarg);
      break;
    case '?':
      break;
    default:
      break ;
    } /* switch c */
  } /* while 1 */

  // Check if sett->nod was set up, if not, exit
  if(!(sett->nod)) { 
    printf("Number of days not set... Exiting\n"); 
    exit(EXIT_FAILURE); 
  } 

  printf("Number of days is %d\n", sett->nod); 

  printf("Input data directory is %s\n", opts->dtaprefix);
  printf("Frame and band numbers are %d and %d\n", opts->ident, opts->band);

  // Starting band frequency:
  // fpo_val is optionally read from the command line
  // Its initial value is set to -1
  if(!(sett->fpo >= 0))

    // The usual definition (multiplying the offset by B=1/(2dt))
    // !!! in RDC_O1 the fstart equals 10, not 100 like in VSR1 !!! 
    // 
    sett->fpo = 10. + 0.96875*opts->band*(0.5/sett->dt);

  printf("The reference frequency fpo is %f\n", sett->fpo);
  printf("The data sampling time dt is %f\n", sett->dt); 

  if (strlen(opts->addsig))
    printf ("Adding signal from '%s'\n", opts->addsig);

} // end of command line options handling: fisher  

