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

int FalseAlarmProb(int, int, double, int*, double *);


int main (int argc, char *argv[]) {

    Command_line_opts opts;
    Search_settings sett;

    size_t i;
    short int noc; 
    int status, band=0, cellf=4, cells=4, celld=4, cella=4, To; 
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
    opts.overlap = -1.;
    
    opts.help_flag=0;
    static int help_flag=0;  
     
    // Default value of noc 
    noc = 2;

    // Reading arguments 
    while (1) {
	static struct option long_options[] = {
	    {"help", no_argument, &help_flag, 1},
	    // frequency band number
	    {"band", required_argument, 0, 'b'},
	    // overlap
	    {"overlap", required_argument, 0, 'o'},
	    // cell size f
	    {"cellf", required_argument, 0, 'i'},
	    // cell size s
	    {"cells", required_argument, 0, 'j'},
	    // cell size d
	    {"celld", required_argument, 0, 'k'},
	    // cell size a
	    {"cella", required_argument, 0, 'l'},
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
	    // number of coincidences  
	    {"noc", required_argument, 0, 'n'},
	    {0, 0, 0, 0}
	};
	
	if (help_flag) {

	    printf("polgraw-allsky periodic GWs FAP of coincidence calculator\n");
	    printf("Usage: ./fap -[switch1] <value1> -[switch2] <value2> ...\n") ;
	    printf("Switches are:\n\n");
	    printf("-band         Band number [-b]\n");
	    printf("-overlap      Band overlap [-o]\n");
	    printf("-cellf        Cell size f (default value: 4) [-i]\n");
	    printf("-cells        Cell size s (default value: 4) [-j]\n");
	    printf("-celld        Cell size d (default value: 4) [-k]\n");
	    printf("-cella        Cell size a (default value: 4) [-l]\n");
	    printf("-data         Coincidence summary file [-d]\n");
	    printf("-grid         Grid matrix directory (default value: .) [-g]\n");
	    printf("-dt           Data sampling time dt (default value: 2) [-s]\n");
	    printf("-threshold    FAP threshold (default value: 0.1) [-t]\n");
	    printf("-nod          Number of days [-y]\n");
	    printf("-vetofrac     Vetoed fraction of the band (default value: 0) [-v]\n");
	    printf("-noc          Number of coincidences [-n]\n\n");
	    
	    printf("Also:\n\n");
	    printf("--help            This help [-h]\n");
	    
	    exit(EXIT_SUCCESS);
	}

	int option_index = 0;
	int c = getopt_long_only(argc, argv, "b:c:d:g:s:t:y:v:n:i:j:k:l", 
				 long_options, &option_index);
	if (c == -1)
	    break;

	switch (c) {
	case 'b': // band 
	    opts.band = atoi(optarg);
	    break;
	case 'o': // overlap
	    opts.overlap = atof(optarg);
	    break;
	case 'i': // cellsize 
	    cellf = atoi(optarg);
	    break;
	case 'j': // cellsize 
	    cells = atoi(optarg);
	    break;
	case 'k': // cellsize 
	    celld = atoi(optarg);
	    break;
	case 'l': // cellsize 
	    cella = atoi(optarg);
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
	case 'n': // Number of coincidences   
	    noc = atoi(optarg);
	    break;
	case '?':
	    break;
	default:
	    break ;
	} /* switch c */
    } /* while 1 */
    
    
    // Check if sett.nod was set up, if not, exit
    if(!(sett.nod)||!(opts.band)|(opts.overlap < 0.)) { 
        printf("Number of days or band number not set... Exiting\n"); 
	exit(EXIT_FAILURE); 
    } 

    printf("Number of days in time segments: %d\n", sett.nod); 

    printf("Input data: %s\n", datafile);
    printf("Grid matrix data directory: %s\n", griddir);

    // Starting band frequency:
    //sett.fpo = 10. + 0.96875*opts.band*(0.5/sett.dt);
    sett.fpo = 10. + (1. - opts.overlap)*opts.band*(0.5/sett.dt);

    printf("Band number, reference frequency fpo: %04d , %f\n", opts.band, sett.fpo);
    printf("Band veto fraction: %f\n", vetofrac);

    printf("The data sampling time dt: %f\n", sett.dt); 
    printf("FAP threshold: %f\n", threshold); 

    printf("Cell sizes: f=%d s=%d d=%d a=%d\n", cellf, cells, celld, cella); 

	
    // Search settings
    //----------------

    search_settings(&sett); 

    // Observation time  
    To = sett.N*sett.dt;

    // Read grid: gamrn the normalized Fisher matrix  
    //sprintf (filename, "%s/grid.bin", griddir);

    if ((data=fopen (griddir, "rb")) != NULL) {

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
	Nc = round((1.0 - vetofrac)*Nc/(cellf*cells*celld*cella));
    else 
	Nc = round(Nc/(cellf*cells*celld*cella));

    // Read the coincidence data 
    //--------------------------  

    if ((data=fopen (datafile, "r")) != NULL) {

	short int i, shift, nof, band, hemi;
	double fpofile;  

	status = fscanf(data, "%hu", &nof); 

	int Nk[nof], Nku[nof], frn[nof], Nkall=0;

	// Frames information: frame number, no. of candidates, no. of unique candidates 
	for(i=0; i<nof; i++) { 
	    status = fscanf(data, "%d %d %d", &frn[i], &Nk[i], &Nku[i]);
	    //fprintf(stderr,"%d %d %d ", frn[i], Nk[i], Nku[i]); 
	    Nkall += Nku[i]; 
	} 


	double *fap;
	fap = (double *)malloc( (nof+1)*sizeof(double)  );

	FalseAlarmProb(noc, nof, Nc, &Nku[0], &fap[0]);

	printf("FAP results:\n");
	fflush(stdout);

	fprintf(stderr, "%d %d ", nof, Nkall);

	for(i=noc; i<=nof; i++) {
	    //#mb hack - fabs because for nonphysical noc FAP equals -inf  
	    //if(fabs(FAP) < threshold) 
	    fprintf(stderr, "%hu %le ", i, fap[i]); 
	}
    
	fprintf(stderr, "\n"); 
 
    } else {
	perror (filename);
	exit(EXIT_FAILURE);
    }
  
    fclose (data);
  
    return 0;
 
}



// Calculates false alarm probablility of noc coincidences out of L
int FalseAlarmProb(
    int noc,      // minimum number of coincidences
    int L,        // number of all frames = maximum number of coincidences
    double Nc,    // number of cells in each frame
    int *Nk,      // array with numbers of unique candidates in each frame
    double *fap   // array: false alarm probablility of n coincidences for noc..L
    ) {
     
    int i, j, k, l;
    double ee[L][5], eeavg[5]={0.};
    double C[L+1][5], pf[L+1][5];

    gsl_combination *cp, *cq;
    size_t *cpd, *cqd;


    for(i=0; i<L; i++){    //#mb length(Nk)
	for(l=0; l<5; l++){
	    ee[i][l] = Nk[i]/(Nc*pow(2,l));
	    eeavg[l] += ee[i][l];
	}
    }
    for(l=0; l<5; l++)
	eeavg[l] /= L;

    // Calculate C[i] - probability of i coincidences in any given call (and l-th correction for shifts)
    for(i=noc; i<=L; i++) {

	printf("n=%3d ", i);
	double P[5], Q[5], Ctmp[5]={0.};

	double ncomb = gsl_sf_choose(L,i);
	if ( ncomb < 1.e9 ) {
	    printf("[ncomb=%.1e][exa]", ncomb);
	    cp = gsl_combination_calloc(L, i);
	    cq = gsl_combination_alloc(L, L-i);
	    gsl_combination_init_last(cq);
     
	    while(1) {
		cpd = gsl_combination_data(cp);
		cqd = gsl_combination_data(cq);
     
		for(l=0; l<5; ++l) {
		    P[l] = 1.;
		    Q[l] = 1.;
		}
     
		for(j=0; j<i; ++j){
		    for(l=0; l<5; ++l)
			P[l] *= ee[cpd[j]][l];
		}
		for(j=0; j<(L-i); ++j){
		    for(l=0; l<5; ++l)
			Q[l] *= (1. - ee[cqd[j]][l]);
		}

		for(l=0; l<5; ++l)
		    Ctmp[l] += P[l]*Q[l];
	    
		if (gsl_combination_next(cp) == GSL_FAILURE || 
		    gsl_combination_prev(cq) == GSL_FAILURE) 
		    break;
	    } 
  
	    gsl_combination_free (cp);
	    gsl_combination_free (cq);
	       
	} else {

	    printf("[ncomb=%8e][avg]", ncomb);
	    for(l=0; l<5; l++)
		Ctmp[l]= ncomb*pow(eeavg[l],i)*pow(1.-eeavg[l],L-i);

	}

	printf(" C[%3d] = [ ", i);
	for(l=0; l<5; l++){
	    C[i][l] = Ctmp[l]; // probability of i coincidences between L frames
	    printf("%12.6e  ", C[i][l]);
	}
	printf("]\n");

    } // i


    // pf - probability of i or more coincidences in any given cell
    for(i=noc; i<=L; i++) {  
	for(l=0; l<5; l++){
#if 0
	    // like in the old version
	    pf[i][l] = C[i][l];
#else
	    // correct
	    pf[i][l] = 0.;
	    for(j=i; j<=L; j++){
		pf[i][l] += C[j][l] ;
	    }
#endif
	}
    }

    // PF0 = 1 - (1 - pfe).^Nc
    for(i=noc; i<=L; i++) {  

	fap[i] = pow(2,4)*pf[i][0]
	    - ( gsl_sf_choose(4,1)*pf[i][1] 
		+ gsl_sf_choose(4,2)*pf[i][2] 
		+ gsl_sf_choose(4,3)*pf[i][3]
		+ gsl_sf_choose(4,4)*pf[i][4] )
	    - ( gsl_sf_choose(4,2)*pf[i][2]
		+ gsl_sf_choose(4,3)*pf[i][3] 
		+ gsl_sf_choose(4,4)*pf[i][4] )
	    - ( gsl_sf_choose(4,3)*pf[i][3] 
		+ gsl_sf_choose(4,4)*pf[i][4] )
	    - pf[i][4];
	  
	fap[i] = 1. - pow(1. - fap[i], Nc);
	  
    }

#if 0
    printf("FAP results:\n");
    for(i=noc; i<=L; i++){
	printf("%hu %le ", i, fap[i]);
    }
    printf("\n");

    exit(1);
#endif
    return 0; 

}
