#include "stdio.h"
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <fftw3.h>

struct line{
double freq, ampl, phase, snr;
};

// funtion to calculate the threshold
double treso(double y, int m){ 
    double t= -log(1. - pow(1. - y, 1./(double)m));
    return t; 
}

// funtion to calculate the mean for an array
double mean(unsigned int lngth, double * a){
    if (lngth <=0){perror("Cant calculate mean value of array, array is empty\n"); abort();}
    double m=0;
    unsigned int i;
    for (i=0; i < lngth; i++){ 
        m +=a[i];
    }
    m /= lngth;
    return m;
}

// size for zero padding (in power of 2)
unsigned int nextpow2(unsigned int x){
    unsigned int w;
    double d; 
    if (x > 1){
        d = log2(x);
        w = (unsigned int)trunc(log2(x));
        if (w == d){ return w;} else {return w+1;}
    }
    else {return 1; }
}


int main(int argc, char *argv[]){

    FILE *xdatc_file;
    char xdatc_name[256];
    int indx, yndx, google;
    int P=1;
    double *x;
    double alpha=0.01;
    static int help_flag=0, fine_flag=0;
    double fpo=0; //initial frequency
    double dt=0;
    int exc=8192;  // 2^13
    int nav=16384; // 2^14

    while (1) {
        static struct option long_options[] = {
            {"help", no_argument, &help_flag, 1},    
            {"infile", required_argument, 0, 'r'}, // triggers file
            {"dt value", required_argument, 0, 'd'}, // dt value
            {"fpo value", required_argument, 0, 'v'}, // fpo value
            {"P value", required_argument, 0, 'p'}, // Zerro padding
            {"alpha value", required_argument, 0, 'a'}, // alpha value (fals alarm probability)
            {"exc value", required_argument, 0, 'e'}, // 
            {"nav value", required_argument, 0, 'n'}, // size of the block of bins of the spectrum
            {"fine", no_argument, &fine_flag, 1},    // perform fine calculation 
            {0, 0, 0, 0}
        };

        if(help_flag) {
            printf("*** Searching  monochromatic signals in data sequence ***\n"); 
            printf("Usage: ./pspm2fx2 [switch1] <value1> [switch2] <value2> ...\n") ;
            printf("Switches are:\n\n"); 
            printf("-r                  input file\n"); 
            printf("--dt    (or -d)     dt  value\n");
            printf("--fpo   (or -v)     fpo (starting frequency) value\n");
            printf("--P     (or -p)     zero padding value\n");
            printf("--nav   (or -n)     nav (default value %d)\n", nav);
            printf("--alpha (or -a)     false alarm probability (default value %2.2f)\n", alpha);
            printf("--exp   (or -e)     exc (default value %d)\n", exc);
            printf("--fine              perform fine frequency value calculation\n\n");
            printf("Also:\n\n"); 
            printf("--help    This help\n");         
            exit (0);
        }
        int option_index = 0;
        int cc = getopt_long (argc, argv, "r:d:v:p:n:a:e:", long_options, &option_index);
        if (cc == -1)     break;

        switch (cc) {
            case 'r':
                strcpy (xdatc_name, optarg);
                break;
            case 'v':
                fpo = atof(optarg); 
                break; 
            case 'a':
                alpha = atof(optarg);
                break;
            case 'd':
                dt = atof(optarg);
                break;
            case 'e':
                exc = atoi(optarg);
                break;
            case 'n':
                nav = atoi(optarg);
                break;
            case 'p':
                P = atoi(optarg);
                break;
            case '?': abort(); 
            default: break;
                 
        }
    }

    if ((dt == 0) || (fpo == 0)) {perror("The \"dt\" or \"fpo\" is not set\n"); abort();}
    
    xdatc_file = fopen(xdatc_name,"r");
    fseek(xdatc_file, 0, SEEK_END); // seek to end of file
    unsigned int file_length = ftell(xdatc_file) / sizeof(double); // get size od the file
    rewind(xdatc_file); // go bet to begining of the file
    
    unsigned int lenX =(P+1)*pow(2, nextpow2(file_length));       //set to power of 2 and zero padding
    x = (double*)malloc(lenX*sizeof(double));
    
    google = 0;
    
    // read xdat file
    while (fread(&x[google], sizeof(double), 1, xdatc_file)){
        if(feof(xdatc_file)) break;
        else{
            google += 1;
        }
    }
    fclose(xdatc_file);

    //Remove the mean
    double m = mean(google, x);
    for (indx=0; indx < google; indx++){ x[indx] -= m; }
    // Fulfill the reast of the array by zeros
    for(indx = google; indx < lenX; indx++){x[indx]=0;}

    /**** FFT ****/
    fftw_complex *fft_in, *fft_out;
    fftw_plan p;
    fft_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * lenX);
    for (indx=0; indx < lenX; indx++){fft_in[indx][0] = x[indx]; fft_in[indx][1] = 0; }
    fft_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * lenX);
    p = fftw_plan_dft_1d(lenX, fft_in, fft_out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    //Take only positive frequencies
    unsigned int lenxp=lenX / 2;
    double * xp = (double*)malloc(sizeof(double)*lenxp);
    for (indx=0; indx < lenxp; indx++){ xp[indx] = fft_out[indx][0]*fft_out[indx][0] + fft_out[indx][1]*fft_out[indx][1];}
    
    fftw_destroy_plan(p);
    fftw_free(fft_in); 
        
    if (exc != 0){ 
        for (indx=0; indx < exc; indx++) { xp[indx] = 0;  xp[lenxp - indx - 1] = 0;}
    }

    // Calculate threshold corresponding to the false alarm probability alpha
    double tr = treso(alpha, file_length);
    
    // Calculate mean values
    unsigned int   part=(lenxp - 2*exc) / nav;
    double * mean_ar=(double*)malloc(sizeof(double)*part);
    double g;
    for (indx=0;  indx < part; indx++){
        g=0;
        for( yndx=0; yndx < nav; yndx++){ g += xp[exc + indx*nav + yndx]; }
        mean_ar[indx]=g/nav;
    }
    
    // Remove mean
    for (indx=0;  indx < part; indx++){
        for (yndx=0; yndx < nav; yndx++){ 
            xp[exc + indx*nav + yndx] /= mean_ar[indx];
        }
    }

    // Find and count local maxima
    google=0;
    for(indx=1; indx < lenxp-1; indx++){if ((xp[indx] > xp[indx-1]) && (xp[indx] > xp[indx+1]) && (xp[indx] > tr)){google++;};}

    // Create the array of local maxima
    int maxarsize=google;
    struct line *maxar=malloc(sizeof(struct line)*maxarsize);
    
    double bin, am;
    google=0;
    if(maxarsize == 0){
        perror("NO SIGNALS FOUND !");
    }else{
        if (fine_flag){
            /* Fine search */
            for(indx=1; indx < lenxp-1; indx++){
                if ((xp[indx] > xp[indx-1]) && (xp[indx] > xp[indx+1]) && (xp[indx] > tr)){
                    /* Parabolic interpolation for frequency */
                    bin=indx + (xp[indx-1] - xp[indx+1]) / (2*(xp[indx-1] - 2*xp[indx] + xp[indx+1])); 
                    maxar[google].freq = fpo + bin/(lenX*dt);
                    /* Liniar interpolation for the phase of signal */
                    if( bin > indx){
                        maxar[google].phase = atan(fft_out[indx][1]/fft_out[indx][0]) + 
                            (atan(fft_out[indx+1][1]/fft_out[indx+1][0]) - atan(fft_out[indx][1]/fft_out[indx][0])) * (bin - indx);
                    }else{
                        maxar[google].phase = atan(fft_out[indx][1]/fft_out[indx][0]) + 
                            (atan(fft_out[indx][1]/fft_out[indx][0]) - atan(fft_out[indx-1][1]/fft_out[indx-1][0])) * (indx - bin);}
                        /* ------ Amplitude of a signal ----- */  
                        am = xp[indx] - (xp[indx+1] - xp[indx-1])*(xp[indx+1] - xp[indx-1])/(8*(xp[indx+1] - 2*xp[indx] + xp[indx-1])); // linear interpolation
                        maxar[google].ampl =(2/(double)nav)*sqrt(am*mean_ar[(maxarsize-2*exc)/nav]);
                        /* ----- Signal to nois ratio of a signal ------ */  
                        maxar[google].snr = sqrt(2*(am - 1));
                        google++;        
                    };
            };
        }else{
            /* Coarse calculation of monochromatic signales */
            for( indx=1; indx<lenxp; indx++){
                if ((xp[indx] > xp[indx-1]) && (xp[indx] > xp[indx+1]) && (xp[indx] > tr)){
                    maxar[google].freq = fpo + indx/(lenX*dt);
                    /*----- Phase of signal ----*/
                    maxar[google].phase = atan(fft_out[indx][1]/fft_out[indx][0]); 
                    am = xp[indx];
                    /* ------ Amplitude of a signal ----- */
                    maxar[google].ampl = (2/(double)nav)*sqrt(am*mean_ar[(maxarsize-2*exc)/nav]);
                    /* ----- Signal to nois ratio of a signal ------ */
                    maxar[google].snr = sqrt(2*(am - 1));
                    google++;                    
                }
            }
        }

        // Print summary
        printf("%% Summary for the file %s\n", xdatc_name);
        if (fine_flag){printf("%% Fine search\n");} else {printf("%% Coars search\n");}
        printf("%% fpo = %f\n", fpo);
        printf("%% dt = %f\n", dt);
        printf("%% alpha = %5.4f\n", alpha);
        printf("%% nav = %d\n", nav);
        printf("%% exc = %d \n", exc);
        printf("%% P = %d\n", P);
        printf("%%   Nr.    Frequency      Amplitude    h  SNR          Phase\n");
        for(indx=0; indx < google; indx++){printf("%5d %15.6f %12.5f %12.5f %12.5f\n", indx+1, maxar[indx].freq, maxar[indx].ampl, maxar[indx].snr, maxar[indx].phase);}
    }

    fftw_free(fft_out);
    return 0;
}
