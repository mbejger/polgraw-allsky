// MSVC macro to include constants, such as M_PI (include before math.h)
#define _USE_MATH_DEFINES

// Polgraw includes
#include <settings.h>
#include <auxi.h>

// Standard C includes
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>


/// <summary>Create directory for disk output.</summary>
///
void setup_output(struct stat* buff,
                  Command_line_opts* opts)
{
    if (stat(opts->prefix, buff) == -1)
    {
        if (errno == ENOENT)
        {
            // Output directory apparently does not exist, try to create one
#ifdef WIN32
            if (_mkdir(opts->prefix) == -1)
#else
            if (mkdir(opts.prefix, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) == -1)
#endif
            {
                perror(opts->prefix);
                exit(EXIT_FAILURE);
            }
        }
        else // can't access output directory
        {
            perror(opts->prefix);
            exit(EXIT_FAILURE);
        }
    }
}

/// <summary>Search settings: FFT lenghts & other details, bandwidth and Earth parameters.</summary>
///
void search_settings(Search_settings* sett)
{
    double dt, B, oms, omr, Smin, Smax;
    int nod, N, nfft, s, nd, interpftpad;

    dt = sett->dt;                      // data sampling time:  
                                        // set in handle_opts() from the command line
                                        // (the default value is dt=0.5)

    B = 0.5 / dt;                       // Bandwidth
    oms = 2.*M_PI*(sett->fpo)*dt;       // Dimensionless angular frequency

    omr = C_OMEGA_R*dt;

    nod = 2;                            // Observation time in days
    N = lround(nod*C_SIDDAY / dt);      // No. of data points

    nfft = 1 << (int)ceil(log(N) / log(2.));    // length of FFT
    s = 1;                                      // No. of spindowns

//  Smin = 1000.*C_YEARSEC;                     // Minimum spindown time
                                                // [sec.]

//  Smax = 2.*M_PI*(sett->fpo + B)*dt*dt / (2.*Smin); // Maximum spindown (1000 years) [angular, dimensionless]

    //#mb ranges of spindown (RDC O1) 
    double fdotmin, fdotmax;
    fdotmin = 0.5e-8;
    fdotmax = 0.5e-9;

    Smax = 2.*M_PI*fdotmin*dt*dt;
    Smin = 2.*M_PI*fdotmax*dt*dt;

    nd = 2;     // Degree of freedom, 
                // (2*nd = deg. no ofrees of freedom for chi^2)

    interpftpad = 2;

    sett->B = B;          	// bandwidth
    sett->oms = oms;      	// dimensionless angular frequency
    sett->omr = omr;      	// C_OMEGA_R * dt
    sett->nod = nod;      	// number of days of observation
    sett->N = N;          	// number of data points
    sett->nfft = nfft;    	// length of fft
    sett->s = s;          	// number of spindowns
    sett->Smin = Smin;    	// minimum spindown
    sett->Smax = Smax;    	// maximum spindown
    sett->nd = nd;        	// degrees of freedom
    sett->interpftpad = interpftpad;

    // Because of frequency-domain filters, we search
    // F-statistic in range (nmin+1, nmax) of data points
    // 
    // The value of sett->fftpad (zero padding - original grids: 2, new grids: 1) 
    // is read from the grid.bin file in read_grid() (see init.c) 

    sett->nmin = sett->fftpad*NAV*sett->B;
    sett->nmax = (sett->nfft / 2 - NAV*sett->B)*sett->fftpad;

    // initial value of number of known instrumental lines in band 
    sett->numlines_band = 0;

} // search settings

#define buf_size 512

/// <summary>Reads the settings of the detectors.</summary>
/// <remarks>Network of detectors' discovery: finds subdirectories in the main input directory, which by convention should be named like V1, L1, H1 and which contain input data and ephemerids; writes appropriate detector-related data into structs.</remarks>
///
void detectors_settings(Search_settings* sett,
                        Command_line_opts *opts)
{
    int i = 0;
    char dirname[buf_size], x[buf_size];

    // Main input directory name 
    int err = sprintf_s(dirname, 512, "%s/%03d", opts->dtaprefix, opts->ident);
    if (err <= 0)
        perror("Directory name assembly failed.");

    DIR *dp;
    struct dirent *ep;

    char **detnames = (char **)malloc(MAX_DETECTORS * sizeof(char *));
    char **xnames = (char **)malloc(MAX_DETECTORS * sizeof(char *));

    dp = opendir(dirname);
    if (dp != NULL)
    {
        while ((ep = readdir(dp)))
        {
            // Subdirectory names checkup: 
            // check if it's a dir
            // name is 2 char long
            // not a directory name of the type "./" or ".."
            // if usedef is not set (length equal 0), or is set and dir name is substring of it 
            if ((ep->d_type == DT_DIR) &&
                (strlen(ep->d_name) == DETNAME_LENGTH) &&
                (strncmp(&ep->d_name[0], ".", 1)) &&
                (!strlen(opts->usedet) || (strlen(opts->usedet) && (strstr(opts->usedet, ep->d_name)))))
            {
                FILE *data;

                // Input time-domain data handling
                // 
                // We assume that in each subdirectory corresponding 
                // to the detector the input data will look as following:
                // sprintf(x, "%s/%03d/%s/xdatc_%03d%s.bin",
                //         opts->dtaprefix, opts->ident, ep->d_name,
                //         opts->ident, opts->label);

                int err = sprintf_s(x,
                    buf_size,
                    "%s/%03d/%s/xdatc_%03d_%04d%s.bin",
                    opts->dtaprefix,
                    opts->ident,
                    ep->d_name,
                    opts->ident,
                    opts->band,
                    opts->label);
                if (err <= 0)
                    perror("Error assembling string.");

                if ((data = fopen(x, "r")) != NULL) {

                    xnames[i] = (char *)calloc(strlen(x) + 1, sizeof(char));
                    detnames[i] = (char *)calloc(DETNAME_LENGTH + 1, sizeof(char));

                    strncpy(xnames[i], x, strlen(x));
                    strncpy(detnames[i], ep->d_name, DETNAME_LENGTH);
                    i++;

                }
                else {
                    printf("Directory %s exists, but no input file found:\n%s missing...\n",
                        ep->d_name, x);
                    //perror (x);
                }

                fclose(data);
                memset(x, 0, sizeof(x));
            }
        }

        (void)closedir(dp);

    }
    else perror("Couldn't open the input directory...");

    sett->nifo = i;      // number of detectors  
    if (sett->nifo)
    {
        printf("Settings - number of detectors: %d\n", sett->nifo);
    }
    else
    {
        printf("No subdirectories with detector data found. Exiting...\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i<sett->nifo; i++)
    {
        // Virgo detector
        if (!strcmp("V1", detnames[i]))
        {
            strncpy(ifo[i].xdatname, xnames[i], strlen(xnames[i]));
            strncpy(ifo[i].name, detnames[i], DETNAME_LENGTH);

            // Geographical latitude phi in radians
            ifo[i].ephi = (43. + 37. / 60. + 53.0880 / 3600.) / RAD_TO_DEG;
            // Geographical longitude in radians
            ifo[i].elam = (10. + 30. / 60. + 16.1885 / 3600.) / RAD_TO_DEG;
            // Height h above the Earth ellipsoid in meters
            ifo[i].eheight = 53.238;
            // Orientation of the detector gamma
            ifo[i].egam = (135. - (19.0 + 25. / 60.0 + 57.96 / 3600.)) / RAD_TO_DEG;

            printf("Using %s IFO as detector #%d... %s as input time series data\n",
                ifo[i].name, i, ifo[i].xdatname);

            // Hanford H1 detector
        }
        else if (!strcmp("H1", detnames[i]))
        {
            strncpy(ifo[i].xdatname, xnames[i], strlen(xnames[i]));
            strncpy(ifo[i].name, detnames[i], DETNAME_LENGTH);

            // Geographical latitude phi in radians
            ifo[i].ephi = (46 + (27 + 18.528 / 60.) / 60.) / RAD_TO_DEG;
            // Geographical longitude in radians
            ifo[i].elam = -(119 + (24 + 27.5657 / 60.) / 60.) / RAD_TO_DEG;
            // Height h above the Earth ellipsoid in meters
            ifo[i].eheight = 142.554;
            // Orientation of the detector gamma
            ifo[i].egam = 170.9994 / RAD_TO_DEG;

            printf("Using %s IFO as detector #%d... %s as input time series data\n",
                   ifo[i].name, i, ifo[i].xdatname);

            // Livingston L1 detector
        }
        else if (!strcmp("L1", detnames[i]))
        {
            strncpy(ifo[i].xdatname, xnames[i], strlen(xnames[i]));
            strncpy(ifo[i].name, detnames[i], DETNAME_LENGTH);

            // Geographical latitude phi in radians
            ifo[i].ephi = (30 + (33 + 46.4196 / 60.) / 60.) / RAD_TO_DEG;
            // Geographical longitude in radians
            ifo[i].elam = -(90 + (46 + 27.2654 / 60.) / 60.) / RAD_TO_DEG;
            // Height h above the Earth ellipsoid in meters
            ifo[i].eheight = -6.574;
            // Orientation of the detector gamma
            ifo[i].egam = 242.7165 / RAD_TO_DEG;

            printf("Using %s IFO as detector #%d... %s as input time series data\n",
                   ifo[i].name, i, ifo[i].xdatname);
        }
        else
        {

            printf("Meh, unknown detector %s (see settings.c) Exiting...\n",
                   detnames[i]);
            exit(EXIT_FAILURE);
        }
    }

    // memory free for detnames and xdatnames
    for (i = 0; i<sett->nifo; i++)
    {
        free(detnames[i]);
        free(xnames[i]);
    }

    free(detnames);
    free(xnames);

} // detectors settings

#undef buf_size

/// <summary>Coefficients of the amplitude modulation functions of the Virgo detector.</summary>
///
void rogcvir(Detector_settings* ifoi)
{
  // In the notation of Phys. Rev. D 58, 063001 (1998):
  // ephi = lambda (geographical latitude phi in radians)
  // egam = gamma (orientation of the detector)
  // 
  // (see modvir function in jobcore.c for Eqs. 12 and 13)

  ifoi->amod.c1 = .25*sin(2.*ifoi->egam)*(1+sqr(sin(ifoi->ephi)));
  ifoi->amod.c2 = -.5*cos(2.*ifoi->egam)*sin(ifoi->ephi);
  ifoi->amod.c3 = .5*sin(2.*ifoi->egam)*sin(2.*ifoi->ephi);
  ifoi->amod.c4 = -cos(2.*ifoi->egam)*cos(ifoi->ephi);
  ifoi->amod.c5 = .75*sin(2.*ifoi->egam)*sqr(cos(ifoi->ephi));
  ifoi->amod.c6 = cos(2.*ifoi->egam)*sin(ifoi->ephi);
  ifoi->amod.c7 = .5*sin(2.*ifoi->egam)*(1.+sqr(sin(ifoi->ephi)));
  ifoi->amod.c8 = cos(2.*ifoi->egam)*cos(ifoi->ephi);
  ifoi->amod.c9 = .5*sin(2.*ifoi->egam)*sin(2.*ifoi->ephi);

} // rogcvir
