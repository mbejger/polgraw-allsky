// MSVC macro to include constants, such as M_PI (include before math.h)
#define _USE_MATH_DEFINES

// Polgraw includes
#include <init.h>       // all function declarations
#include <struct.h>     // Search_settings, Command_line_opts, OpenCL_handles, ...
#include <settings.h>
#include <auxi.h>
#include <spline_z.h>

// OpenCL includes
#include <CL/cl.h>

// Standard C includes
#include <stdio.h>
#include <stdlib.h>     // EXIT_FAILURE
//#include <unistd.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <getopt.h>
#include <time.h>

#include <direct.h>


/// <summary>Command line options handling: search</summary>
///
void handle_opts(Search_settings* sett,
                 OpenCL_settings* cl_sett,
		         Command_line_opts* opts,
		         int argc, 
		         char* argv[]) {
  
  opts->hemi=0;
  opts->wd=NULL;

  // Default F-statistic threshold 
  opts->trl=20;
	
  strcpy (opts->prefix, TOSTR(PREFIX));
  strcpy (opts->dtaprefix, TOSTR(DTAPREFIX));

  opts->label[0]  = '\0';
  opts->range[0]  = '\0';
  opts->getrange[0] = '\0';
  opts->usedet[0]   = '\0';
  opts->addsig[0] = '\0';
	
  // Initial value of starting frequency set to a negative quantity. 
  // If this is not changed by the command line value, fpo is calculated 
  // from the band number b (fpo = fpo = fstart + 0.96875*b/(2dt))
  sett->fpo = -1;

  // Default initial value of the data sampling time 
  sett->dt = 0.5; 

  opts->help_flag=0;
  opts->white_flag=0;
  opts->s0_flag=0;
  opts->checkp_flag=0;

  static int help_flag=0, white_flag=0, s0_flag=0, checkp_flag=1;

  // Reading arguments 

  while (1) {
    static struct option long_options[] = {
      {"help", no_argument, &help_flag, 1},
      {"whitenoise", no_argument, &white_flag, 1},
      {"nospindown", no_argument, &s0_flag, 1},
      {"nocheckpoint", no_argument, &checkp_flag, 0},
      // frame number
      {"ident", required_argument, 0, 'i'},
      // frequency band number
      {"band", required_argument, 0, 'b'},
      // output directory
      {"output", required_argument, 0, 'o'},
      // input data directory
      {"data", required_argument, 0, 'd'},
      // non-standard label for naming files
      {"label", required_argument, 0, 'l'},
      // narrower grid range parameter file
      {"range", required_argument, 0, 'r'},
      // write full grid range to file
      {"getrange", required_argument, 0, 'g'},
      // change directory parameter
      {"cwd", required_argument, 0, 'c'},
      // interpolation method
      {"threshold", required_argument, 0, 't'},
      // hemisphere
      {"hemisphere", required_argument, 0, 'h'},
      // fpo value
      {"fpo", required_argument, 0, 'p'},
      // add signal parameters
      {"addsig", required_argument, 0, 'x'},
      // which detectors to use
      {"usedet", required_argument, 0, 'u'}, 
      // data sampling time 
      {"dt", required_argument, 0, 's'},
      {0, 0, 0, 0}
    };

    if (help_flag) {

      printf("polgraw-allsky periodic GWs: search for candidate signals with the F-statistic\n");
      printf("Usage: ./search -[switch1] <value1> -[switch2] <value2> ...\n") ;
      printf("Switches are:\n\n");
      printf("-d, -data         Data directory (default is .)\n");
      printf("-o, -output       Output directory (default is ./candidates)\n");
      printf("-i, -ident        Frame number\n");
      printf("-b, -band         Band number\n");
      printf("-l, -label        Custom label for the input and output files\n");
      printf("-r, -range        Use file with grid range or pulsar position\n");
      printf("-g, -getrange     Write grid ranges & exit (ignore -r)\n");
      printf("-c, -cwd          Change to directory <dir>\n");
      printf("-t, -threshold    Threshold for the F-statistic (default is 20)\n");
      printf("-h, -hemisphere   Hemisphere (default is 0 - does both)\n");
      printf("-p, -fpo          Reference band frequency fpo value\n");
      printf("-s, -dt           data sampling time dt (default value: 0.5)\n");
      printf("-u, -usedet       Use only detectors from string (default is use all available)\n");
      printf("-x, -addsig       Add signal with parameters from <file>\n\n");


      printf("Also:\n\n");
      printf("--whitenoise      White Gaussian noise assumed\n");
      printf("--nospindown      Spindowns neglected\n");
      printf("--nocheckpoint    State file won't be created (no checkpointing)\n");
      printf("--help            This help\n");

      exit(EXIT_SUCCESS);
    }

    int option_index = 0;
    int c = getopt_long_only(argc, argv, "i:b:o:d:l:r:g:c:t:h:p:x:s:u:", 
			     long_options, &option_index);
    if (c == -1)
      break;

    switch (c) {
    case 'i':
      opts->ident = atoi (optarg);
      break;
    case 't':
      opts->trl = atof(optarg);
      break;
    case 'h':
      opts->hemi = atof(optarg);
      break;
    case 'b':
      opts->band = atoi(optarg);
      break;
    case 'o':
      strcpy(opts->prefix, optarg);
      break;
    case 'd':
      strcpy(opts->dtaprefix, optarg);
      break;
    case 'l':
      opts->label[0] = '_';
      strcpy(1+opts->label, optarg);
      break;
    case 'r':
      strcpy(opts->range, optarg);
      break;
    case 'g':
      strcpy(opts->getrange, optarg);
      break;
    case 'c':
      opts->wd = (char *) malloc (1+strlen(optarg));
      strcpy(opts->wd, optarg);
      break;
    case 'p':
      sett->fpo = atof(optarg);
      break;
    case 'x':
      strcpy(opts->addsig, optarg);
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

  opts->white_flag = white_flag;
  opts->s0_flag = s0_flag;
  opts->checkp_flag = checkp_flag;	
	
  printf("Input data directory is %s\n", opts->dtaprefix);
  printf("Output directory is %s\n", opts->prefix);
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
  printf("The data sampling time dt is  %f\n", sett->dt); 

  if (opts->white_flag)
    printf ("Assuming white Gaussian noise\n");

  // For legacy: FFT is now the only option 
  printf ("Using fftinterp=FFT (FFT interpolation by zero-padding)\n");

  if(opts->trl!=20)
    printf ("Threshold for the F-statistic is %lf\n", opts->trl);
  if(opts->hemi)
    printf ("Search for hemisphere %d\n", opts->hemi);
  if (opts->s0_flag)
    printf ("Assuming s_1 = 0.\n");
  if (strlen(opts->label))
    printf ("Using '%s' as data label\n", opts->label);

  if(strlen(opts->getrange)){
    printf ("Writing full grid ranges to '%s'\n", opts->getrange);
    if(strlen(opts->range)) {
      opts->range[0] = '\0';
      printf ("     WARNING! -r option will be ignored...\n");
    }
  }

  if (strlen(opts->range))
    printf ("Obtaining grid range from '%s'\n", opts->range);

  if (strlen(opts->addsig))
    printf ("Adding signal from '%s'\n", opts->addsig);
  if (opts->wd) {
    printf ("Changing working directory to %s\n", opts->wd);
#ifdef WIN32
    if (_chdir(opts->wd)) { perror(opts->wd); abort(); }
#else
    if (chdir(opts->wd)) { perror (opts->wd); abort (); }
#endif
  }

} // end of command line options handling 

/// <summary>OpenCL error handling function.</summary>
/// <remarks>If an error occurs, prints it to standard error and exits.</remarks>
///
void checkErr(cl_int err,
              const char * name)
{
    if (err != CL_SUCCESS)
    {
        perror("ERROR: %s (%i)\n", name, err);
        exit(err);
    }
}

/// <summary>Initialize OpenCL devices based on user preference.</summary>
/// <remarks>Currently, only a sinle platform can be selected.</remarks>
///
void init_opencl(OpenCL_handles* cl_handles,
                 OpenCL_settings* cl_sett)
{
    cl_handles->plat = select_platform(cl_sett->plat_id);

    cl_handles->devs = select_devices(cl_handles->plat,
                                      cl_sett->dev_type,
                                      &cl_handles->dev_count);

    cl_handles->ctx = create_standard_context(cl_handles->devs,
                                              cl_handles->dev_count);

    cl_handles->write_queues = create_command_queue_set(cl_handles->ctx);
    cl_handles->exec_queues  = create_command_queue_set(cl_handles->ctx);
    cl_handles->read_queues  = create_command_queue_set(cl_handles->ctx);

    const char* source = load_program_file(NULL);

    cl_handles->prog = build_program_source(cl_handles->ctx, source);

    cl_handles->kernels = create_kernels(cl_handles->prog);

    free(source);
}

/// <summary>Tries selecting the platform with the specified index.</summary>
///
cl_platform_id select_platform(cl_uint plat_id)
{
    cl_int CL_err = CL_SUCCESS;
    cl_platform_id result = NULL;
    cl_uint numPlatforms = 0;

    CL_err = clGetPlatformIDs(0, NULL, &numPlatforms);
    checkErr(CL_err, "clGetPlatformIDs(numPlatforms)");

    if (plat_id > (numPlatforms - 1))
    {
        perror("Platform of the specified index does not exist.");
        exit(-1);
    }
    else
    {
        cl_platform_id* platforms = (cl_platform_id*)malloc(numPlatforms * sizeof(cl_platform_id));
        CL_err = clGetPlatformIDs(numPlatforms, platforms, NULL);
        checkErr(CL_err, "clGetPlatformIDs(platforms)");

        result = platforms[plat_id];
        
        free(platforms);
    }

    return result;
}

/// <summary>Selects all devices of the specified type.</summary>
///
cl_device_id* select_devices(cl_platform_id platform,
                             cl_device_type dev_type,
                             cl_uint* count)
{
    cl_int CL_err = CL_SUCCESS;
    cl_device_id* result = NULL;

    CL_err = clGetDeviceIDs(platform, dev_type, 0, 0, count);
    checkErr(CL_err, "clGetDeviceIDs(numDevices)");

    if (*count == 0u)
    {
        perror("No devices of the specified type are found on the specified platform.");
        exit(-1);
    }

    result = (cl_device_id*)malloc(*count * sizeof(cl_device_id));
    CL_err = clGetDeviceIDs(platform, dev_type, *count, result, 0);
    checkErr(CL_err, "clGetDeviceIDs(devices)");

    return result;
}

/// <summary>Create a contet that holds all specified devices.</summary>
///
cl_context create_standard_context(cl_device_id* devices, cl_uint count)
{
    cl_int CL_err = CL_SUCCESS;
    cl_uint count = 0;
    cl_context result = NULL;
    cl_platform_id platform = NULL;

    CL_err = clGetDeviceInfo(devices[0],
                             CL_DEVICE_PLATFORM,
                             sizeof(cl_platform_id),
                             &platform,
                             NULL);
    checkErr(CL_err, "clGetDeviceInfo(CL_DEVICE_PLATFORM)");

    cl_context_properties cps[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties)platform, 0 };

    result = clCreateContext(cps, count, devices, NULL, NULL, &CL_err);
    checkErr(CL_err, "clCreateContext()");

    return result;
}

/// <summary>Create a set of command queues to all the devices in the context.</summary>
///
cl_command_queue* create_command_queue_set(cl_context context)
{
    cl_int CL_err = CL_SUCCESS;
    cl_uint count = 0;
    cl_device_id* devices = NULL;
    cl_command_queue* result = NULL;

    CL_err = clGetContextInfo(context, CL_CONTEXT_NUM_DEVICES, sizeof(cl_uint), &count, NULL);
    checkErr(CL_err, "clGetContextInfo(CL_CONTEXT_NUM_DEVICES)");

    CL_err = clGetContextInfo(context, CL_CONTEXT_DEVICES, sizeof(cl_device_id*), devices, NULL);
    checkErr(CL_err, "clGetContextInfo(CL_CONTEXT_DEVICES)");

    result = (cl_command_queue*)malloc(count * sizeof(cl_command_queue));

    for (cl_uint i = 0; i < count; ++i)
    {
        result[i] = clCreateCommandQueue(context, devices[i], CL_QUEUE_PROFILING_ENABLE, &CL_err);
        checkErr(CL_err, "clCreateCommandQueue()");

        CL_err = clReleaseDevice(devices[i]);
        checkErr(CL_err, "clReleaseDevice()");
    }

    free(devices);

    return result;
}

/// <summary>Load kernel file from disk.</summary>
///
char* load_program_file(const char* filename)
{
    long int size = 0;
    size_t res = 0;
    char* src = NULL;
    FILE* file = NULL;
    errno_t err = 0;

    err = fopen_s(&file, filename, "rb");

    if (!file)
    {
        printf_s("Failed to open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    if (fseek(file, 0, SEEK_END))
    {
        fclose(file);
        return NULL;
    }

    size = ftell(file);
    if (size == 0)
    {
        fclose(file);
        return NULL;
    }

    rewind(file);

    src = (char *)calloc(size + 1, sizeof(char));
    if (!src)
    {
        src = NULL;
        fclose(file);
        return src;
    }

    res = fread(src, 1u, sizeof(char) * size, file);
    if (res != sizeof(char) * size)
    {
        fclose(file);
        free(src);

        return src;
    }

    src[size] = '\0'; /* NULL terminated */
    fclose(file);

    return src;
}

/// <summary>Build program file.</summary>
///
cl_program build_program_source(cl_context context,
                                const char* source)
{
    cl_int CL_err = CL_SUCCESS;
    cl_program result = NULL;

    cl_uint numDevices = 0;
    cl_device_id* devices = NULL;

    const size_t length = strnlen_s(source, UINT_MAX);

    result = clCreateProgramWithSource(context, 1, &source, &length, &CL_err);
    checkErr(CL_err, "clCreateProgramWithSource()");

    CL_err = clGetContextInfo(context, CL_CONTEXT_NUM_DEVICES, sizeof(cl_uint), &numDevices, NULL);
    checkErr(CL_err, "clGetContextInfo(CL_CONTEXT_NUM_DEVICES)");
    devices = (cl_device_id*)malloc(numDevices * sizeof(cl_device_id));
    CL_err = clGetContextInfo(context, CL_CONTEXT_DEVICES, numDevices * sizeof(cl_device_id), devices, NULL);
    checkErr(CL_err, "clGetContextInfo(CL_CONTEXT_DEVICES)");

    // Warnings will be treated like errors, this is useful for debug
    char build_params[] = { "-Werror" };
    CL_err = clBuildProgram(result, numDevices, devices, build_params, NULL, NULL);

    if (CL_err != CL_SUCCESS)
    {
        size_t len = 0;
        char *buffer;

        CL_err = clGetProgramBuildInfo(result, devices[0], CL_PROGRAM_BUILD_LOG, 0, NULL, &len);
        checkErr(CL_err, "clGetProgramBuildInfo(CL_PROGRAM_BUILD_LOG)");

        buffer = calloc(len, sizeof(char));

        clGetProgramBuildInfo(result, devices[0], CL_PROGRAM_BUILD_LOG, len, buffer, NULL);

        fprintf(stderr, "%s\n", buffer);

        free(buffer);

        exit(CL_err);
    }

    return result;
}

/// <summary>Create a kernel for all entry points in the program.</summary>
///
cl_kernel* create_kernels(cl_program program)
{
    cl_int CL_err = CL_SUCCESS;
    cl_uint kernel_count = 1;
    cl_kernel* result = (cl_kernel*)malloc(kernel_count * sizeof(cl_kernel));

    for (cl_uint i = 0; i < kernel_count; ++i)
        result = obtain_kernel(program, i);

    return result;
}

/// <summary>Obtain the name of the kernel of a given index.</summary>
///
const char* obtain_kernel_name(cl_uint i)
{
    const char* result = NULL;

    switch (i)
    {
    case ComputeSinCosModF:
        result = "compute_sincosmodf";
        break;
    default:
        perror("Unkown kernel index");
        exit(-1);
        break;
    }

    return result;
}

/// <summary>Obtain kernel with the specified index.</summary>
///
cl_kernel obtain_kernel(cl_program program, cl_uint i)
{
    cl_int CL_err = CL_SUCCESS;
    cl_kernel result = NULL;

    result = clCreateKernel(program, obtain_kernel_name(i), &CL_err);
    checkErr(CL_err, "clCreateKernel()");

    return result;
}

/// <summary>Generate grid from the M matrix.</summary>
/// <remarks>Processes the file 'grid.bin'</remarks>
///
void read_grid(Search_settings *sett,
               Command_line_opts *opts)
{
    sett->M = (double *)calloc(16, sizeof(double));

    FILE *data;
    char filename[512];

    // In case when -usedet option is used for one detector
    // i.e. opts->usedet has a length of 2 (e.g. H1 or V1), 
    // read grid.bin from this detector subdirectory 
    // (see detectors_settings() in settings.c for details) 
    if (strlen(opts->usedet) == 2)
        sprintf(filename, "%s/%03d/%s/grid.bin", opts->dtaprefix, opts->ident, opts->usedet);
    else
        sprintf(filename, "%s/%03d/grid.bin", opts->dtaprefix, opts->ident);

    if ((data = fopen(filename, "r")) != NULL)
    {
        printf("Using grid file from %s\n", filename);
        fread((void *)&sett->fftpad, sizeof(int), 1, data);
        printf("Using fftpad from the grid file: %d\n", sett->fftpad);

        // M: vector of 16 components consisting of 4 rows
        // of 4x4 grid-generating matrix
        fread((void *)sett->M, sizeof(double), 16, data);
        fclose(data);
    }
    else
    {
        perror(filename);
        exit(EXIT_FAILURE);
    }

} // end of read grid 

/// <summary>Initialize auxiliary and F-statistics arrays.</summary>
///
void init_arrays(Search_settings* sett,
                 OpenCL_handles* cl_handles,
                 Command_line_opts* opts,
                 Aux_arrays *aux_arr,
                 cl_mem* F_d)
{
    cl_int CL_err = CL_SUCCESS;
    int i, status;

    // Allocates and initializes to zero the data, detector ephemeris
    // and the F-statistic arrays

    FILE *data;

    sett->Ninterp = sett->interpftpad*sett->nfft;
    sett->nfftf = sett->fftpad*sett->nfft;

    for (i = 0; i<sett->nifo; i++)
    {
        ifo[i].sig.xDat = (real_t*)calloc(sett->N, sizeof(real_t));

        ifo[i].sig.xDat_d = clCreateBuffer(cl_handles->ctx,
                                           CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                           sett->N * sizeof(real_t),
                                           ifo[i].sig.xDat,
                                           &CL_err);
        checkErr(CL_err, "clCreateBuffer(ifo[i].sig.xDat_d)");

        // Input time-domain data handling
        // 
        // The file name ifo[i].xdatname is constructed 
        // in settings.c, while looking for the detector 
        // subdirectories
        if ((data = fopen(ifo[i].xdatname, "r")) != NULL)
        {
            status = fread((void *)(ifo[i].sig.xDat),
                           sizeof(real_t),
                           sett->N,
                           data);
            fclose(data);
        }
        else
        {
            perror(ifo[i].xdatname);
            exit(EXIT_FAILURE);
        }

        int j, Nzeros = 0;
        // Checking for null values in the data
        for (j = 0; j < sett->N; j++)
            if (!ifo[i].sig.xDat[j]) Nzeros++;

        ifo[i].sig.Nzeros = Nzeros;

        // factor N/(N - Nzeros) to account for null values in the data
        ifo[i].sig.crf0 = (real_t)sett->N / (sett->N - ifo[i].sig.Nzeros);

        // Estimation of the variance for each detector 
        ifo[i].sig.sig2 = (ifo[i].sig.crf0)*var(ifo[i].sig.xDat, sett->N);

        ifo[i].sig.DetSSB = (real_t*)calloc(3 * sett->N, sizeof(real_t));

        ifo[i].sig.DetSSB_d = clCreateBuffer(cl_handles->ctx,
                                             CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                             3 * sett->N * sizeof(real_t),
                                             ifo[i].sig.DetSSB,
                                             &CL_err);
        checkErr(CL_err, "clCreateBuffer(ifo[i].sig.DetSSB_d)");

        // Ephemeris file handling
        char filename[512];

        sprintf(filename,
                "%s/%03d/%s/DetSSB.bin",
                opts->dtaprefix,
                opts->ident,
                ifo[i].name);

        if ((data = fopen(filename, "r")) != NULL)
        {
            // Detector position w.r.t Solar System Baricenter
            // for every datapoint
            status = fread((void *)(ifo[i].sig.DetSSB),
                           sizeof(real_t),
                           3 * sett->N,
                           data);

            // Deterministic phase defining the position of the Earth
            // in its diurnal motion at t=0 
            status = fread((void *)(&ifo[i].sig.phir),
                           sizeof(real_t),
                           1,
                           data);

            // Earth's axis inclination to the ecliptic at t=0
            status = fread((void *)(&ifo[i].sig.epsm),
                           sizeof(real_t),
                           1,
                           data);
            fclose(data);

            printf("Using %s as detector %s ephemerids...\n", filename, ifo[i].name);

        }
        else
        {
            perror(filename);
            return;
        }

        // sincos 
        ifo[i].sig.sphir = sin(ifo[i].sig.phir);
        ifo[i].sig.cphir = cos(ifo[i].sig.phir);
        ifo[i].sig.sepsm = sin(ifo[i].sig.epsm);
        ifo[i].sig.cepsm = cos(ifo[i].sig.epsm);

        sett->sepsm = ifo[i].sig.sepsm;
        sett->cepsm = ifo[i].sig.cepsm;

        ifo[i].sig.xDatma_d = clCreateBuffer(cl_handles->ctx,
                                             CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR,
                                             sett->N * sizeof(complex_devt),
                                             NULL,
                                             &CL_err);
        checkErr(CL_err, "clCreateBuffer(ifo[i].sig.xDatma_d)");

        ifo[i].sig.xDatmb_d = clCreateBuffer(cl_handles->ctx,
                                             CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR,
                                             sett->N * sizeof(complex_devt),
                                             NULL,
                                             &CL_err);
        checkErr(CL_err, "clCreateBuffer(ifo[i].sig.xDatmb_d)");

        ifo[i].sig.aa_d = clCreateBuffer(cl_handles->ctx,
                                         CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR,
                                         sett->N * sizeof(real_t),
                                         NULL,
                                         &CL_err);
        checkErr(CL_err, "clCreateBuffer(ifo[i].sig.aa_d)");

        ifo[i].sig.bb_d = clCreateBuffer(cl_handles->ctx,
                                         CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR,
                                         sett->N * sizeof(real_t),
                                         NULL,
                                         &CL_err);
        checkErr(CL_err, "clCreateBuffer(ifo[i].sig.bb_d)");

        ifo[i].sig.shft_d = clCreateBuffer(cl_handles->ctx,
                                           CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR,
                                           sett->N * sizeof(real_t),
                                           NULL,
                                           &CL_err);
        checkErr(CL_err, "clCreateBuffer(ifo[i].sig.shft_d)");

        ifo[i].sig.shftf_d = clCreateBuffer(cl_handles->ctx,
                                            CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR,
                                            sett->N * sizeof(real_t),
                                            NULL,
                                            &CL_err);
        checkErr(CL_err, "clCreateBuffer(ifo[i].sig.shftf_d)");

    } // end loop for detectors 

      // Check if the ephemerids have the same epsm parameter
    for (i = 1; i<sett->nifo; i++)
    {
        if (!(ifo[i - 1].sig.sepsm == ifo[i].sig.sepsm))
        {
            printf("The parameter epsm (DetSSB.bin) differs for detectors %s and %s. Aborting...\n",
                   ifo[i - 1].name,
                   ifo[i].name);
            exit(EXIT_FAILURE);
        }

    }

    // if all is well with epsm, take the first value 
    sett->sepsm = ifo[0].sig.sepsm;
    sett->cepsm = ifo[0].sig.cepsm;

    *F_d = clCreateBuffer(cl_handles->ctx,
                          CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR,
                          2 * sett->nfft * sizeof(real_t),
                          NULL,
                          &CL_err);
    checkErr(CL_err, "clCreateBuffer(F_d)");

    // Auxiliary arrays, Earth's rotation

    aux_arr->t2_d = clCreateBuffer(cl_handles->ctx,
                                   CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR,
                                   sett->N * sizeof(real_t),
                                   NULL,
                                   &CL_err);
    checkErr(CL_err, "clCreateBuffer(aux_arr->t2_d)");

    aux_arr->cosmodf_d = clCreateBuffer(cl_handles->ctx,
                                        CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR,
                                        sett->N * sizeof(real_t),
                                        NULL,
                                        &CL_err);
    checkErr(CL_err, "clCreateBuffer(aux_arr->cosmodf_d)");

    aux_arr->sinmodf_d = clCreateBuffer(cl_handles->ctx,
                                        CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR,
                                        sett->N * sizeof(real_t),
                                        NULL,
                                        &CL_err);
    checkErr(CL_err, "clCreateBuffer(aux_arr->sinmodf_d)");

    aux_arr->tshift_d = clCreateBuffer(cl_handles->ctx,
                                       CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR,
                                       sett->N * sizeof(real_t),
                                       NULL,
                                       &CL_err);
    checkErr(CL_err, "clCreateBuffer(aux_arr->tshift_d)");

    aux_arr->ifo_amod_d = clCreateBuffer(cl_handles->ctx,
                                         CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR,
                                         sett->nifo * sizeof(Ampl_mod_coeff),
                                         NULL,
                                         &CL_err);
    checkErr(CL_err, "clCreateBuffer(aux_arr->ifo_amod_d)");

    init_spline_matrices(cl_handles,
                         &aux_arr->diag_d,
                         &aux_arr->ldiag_d,
                         &aux_arr->udiag_d,
                         &aux_arr->B_d,
                         sett->Ninterp);

    CL_err = clSetKernelArg(cl_handles->kernels[ComputeSinCosModF], 0, sizeof(cl_mem), &aux_arr->sinmodf_d);
    checkErr(CL_err, "clSetKernelArg(0)");
    CL_err = clSetKernelArg(cl_handles->kernels[ComputeSinCosModF], 1, sizeof(cl_mem), &aux_arr->cosmodf_d);
    checkErr(CL_err, "clSetKernelArg(1)");
    CL_err = clSetKernelArg(cl_handles->kernels[ComputeSinCosModF], 2, sizeof(real_t), &sett->omr);
    checkErr(CL_err, "clSetKernelArg(2)");
    CL_err = clSetKernelArg(cl_handles->kernels[ComputeSinCosModF], 3, sizeof(cl_int), &sett->N);
    checkErr(CL_err, "clSetKernelArg(3)");

    cl_event exec;
    CL_err = clEnqueueNDRangeKernel(cl_handles->exec_queues[0],
                                    cl_handles->kernels[ComputeSinCosModF],
                                    1,
                                    NULL,
                                    &sett->N,
                                    NULL,
                                    0,
                                    NULL,
                                    &exec);
    checkErr(CL_err, "clEnqueueNDRangeKernel(ComputeSinCosModF)");

    CL_err = clWaitForEvents(1, &exec);
    checkErr(CL_err, "clWaitForEvents(exec)");

    // OpenCL cleanup
    clReleaseEvent(exec);

} // end of init arrays 

/// <summary>Set search ranges based on user preference.</summary>
///
void set_search_range(Search_settings *sett,
                      Command_line_opts *opts,
                      Search_range *s_range)
{
    // Hemispheres (with respect to the ecliptic)
    if (opts->hemi)
    {
        s_range->pmr[0] = opts->hemi;
        s_range->pmr[1] = opts->hemi;
    }
    else
    {
        s_range->pmr[0] = 1;
        s_range->pmr[1] = 2;
    }

    // If the parameter range is invoked, the search is performed
    // within the range of grid parameters from an ascii file
    // ("-r range_file" from the command line)
    FILE *data;

    if (strlen(opts->range))
    {
        if ((data = fopen(opts->range, "r")) != NULL)
        {
            int aqq = fscanf(data, "%d %d %d %d %d %d %d %d",
                             s_range->spndr, 1 + s_range->spndr, s_range->nr,
                             1 + s_range->nr, s_range->mr, 1 + s_range->mr,
                             s_range->pmr, 1 + s_range->pmr);

            if (aqq != 8)
            {
                printf("Error when reading range file!\n");
                exit(EXIT_FAILURE);
            }

            fclose(data);
        }
        else
        {
            perror(opts->range);
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        // Establish the grid range in which the search will be performed
        // with the use of the M matrix from grid.bin
        gridr(sett->M,
              s_range->spndr,
              s_range->nr,
              s_range->mr,
              sett->oms,
              sett->Smax);

        if (strlen(opts->getrange))
        {
            FILE *data;
            if ((data = fopen(opts->getrange, "w")) != NULL)
            {
                fprintf(data, "%d %d\n%d %d\n%d %d\n%d %d\n",
                        s_range->spndr[0], s_range->spndr[1],
                        s_range->nr[0], s_range->nr[1],
                        s_range->mr[0], s_range->mr[1],
                        s_range->pmr[0], s_range->pmr[1]);

                printf("Wrote input data grid ranges to %s\n", opts->getrange);
                fclose(data);
                //exit(EXIT_SUCCESS);
            }
            else
            {
                printf("Can't open %s file for writing\n", opts->getrange);
                exit(EXIT_FAILURE);
            }
        }
    }

    printf("set_search_range() - the grid ranges are maximally this:\n");
    printf("(spndr, nr, mr, pmr pairs): %d %d %d %d %d %d %d %d\n",
           s_range->spndr[0], s_range->spndr[1], s_range->nr[0], s_range->nr[1],
           s_range->mr[0], s_range->mr[1], s_range->pmr[0], s_range->pmr[1]);

    printf("Smin: %le, -Smax: %le\n", sett->Smin, sett->Smax);

} // end of set search range

/// <summary>Sets up BLAS plans.</summary>
///
void init_blas(Search_settings* sett,
               OpenCL_handles* cl_handles,
               BLAS_handles* blas_handles)
{
    blas_handles->BLAS_err = clblasSetup();
    checkErr(blas_handles->BLAS_err, "clblasSetup()");
}

/// <summary>Sets up FFT plans.</summary>
///
void plan_fft(Search_settings* sett,
              OpenCL_handles* cl_handles,
              FFT_plans* plans,
              FFT_arrays* fft_arr)
{
    cl_int CL_err = CL_SUCCESS;

    fft_arr->arr_len = (sett->fftpad*sett->nfft > sett->Ninterp ?
        sett->fftpad*sett->nfft :
        sett->Ninterp);

    fft_arr->xa_d = clCreateBuffer(cl_handles->ctx,
                                   CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR,
                                   fft_arr->arr_len * sizeof(complex_devt),
                                   NULL,
                                   &CL_err);
    checkErr(CL_err, "clCreateBuffer(fft_arr->xa_d)");

    fft_arr->xb_d = clCreateBuffer(cl_handles->ctx,
                                   CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR,
                                   fft_arr->arr_len * sizeof(complex_devt),
                                   NULL,
                                   &CL_err);
    checkErr(CL_err, "clCreateBuffer(fft_arr->xb_d)");

    CL_err = clfftSetPlanPrecision(plans->plan, CLFFT_TRANSFORM_PRECISION);
    checkErr(CL_err, "clfftSetPlanPrecision(CLFFT_SINGLE)");
    CL_err = clfftSetLayout(plans->plan, CLFFT_TRANSFORM_LAYOUT, CLFFT_TRANSFORM_LAYOUT);
    checkErr(CL_err, "clfftSetLayout(CLFFT_COMPLEX_INTERLEAVED, CLFFT_COMPLEX_INTERLEAVED)");
    CL_err = clfftSetResultLocation(plans->plan, CLFFT_INPLACE);
    checkErr(CL_err, "clfftSetResultLocation(CLFFT_INPLACE)");

    CL_err = clfftBakePlan(plans->plan,
                           cl_handles->dev_count,
                           cl_handles->exec_queues,
                           NULL,
                           NULL);
    checkErr(CL_err, "clfftBakePlan()");

} // plan_fft


/* Checkpointing */

void read_checkpoints(Command_line_opts *opts, 
		      Search_range *s_range, 
		      int *FNum) {

  if(opts->checkp_flag) {

    // filename of checkpoint state file, depending on the hemisphere
    if(opts->hemi)
      sprintf(opts->qname, "state_%03d_%04d%s_%d.dat",  
	      opts->ident, opts->band, opts->label, opts->hemi);
    else
      sprintf(opts->qname, "state_%03d_%04d%s.dat", 
	      opts->ident, opts->band, opts->label);

    FILE *state;
    if((state = fopen(opts->qname, "r")) != NULL) {

      // Scan the state file to get last recorded parameters
      if((fscanf(state, "%d %d %d %d %d", &s_range->pst, &s_range->mst,
		 &s_range->nst, &s_range->sst, FNum)) == EOF) {

	// This means that state file is empty (=end of the calculations)
	fprintf (stderr, "State file empty: nothing to do...\n");
	fclose (state);
	return;

      }

      fclose (state);

      // No state file - start from the beginning
    } else {
      s_range->pst = s_range->pmr[0];
      s_range->mst = s_range->mr[0];
      s_range->nst = s_range->nr[0];
      s_range->sst = s_range->spndr[0];
      *FNum = 0;
    } // if state

  } else {
    s_range->pst = s_range->pmr[0];
    s_range->mst = s_range->mr[0];
    s_range->nst = s_range->nr[0];
    s_range->sst = s_range->spndr[0];
    *FNum = 0;
  } // if checkp_flag

} // end reading checkpoints


   /* Cleanup & memory free */

void cleanup(
	     Search_settings *sett,
	     Command_line_opts *opts,
	     Search_range *s_range,
	     FFT_plans *plans,
	     FFT_arrays *fft_arr,
	     Aux_arrays *aux,
	     double *F_d) {

  int i; 
  
  for(i=0; i<sett->nifo; i++) {
    //CudaSafeCall( cudaFreeHost(ifo[i].sig.xDat) );
    //CudaSafeCall( cudaFreeHost(ifo[i].sig.DetSSB) );
    //
    //CudaSafeCall( cudaFree(ifo[i].sig.xDatma_d) );
    //CudaSafeCall( cudaFree(ifo[i].sig.xDatmb_d) );
    //
    //CudaSafeCall( cudaFree(ifo[i].sig.aa_d) );
    //CudaSafeCall( cudaFree(ifo[i].sig.bb_d) );
    //
    //CudaSafeCall( cudaFree(ifo[i].sig.shft_d) );
    //CudaSafeCall( cudaFree(ifo[i].sig.shftf_d) );
  } 

  //CudaSafeCall( cudaFree(aux->cosmodf_d) );
  //CudaSafeCall( cudaFree(aux->sinmodf_d) );
  //CudaSafeCall( cudaFree(aux->t2_d) );
  //
  //CudaSafeCall( cudaFree(F_d) );
  //
  //CudaSafeCall( cudaFree(fft_arr->xa_d) );

  free(sett->M);

  //cufftDestroy(plans->plan);
  //cufftDestroy(plans->pl_int);
  //cufftDestroy(plans->pl_inv);

} // end of cleanup & memory free 



 /* Command line options handling: coincidences  */ 

void handle_opts_coinc(Search_settings *sett,
                       Command_line_opts_coinc *opts,
                       int argc, 
                       char* argv[])
{

  opts->wd=NULL;

  strcpy (opts->prefix, TOSTR(PREFIX));
  strcpy (opts->dtaprefix, TOSTR(DTAPREFIX));

  // Default initial value of the data sampling time 
  sett->dt = 0.5;

  opts->help_flag=0;
  static int help_flag=0;

  // Default value of the minimal number of coincidences 
  opts->mincoin=3; 

  // Default value of the narrow-down parameter 
  opts->narrowdown=0.5; 

  // Default value of the cell shift: 0000 (no shifts)
  opts->shift=0;

  // Default value of the cell scaling: 1111 (no scaling)
  opts->scale=1111;

  // Default signal-to-noise threshold cutoff
  opts->snrcutoff=6;

  // Reading arguments 

  while (1) {
    static struct option long_options[] = {
      {"help", no_argument, &help_flag, 1},
      // Cell shifts  
      {"shift", required_argument, 0, 's'},
      // Cell scaling 
      {"scale", required_argument, 0, 'z'},
      // Reference frame number 
      {"refr", required_argument, 0, 'r'},
      // output directory
      {"output", required_argument, 0, 'o'},
      // input data directory
      {"data", required_argument, 0, 'd'},
      // fpo value
      {"fpo", required_argument, 0, 'p'},
      // data sampling time 
      {"dt", required_argument, 0, 't'},
      // triggers' name prefactor 
      {"trigname", required_argument, 0, 'e'},
      // Location of the reference grid.bin and starting_date files  
      {"refloc", required_argument, 0, 'g'},
      // Minimal number of coincidences recorded in the output  
      {"mincoin", required_argument, 0, 'm'},
      // Narrow down the frequency band (+- the center of band) 
      {"narrowdown", required_argument, 0, 'n'},
      // Signal-to-noise threshold cutoff  
      {"snrcutoff", required_argument, 0, 'c'},
      {0, 0, 0, 0}
    };

    if (help_flag) {

      printf("polgraw-allsky periodic GWs: search for concidences among candidates\n");
      printf("Usage: ./coincidences -[switch1] <value1> -[switch2] <value2> ...\n") ;
      printf("Switches are:\n\n");
      printf("-data         Data directory (default is ./candidates)\n");
      printf("-output       Output directory (default is ./coinc-results)\n");
      printf("-shift        Cell shifts in fsda directions (4 digit number, e.g. 0101, default 0000)\n");
      printf("-scale        Cell scaling in fsda directions (4 digit number, e.g. 4824, default 1111)\n");
      printf("-refr         Reference frame number\n");
      printf("-fpo          Reference band frequency fpo value\n");
      printf("-dt           Data sampling time dt (default value: 0.5)\n");
      printf("-trigname     Part of triggers' name (for identifying files)\n");
      printf("-refloc       Location of the reference grid.bin and starting_date files\n");
      printf("-mincoin      Minimal number of coincidences recorded\n");
      printf("-narrowdown   Narrow-down the frequency band (range [0, 0.5] +- around center)\n");
      printf("-snrcutoff    Signal-to-noise threshold cutoff (default value: 6)\n\n");

      printf("Also:\n\n");
      printf("--help		This help\n");

      exit (0);
    }

    int option_index = 0;
    int c = getopt_long_only (argc, argv, "p:o:d:s:z:r:t:e:g:m:n:c:", long_options, &option_index);
    if (c == -1)
      break;

    switch (c) {
    case 'p':
      sett->fpo = atof(optarg);
      break;
    case 's': // Cell shifts 
      opts->shift = atof(optarg);
      break;
    case 'z': // Cell scaling   
      opts->scale = atoi(optarg);
      break;
    case 'r':
      opts->refr = atoi(optarg);
      break;
    case 'o':
      strcpy(opts->prefix, optarg);
      break;
    case 'd':
      strcpy(opts->dtaprefix, optarg);
      break;
    case 't':
      sett->dt = atof(optarg);
      break;
    case 'e':
      strcpy(opts->trigname, optarg);
      break;
    case 'g':
      strcpy(opts->refloc, optarg);
      break;
    case 'm':
      opts->mincoin = atoi(optarg);
      break;
    case 'n':
      opts->narrowdown = atof(optarg);
      break;
    case 'c':
      opts->snrcutoff = atof(optarg);
      break;
    case '?':
      break;
    default:
      break ;
    } /* switch c */
  } /* while 1 */

  // Putting the parameter in triggers' frequency range [0, pi] 
  opts->narrowdown *= M_PI; 

  printf("#mb add info at the beginning...\n"); 
  printf("The SNR threshold cutoff is %.12f, ", opts->snrcutoff); 
  printf("corresponding to F-statistic value of %.12f\n", 
    pow(opts->snrcutoff, 2)/2. + 2); 

} // end of command line options handling: coincidences  



#if 0
/* Manage grid matrix (read from grid.bin, find eigenvalues 
 * and eigenvectors) and reference GPS time from starting_time
 * (expected to be in the same directory)    
 */ 

void manage_grid_matrix(
			Search_settings *sett, 
			Command_line_opts_coinc *opts) {

  sett->M = (double *)calloc(16, sizeof (double));

  FILE *data;
  char filename[512];
  sprintf (filename, "%s/grid.bin", opts->refloc);

  if ((data=fopen (filename, "r")) != NULL) {

    printf("Reading the reference grid.bin at %s\n", opts->refloc);

    fread ((void *)&sett->fftpad, sizeof (int), 1, data);

    printf("fftpad from the grid file: %d\n", sett->fftpad); 

    fread ((void *)sett->M, sizeof(double), 16, data);
    // We actually need the second (Fisher) matrix from grid.bin, 
    // hence the second fread: 
    fread ((void *)sett->M, sizeof(double), 16, data);
    fclose (data);
  } else {
    perror (filename);
    exit(EXIT_FAILURE);
  }

  /* //#mb seems not needed at the moment 
     sprintf (filename, "%s/starting_date", opts->refloc);

     if ((data=fopen (filename, "r")) != NULL) {
     fscanf(data, "%le", &opts->refgps);

     printf("Reading the reference starting_date file at %s The GPS time is %12f\n", opts->refloc, opts->refgps);
     fclose (data);
     } else {
     perror (filename);
     exit(EXIT_FAILURE);
     }
  */ 

  // Calculating the eigenvectors and eigenvalues 
  gsl_matrix_view m = gsl_matrix_view_array(sett->M, 4, 4);

  gsl_vector *eval = gsl_vector_alloc(4);
  gsl_matrix *evec = gsl_matrix_alloc(4, 4);

  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(4); 
  gsl_eigen_symmv(&m.matrix, eval, evec, w);
  gsl_eigen_symmv_free(w);

  double eigval[4], eigvec[4][4]; 
  // Saving the results to the settings struct sett->vedva[][]
  { int i, j;
    for(i=0; i<4; i++) { 
      eigval[i] = gsl_vector_get(eval, i); 
      gsl_vector_view evec_i = gsl_matrix_column(evec, i);

      for(j=0; j<4; j++)   
	eigvec[j][i] = gsl_vector_get(&evec_i.vector, j);               
    } 

    // This is an auxiliary matrix composed of the eigenvector 
    // columns multiplied by a matrix with sqrt(eigenvalues) on diagonal  
    for(i=0; i<4; i++) { 
      for(j=0; j<4; j++) { 
	sett->vedva[i][j]  = eigvec[i][j]*sqrt(eigval[j]); 
	//        printf("%.12le ", sett->vedva[i][j]); 
      } 
      //      printf("\n"); 
    }

  } 

  /* 
  //#mb matrix generated in matlab, for tests 
  double _tmp[4][4] = { 
  {-2.8622034614137332e-001, -3.7566564762376159e-002, -4.4001551065376701e-012, -3.4516253934827171e-012}, 
  {-2.9591999145463371e-001, 3.6335210834374479e-002, 8.1252443441098394e-014, -6.8170555119669981e-014}, 
  {1.5497867603229576e-005, 1.9167007413107127e-006, 1.0599051611325639e-008, -5.0379548388381567e-008}, 
  {2.4410008440913992e-005, 3.2886518554938671e-006, -5.7338464150027107e-008, -9.3126913365595100e-009},
  };

  { int i,j; 
  for(i=0; i<4; i++) 
  for(j=0; j<4; j++) 
  sett->vedva[i][j]  = _tmp[i][j]; 
  }

  printf("\n"); 

  { int i, j; 
  for(i=0; i<4; i++) { 
  for(j=0; j<4; j++) {
  printf("%.12le ", sett->vedva[i][j]);
  }
  printf("\n"); 
  } 

  } 
  */ 

  gsl_vector_free (eval);
  gsl_matrix_free (evec);

} // end of manage grid matrix  

#endif

/*---------------------------------------------------------------------------*/

/*
  Initialize CUDA: cuinit
  - sets cuda device to (in priority order): cdev, 0 
  - returns: device id or -1 on error
*/
int cuinit(int cdev)
{
  //int dev, deviceCount = 0;
  //cudaDeviceProp deviceProp;
  
//  if (cudaGetDeviceCount(&deviceCount) != cudaSuccess) {
//    printf("ERROR: cudaGetDeviceCount FAILED CUDA Driver and Runtime version may be mismatched.\n");
//    return(-1);
//  }
//  if (deviceCount == 0) {
//    printf("ERROR: There is no device supporting CUDA\n");
//    return(-1);
//  }
//  if (cdev < 0 && cdev >= deviceCount) {
//    printf("\nWARNING: Device %d is not available! Trying device 0\n", cdev);
//    cdev = 0;
//  }
//
//  printf("__________________________________CUDA devices___________________________________\n");
//  printf("Set | ID |        Name        |   Gmem(B)   | Smem(B) | Cmem(B) | C.Cap. | Thr/bl |\n");
//  
//  for (dev = 0; dev < deviceCount; ++dev) {
//    cudaGetDeviceProperties(&deviceProp, dev);
//    if (deviceProp.major == 9999 && deviceProp.minor == 9999) {
//      printf("- | %1d | %16s | Error | Error | Error | Error | Error |\n", dev, deviceProp.name );
//      if ( dev==cdev ) {
//	printf("ERROR: Can't set device %d\n", cdev);
//	return(-1);
//      }
//    }
//    if (dev==cdev) {
//      printf(" *  |");
//      cudaSetDevice(cdev);
//    } else {
//      printf("    |");
//    }
//    printf(" %1d  | %18.18s | %11Zu | %7Zu | %7Zu |   %d.%d  | %6d |\n", 
//	   dev, deviceProp.name, deviceProp.totalGlobalMem, deviceProp.sharedMemPerBlock, 
//	   deviceProp.totalConstMem, deviceProp.major, deviceProp.minor, deviceProp.maxThreadsPerBlock );
//  }
//  printf("---------------------------------------------------------------------------------\n");
//  
//  /* enable mapped memory */
//  cudaSetDeviceFlags(cudaDeviceMapHost);
//
//  /* force initialization */
//  cudaThreadSynchronize();
  return(cdev);
}
