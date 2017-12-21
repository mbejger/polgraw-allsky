# Grid generation 

## Implementation of optimal grid of templates 

The implementation of the [banks of templates for all-sky narrow-band searches of gravitational waves from spinning neutron stars](http://iopscience.iop.org/article/10.1088/0264-9381/32/14/145014) paper can be found [here](https://github.com/mbejger/polgraw-allsky/tree/master/gridgen).

### Prerequisites 

C++ compiler.

### Example   

Using the [test Gaussian data](https://polgraw.camk.edu.pl/H1L1_2d_0.25.tar.gz):  
```
% ./gridgen -m 0.5 -p dfg -d ../testdata/2d_0.25/001/H1/ -n 17 
```
where   

* `-m` is the minimal-match parameter
* `-d` is the directory with the ephemerids data 
* `-n` is the $\log_2$ of the number of data points 
* `-p` denotes the full output (density of covering, Fisher and grid matrices)

The output is 

```bash 
grid.bin will be saved in ../testdata/2d_0.25/001/H1/grid.bin
Covariance, Density of covering
Density of covering: 0.25 1.76768559133
fftpad: 1
Normalized grid matrix:
4.793689962143e-05      0.000000000000e+00      0.000000000000e+00      0.000000000000e+00
-1.742438091399e-04     2.161320668089e-09      0.000000000000e+00      0.000000000000e+00
-2.479053689156e-02     -1.684378936606e-09     2.558634258126e+02      0.000000000000e+00
9.606645586405e-03      3.102095432504e-09      -1.321384079961e+02     1.941709922254e+02

Fisher matrix:
8.333333333333e-02      8.333333333333e-02      8.125315899793e-06      1.296986289673e-06
8.333333333333e-02      8.888888888889e-02      8.127119310149e-06      1.288789877048e-06
8.125315899793e-06      8.127119310149e-06      7.922523031349e-10      1.264584758632e-10
1.296986289673e-06      1.288789877048e-06      1.264584758632e-10      2.020158195389e-11

Grid matrix:
4.130435018981e+00      0.000000000000e+00      0.000000000000e+00      0.000000000000e+00
-1.501354357073e+01     1.604615232547e+01      0.000000000000e+00      0.000000000000e+00
-2.136051820725e+03     -1.250522487924e+01     2.204621622171e+07      0.000000000000e+00
8.277470103070e+02      2.303068516072e+01      -1.138557378658e+07     1.673054937411e+07
```

### Full description of options  

```bash 
DESCRIPTION
         GridsGenerator (GG) is designated to be used in all-sky narrow-band
         searches of continuous gravitational waves. Program allow to:
         - generate efficient grid(s) for chosen initial time of observation (1).
         - generate reduced Fisher matrix for chosen initial time of observation (1),
         - generate density of covering (2).
         (1) To get result, ephemeris must be provided.
         (2) Result can be get without ephemeris (for more information see flags: -nd
         (--ndata)).

FLAGS
   Flags with (optional) argument(s):
         -c or --covariance <min> <max> <step>
                 Covariance. Flag -c is required (even without argument) to get result.
                 <min> - minimum value of covariance but not less than 0;
                 default set: 0.75.
                 <max> - optional maximum value of covariance but less than 1.
                 <step> - optional step value of covariance;
                 default set: 0.01.

                 # Calculation are preform only in two cases:
                 # 1. No flag are provided. Sets are read from file 'gg.ini'.
                 # Result(s) is(are) printed to file(s) - see configuration file: 'gg.ini'.
                 # 2. Flag -c (--covariance) or -m (--match) is used.
                 # Result(s) is(are) printed to tty.

         -m or --match <min> <max> <step>
                 Minimal match (MM^2 == covariance).
                 <min> - minimum value of minimal match but not less than 0;
                 default set: MM^2 = 0.75.
                 <max> - optional maximum value of minimal match but less than 1.
                 <step> - optional step value of minimal match;
                 default set: 0.1.
                 ## If flags -c (--covariance) and -m (--match) provided simultaneously,
                 ## program will read options from -c flag only.


         -d or --directory <path>
                 Directory containing ephemeris (need to contain binary files:
                 'rSSB.bin', 'rDet.bin', 'rDetSSB.bin').
                 <path> - take path to directory.
                 E.g. -d 001: directory '001' located inside folder with GridsGenerator.
                 If directory is not provided, program will try to find ephemeris
                 in directory with GridsGenerator.

         -i or --initial <min> <max> <step>
                 Initial time of observation.
                 <min> - minimum value of initial time;
                 default set: 0.5.
                 <max> - optional maximum value of initial time.
                 <step> - optional step value of minimal match;
                 if not provided step will be set on step = max - min.


         -a or --algorithm <type>
                 Algorithm type to choose.
                 <type> take option: s1, s2, a. Algorithms:
                 s1 - based on fully analytic formula,
                 s2 - partially numeric formula.
                 Accuracy for algorithm 's2' depended on -na (--nalpha) and -nr (--nradius)
                 flags.
                 a - automatically choose algorithm 's1' or 's2'. Use this argument to allow
                 GridsGenerator to decide which algorithm (for given parameters) should be
                 used to get grid with better density of covering.
                 Information about implemented algorithms can be found in article:
                 http://dx.doi.org/10.1088/0264-9381/32/14/145014
                 Default set: a.

         -n or --nfft <int>
                 Number of Fourier bins.
                 <int> take positive integer number without zero (exponent).
                 E.g. to get grid destined to work with discreet
                 Fourier transform (DFT, FFT) with length 1024 put same exponent of 2: 10
                 (1024 == 2^10).
                 Default set: 20 (1048576 = 2^20).

         -nd or --ndata <int>
                 Number of data (data length collected by detector, equal to length
                 of ephemeris).
                 <int> take positive integer number including zero.
                 If <int> is set to zero data length will be read from ephemeris*.
                 If <int> is bigger than zero program will able to obtain density of
                 covering only**.
                 Default set: 0.
                 * With this set (-nd 0 or -ndata 0) density is obtained without using
                 information stored in file 'dataS2.txt'.
                 ### With this set density of covering for algorithm 's2' is always
                 ### obtaining with maximal accuracy, depending only from flags: -na, -nr.
                 ### Density of covering for algorithm 's1' is always obtained with maximal
                 ### accuracy regardless to sets of -nd flag.
                 ** Ephemeris are not used and because of that grid(s) and reduced Fisher
                 matrix cannot be obtained. Density of covering for algorithm 's2' will be
                 be obtained in approximated way based on information*** collected (in 
                 'dataS2.txt' file) during previous runs with ephemeris. File 'dataS2.txt'
                 collect information about coverings only for algorithm 's2'.
                 Program inform which algorithm has been used only when work without
                 ephemeris (<int> bigger that zero).
                 E.g. (without ephemeris) -nd 344656:
                 Covariance, Density, Algorithm
                 0.86      1.82      s2
                 0.87      1.8494      s1
                 0.88      1.77685      s1
                 E.g. (data length taken from ephemeris) -nd 0:
                 Covariance, Density of covering
                 0.86 1.82299890516
                 0.87 1.84940429814
                 0.88 1.77685017541
                 *** Information stored in file 'dataS2.txt' are independent of ephemeris.
                 Based on this information density of covering can be obtained fast but
                 in approximated way (generally speaking more collected data allows
                 to get more accurate results).

         -na or --nalpha <int>
                 Number of loops (incrementation deep) in root finding algorithm.
                 <int> take positive integer number without zero.
                 This flag affect only on 's2' algorithm.
                 Default set: 35.

         -nr or --nradius <int>
                 Number of loops in covering radius (deep hole) finding algorithm.
                 <int> take positive integer number without zero.
                 This flag affect only on 's2' algorithm.
                 Default set: 20.

         -cv or --convert <bool>
                 Convert grid from space with hyper-sphere to hyper-ellipsoid space.
                 <bool> take argument: t (true), f (false).
                 Default set: t.

         -ch or --chop <bool>
                 Chop result to zero if value of result is smaller than |10^-12|.
                 <bool> take argument: t (true), f (false).
                 Default set: f.

         -p or --print <string>
                 Print result(s).
                 <string> take argument(s): d, f, g.
                 d - density of covering,
                 f - Fisher reduced matrix,
                 g - grid.
                 Option can be joined (order is not important), e.g. df, dg, fg, dfg.
                 Default set: g.
                 #### If argument of flag -nd is set to be bigger than zero, flag -p
                 #### can accept only 'd' argument.

   Flags with no arguments:
         -h or --help
                 Display this help.

         -ht or --helptxt
                 Print this help to text file.

         -v or --version
                 Display information about version of GridsGenerator.

         -aa or --author
                 About author(s). Display information about developer(s).

```
