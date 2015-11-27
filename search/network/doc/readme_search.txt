search - narrow band all-sky search for periodic signals of GW

1. Prerequisites

   Gnu  Scientific Library  (GSL)  and FFTW  library  (version 3.0  or
   later) are needed to compile and run the search program.

2. How to compile the program?

   Run  "make  search"  or  "make".  The resulting  binary  is  called
   "search". Modify CFLAGS and LDFLAGS in the Makefile if needed.

3. How to run the program?

   "./search <-d  data_dir> <-o output_dir>  -i frame -b  band", where
   data_dir is the  base directory of input data  files and output_dir
   is a directory  to write output files. Both -d  and -o are optional
   and set to "." and "./candidates" by default. "frame" is the number
   of 2-day  long time frame to  be analyzed. "band" is  the number of
   frequency  band  (see  "Structure  of  VSR1  data"  below  for  the
   discussion). The following  command should  work with a testfile
   provided (for example, from a directory containing the data):

   ./search -d ./data -o ./candidates -i 42 -b 271 --whitenoise
 
   The program will proceed assuming that data directory ./data/042 
   contain subdirectories (H1, L1, and/or V1) with input time series 
   (xdat_271_042.bin in each subdirectory, corresponding to each 
   detector, as well as corresponding ephemerids DetSSB.bin. 
   The grid.bin file should be located in ./data/042/grid.bin 

   Additional switches:

  -d  Data directory (default is .)
  -o  Output directory (default is ./candidates)
  -i  Frame number
  -b  Band number
  -l  Custom label for the input and output files
  -r  File with grid range or pulsar position
  -c  Change to directory <dir>
  -t  Threshold for the F-statistic (default is 20)
  -h  Hemisphere (default is 0 - does both)
  -p  fpo (starting frequency) value
  -x  Add signal with parameters from <file>

  Also:

  --whitenoise  white Gaussian noise assumed
  --nospindown  spindowns neglected
  --nocheckpoint  state file won't be created (no checkpointing)
  --help    This help

	The range file may either contain 8 integers that denote the ranges 
	on the spindown and sky positions grid and the range of hemispheres 
	to search in: spndr1 spndr2 nr1 nr2 mr1 mr2 hemi1 hemi2, for example
	131 137
	8 14
	41 47
	2 2
	
	The same switch -r may also be used to feed the code with an
	astronomical position and parameters of the object - these
	parameters will be translated to grid positions. In this case 
	the order of parameters is: 
	EPOCH (in MJD)
	alpha
	delta (both in radians)
	FO
	F1
	F2	(spin frequency, first and second derivative of a star)
	gsize (integer, how large the grid should be). 

4. Structure of VSR1 data

   The  VSR1 data  is divided  into  68 time  frames, each  of them  2
   sideral  days  long, rounded  to  half  second (172328.0  seconds).
   Beginning time of the  first frame is MJD=54239.00 (2007/05/18, UTC
   12:00).   Frames are  labelled with  two-digit  consecutive number,
   from 01 to  68. Frame label is the name  of a subfolder, containing
   all narrowband sequencies corresponding to that frame, e.g. "./42".

   Beginning MJD of  each time frame is saved  in the nn/starting_date
   file, for example:

   % cat 42/starting_date
   54320.7760185185

   The database contains  929 narrow (1 Hz) bands,  covering the range
   100 - 1000 Hz.  Neighboring bands overlap by 0.03125 Hz.  Beginning
   frequency f_{po} of each band can be computed from

   fpo = 100.0 + df * bbb [Hz],

   where df = 1-2^{-5} = 0.96875 Hz and bbb is the band number, from 0
   to 928.

5. Input data files

   There are 3 data files, stored in data_dir/nn (nn is the time frame
   number):

   * DetSSB.bin  -  location  of  the  detector  w.r.t.  the  SSB,  in
     cartesian coordinates,  sampled at half second.  Last two records
     of  this file  are the  angle phir,  determining the  position of
     Earth in  its diurnal motion,  and the obliquity of  the ecliptic
     (epsm), both calculated at the first sample of the data.

   * grid.bin  - generator  matrix of  an optimal  grid  of templates,
     described in readme_gridopt.txt.

   * xdat_nn_fff.bin - time-domain  narrow-band data sequence, sampled
     at  half second.   nn is  the number  of time  frame, fff  is the
     number of frequency band.

6. Output files

   Output  files,   containing  trigger  events   above  an  arbitrary
   threshold, are written to  the output_dir directory.  There are two
   output files  per every input  data sequence: triggers_nn_fff_1.bin
   and  triggers_nn_fff_2.bin,  where  1  and  2  corresponds  to  the
   northern  and  southern  ecliptic  hemisphere.   Output  files  are
   written  in  a  binary  form.   Every trigger  event  occupies  40
   consecutive bytes (5 double numbers), with the following meaning:

   record no.   meaning
   -----------------------------------------------------
   1		frequency [radians, between 0 and \pi] above f_{po}
   2		spindown [Hz s^-1]
   3		declination [radians]
   4		right ascension [radians]
   5		signal to noise ratio

7. Temporary files

   There are several temporary files created by the search program in
   the working directory:

   * wisdom-hhh.dat  - "wisdom"  file  created by  fftw.   hhh is  the
     hostname, determined by a call to gethostname().

   * state_nn_fff.dat  - checkpoint  file.  the  search can  be safely
     restarted, calculations will continue  from the position saved to
     this file.  After successful termination, checkpoint file is left
     empty.

8. The algorithm

   For an overview of the search algorithm, consult the Flowchart.pdf
   file. 

9. Additional documentation 
 
   Doxygen documentation may be generated by typing "make doc" from 
   the src directory - it will reside in doc/html. 
