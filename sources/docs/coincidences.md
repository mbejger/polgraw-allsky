# Coincidences between candidates 

In order to establish the probability of detection of a real gravitational wave, after finding the candidate signals with the [F-statistic candidate signal search](../polgraw-allsky/candidate_search/), the pipeline searches for coincidences between candidates found in different time frames. 

The coincidences code is available at [github](https://github.com/mbejger/polgraw-allsky/tree/master). To get it, run `git clone https://github.com/mbejger/polgraw-allsky.git`.

#### Prerequisites 

The code is written in standard `C`. The only dependency is [GNU Scientific Library (GSL)](http://www.gnu.org/software/gsl/), used to manipulate the Fisher matrix (calculate the eigenvectors and eigenvalues). [GNU struct dirent](http://www.gnu.org/software/libc/manual/html_node/Accessing-Directories.html#Accessing-Directories) objects are used to read the directories. 
 
### The idea behind coincidences 

After finding the candidate signals in different time frames (`search`), we want to confirm the existence of signals with the same parameters along the span of time segments. to further perform a validation search for high-coincidence, or otherwise interesting candidates (the [followup](https://github.com/mbejger/polgraw-allsky/tree/master/followup), currently under construction). To do this, the candidates from a list of trigger files (time frames) are read, and for each trigger file

* a candidate is transformed to a well-defined time moment (frequency shifted to a reference frame), 
* translated into linear (integer) coordinates, 
* duplicates within each frame are removed, 
* list of unique candidates from all the frames is created and sorted, 
* duplicates are counted - these are the coincidences. 

*TODO: describe cell shifts (16 different shifts in total: 0101, 0110, 0010 etc. in f, s, d, a directions) and scaling of cells (used to define the linear parameters for a given cell to subsequently compare the candidate values)* 

### Compilation

Run `make coincidences`; resulting binary is called `coincidences`. Modify the `Makefile` to fit your system.

#### Full list of switches 
Type 
```
% ./coincidences --help 
```
to obtain the following description: 

| Switch          | Description       |
|-----------------|:------------------|
|-data            | Data directory (default is `./candidates`)
|-output          | Output directory (default is `./coinc-results`)
|-shift           | Cell shifts in `fsda` directions (4 digit number, e.g. `0101`, default `0000`)
|-scalef          | Cell scaling in f direction (a number, e.g. 32, default 1)
|-scales          | Cell scaling in s direction (a number, e.g. 8, default 1)
|-scaled          | Cell scaling in d direction (a number, e.g. 4, default 1)
|-scalea          | Cell scaling in a direction (a number, e.g. 4, default 1)
|-refr            | Reference frame number
|-fpo             | Reference band frequency `fpo` value
|-dt              | Data sampling time dt (default value: `0.5`)
|-trigname        | Part of triggers' name (for identifying files)
|-refloc          | Location of the reference grid.bin and starting_date files
|-mincoin         | Minimal number of coincidences recorded
|-narrowdown      | Narrow-down the frequency band (range [0, 0.5] +- around center)
|-snrcutoff       | Signal-to-noise threshold cutoff (default value: 6)

Also:

|                 |             | 
|-----------------|:------------|
| --help          |This help    |


### Example 

Using the software injection added to 2-day Gaussian noise data segments (see [minimal example of the pipeline](../polgraw-allsky/pipeline_script)):
```
% for s in {0..1}{0..1}{0..1}{0..1}; do ./coincidences -data ../../search/network/src-cpu -output . -shift $s -scalef 4 -scales 4 -scaled 4 -scalea 4 -refr 4 -dt 2 -trigname 1234_2 -refloc ../../testdata/2d_0.25/004 -nod 2 -fpo 308.859375 -snrcutoff 5; done 2>> summary
```
This assumes that for the band $bbbb=1234$ and the sampling time $dt=2\ \mathrm{s}$ the band frequency $fpo=308.859375\ \mathrm{Hz}$, because 
$$
fpo = 10 + (1 - 2^{-5})\cdot bbbb\cdot \frac{1}{2dt}\ \mathrm{[Hz]}.
$$
The reference grid file `grid.bin`, corresponding to the reference frame `004` (`-refr 4`) is located at the `-refloc` location. Signal-to-noise ratio cutoff is set to 5 (`-snrcutoff 5`). Some output is directed to `stdin`. For example, the output for shift `0000` is 
```
Number of days is 2
The SNR threshold cutoff is 5.000000000000, corresponding to F-statistic value of 14.500000000000
Reading the reference grid.bin at ../../testdata/2d_0.25/004
fftpad from the grid file: 1
Settings dt: 2.000000, oms: 3881.241374
Reference frame number: 4
Cell shifts  in f, s, d, a directions: 0 0 0 0 
Cell scaling in f, s, d, a directions: 4 4 4 4 
Reading triggers... Frame 5: 1966/2040
Reading triggers... Frame 1: 2409/2483
Reading triggers... Frame 4: 2176/2384
Reading triggers... Frame 3: 2132/2247
Reading triggers... Frame 8: 2372/2408
Reading triggers... Frame 2: 2197/2249
Reading triggers... Frame 6: 2225/2305
Reading triggers... Frame 7: 2175/2226
Total number of candidates from all frames: 17652
```

The highest multiplicity coincidence for each shift is streamed to the `summary` file via the `stderr` (`2>> summary`). From this list of 16 lines, the highest multiplicity coincidence with the highest signal-to-noise is selected (`sort -gk5 -gk10 summary | tail -1`):
```
1234_2 1111 308.859375     8     5  9.95663703e-01 -1.10830358e-09 -1.12585347e-01 1.97463002e+00 1.246469e+01 5 2040 1987 1 2483 2419 4 2384 2193 3 2247 2137 8 2408 2363 2 2249 2172 6 2305 2220 7 2226 2191 6 2 8 3 5
```  
This output contains the 

* `band_hemisphere` identifier ($1234\_2$), 
* the `shift` value ($1111$), 
* the `fpo` reference band frequency ($308.859375\ \mathrm{Hz}$), 
* the number of triggers files read ($8$), 
* the multiplicity of coincidence found ($N_{coin}=5$), 
* arithmetic mean values of the frequency $\bar{f}$ (range $[0:2\pi]$, corresponding to the width of the band), mean frequency derivative $\bar{s}$ (spindown, in $Hz/s$), equatorial coordinate system sky positions $\bar{\delta}$ ($[\pi:-\pi]$) and $\bar{\alpha}$ ($[0:2\pi]$), and the mean signal-to-noise ratio, $\widetilde{\mathrm{snr}}=\sqrt{\sum_i \mathrm{snr}_i^2}$ ($5$ floating-point numbers), 
* 8 triplets of the frame number, number of all candidates in the corresponding triggers file, and the number of unique candidates after sorting to unique cells ($8\times 3$ integers), 
* frame numbers participating in the coincidence ($5$ integers). 

Coincidences larger or equal `mincoin` (default value `3`) are recorded in binary files `.coi`, separately for each shift, in the `-output` directory. Each coincidence is a set of following numbers: 
$$
N_{coin},\quad\bar{f},\quad\bar{s},\quad\bar{\delta},\quad\bar{\alpha},\quad\widetilde{\mathrm{snr}},\quad\mathrm{fr}_{1},\,\dots\,\mathrm{fr}_{N_{coin}},\quad\mathrm{p}_{1},\,\dots\,\mathrm{p}_{N_{coin}}
$$
where 

* $N_{coin}$ is the multiplicity of coincidence (written as one `unsigned short int`), 
* $\bar{f}$, $\bar{s}$, $\bar{\delta}$, $\bar{\alpha}$ and $\widetilde{\mathrm{snr}}$ are the mean parameters of the signal ($5\times$`float`),
* $\mathrm{fr}_{1},\,\dots\,\mathrm{fr}_{N_{coin}}$ are the frame numbers ($N_{coin}\times$`unsigned short int`), 
* $\mathrm{p}_{1},\,\dots\,\mathrm{p}_{N_{coin}}$ are the positions of candidate signals that took part in the coincidences, in their corresponding trigger files ($N_{coin}\times$`int`) 

in order to recover the information about the original triggers for further studies. `.coi` files can be read with the auxilary `read_coi` code: 
```
% gcc -o read_coi read_coi.c -lm 
% ./read_coi 
# num_of_coincidences    mean_val_of_pars (f, s, d, a), snr    frame_num:trigger_num_in_trigger_file
# (see http://mbejger.github.io/polgraw-allsky/coincidences for details)
5 9.956637e-01 -1.108304e-09 -1.125853e-01 1.974630e+00 1.246469e+01 6:777 2:708 8:968 3:701 5:603 
3 9.972902e-01 -1.174760e-10 -1.235594e-01 1.968819e+00 9.154494e+00 8:983 3:693 1:669 
```
where the pairs `frame-number:candidate-position-in-trigger-file` have been arranged for readability.   
