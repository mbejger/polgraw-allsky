# Coincidences between candidates 

In order to establish the probability of detection of a real gravitational wave, after finding the candidate signals with the [F-statistic candidate signal search](../polgraw-allsky/candidate_search/), the pipeline searches for coincidences between candidates found in different time frames. 

The coincidences code is available at [github](https://github.com/mbejger/polgraw-allsky/tree/master). To get it, run `git clone https://github.com/mbejger/polgraw-allsky.git`.

#### Prerequisites 

The code is written in standard `C`. The only dependency is [GNU Scientific Library (GSL)](http://www.gnu.org/software/gsl/), used to manipulate the Fisher matrix (calculate the eigenvectors and eigenvalues). [GNU struct dirent](http://www.gnu.org/software/libc/manual/html_node/Accessing-Directories.html#Accessing-Directories) objects are used to read the directories. 
 
#### The idea behind coincidences 

After finding the candidate signals in different time frames (`search`), we want to confirm the existence of signals with the same parameters along the span of time segments. to further perform a validation search for high-coincidence, or otherwise interesting candidates (the [followup](https://github.com/mbejger/polgraw-allsky/tree/master/followup), currently under construction). To do this, the candidates from a list of trigger files (time frames) are read, and for each trigger file

* a candidate is transformed to a well-defined time moment (frequency shifted to a reference frame), 
* translated into linear (integer) coordinates, 
* duplicates within each frame are removed, 
* list of unique candidates from all the frames is created and sorted, 
* duplicates are counted - these are the coincidences. 
```
[-x----][-x----][-x----][-x----][-----y][-x----][-x----][---z--]
6 coincidences of x in 8 data segments, not bad...
```
*TODO: describe cell shifts (16 different shifts in total: 0101, 0110, 0010 etc. in f, s, d, a directions) and scaling of cells (used to define the linear parameters for a given cell to subsequently compare the candidate values)* 

### 7. Compilation

Run `make coincidences`; resulting binary is called `coincidences`. Modify the `Makefile` to fit your system.

### 8. How to run the program?

Sample call of `coincidences` is as follows:
```
> ./coincidences -data . -output 0666-coincidences -mincoin 3 -snrcutoff 5 -trigname triggers_ -refloc 6d_xdat_0.25/010 -refr 010 -fpo 171.296875 -shift 0110 -scale 4444 -dt 2 2>> summary
```
This assumes that the 0666 band frequency `fpo` is 171.296875, since we define 
```
fpo = 10. + (1-2^(-5))*band*(0.5/dt) 
```
The reference grid, corresponding to the reference frame 010 is located at `-refloc` location. Some output is directed to `stdin`. The highest-coincidence is outputed to `stderr`, redirected to a `summary` file (`2>> summary`). In principle one has to run the code 16 times for all the $2^4$ shifts `0000--1111`.   

#### 8.1. Full list of switches 
Type 
```
> ./coincidences --help 
```
to obtain the following description: 

| Switch          | Description       |
|-----------------|:------------------|
|-data            | Data directory (default is `./candidates`)
|-output          | Output directory (default is `./coinc-results`)
|-shift           | Cell shifts in `fsda` directions (4 digit number, e.g. `0101`, default `0000`)
|-scale           | Cell scaling in `fsda` directions (4 digit number, e.g. `4824`, default `1111`)
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

### 9. Output
Output to the screen (`stdout`) in the case of 
```
> ./coincidences -data . -output 0666-coincidences -mincoin 3 -snrcutoff 5 -trigname triggers_ -refloc 6d_xdat_0.25/010 -refr 010 -fpo 171.296875 -shift 0110 -scale 4444 -dt 2 2>> summary
```
is
```
The SNR threshold cutoff is 5.000000000000, corresponding to F-statistic value of 14.500000000000
Reading the reference grid.bin at 6d_xdat_0.25/010
fftpad from the grid file: 1
Settings dt: 2.000000, oms: 2152.580016
Reference frame number: 10
Cell shifts  in f, s, d, a directions: 0 1 1 0 
Cell scaling in f, s, d, a directions: 4 4 4 4 
Reading triggers_005_0666_2.bin... Frame 5: 941/1879
Reading triggers_019_0666_2.bin... Frame 19: 1320/1794
Reading triggers_003_0666_2.bin... Frame 3: 1512/1920
Reading triggers_008_0666_2.bin... Frame 8: 3090/6849
Reading triggers_002_0666_2.bin... Frame 2: 1014/1173
Reading triggers_014_0666_2.bin... Frame 14: 779/1542
Reading triggers_006_0666_2.bin... Frame 6: 2884/8822
Reading triggers_011_0666_2.bin... Frame 11: 918/2896
Total number of candidates from all frames: 12458
```
The highest-coincidence is streamed to a `summary` file via `stderr`:  
```
> cat summary 
triggers_ 0110 171.296875 8 8 1.05795986e+00 -7.69412393e-09 -7.42581733e-01 5.39896106e+00 2.417267e+01 5 1879 941 19 1794 1320 3 1920 1512 8 6849 3090 2 1173 1014 14 1542 779 6 8822 2884 11 2896 918
```
It contains the `trigname` identifier, the `shift` value, the `fpo` band frequency, the number of files read (8), and the highest coincidence found (8). Next 5 numbers are arithmetic mean values of the frequency $\,f$ (in radians), frequency derivative $\,\dot{f}$ (spindown), sky positions $\delta$ and $\alpha$, and the mean signal-to-noise ratio, $\widetilde{\mathrm{snr}}=\sqrt{\sum_i \mathrm{snr}_i^2}$. The following integers are grouped in threes and denote the triggers/time frame number, number of all candidates in that triggers' file, and the number of unique candidates. 

Coincidences above `mincoin` are recorded in a binary file `.coi`, separately for each shift, in the `-output` directory. Each coincidence is a set of following numbers: 
$$
N,\quad\bar{f},\quad\bar{s},\quad\bar{d},\quad\bar{a},\quad\widetilde{\mathrm{snr}},\quad\mathrm{fr}_{1},\,\dots\,\mathrm{fr}_{N},\quad\mathrm{p}_{1},\,\dots\,\mathrm{p}_{N}
$$
where 

* $N$ is the size of coincidence (written as one `unsigned short int`), 
* $\bar{f}$, $\bar{s}$, $\bar{d}$, $\bar{a}$ and $\widetilde{\mathrm{snr}}$ are the mean parameters of the signal ($5\times$`float`),
* $\mathrm{fr}_{1},\,\dots\,\mathrm{fr}_{N}$ are the frame numbers ($N\times$`unsigned short int`), 
* $\mathrm{p}_{1},\,\dots\,\mathrm{p}_{N}$ are the positions of candidate signals that took part in the coincidences, in their corresponding trigger files ($N\times$`int`).

