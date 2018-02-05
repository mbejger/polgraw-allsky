# Followup 

Currently work in progress; for the present state of the code see [here](https://github.com/mbejger/polgraw-allsky/tree/master/followup). 

## Introduction

On the last stage we can precisely estimate signal parameters. At this stage we focus only on a very narrow space around a candidate. Our goal is to find the maximum of the F-statistic and the corresponding signal parameters. To increase the efficiency of the code, we use special libraries ([GNU Scientific Library (GSL)](http://www.gnu.org/software/gsl/), [YEPPP!](http://www.yeppp.info)) and parallelisation ([OpenMP](http://www.openmp.org/)). We already implemented few algorithms to do this: Simplex algorithm, Mesh Adaptive Direct Search - MADS and modified MADS (called 'inverted MADS').

## Compilation
To compile a code, go to the `followup/src` and run: `make followup`

## Switches
Run `./followup --help`

| Switch      | Description       |
|-------------|:------------------| 
|-d, -data         |Data directory (default is .)
|-o, -output       |Output directory (default is ./candidates)
|-i, -ident        |Frame number
|-b, -band         |Band number
|-l, -label        |Custom label for the input and output files
|-c, -cwd          |Change to directory [dir]
|-t, -threshold    |Threshold for the F-statistic (default is 20)
|-p, -fpo          |Reference band frequency fpo value
|-s, -dt           |Data sampling time dt (default value: 0.5)
|-u, -usedet       |Use only detectors from string (default is use all available)
|-y, -nod          |Number of days
|-x, -addsig       |Add signal with parameters from [file]
|-a, -candidates   |As a starting point in followup use parameters from [file]
|-r, -refr         |Reference frame number

Also:

|             |                   | 
|-------------|:------------------| 
|--vetolines     |Veto known lines from files in data directory
|--simplex       |Direct search of maximum using Nelder-Mead (simplex) algorithm
|--mads          |Direct search of maximum using MADS algorithm
|--gauss         |Generate Gaussian noise instead of reading data. Amplitude and sigma of the noise declared in init.c
|--neigh         |Function neigh() generate area as % from initial value instead of taking it from grid.bin
|--naive         |Function naive() generate area as +/- points taking it from grid.bin and divide it into smaller grid.
|--onepoint      |Calculate Fstatistic only in one point taken from file with candidates (without generating any grid).
|--help          |This help


By default code is searching maximum of the F-statistic on the optimal, 4-dimensional grid (parameters of the grid, like e.g. minimal match, ale defined inside the code, in `followup.c`). Main idea was to search on a denser grid than in `search` code, but to focus only on a few points around candidate. Additionally, when the point (on the optimal grid) with the maximal value of the F-statistic is established, one can run one of the direct maximum search algorithm: Nelder-Mead (`--simplex` switch) or MADS (`--mads` switch), to determine parameters of the maximum more precisely. In the current version of the code, by using `--mads` switch, one will use modified, 'inverted' MADS. To skip calculations on the optimal grid and run one of the algorithms directly from the initial point, use `--onepoint` switch.

## Input data

Files required to run the code: `xdatc*`, `DetSSB.bin` and `grid.bin`. Path to the directory is taken from `-data` switch; there directories named as `iii/H1`,`iii/L1` etc. are expected, where `iii` is a frame number taken from `-ident` switch.

Initial point of the calculations (switch `-candidates`) is required. Usually the best result from coincidences (previous step of the pipeline) is taken. Parameters of the candidate need to be written into the file, as:

`frequency spin-down declination right_ascension`

## Minimal example how to run a code
Run: 

```bash
LD_LIBRARY_PATH=../../search/network/src-cpu/lib/yeppp-1.0.0/binaries/linux/x86_64 ./followup -data data_path/ -band 0666 -dt 16 -candidates path_to_candidate/candidate.txt -ident 001 -nod 6 --mads> output.txt 
```

As a result (last line in the output file) one can get parameters of the found maximum:

`frequency spin-down declination right_ascension F-statistic_value SNR`


