# Pipeline: a minimal example 

This is a demonstration of the pipeline using Gaussian noise test data with randomly-selected signal added to the data (software injection).  

### Structure of the input data 

The generic directory structure of the input data is

```bash 
001
├── grid.bin
├── H1
│   ├── DetSSB.bin
│   ├── grid.bin
│   ├── rDet.bin
│   ├── rSSB.bin
│   ├── starting_date
│   └── xdatc_001_1234.bin
└── L1
    ├── DetSSB.bin
    ├── grid.bin
    ├── rDet.bin
    ├── rSSB.bin
    ├── starting_date
    └── xdatc_001_1234.bin
```
(here for two LIGO detectors H1 and L1, and frame `001`). Test data frames $nnn=001-008$ with pure Gaussian noise 2-day time segments with sampling time equal to 2s (`xdatc_nnn_1234.bin`) are [available here](https://polgraw.camk.edu.pl/H1L1_2d_0.25.tar.gz). 

In principle, given the ephemerides (`DetSSB.bin`, `rDet.bin` and `rSSB.bin`) for each detector and frame `001-008`, one can create the grid matrices using the [gridgen](../polgraw-allsky/grid_generation) implementation (see the [code](https://github.com/mbejger/polgraw-allsky/tree/master/gridgen) for details): 

```bash 
# grid generation
cd gridgen
make
for ifo in H1 L1; do for d in $(seq -f %03g 1 8); do ./gridgen -m 0.5 -p dfg -d ../testdata/2d_0.25/${d}/${ifo}/ -n 17; done; done

# copying the H1 grid file one level up for the case of the network search 
for d in $(seq -f %03g 1 8); do cp -v ../testdata/2d_0.25/${d}/H1/grid.bin ../testdata/2d_0.25/${d}; done
```
Test Gaussian noise time series data were created as follows:  

```bash 
#!/bin/bash 

band=1234

# Gaussian data generation
cd ../search/network/scr-cpu
gcc gauss-xdat.c -o gauss-xdat -lm -lgsl -lgslcblas
# 86164: number of points in 2-day segment with 2s sampling time 
for ifo in H1 L1; do for d in $(seq -f %03g 1 8); do echo $d $ifo; ./gauss-xdat 86164 1 1 ../../../testdata/2d_0.25/${d}/${ifo}/xdatc_${d}_${band}.bin; done; done
```
Given the complete input data, this pipeline minimal example consists of 

* Adding an artificial signal to the data (random parameters of the signal generated with `sigen`), 
* Performing a search around the injection for each time segment (`gwsearch-cpu`), 
* Looking for coincidences between the candidate signals from different time frames (`coincidences`), 
* Establishing the false alarm probability of the best coincidence (`fap`).   

### 
### Generating random parameters for the signal 

Random parameters of the signal are chosen using the [sigen](https://github.com/mbejger/polgraw-allsky/blob/master/search/network/src-cpu/sigen.c) and added to the data time series with the `add_signal()` function in ([init](https://github.com/mbejger/polgraw-allsky/blob/master/search/network/src-cpu/init.c)). 

```bash 
# Create random parameters of a signal signal
cd search/network/src-cpu
make sigen
band=1234; dt=2; nod=2; ./sigen -amp 4.e-2 -band $band -dt $dt -gsize 10 -reffr 4 -nod $nod 1> sig1 
```
Signal parameters used in this example:

```bash 
% cat sig1 
amp 4.000000e-02
10
4
9.9791082090028898e-01
-1.6533871297433800e-09
-1.1821269273133420e-01
1.9839903273071489e+00
4.7717937494571394e-01
7.5715524886052021e-01
7.5154297884129850e-01
-4.7938541489358644e-01
``` 

### 
### Adding signal to the Gaussian data and searching for candidates

```bash 
band=1234; dt=2; nod=2; for d in $(seq -f %03g 1 8); do 

  LD_LIBRARY_PATH=lib/yeppp-1.0.0/binaries/linux/x86_64 ./gwsearch-cpu \
  -data ../../../testdata/2d_0.25/ \
  -ident ${d} \
  -band $band \
  -dt $dt \ 
  -nod $nod \ 
  -addsig sig1 \  
  -output . \
  -threshold 14.5 \
  --nocheckpoint \ 

done
``` 

This produces trigger files for each frame (size in bytes also listed): 
```
99320 triggers_001_1234_2.bin
89960 triggers_002_1234_2.bin
89880 triggers_003_1234_2.bin
95360 triggers_004_1234_2.bin
81600 triggers_005_1234_2.bin
92200 triggers_006_1234_2.bin
89040 triggers_007_1234_2.bin
96320 triggers_008_1234_2.bin
```

First 10 triggers from `triggers_001_1234_2.bin` are 

```bash
3.05617018e+00 -3.42376198e-08 -7.68007347e-02 2.59248668e+00 5.06667333e+00 
1.18243015e+00 -3.20762991e-08 -7.68007347e-02 2.59248668e+00 5.05528873e+00 
1.08103361e-01 -2.77536578e-08 -7.68007347e-02 2.59248668e+00 5.07085254e+00 
1.90022435e+00 -2.77536578e-08 -7.68007347e-02 2.59248668e+00 5.15191593e+00 
1.90000217e+00 -2.55923371e-08 -7.68007347e-02 2.59248668e+00 5.42638039e+00 
2.09224664e+00 -2.34310165e-08 -7.68007347e-02 2.59248668e+00 5.20879551e+00 
2.38731576e+00 -2.12696958e-08 -7.68007347e-02 2.59248668e+00 5.31983396e+00 
3.00543165e+00 -1.91083751e-08 -7.68007347e-02 2.59248668e+00 5.29454616e+00 
7.49333983e-01 -1.26244131e-08 -7.68007347e-02 2.59248668e+00 5.08724856e+00 
2.08710778e-01  3.43510887e-10 -7.68007347e-02 2.59248668e+00 5.17537018e+00 
```
### 
### Coincidences among these trigger files 

```bash 
cd ../../../coincidences/src
make
band=1234; dt=2; nod=2; fpo=$(echo $band $dt |awk '{printf("%.6f", 10 + 0.96875*$1/(2.0*$2))}'); for s in {0..1}{0..1}{0..1}{0..1}; do 

  ./coincidences \ 
  -data ../../search/network/src-cpu \ 
  -output . \ 
  -shift $s \ 
  -scalef 4 \ 
  -scales 4 \ 
  -scaled 4 \ 
  -scalea 4 \ 
  -refr 4 \ 
  -dt $dt \ 
  -trigname ${band}_2 \ 
  -refloc ../../testdata/2d_0.25/004 \ 
  -nod $nod \ 
  -fpo $fpo \ 
  -snrcutoff 5 \ 
 
  done 2>> summary

# best shift (highest multiplicity with largest snr)
sort -gk5 -gk10 summary | tail -1
```
The highest coincidence with the largest signal-to-noise ratio is  
```
1234_2 1111 308.859375     8     5  9.95663703e-01 -1.10830358e-09 -1.12585347e-01 1.97463002e+00 1.246469e+01 5 2040 1987 1 2483 2419 4 2384 2193 3 2247 2137 8 2408 2363 2 2249 2172 6 2305 2220 7 2226 2191 6 2 8 3 5
```

###
### False alarm probability 

```
make fap 
fap.sh <(sort -gk5 -gk10 summary | tail -1) <(echo $band 0.0) ../../testdata/2d_0.25/004
```
resulting in 
```
Number of days in time segments: 2
Input data: /dev/fd/63
Grid matrix data directory: ../../testdata/2d_0.25/004
Band number: 1234 (veto fraction: 0.000000)
The reference frequency fpo: 308.859375
The data sampling time dt: 2.000000
FAP threshold: 1.000000
Cell size: 4
1234 3.088594e+02 3.091094e+02 7.665713e-08 5 17682 9.956637e-01 -1.108304e-09 -1.125853e-01 1.974630e+00 1.246469e+01 2
```
The false alarm probability in this case is `7.665713e-08`. It's low enough to be an interesting outlier for a [followup](../polgraw-allsky/followup) procedure. 

