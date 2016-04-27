# Input data generation

## `extract_band`

This program converts Short Fourier Transformation series to time series. 
Written by Pia Astone (INFN, Physics Department of University of Rome "La Sapienza").

### Prerequisites 

C compiller. Uses standard C libraries, `libm`. Links to the PSS library (created by Pia Astone).

### How to run it? 

```
> extract_band < input_file
```
where `input_file` is an ASCII file contaning the following rows:  

* Maximal number of SFT
* The name of the output file
* The list of SFT files
* The frequency band in Hz
* The width of frequency band in Hz

for example: 
```
100000
J0034+1612_2010-10-10.out
J0034+1612_2010-10-10.list
718.2480
1
```

### Output

Example output: 
```
% Beginning freq- Band- Samples in one stretch- Subsampling factor- inter (overlapping, 2 if data were overlapped)- Frequency step- Scaling factor- ***The data are real and imag of the FFT
% 908.152344 0.250000 256 8192.000000 2  0.0009766 1.000000e-20
% FFT number in the file; Beginning mjd days; Gps s; Gps ns;
% 100 55099.5879745370 937922816 0
 4.59662571e+02  2.27630825e+01
-3.50387007e+02 -2.20005558e+02
 3.57587904e+02  1.01217077e+02
 1.74400486e+02  2.62086552e+02
 2.21804800e+02 -5.20278366e+02
-3.87826732e+02 -1.55758978e+02
```


## `gen2day` description 

... 

