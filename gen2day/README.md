gen2day
==============

Time series input data generator. 

This code uses the Short Fourier Transform database to generate input data in time domain for the search code. 
Additionally, linked with the LAL in can generate the emphemerids for the detectors. 

## Example call for O3

`genseg 0044_6d.g2d`

## Major changes in O3

* major code reorganization
* outliers removed from the final time series (before - from smaller chunks, out files)
* use only data with science flag
* when large outliers are found at the edges of science region - narrow it down
* Tukey window is applied to science segments (window taper size = 600 s for O3)
