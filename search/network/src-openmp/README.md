# Search code

## Example call for O3
`gwsearch-omp-haswell-dev7 -o . -nod 6 -data /net/archive/groups/plggpolgraw/O3/allsky_o3_c01 -ident 001 -b 0340 -v 0.1 --threshold 14.5 -dt 2 --narrowdown 0.45 --nocheckpoint --whitenoise`

band (-b) nad overlap (-v) are required now; base frequency fpo = 10 + (1 - overlap) * band / (2 dt)

## Major changes in O3

* Use grids for a network of detectors
* Bug fix: modulation factors are set to 0 when data are == 0
* Use `--whitenoise` (do not normalize F-statistics)
* Observationaly motivated spindown range
* Bug fix: line widths of combs calculated properly

