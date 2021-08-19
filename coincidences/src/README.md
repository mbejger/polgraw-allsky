# coincidences code

## Example call for O3
`./coincidences-haswell-dev7-grid -nod 6 -dt 2 -v 0.1 -mincoin 3 -snrcutoff 5.1961 -narrowdown 0.45 -refr 029 -b 0340 -scale 16,8,2,2 -refgrid /work/chuck/virgo/O3/allsky_o3_c01/029/grids/grid_029_0340_H1L1c.bin -output . -shift 0000 -infile coinc_0340_hemi1.in`

Band (-b) and overlap (-v) are required.

The coinc_0340_hemi1.in file contains a list of input trigger files with corresponding veto files.
```
/work/chuck/virgo/O3/allsky_o3_c01-results2/0340/triggers_001_0340_1.bin /work/chuck/virgo/O3/allsky_o3_c01-results2/0340/triggers_001_0340.vlines
/work/chuck/virgo/O3/allsky_o3_c01-results2/0340/triggers_002_0340_1.bin /work/chuck/virgo/O3/allsky_o3_c01-results2/0340/triggers_002_0340.vlines
/work/chuck/virgo/O3/allsky_o3_c01-results2/0340/triggers_004_0340_1.bin /work/chuck/virgo/O3/allsky_o3_c01-results2/0340/triggers_004_0340.vlines
```

Hemisphere in trigger file names can be 0,1 or *. In the last case both hemispheres are read.

## major changes in O3

* read trigger file names from file (before we used all files in the current directory)
* use grids for network of detectors (consistent with search)
* apply veto lines in coincidences (before we did this during the search step)
* coincidence cells are based on the same grid matrix as in search

# fap code

## Example call for O3

`./fap2-grid -nod 6 -band 44 -o 0.1 -vetofrac 0.0 -cellf 16 -cells 8 -celld 2 -cella 2 -threshold 15.5 -noc 3 -mcomb 1e7 -grid /work/chuck/virgo/O3/allsky_o3_c01/029/grids/grid_029_0044_H1L1c.bin -data coinc-summary-grid-th15.5-16-8-2-2.txt`

## major changes in O3

* major cleaning and simplification of the code
* approximate probabilities with median if number of combinations is too large (-mcomb parameter)
* parameter space volume consisten with coincidences (new grid)


# --- tools ---

## coi_digger.py 

Script to read the `.coi` files, calculate the False Alarm Probability for the candidates 
with a given coincidence multiplicity, extract information (also from the trigger files) 
for followup. 

`coi_digger.py` calls the `fap` code (run `make fap-many`) to compile. 

### Sample call: 

```
python coi_digger.py config.ini 0081 1
```

