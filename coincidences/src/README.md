# coincidences code

## Sample call

`./coincidences -nod 6 -dt 2 -v 0.1 -mincoin 10 -snrcutoff 4. -narrowdown 0.45 -refr 029 -b 0340 -scale 16,16,1,1 -refgrid <input data path>/029/grids/grid_029_0340_H1L1c.bin -output . -shift 1010 -infile coinc.in`

Band (-b) and overlap (-v) are required.

The coinc.in file contains a list of input trigger files with corresponding veto files.

```
/work/chuck/virgo/O3/allsky_o3_c01-results2/0340/triggers_001_0340_1.bin /work/chuck/virgo/O3/allsky_o3_c01-results2/0340/triggers_001_0340.vlines
/work/chuck/virgo/O3/allsky_o3_c01-results2/0340/triggers_002_0340_1.bin /work/chuck/virgo/O3/allsky_o3_c01-results2/0340/triggers_002_0340.vlines
/work/chuck/virgo/O3/allsky_o3_c01-results2/0340/triggers_004_0340_1.bin /work/chuck/virgo/O3/allsky_o3_c01-results2/0340/triggers_004_0340.vlines
...
```

Hemisphere in trigger file names can be 0,1 or *. In the last case both hemispheres are read.


# coi_digger.py 

Script to read the `.coi` files, calculate the False Alarm Probability for the candidates 
with a given coincidence multiplicity, extract information (also from the trigger files) 
for followup. 

`coi_digger.py` calls the `fap` code (run `make fap-many`) to compile. 

## Sample call: 

```
python coi_digger.py config.ini 0081 1
```

