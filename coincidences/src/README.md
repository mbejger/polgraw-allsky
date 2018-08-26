# coi_digger.py 

Script to read the `.coi` files, calculate the False Alarm Probability for the candidates 
with a given coincidence multiplicity, extract information (also from the trigger files) 
for followup. 

`coi_digger.py` calls the `fap` code (run `make fap-many`) to compile. 

## Sample call: 

```
python coi_digger.py config.ini 0081 1
```

