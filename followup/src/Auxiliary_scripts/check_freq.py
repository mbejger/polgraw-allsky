import numpy as np
#from astropy.time import Time 
import sys 

# read in the pulsar .par file  
pdata = {}
with open(sys.argv[1]) as f:
    for line in f:
        l = line.split()
        pdata[l[0]] = l[1]

# convert the MJD time to GPS 
reftime = 930582151 #Time(float(pdata['PEPOCH']), format='mjd').gps 

# read the gps starting time of frame from starting_date file 
with open(sys.argv[2]) as f:
    for line in f:
        l = line.split()
        # expecting only one number, the gps time in seconds 
        stime = float(l[0]) 

# GW frequency and frequency derivative (2 \times F0 and F1) 
freq = 2.*float(pdata['F0']) 
f1dot = 2.*float(pdata['F1']) 

# Our setup: sampling time dt and band number 
dt = 2
#band = 67
#check every band to find in which pulsar will be visible
for band in range(41,3999+1): 
# band reference frequency 
	fpo = 10. + 0.96875*band/(2.0*dt)  

# Frequency of the pulsar at stime 
	fnow = freq + (stime - reftime)*f1dot

# Frequency in the band of width 2*dt in [0, pi]
	finband = (fnow - fpo)*(2*dt)*np.pi 
	if (finband>=0.0 and finband<=np.pi):
		print(sys.argv[1], band, finband) 





 
