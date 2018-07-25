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
band = float(sys.argv[3])

# band reference frequency 
fpo = 10. + 0.96875*band/(2.0*dt) 

# Frequency of the pulsar at stime 
fnow = freq + (stime - reftime)*f1dot

# Frequency in the band of width 2*dt in [0, pi]
finband = (fnow - fpo)*(2*dt)*np.pi 

# sky position
delta = pdata['DEC'] 
for i in range (0, len(delta)):
	if (delta[i] == ':'):
		colon1=i
		break
delta_hours = float(delta[0:colon1])
delta_minutes = float(delta[colon1+1:colon1+3])*60
delta_seconds = float(delta[colon1+4:])
delta_part = (delta_minutes+delta_seconds)/3600
if (delta[0] == '-'):
	delta_fin = (delta_hours-delta_part)*np.pi/180
else:
	delta_fin = (delta_hours+delta_part)*np.pi/180

alpha = pdata['RA'] 
for i in range (0, len(alpha)):
	if (alpha[i] == ':'):
		colon2=i
		break
alpha_hours = float(alpha[0:colon2])
alpha_minutes = float(alpha[colon2+1:colon2+3])*60
alpha_seconds = float(alpha[colon2+4:])
alpha_part = (alpha_minutes+alpha_seconds)/3600
alpha_part= alpha_part
alpha_fin = (alpha_hours+alpha_part)*2*np.pi/24


#  sgnlo[4] =  cos(2.*psik)*hop*cos(ph_o) - sin(2.*psik)*hoc*sin(ph_o) ;
#  sgnlo[5] =  sin(2.*psik)*hop*cos(ph_o) + cos(2.*psik)*hoc*sin(ph_o) ;
#  sgnlo[6] = -cos(2.*psik)*hop*sin(ph_o) - sin(2.*psik)*hoc*cos(ph_o) ;
#  sgnlo[7] = -sin(2.*psik)*hop*sin(ph_o) + cos(2.*psik)*hoc*cos(ph_o) ;

h0 = float(pdata['H0']) 
cosi = float(pdata['COSIOTA']) 
psi = float(pdata['PSI'])*np.pi 
phi = float(pdata['PHI0'])*np.pi 

h0p = (1.0+cosi*cosi)/2.0
h0c = cosi

a1 = np.cos(2.0*psi)*h0p*np.cos(phi) - np.sin(2.0*psi)*h0c*np.sin(phi)
a2 = np.sin(2.0*psi)*h0p*np.cos(phi) + np.cos(2.0*psi)*h0c*np.sin(phi)
a3 = -np.cos(2.0*psi)*h0p*np.sin(phi) - np.sin(2.0*psi)*h0c*np.cos(phi)
a4 = -np.sin(2.0*psi)*h0p*np.sin(phi) + np.cos(2.0*psi)*h0c*np.cos(phi)

# 3 first params are not needed in fisher.c
# put something random

print 'amp  1.00 \n5 \n003 \n', finband, '\n', f1dot, '\n', delta_fin, '\n', alpha_fin, '\n', a1, '\n', a2, '\n', a3, '\n', a4
