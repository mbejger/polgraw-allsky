import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys

import matplotlib.style
import matplotlib as mpl
mpl.style.use('classic')

dt= np.dtype([('frequency','f8'), ('spindown', 'f8'), ('delta', 'f8'), ('alpha', 'f8'), ('snr', 'f8')])

#Files with triggers
fp_v1 = open(sys.argv[1],"rb") #results with veto; 1st hemisphere
fp1 = open(sys.argv[2],"rb") #results without veto; 1st hemisphere
data_v1 = np.fromfile(fp_v1, dtype=dt)
data1 = np.fromfile(fp1, dtype=dt)
fp_v2 = open(sys.argv[3],"rb") #results with veto; 2nd hemisphere
fp2 = open(sys.argv[4],"rb") #results without veto; 2nd hemisphere
data_v2 = np.fromfile(fp_v2, dtype=dt)
data2 = np.fromfile(fp2, dtype=dt)

#calculated lines in band
l1_003=[3.136128, 3.141593]
l2_003=[0.000000, 0.593648]
l3_003=[1.879523, 2.479437]
l4_003=[3.136008, 3.141593]

#czestosc odnosna:
fpo=26.226562

#linia H1:
#Maybe OMC length dither
line1 = 26.2499612
line1_band=(line1-fpo)*np.pi/0.25
#linia L1:
#First comb
line2=26.40000000000
line2_band=(line2-fpo)*np.pi/0.25

fig = plt.figure()

ax1=fig.add_subplot(111) 
ax1.grid(True)
ax1.set_xlim(0,np.pi)
ax1.set_ylim(0.,50.)

ax1.set_title(r'Both hemispheres')

ax1.set_ylabel(r'SNR',fontsize=16)
ax1.set_xlabel(r'Frequency',fontsize=16)
ax1.fill_between(l1_003, 50., color = 'grey', alpha = 0.5)
ax1.fill_between(l2_003, 0, 50., color = 'grey', alpha = 0.5)
ax1.fill_between(l3_003, 0, 50., color = 'grey', alpha = 0.5)
ax1.fill_between(l4_003, 0, 50., color = 'grey', alpha = 0.5)
f_v1=data_v1['frequency']
snr_v1=data_v1['snr']
f1=data1['frequency']
snr1=data1['snr']

f_v2=data_v2['frequency']
snr_v2=data_v2['snr']
f2=data2['frequency']
snr2=data2['snr']

ax1.plot(f1[::200], snr1[::200], linestyle='None', marker = 'o', color='darkturquoise', markeredgecolor='none', label='without vetos')
ax1.plot(f_v1[::200], snr_v1[::200], linestyle='None', marker = 'o', color='deeppink', markeredgecolor='none', label='with vetos')
ax1.plot(f2[::200], snr2[::200], linestyle='None', marker = 'o', color='darkturquoise', markeredgecolor='none')
ax1.plot(f_v2[::200], snr_v2[::200], linestyle='None', marker = 'o', color='deeppink', markeredgecolor='none')
#ta wartosc jest inna dla kazdego pasma
pulsar10=1.40878132563
ax1.axvline(x=pulsar10, ymin=0.0, ymax=50.0, c = 'gold', linestyle = '-', linewidth = 2.5, label='pulsar10')

ax1.axvline(x=line1_band, ymin=0.0, ymax=50.0, c = 'black', linestyle = '-', linewidth = 1.5, label= 'unknown comb (H1)')
ax1.axvline(x=line2_band, ymin=0.0, ymax=50.0, c = 'black', linestyle = '--', linewidth = 1.5, label= 'first 0.0 comb (L1)')
ax1.annotate('line 26.2499612 Hz', xy=(line1_band+0.01, 30.01), xytext=(line1_band+0.05, 31))
ax1.annotate('line 26.4000000 Hz', xy=(line2_band+0.01, 30.01), xytext=(line2_band+0.05, 31))

#koincydencje 1.40253345e+00 -3.22116563e-09 5.66966455e-01 4.53142908e+00 2.970208e+01
f_coi=1.40253345e+00
ax1.axvline(x=f_coi, ymin=0.0, ymax=50.0, c = 'purple', linestyle = ':', linewidth = 2.5, label='coincidences')

#followup 1.414145e+00 -4.466309e-09 5.904958e-01 4.498395e+00 -3.489938e+02 2.634364e+01
f_followup=1.414145e+00
snr_followup=2.634364e+01
ax1.plot(f_followup,snr_followup, marker='*', color='lime', markersize=10)
ax1.annotate('Followup', xy=(f_followup+0.01, snr_followup+0.01), xytext=(f_followup+0.05, snr_followup+0.05))

ax1.legend(loc='upper right', fontsize = '10')

#plt.show()
plt.savefig("summary_plot.png", format="png", bbox_inches="tight")
