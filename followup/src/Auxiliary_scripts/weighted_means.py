import numpy as np
import sys

data = np.genfromtxt(sys.argv[1],usecols=np.arange(0,6))
N = int(sys.argv[2])
lf = int(sys.argv[3])
ref = int(sys.argv[4])
shiftedf=np.zeros(lf)
start_frame = int(sys.argv[5])
lf2=lf
m=0
k=start_frame
data2={}
for i in range(0, lf):
	if data[i,0] == -1000 and data[i,1] == -1000 and data[i,2] == -1000 and data[i,3] == -1000:
		lf2=lf2-1
	else:
		shiftedf[m] = data[i,0] + 2.0*data[i,1]*N*(ref-k)
		data2[m,0] = data[i,0]
		data2[m,1] = data[i,1]
		data2[m,2] = data[i,2]
		data2[m,3] = data[i,3]
		data2[m,4] = data[i,4]
		data2[m,5] = data[i,5]
		m=m+1
		k=start_frame+1

meanf=means=meand=meana=0
x=0
for j in range(0, lf2):
	meanf += shiftedf[j]*data2[j,5]
	means += data2[j,1]*data2[j,5]
	meand += data2[j,2]*data2[j,5]
	meana += data2[j,3]*data2[j,5]
	x=x+data2[j,5]

meanf = meanf/x
means = means/x
meand = meand/x
meana = meana/x

print meanf, means, meand, meana

#python weighted_means.py /work/ns/msieniawska/CGW_dt16_stattest/output_densegrid/8/followup/test.txt 32312 8 004