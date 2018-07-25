import numpy as np
import sys

# For fixed-width combs, veto the band:
#     [offset+index*spacing-leftwidth, offset+index*spacing+rightwidth]
# For scaling-width combs, veto the band:
#     [offset+index*spacing-index*leftwidth, offset+index*spacing+index*rightwidth]

b=26.226562 
e=26.476562

#cleaned files with lines (without header and comments)
lh1=np.genfromtxt("linie_H1.txt")
ll1=np.genfromtxt("linie_L1.txt")

print 'H1:'
for i in range(0, len(lh1)):
	if lh1[i][1]==0:
		line1=lh1[i][0]
		if e>line1>b or b<line1<e:
			print i,line1,"THIS!"
		else:
			print i,line1
	else:
		beg=lh1[i][3]
		end=lh1[i][4]
		print beg, end
		for j in range(int(beg), int(end)+1):
			line1=lh1[i][2]+j*lh1[i][0]
			if e>line1>b or b<line1<e:
				print i,line1,"THIS!"
			else:
				print i,line1

print 'L1:'
for i in range(0, len(ll1)):
	if ll1[i][1]==0:
		line1=ll1[i][0]
		if e>line1>b or b<line1<e:
			print i,line1,"THIS!"
		else:
			print i,line1
	else:
		beg=ll1[i][3]
		end=ll1[i][4]
		print beg, end
		for j in range(int(beg), int(end)+1):
			line1=ll1[i][2]+j*ll1[i][0]
			if e>line1>b or b<line1<e:
				print i,line1,"THIS!"
			else:
				print i,line1