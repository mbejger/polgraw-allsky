from glob import iglob
import shutil, os, sys

datapath = sys.argv[3]
newpath  = sys.argv[4]

ifo = sys.argv[1] #L1 or H1
band = sys.argv[2] #e.g. 0067

frame = '001'
start = int(sys.argv[5])
end =  int(sys.argv[6])

N = 1033969 #int(sys.argv[7])
newsize = 258492 #int(sys.argv[8])

# Open big file

xdat = datapath + '/' + frame + '/' + ifo + '/xdatc_' + frame + '_' + band + '.bin'

nbytes1=newsize*8
nbytes2=2*nbytes1
nbytes3=3*nbytes1
nbytes4=4*nbytes1

f = file(xdat,'rb')
data1 = f.read()[:nbytes1]
f.close()
f = file(xdat,'rb')
data2 = f.read()[nbytes1:nbytes2]
f.close()
f = file(xdat,'rb')
data3 = f.read()[nbytes2:nbytes3]
f.close()
f = file(xdat,'rb')
data4 = f.read()[nbytes3:nbytes4]
f.close()

# Create small files
for i in range(start, end + 1):
    # add leading zeros 
    ii = str(i).zfill(3)
    xname = newpath + '/' + ii + '/' + ifo + '/xdatc_' + ii + '_' + band + '.bin'
    x = file(xname, 'wb') 
    if i==start:
    	x.write(data1)
    elif i==start+1:
    	x.write(data2)
    elif i==start+2:
    	x.write(data3)
    elif i==start+3:
    	x.write(data4)
    x.close() 



