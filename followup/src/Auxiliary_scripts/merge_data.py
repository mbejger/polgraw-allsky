from glob import iglob
import shutil, os, sys

datapath = sys.argv[3]
newpath  = sys.argv[4]

ifo = sys.argv[1] 
band = sys.argv[2]

start = 9 #int(sys.argv[5])
end =  12 #int(sys.argv[6])
multi = 2
#howmany = end/multi

#-------------------
# Merge small files 
#-------------------

# N: 6d xdat files contain this number of data points 
# Nx3 for the ephemerids 
N = 258492 #int(sys.argv[7])
newsize = 1033969 #int(sys.argv[8])

# Destination files
xdatdname = 'xdatc_' + ifo + '_' + band + '.bin' 

xdatd = open(xdatdname, 'wb')

for i in range(start, end + 1):

    # add leading zeros 
    ii = str(i).zfill(3)

    xdat = datapath + '/' + ii + '/' + ifo + '/xdatc_' + ii + '_' + band + '.bin'
    shutil.copyfileobj(open(xdat, 'rb'), xdatd)

xdatd.close()

#-----------------
# Split big files 
#-----------------

# xdat size 
M = newsize*8

# xdat 
f = open(xdatdname, 'rb')

#MS change

i= 003 #int(end/2)
ii = str(i).zfill(3) 
xname = newpath + '/' + ii + '/' + ifo + '/xdatc_' + ii + '_' + band + '.bin'
x = file(xname, 'wb') 
x.write(f.read(M))
x.close() 

# check if file is complete, if not add zeros at the end 
b = os.path.getsize(xname)
if M-b > 0: 
    x = file(xname, 'ab')
    zeros = [0] * (M-b)
    x.write(bytearray(zeros))
    x.close() 

f.close() 

