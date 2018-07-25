from glob import iglob
import shutil, os, sys

datapath = sys.argv[3]
newpath  = sys.argv[4]

ifo = sys.argv[1] 
band = sys.argv[2]

start = 1
end = 2
multi = 1
howmany = 1
new_frame = 1

#-------------------
# Merge small files 
#-------------------

# N: 6d xdat files contain this number of data points 
# Nx3 for the ephemerids 
N = 129246
newsize = 258492

# Destination files
detssbdname= 'DetSSB_' + ifo + '.bin'
rdetdname = 'rDet_' + ifo + '.bin'
rssbdname = 'rSSB_' + ifo + '.bin'

detssbd = file(detssbdname, 'w+b')
rdetd = open(rdetdname, 'wb')
rssbd = open(rssbdname, 'wb')

for i in range(start, end + 1):

    # add leading zeros 
    ii = str(i).zfill(3)

    detssb = datapath + '/' + ii + '/' + ifo + '/DetSSB.bin'
    f = file(detssb,'rb')
    # read all except last 16 bytes 
    data = f.read()[:-16]
    detssbd.write(data)
    f.close()

    rdet = datapath + '/' + ii + '/' + ifo + '/rDet.bin'
    shutil.copyfileobj(open(rdet, 'rb'), rdetd)

    rssb = datapath + '/' + ii + '/' + ifo + '/rSSB.bin'
    shutil.copyfileobj(open(rssb, 'rb'), rssbd)

detssbd.close()
rdetd.close() 
rssbd.close()

#-----------------
# Split big files 
#-----------------

# xdat size 
M = newsize*8

#MS change
start = new_frame

# The ephemerids size 
M = newsize*8*3 

# DetSSB 
f = open(detssbdname, 'rb')

for i in range(start, howmany + 1):

    ii = str(i).zfill(3) 
    xname = newpath + '/' + ii + '/' + ifo + '/DetSSB.bin'
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

    # add two additional doubles at the end 
    # from the original DetSSB.bin, first in the new sequence 
    # hence new i -> old multi*(i-1)+1 
    detssb = datapath + '/' + str(multi*(i-1)+1).zfill(3) + '/' + ifo + '/DetSSB.bin'
    d = file(detssb,'rb')
    # read last 16 bytes 
    data = d.read()[-16:]
    d.close() 

    x = file(xname, 'ab') 
    x.write(data)
    x.close() 

f.close()

# rSSB
f = open(rssbdname, 'rb')

for i in range(start, howmany + 1):

    ii = str(i).zfill(3) 
    xname = newpath + '/' + ii + '/' + ifo + '/rSSB.bin'
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

# rDet
f = open(rdetdname, 'rb')

for i in range(start, howmany + 1):

    ii = str(i).zfill(3) 
    xname = newpath + '/' + ii + '/' + ifo + '/rDet.bin'
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

