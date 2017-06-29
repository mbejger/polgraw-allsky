import os, sys
import numpy as np

band   = sys.argv[1]
thresh = float(sys.argv[2]) 
howmany= float(sys.argv[3])

# list directories ending with band
file_names = [d for d in os.listdir('.') if os.path.isdir(d) and d.endswith(band)]

for d in file_names:  
    m = 0 

    for f in os.listdir(d):

        # look for .sum file in a given subdirectory 
        if f.endswith('.sum'):  
            s = os.path.join(d, f) 

            # read the file line by line 
            with open (s, 'r') as sf:
                for line in sf:
                    l = line.split(' ')

                    if len(l)>15: 
                        # check how many crossed the threshold 
                        # col. 16: number of frames 
                        # col. 17: number of coincidences 
                        if (int(l[17]) > thresh*int(l[16])):
                            m += 1
    
    # output the results: band number, h0 amplitude, fraction of significant coincidences 
    print ('{} {:.3f} {:.2f}'.format(band, float(d.replace('_'+band,'')), m/howmany))         
