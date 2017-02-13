#!/env/python

import numpy as np 
import pandas as pd
import sys

df = pd.DataFrame()

print "band h   ul" 

with open(sys.argv[1], 'r') as f:
    for line in f:
        l = pd.DataFrame([tuple(line.strip().split(' '))])
        d = l.values

        # first column: band number 
        band = d.item(0) 

        # deleting the first element of the array
        d = np.delete(d, 0, axis=None)

        # spliting the array into two 
        # (first: amplitudes, second: % upper limits) 
        d = np.split(d, 2)

        # converting the elements to float values 
        h = d[0].astype(float)
        ul = d[1].astype(float)  

        i = 0 
        while i < len(h): 
            print "%s %4.3f %4.2f" % (band, h[i], ul[i]) 
            i += 1 

