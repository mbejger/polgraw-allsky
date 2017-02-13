#!/env/python

### python ul.py 0666_results 0666 0.01 test.pdf 

import numpy as np 
import pandas as pd
import sys
from scipy.optimize import curve_fit

def sigmoid(x, x0, k):
     y = 1.0 / (1.0 + np.exp(k*(x0-x)))
     return y

# Reading the data into pandas DataFrame
d = pd.read_table(sys.argv[1], sep= "\s+");

# This selects the entries with 'band' column equal to sys.argv[2]
data = d[d.band == int(sys.argv[2])] 
# casting the other columns into arrays
xdata = data['h'].as_matrix()
ydata = data['ul'].as_matrix()
# print xdata, ydata 

# fitting with initial quess values of x0 and k (=p0[0], p0[1])
popt, pcov = curve_fit(sigmoid, xdata, ydata, p0=[float(sys.argv[3]), 50])
#print sys.argv[2], popt

x = np.linspace(0.0, 3., 500)
y = sigmoid(x, *popt)

# our desired upper limit value equal to 0.95  
ul95 = 0.95
Ax = - np.log(1./ul95 - 1)/(popt[1]) + popt[0]

print "%04d %.4e" % (int(sys.argv[2]), Ax)

if len(sys.argv) > 4: 

    import matplotlib.pyplot as plt
    import matplotlib.ticker as plticker  
    from time import sleep 

    fig, ax = plt.subplots() 
    plt.plot(xdata, ydata, 'o', label='data')
    plt.plot(x,y, label='fit')
    plt.plot(Ax, ul95, 'ro', markersize=12, alpha=0.5, label='res')
    plt.xlim(0, 1.75*Ax)
    plt.ylim(0, 1.05)
    plt.grid(True)
    plt.xlabel("h", fontsize=16)
    plt.ylabel("ul [%]", fontsize=16) 
    plt.title("band #" + sys.argv[2], fontsize=16) 

#    ax.xaxis.set_minor_locator(plticker.AutoMinorLocator(4))
#    ax.yaxis.set_minor_locator(plticker.AutoMinorLocator(4))
#    ax.xaxis.set_major_locator(plticker.MultipleLocator(base=0.05))

    plt.savefig(sys.argv[4], filetype="pdf")

