import numpy as np 
import sys
import os
import subprocess 
import re
import itertools  
import tempfile

from configparser import ConfigParser, ExtendedInterpolation

# Parser initialisation 
p = ConfigParser(interpolation=ExtendedInterpolation())

# Parse the ini file 
p.read(sys.argv[1])

band = sys.argv[2]
hemi = sys.argv[3]

coi_res = p.get('paths', 'coi_res')
veto_frac_file = p.get('paths', 'veto_frac_file')
griddir = p.get('paths', 'griddir')

coi_out_prefix = p.get('settings', 'coi_out_prefix')
coi_out_suffix = p.get('settings', 'coi_out_suffix')
coi_bin_suffix = p.get('settings', 'coi_bin_suffix')

mincoin = int(p.get('settings', 'mincoin')) 
nod = p.get('settings', 'nod') 
cellsize = p.get('settings', 'cellsize') 
threshold = p.get('settings', 'threshold') 


# read in the veto fraction file 
try: 
    fin = open(veto_frac_file, 'r') 
    print('Opening {}...'.format(veto_frac_file))
except IOError: 
    print('Problem with {} file. Exiting...'.format(veto_frac_file))
    sys.exit(1)    
else: 
    with fin:
        veto_frac = [line.rstrip('\n').split(' ') for line in fin] 

fin.close() 


# select the veto fraction corresponding to the band (one element list) 
vetofrac = [x[1] for x in veto_frac if x[0] == band][0] 

# FAP command evaluated for informations taken from the .coi files 
fap_cmd = './fap_v2 -nod ' + nod + ' -band ' + band + ' -vetofrac ' + vetofrac + ' -cellsize ' + cellsize + ' -threshold ' + threshold + ' -grid ' + griddir 


#list_of_shifts = list(itertools.product([0,1], repeat=4)) 
#mb testing  
list_of_shifts = [(1, 0, 0, 1)] 

for index, item in enumerate(list_of_shifts): 

    shift = ''.join(map(str, item)) 

    # reading .coi file (binary output from the coincidences code) 
    coi_bin = coi_res + band + '/' + shift + '_' + band + '_' + hemi + coi_bin_suffix  

    try: 
        fin = open(coi_bin, "rb") 
        print("Opening {}...".format(coi_bin)) 
    except IOError: 
        print("Problem with {} file. Exiting...".format(coi_bin))
        break  

    # number of all frames 
    frcount = int(np.fromfile(fin, dtype=np.uint16, count=1)) 

    triginfo = [] 
    for _ in range(frcount):  
        triginfo.append(np.fromfile(fin, dtype=np.int32, count=4))

    fin.close() 

    frame_cand_info_flatten = ' '.join([str(item) for sublist in triginfo for item in sublist[:-1]])    
    frdata = str(frcount) + ' ' + frame_cand_info_flatten

    #mb testing 
    frdata = '8 10 10764733 10750313 4 51309067 30399820 6 17278742 11576680 11 10033267 10020382 7 14104141 12409996 3 67469909 38261008 5 51557645 28246615 9 13038132 12858232'

    # evaluate the False Alarm Probability
    cmd = fap_cmd + ' -data <( echo "' + frdata + '")'

    x = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
    x.communicate()

#$ python fap_from_coi.py config.ini 0081 1 
#Opening vf_0041-0107_lines_v2...
#Opening /work/chuck/virgo/O2/O2-lf-results/0081/1001_0081_1.coi...
#Number of days in time segments: 24
#Input data: /dev/fd/63
#Grid matrix data directory: /work/chuck/virgo/O2/C02_24d_clean_noscience/006
#Band number: 0081 (veto fraction: 0.381069)
#The reference frequency fpo: 29.617188
#The data sampling time dt: 2.000000
#FAP threshold: 1.000000
#Cell size: 4
#3.869966e-03 8 8 154523046
#6.984099e-01 8 7 154523046


