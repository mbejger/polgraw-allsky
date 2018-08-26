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

fap_cmd = './fap-many -nod ' + nod + ' -band ' + band + ' -vetofrac ' + vetofrac + ' -cellsize ' + cellsize + ' -threshold ' + threshold + ' -grid ' + griddir 


list_of_shifts = list(itertools.product([0,1], repeat=4)) 
#list_of_shifts = [(1, 0, 0, 1)] 

for index, item in enumerate(list_of_shifts): 

    shift = ''.join(map(str, item)) 

    # reading coi_out file (stdout from the coincidences code)
    #mb this part should be replaced by reading this data from 
    # the .coi file header  
    coi_out = coi_res + band + '/' + coi_out_prefix + band + '_' + hemi + '-' + shift + coi_out_suffix

    frame_cand_info = [] 
 
    try: 
        fin = open(coi_out, "r") 
        print("Opening {}...".format(coi_out))
    except IOError: 
        print("Problem with {} file. Exiting...".format(coi_out))
        break  
    else: 
        with fin:
            lines = [line.rstrip('\n') for line in fin] 

    fin.close()         
   
    for ind, line in enumerate(lines):

        # first line - info on best coincidence 
        if ind == 0: 
            best_coincidence = list(filter(None, line.split(' '))) 
            bc = ' '.join(best_coincidence[0:4])   

        res = re.search('.*Frame (.*): (.*)/(.*)', line)
    
        if res: 
            frame_number = int(res.group(1))
            all_cands    = int(res.group(3))
            unique_cands = int(res.group(2))
     
            frame_cand_info.append([frame_number, all_cands, unique_cands]) 

    frame_cand_info_flatten = [item for sublist in frame_cand_info for item in sublist]


    # reading .coi file (binary output from the coincidences code) 
    coi_bin = coi_res + band + '/' + shift + '_' + band + '_' + hemi + coi_bin_suffix  

    # starting value of coin, set to maximal for this shift  
    coin = int(best_coincidence[3])  

    try: 
        fin = open(coi_bin, "rb") 
        print("Opening {}...".format(coi_bin)) 
    except IOError: 
        print("Problem with {} file. Exiting...".format(coi_bin))
        break  


    # temporary file for the coincidences info - input for the fap code  
    fp = tempfile.NamedTemporaryFile(delete=False)

    # dtype for first 6 numbers: coincidence multiplicity, 
    # mean parameters of the coincidence (f, s, d, a) and snr 
    dtype = [('coin', 'uint16'), 
            ('f', 'float32'),
            ('s', 'float32'),
            ('d', 'float32'),
            ('a', 'float32'),
            ('snr', 'float32'),]

    with fin:

        while coin >= mincoin: 

            coin_mean_vals_of_pars = np.fromfile(fin, dtype=dtype, count=1)
             
            coin = int(coin_mean_vals_of_pars['coin']) 
            cmvop = ' '.join(map(str, coin_mean_vals_of_pars[0,])) 
   
            # reading frames info  
            frames = []      
            for _ in range(coin):  
                frames.append(int(np.fromfile(fin, dtype=np.uint16, count=1)))

            # positions of triggers in the trigger files (not used right now) 
            tripos = np.fromfile(fin, dtype=np.int32, count=coin) 

            # all frames information (number, all cands, unique cands) 
            fci = ' '.join([str(x) for x in frame_cand_info_flatten]) 
            # numbers of frames in coincidence 
            fic  = ' '.join([str(x) for x in frames]) 

            coin_data = bc + ' ' + cmvop + ' ' + fci + ' ' + fic 

            try: 
                fp.write(coin_data.encode())
                fp.write('\n'.encode()) 
            except IOError: 
                print('Problem with the temporary file {}. Exiting...'.format(fp.name))
                break            
 
    fin.close()
    fp.close() 

    # FAP results files 
    fin = open(shift + '_' + band + '_' + hemi + '.fap', 'w') 

    # evaluate the False Alarm Probability 
    cmd = fap_cmd + ' -data ' + fp.name
    subprocess.call([cmd], shell=True, stderr=fin)
    os.unlink(fp.name) 
    
    fin.close() 

    #mb better output handling needed
    # do not write if output empty (use Popen instead of call) 
 
