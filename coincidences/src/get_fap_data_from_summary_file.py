import os
import sys
import subprocess
import numpy as np 
from configparser import ConfigParser, ExtendedInterpolation

# Parser initialisation 
p = ConfigParser(interpolation=ExtendedInterpolation())

# Parse the ini file 
p.read(sys.argv[1])

# data from the config file 
summaryfile  = p.get('paths', 'summaryfile')
vetofracfile = p.get('paths', 'vetofracfile')
griddir = p.get('paths', 'griddir')
 
noc   = p.get('settings', 'noc') 
nod   = p.get('settings', 'nod') 
cellf = int(p.get('settings', 'cellf'))
cells = p.get('settings', 'cells')
celld = p.get('settings', 'celld')
cella = p.get('settings', 'cella')
threshold = p.get('settings', 'threshold') 


with open(summaryfile) as f:

    lines=f.readlines()

    for line in lines:
        
        str_list = filter(None, line.split(" ")) 

        # Band number and hemisphere 
        band = str_list[0].split("_")[0] 
        hemi = str_list[0].split("_")[1]
        
        # Shift value 
        shift = str_list[1] 

        # total number of frames 
        nof = int(str_list[3])
    
        # total number of frames, frame numbers, all candidates and unique candidates for each frame  
        frdata = str(nof) + " " + " ".join(str_list[10:10+3*nof])

        # Veto fraction from vetofracfile 
        try: 
            fin = open(vetofracfile, 'r') 
            #print('Opening {}...'.format(vetofracfile))
        except IOError: 
            print('Problem with {} file. Exiting...'.format(vetofracfile))
            sys.exit(1)    
        else: 
            with fin:
                veto_frac = [lv.rstrip('\n').split(' ') for lv in fin] 

        fin.close() 

        # select the veto fraction corresponding to the band (one element list) 
        vetofrac = [x[1] for x in veto_frac if x[0] == band][0] 
        #print(band, vetofrac) 

        # FAP command 
        fap_cmd = './fap_new_coi_header_hyper_improved -nod ' + nod + ' -band ' + band + ' -vetofrac ' + vetofrac + ' -cellf ' + str(cellf) + ' -cells ' + cells + ' -celld ' + celld +  ' -cella ' + cella + '-noc ' + noc + ' -threshold ' + threshold + ' -grid ' + griddir + ' -data <( echo "' + frdata + '")' 

        p = subprocess.Popen(fap_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, executable='/bin/bash')
        out, err = p.communicate()

        #print(out) 
        print("{} {} {} {}".format(band, hemi, shift, err)) 

