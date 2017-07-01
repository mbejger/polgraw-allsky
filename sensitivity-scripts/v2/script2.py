from configparser import ConfigParser, ExtendedInterpolation
import sys, os, re 
import numpy as np 

# Parser initialisation 
p = ConfigParser(interpolation=ExtendedInterpolation())

# Parse the ini file 
p.read(sys.argv[1])

# replacement dictionary 
rd = {'GSIZE':  p.get('settings', 'gsize'), 
      'DT':     p.get('settings', 'dt'), 
      'NOD':    p.get('settings', 'nod'), 
      'THRESH': p.get('settings', 'thresh'),

      'LDLP':   p.get('paths', 'ldlp'), 
      'DATA':   p.get('paths', 'data'),
      'LOF':    p.get('paths', 'list_of_frames'), 
      'SEARCH': p.get('paths', 'sear_path'),
      'SIGEN':  p.get('paths', 'sige_path'),
      'COINCID':p.get('paths', 'coin_path'), 

      'CELL':   p.get('coincidences', 'cell'), 
      'REFFR':  p.get('coincidences', 'reffr'), 
      'USEDET': p.get('coincidences', 'usedet'), 
      'SNRCUT': p.get('coincidences', 'snrcut'),
      'MINCOIN':p.get('coincidences', 'mincoin'),} 

def replfunc(match):
    return rd[match.group(0)]

# read band-amplitudes file
with open(sys.argv[2], 'r') as f:

    for line in f:
        l = np.array(tuple(line.strip().split(' ')))

        # first column: band number 
        band = l.item(0) 

        # amplitudes: removing band from the array 
        h = np.delete(l, 0, axis=None) 

        # remove job_BAND.sub, if present 
        if os.path.exists('job_'+band+'.sub'): 
            os.remove('job_'+band+'.sub')

        # PBS jobs for one band and different amplitudes
        # stacked one after another (to have one longer job)
        js = open('job_'+band+'.sub', 'a')

        js.write(p.get('pbs', 'header')+'\n\n') 

        for a in np.nditer(h, flags=["refs_ok"]):

            # subdirectory name 
            di = str(a) + "_" + band

            # starting number of simulations
            sums=1

            # creating subdirectories        
            if not os.path.isdir(di): 
                print (di + ' does not exist! Creating...')
                os.mkdir(di)
 
                rd.update({'H0': a.item(0), 'BAND': band})
                regex = re.compile('|'.join(re.escape(i) for i in rd))

                with open(p.get('paths','scri_path')+'/dummy.sh','r') as fin, open(di+'/script.sh','w') as fout:
                    for line in fin:
                        fout.write(regex.sub(replfunc,line))

            else: 
                import os.path
                for file in os.listdir(di):
                    if file.endswith('.sum'):
                        sums = sum(1 for line in open(di+'/'+file))+1
                print (di + ' exists! Continuing with sim. #' + str(sums))  

            # At each di, check for job.sub 
            # if it doesn't exist, create it 
            if not os.path.exists(di+'/job.sub'): 
                j = open(di+'/job.sub', 'a') 
                j.write(p.get('pbs', 'header')+'\n\n') 
                j.write('cd $PBS_O_WORKDIR\n')
                j.write(p.get('pbs', 'forloop')+'\n')
                j.close() 

            # This script enters each subdirectory one after another 
            # and does the loops 
            js.write('cd %s/%s\n' % (os.getcwd(), di))
            # the actual bash for loop is defined in the config file 
            # - variable '$start' is replaced by sums 
            js.write(p.get('pbs', 'forloop').replace('$start', str(sums)) +'\n\n')
            
        js.write('exit 0\n')
        js.close()
