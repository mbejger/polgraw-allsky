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

# remove run.sh, if present
if os.path.exists('run.sh'): 
    os.remove('run.sh')

runsh = open('run.sh', 'a')

# read band-amplitudes file
with open(sys.argv[2], 'r') as f:

    for line in f:
        l = np.array(tuple(line.strip().split(' ')))

        # first column: band number 
        band = l.item(0) 

        # amplitudes: removing band from the array 
        h = np.delete(l, 0, axis=None) 

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
                print (di + ' exists ('+ str(sums-1)+'/'+ p.get('settings', 'howmany')+' completed)')  

            if sums < int(p.get('settings', 'howmany')):
                runsh.write('cd %s/%s; qsub -N %s -v start=%d,howmany=%s job.sub\n' % (os.getcwd(), di, di, sums, p.get('settings', 'howmany')))

            from shutil import copyfile
            copyfile(p.get('paths','scri_path')+'/job.sub', di+'/job.sub')


runsh.close()
#os.system('bash run.sh')

