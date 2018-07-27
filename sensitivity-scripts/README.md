# Monte-Carlo simulations for the sensitivity upper limits 
This directory contains a set of scripts that prepare and run a pipeline search in a small area in the parameter space around the injected signal. 

### Prerequisites 

The scripts require `python3`. A working solution is to install a `python` virtual environment (`python3` comes with a built-in `pyvenv` virtual environment software).  

#### Install `python3.4.5` locally

The installation directory is, e.g., 
```bash 
installdir=/path/to/installdir
```
then 
```bash 
mkdir -p ${installdir}; cd ${installdir} 
wget https://www.python.org/ftp/python/3.4.5/Python-3.4.5.tgz
tar zxvf Python-3.4.5.tgz
cd Python-3.4.5
make clean
./configure --prefix=$(dirname "${PWD}") 
make -j4
make install
cd ../; rm -fr Python-3.4.5*
```

#### Create virtual environment 

In a selected location (e.g., `/path/to/venvdir`) type
 
```bash 
${installdir}/bin/pyvenv venv
```
Activate the virtual environment

```bash
. /path/to/venvdir/bin/activate
```

(to leave the environment, type `deactivate`). You can now install specific packages using the `pip` installer: 

```bash
pip install numpy
pip install scipy
pip install matplotlib
pip install pandas
```

### Running the scripts 

The steps of the procedure is as follows:

1. Chose the GW strain amplitude $h_0$,
2. Randomly chose other signal parameters (with signal generator `sigen`)
3. Add signal to the data (with the `gwdetect-cpu --addsig` feature) to selected time segments, and perform the search for candidates in each of them,
4. Perform the search for coincidences (`coincidences/src`)
5. Find if the signal was detected (find the highest coincidences for a given band and compare them with the number of time segments analyzed).

The script `script.py` creates a subdirectory in which the pipeline will be launched based on the following input files:
1. `config.ini` which contains the paths to codes and the input data, and the parameters of the search: 
    * F-statistic threshold, 
    * how many simulations, 
    * which detectors to use, 
    * size of the region to search, 
    * how to perform the search for coincidences etc. 

2. `bandlist` which is a list of bands with strain amplitudes, e.g. 
```bash
0164 2.25e-1 
0165 1.5e-1 2e-1 2.5e-1 3e-1
0166 2e-1 4e-1
```
The call is
```bash
$ python script.py config.ini bandlist
```
Two other auxiliary files are:
1. Dummy `bash` script `dummy.sh` with the actual pipeline calls (variables replaced with actual values by `script.py`),
2. `PBS/Torque` script `job.sub`, launched into the cluster queue, which contains the call to `script.sh`.

Script `script.py` creates a `run.sh` file which contains commands to send the jobs into the queue. The results are summary files (`.sum`) for the requested number of simulations. In order to process them, call the `summary.py` script
```bash
$ python summary.py band coincidence_threshold number_of_simulations
```
e.g.
```bash
$ python summary.py 0165 0.7 100
```
The result will be something as follows (columns are `band` number, amplitude `h`, upper limit `ul`): 
```bash
band h   ul 
0165 0.150 0.61
0165 0.200 0.78
0165 0.250 0.95
0165 0.300 0.99
```
#### Serial (stacked) version for longer jobs 

`script2.py` creates subdirectories and a `job_BAND.sub` file for a list of amplitudes for BAND from `bandlist`, stacked one after another (can be handy to send one band as one job to the queue). Call: 

```bash
$ python script2.py config.ini bandlist
```
and then (for e.g., band 0165) send it to the queue 
```bash 
$ qsub -N 0165 -v howmany=100 job_0165.sub
```

#### Slurm scripts

Versions of scripts using Slurm instead of PBS are denoted with the `-slurm` suffix.  
