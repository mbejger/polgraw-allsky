# Pipeline: a minimal example 

### 1. Sample PBS/Torque script

This is a sample `pipeline.sh` script for the spotlight version of the `search` code, designed to look for signal around a specific location given by the spotlight range data file for 2-day narrow-band data segments. After the candidate signals are found, the search for coincidences is performed. It is currently used in the Mock Data Challenge with injected signals.  

This script uses a lower threshold (option `-threshold 10`) while searching for candidate signals, looks for coincidences in a `[4 4 4 4]` cell (see options to `coiproc`) and prepares the summary of the estimation (`summary` file) at the end of execution. The binary files `*.coi` - results of coincidences - are archived with the `lz4.gz` ([lz4](https://code.google.com/p/lz4)) compression to save space, and deleted after the coincidences procedure is done.  

```bash 
#!/bin/bash 
#PBS -m n 
#PBS -j oe
#PBS -q medium
#PBS -l mem=2GB
#PBS -l walltime=04:00:00

# Data directory 
data=/work/psk/bejger/mdc/pulsar-data
# Candidate (triggers) files output directory 
candout=/work/psk/bejger/mdc/candidates-mdc-thr10

# Directory with coincidence codes 
coisrc=/work/psk/bejger/mdc/ccoin
# Coincidences output directory 
coiout=/work/psk/bejger/mdc/coin-mdcthr10

# Path to the libraries (in the src/lib directory of the search code)
ldlp=/work/psk/bejger/mdc/spotlight/src/lib/yeppp-1.0.0/binaries/linux/x86_64

# Path to lz4 compressor
lz4path=/work/psk/bejger/codes/lz4-r131/programs

cd $PBS_O_WORKDIR

while read line; do

  # pulsar name and frequency, and band frequency (from pulsar_file_$FILENUM)
  # the master file is mdc_1561pulsars (name fpo f_gw) 
  name=$(echo $line | awk '{print $1}')
  fgw=$(echo $line | awk '{print $3}')
  fpo=$(echo $line | awk '{print $2}')

  #-----------------------
  # Search for candidates 
  #-----------------------
  echo $name $(date +%H:%M:%S) "candidate search" >> $PBS_O_WORKDIR/timing_$FILENUM

  # search in frames xxx -- yyy  
  for i in $(seq -f %03g 1 210); do

    # spotlight file range location 
    spoth1=${data}/${i}/H1/spot_${name}.dat
    spotl1=${data}/${i}/L1/spot_${name}.dat

    # Here we assume that both H1 and L1 data is present for segment $i, 
    # and take the H1 spotlight file to analyze it. If $spoth1 is not present, 
    # it means there is no good H1 data. We then try $spotl1 as the spotlight file.  
    # If both $spoth1 and $spotl1 are missing, it means there is no good data 
    # and the code exits (since there is nothing to analyze).

    spotfile=''
    if [ -f $spoth1 ]; then spotfile=$spoth1 
    else spotfile=$spotl1
    fi 

    # threshold is set to -threshold 10
    LD_LIBRARY_PATH=$ldlp ./search -data $data -ident $i -band 000 --nocheckpoint -output $candout -spotlight $spotfile -fpo $fpo -label $name -dt 2 -threshold 10 1>> $name.$PBS_JOBID.out 2>> $name.$PBS_JOBID.err

  done 

  #------------------------------------------ 
  # Search for coincidences among candidates 
  #------------------------------------------
  echo $name $(date +%H:%M:%S) "coinc: start" >> $PBS_O_WORKDIR/timing_$FILENUM

  # Creating a directory for each pulsar coincidence files 
  mkdir -p $coiout/$name

  # Searching for coincidences 
  for shi in {0..1}{0..1}{0..1}{0..1}; do
    ./coincidences -data $candout -trigname $name -refloc ${coisrc}/coinc-testdata -fpo $fpo -shift $shi -scale 4444 -refr 100 -dt 2 -mincoin 15 -narrowdown 0.2 -output ${coiout}/${name} 2>> ${coiout}/${name}/summary 1> /dev/null  
  done

  # Selecting maximal coincidence among all 16 shifts (col. 5) 
  # If many, select the one with highest SNR (col. 10)  
  sort -rgk5 -gk10 ${coiout}/${name}/summary | head -1 >> ${coiout}/summary 

  #------------------------------------------ 
  # Cleanup 
  #------------------------------------------

  echo $name $(date +%H:%M:%S) "coinc: archiving coi files" >> $PBS_O_WORKDIR/timing_$FILENUM

  for r in $(ls ${coiout}/${name}/*.coi); do 
    ${lz4path}/lz4 -6 ${r} ${r}.lz4 && gzip ${r}.lz4
    if [ -f ${r}".lz4.gz" ]; then 
      rm $r
    fi 
  done 

  echo $name $(date +%H:%M:%S) "coinc: end" >> $PBS_O_WORKDIR/timing_$FILENUM

done < pulsar_file_$FILENUM

exit 0
```
