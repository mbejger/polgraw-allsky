# Pipeline: search for candidate signals and coincidences among them 

### 1. Sample PBS/Torque script

This is a sample `pipeline.sh` script for the spotlight version of the `search` code, designed to look for signal around a specific location given by the spotlight range data file for 2-day narrow-band data segments. After the candidate signals are found, the search for coincidences is performed. It is currently used in the Mock Data Challenge with injected signals.  

This script uses a lower threshold (option `-t 10`) while searching for candidate signals, looks for coincidences in a `[4 4 4 4]` cell (see options to `coiproc`) and prepares the summary of the estimation (`summary` file) at the end of execution. The binary files `*.coi` - results of coincidences - are archived with the `lz4.gz` ([lz4](https://code.google.com/p/lz4)) compression to save space, and deleted after the coincidences procedure is done.  

```bash
#!/bin/bash 
#PBS -m n 
#PBS -j oe
#PBS -q medium
#PBS -l walltime=48:00:00

# Data directory 
data=/work/psk/bejger/mdc/pulsar-data
# Candidate (triggers) files output directory 
candout=/work/psk/bejger/mdc/cand-mdcthr10_$FILENUM

# Directory with coincidence codes 
coisrc=/work/psk/bejger/mdc/coincidences-codes-binary
# Coincidences output directory 
coiout=/work/psk/bejger/mdc/coin-mdcthr10_$FILENUM

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
  echo $name $(date +%H:%M:%S) "coincidences" >> $PBS_O_WORKDIR/timing_$FILENUM


  mkdir -p $coiout/$name; cd $coiout/$name
 
  # Copying the trigger files to the working directory 
  cp $candout/triggers*${name}*.bin .

  # Name for the list of files 
  st=$(ls triggers*${name}*.bin | tail -1); listfile=${st: -16:10}'_000'${st: -6:2}'.list';
  echo "Name for the list of files: " $listfile 

  echo "Removing the name of the pulsar from the triggers' files..."
  for f in $(ls triggers_*.bin); do mv $f ${f:0:16}${f: -6:6}; done

  echo "Vetoing: processing the triggers' files..." 
  for f in $(ls triggers_*.bin); do 
    $coisrc/trigveto -noveto -fpo $fpo $f  
    # Removing original triggers files 
    rm $f 
  done 

  echo "Removing candidates outside a narrow frequency band (+- 0.05)..."
  for f in $(ls pvc_*.bin); do
    $coisrc/decym $f $fgw 0.05
  done

  echo "Writing list of pvc files to " $listfile 
  ls pvc_*.bin > $listfile

  echo "Searching for coincidences..." 
  for i in {0..1}{0..1}{0..1}{0..1}; do 
    $coisrc/coiproc -binary -refr 100 -scale_f 4 -scale_s 4 -scale_a 4 -scale_d 4 -shift $i -fpo $fpo $listfile 
  done

  echo "Generating the stat file..." 
  st=$(grep "Max value" *.resf | sort -gk 5 | tail -1)
  mcoi=$(echo $st | awk '{print $5}')
  resf=$(echo $st | awk '{print $1}')
  resf=${resf:0:10} 
  $coisrc/resf2stat -threshold $mcoi -binary $resf 
  
  echo "Estimating the parameters..."
  $coisrc/stat2data -threshold $mcoi -refr 100 . *.stat

  echo "Preparing the summary..." 
  st=$(find . -name "*.stat.dat" -printf '%f\n')
  if [ -z $st ]; then
    echo $name $fpo "no coincidences" >> ../summary
  else
    echo $name $fpo $(wc -l < $listfile) $(awk 'END{printf "%20s %.9e %.9e %.9e %.9e %.9e",FILENAME,$3,$4,$5,$6,$7}' $st) >> ../summary
  fi 

  # Cleanup (archiving and deleting the coi files) 
  echo $name $(date +%H:%M:%S) "archiving coi files" >> $PBS_O_WORKDIR/timing_$FILENUM

  for r in $(ls *.coi); do 
    ${lz4path}/lz4 -6 ${r} ${r}.lz4 && gzip ${r}.lz4
    if [ -f ${r}".lz4.gz" ]; then 
      rm $r
    fi 
  done 
  
  cd $PBS_O_WORKDIR

  echo $name $(date +%H:%M:%S) "end" >> $PBS_O_WORKDIR/timing_$FILENUM

done < pulsar_file_$FILENUM

exit 0
```
where the `pulsar_file_$FILENUM` contains 3 columns: pulsar `name`, `fpo` reference frequency and gravitational wave frequency of the pulsar `fgw`. For an exemplary `pulsar_file_666`, contaning 8 pulsars,  
```
J0440+6416 516.601562 516.721000
J0714-6355 538.630859 538.755400
J0424+3800 561.185547 561.302200
J0303-5616 1089.595703 1089.719400
J0612+3324 390.177734 390.297600
J0319-7629 466.138672 466.263800
J0803-6214 1700.154297 1700.281800
J0534+0935 258.660156 258.781000
```
the job can be sent to the PBS queue in the following way:  
```
qsub -v FILENUM=666 -N mdc-666 pipeline.sh
```

### 2. Sample PBS/Torque script using the `C` coincidences code 

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
