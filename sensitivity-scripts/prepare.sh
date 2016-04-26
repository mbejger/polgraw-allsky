#!/bin/bash

#--- General settings (please configure in mcsim.sh) -------------------
band=${1}               # band number 
h0=${2}                 # GW amplitude 
#-----------------------------------------------------------------------
code_home=${3}          # search codes 
coin_home=${4}          # coincidence codes   
script_home=${5}        # scripts 
data=${6}               # input data
list_of_frames=${7}     # file with frame list 
lz4path=${8}            # archiver (coincidences) 
ldlp=${9}               # yeppp library path 
dt=${10}                # sampling time 
howmany=${11}           # how many simulations
gsize=${12}             # +- gsize around the grid
#--- Coincidences ------------------------------------------------------
reffr=${13}             # the reference frame (coincidences)
cell=${14}              # Cell size (coincidences) 
snrcut=${15}            # Signal-to-noise cutoff 
mincoin=${16}           # Minimal no. of coincidences to register 
#-----------------------------------------------------------------------

cp ${script_home}/script_manyframes.sh script.sh
#cp ${script_home}/clean.sh .
#cp ${script_home}/get_summary.sh .

# for the case when signal enters another band
next_band=$((${band##+(0)}+1))
next_band=$(printf "%03d\n" ${next_band})

# replacements in script.sh 
sed -i 's|BAND|'$band'|g;s|GSIZE|'$gsize'|g;s|REFFR|'$reffr'|g' script.sh
sed -i 's|H0|'$h0'|g;s|CELL|'$cell'|g;s|DT|'$dt'|g' script.sh
sed -i 's|SNRCUT|'$snrcut'|g;s|MINCOIN|'$mincoin'|g' script.sh 
sed -i 's|DATA|'$data'|g;s|LOF|'$list_of_frames'|g' script.sh
#sed -i 's|NEXTB|'$next_band'|g;s|LONF|'$list_of_frames'/'$next_band'.d|g' script.sh
sed -i 's|LDLP|'$ldlp'|g;s|LZ4PATH|'$lz4path'|g' script.sh 

# loop below prepares a number of catalogues: 01, 02,..., $howmany
# makes symbolic links to programs
# writes run.sh script for the queuing system 

# removing the contents of an old run.sh
if [ -e run.sh ]; then >run.sh; fi

for dirs in $(seq -f "%03g" 1 $howmany); do

	if [ ! -d "$dirs" ]; then 
		mkdir $dirs
	fi

    # symbolic links to the codes 
	ln -sf ${code_home}/gwsearch-cpu $dirs/search
	ln -sf ${code_home}/sigen $dirs/sigen
    ln -sf ${coin_home}/coincidences $dirs/coincidences 

 	ln -sf ../script.sh $dirs/script.sh
	ln -sf ${script_home}/job.sub 

done

echo "cd "${PWD}"; qsub -N "${h0}_${band} " job.sub" >> run.sh

exit 0
