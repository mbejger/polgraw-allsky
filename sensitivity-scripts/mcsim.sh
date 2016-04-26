#!/bin/bash 

#--- General settings --------------------------------------------------
code_home=/home/michal/polgraw-allsky/search/network/src-cpu    # search codes 
coin_home=/home/michal/polgraw-allsky/coincidences/src          # coincidence codes   
script_home=/home/michal/polgraw-allsky/sensitivity-scripts     # scripts 
data=/storage/mdc/RDC_O1_0.25                                   # input data
list_of_frames=${data}/good_frames_H1L1                         # file with frame list 
lz4path=/scratch2/lz4-r131/programs                             # archiver (coincidences) 
ldlp=${code_home}/lib/yeppp-1.0.0/binaries/linux/x86_64         # yeppp library path 
dt=2                      # sampling time 
howmany=2                 # how many simulations
gsize=2                   # +- gsize around the grid
#--- Coincidences ------------------------------------------------------
reffr=031                 # the reference frame (coincidences)
cell=4444                 # Cell size (coincidences) 
snrcut=6                  # Signal-to-noise cutoff 
mincoin=5                 # Minimal no. of coincidences to register 
#-----------------------------------------------------------------------

# Band number and h0 amplitudes are read from a file 
# (command-line parameter to this script) 
while read line; do

	# line contains N+1 columns: band number + N h0 amplitudes 
	h=($line)           # whole array 
  b=$h                # band number (first column)  
  N=$((${#h[@]}-1))   # number of amplitudes 

  echo $N 

	for i in $(seq 1 $N); do

		diri=${h[$i]}_${b}

		if [ ! -d "$diri" ]; then 
			echo "Dir "$diri" does not exist! Creating..."
			mkdir $diri
		fi

        # Copy the helper script prepare.sh 
        cp ${script_home}/prepare.sh ${diri}

		cd $diri
		# Prepare subdirectories and links
		bash prepare.sh ${b} ${h[$i]} ${code_home} ${coin_home} ${script_home} ${data} ${list_of_frames} ${lz4path} ${ldlp} ${dt} ${howmany} ${gsize} ${reffr} ${cell} ${snrcut} ${mincoin}  
		# Execute jobs
		echo "Sending "$diri" jobs into the queue..."
    bash run.sh
    cd ../

	done  

done < $1

exit 0
