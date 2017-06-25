#!/bin/bash 

nod=6
cellsize=4
threshold=0.2

summaryfile=$1
vetofracfile=$2
griddir=$3

# e.g.: 
# $ bash fap.sh ../../matlabO1/files/summary_0600_0699.txt ../../matlabO1/files/veto_fraction_frame11.txt path/to/grid 

while read line 
do 

  band=$(echo $line | cut -f1 -d' ') 
  # stripping off leading 0s  
  band=$((10#${band%_?})) 

  # veto fraction corresponding to $band  
  vetofrac=$(awk '{if($1=='$band') print $2}' $vetofracfile)

  # call the code  
  ./fap -nod $nod -band $band -data <(echo $line) -grid $griddir \
   -vetofrac $vetofrac -cellsize $cellsize -threshold $threshold

done < $summaryfile 

exit 0
