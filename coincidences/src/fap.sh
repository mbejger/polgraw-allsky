#!/bin/bash 

nod=2
cellsize=4
threshold=1.0

summaryfile=$1
vetofracfile=$2
griddir=$3

# bash fap.sh path/to/summary_file path/to/veto_fraction_file location/of/grid.bin
# e.g.: 
# $ bash fap.sh fap-test/summary_6666_1.txt fap-test/veto_fraction.txt fap-test 

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
