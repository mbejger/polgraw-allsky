#!/bin/bash 

thresh=$2
howmuch=$3

# band number is read from a list given from 
# command line 
for b in $(cat $1); do

	if [ n ]; then unset n; fi 

    # list file contains the following columns:
    # 1. band number, 2-..., h0 GW amplitude
    h=$(ls -d *_${b}|sed -e 's/_'${b}'//g'|sort -g)
    h=($b $h)
    # length of h array
    hlen=$((${#h[@]}-1))

    for i in $(seq 1 ${hlen}); do

        diri=${h[$i]}_${b}
		# Summary for each band--amplitude simulation
        if [ -d "$diri" ]; then
			cd $diri
			find . -name "${diri}*.sum" -exec cat {} > ${diri}.summary \;
			n[$i]=$(awk '$20>='$thresh'*$19 {n++;} END {if(n) print n/'$howmuch'; else print 0}' ${diri}.summary)
			cd ../
        fi

    done

	echo ${h[@]} ${n[@]}

done

exit 0

