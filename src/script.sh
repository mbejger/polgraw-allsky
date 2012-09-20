#!/bin/bash

sig_pars=''
until [ -n "$sig_pars" ]; do  
	sig_pars=$(./sigen -d DATA -b BAND -s GSIZE -a H0)
done 

echo "$sig_pars" > sig_BAND

num_of_frames=0 

for frame in $(cat LOF|sort -r); do 

	./search -d DATA -i ${frame} -b BAND --whitenoise -a <(echo $sig_pars)
	let "num_of_frames += 1"
	if [[ $? -eq 137 ]]; then
		last_frame=${frame}
		break
	fi	

done 

if [[ $last_frame ]]; then 

	for frame in $(awk '{if($1<='${last_frame}') print $1}' LONF|sort -r); do

	        ./search -d /home/mbejger/xdat -i ${frame} -b NEXTB --whitenoise -a <(echo $sig_pars)	
		let "num_of_frames += 1"

	done 
fi 

# signal parameters and the SNR estimate
#num_of_frames=$(awk 'END{print NR}' LOF)
#sig_pars=$(awk '{if(NR>3) printf("%.16f ",$0)}' sig_BAND)

snr=$(grep "SNR:" *.err | awk '{sum +=$2*$2} END {print sqrt(sum/'$num_of_frames')}')
sim_num=${PWD##*/} 

# vetoing candidates
cd ./candidates/
for f in $(ls tri*.bin) 
	do 
		if [[ $f =~ triggers_([0-9]{2})_* ]] 
		then 
			/home/orest/trigveto/trigveto $f /home/orest/sum_0.05/sum_${BASH_REMATCH[1]}.hum /storage/VSR1-allsky-data/xdat/${BASH_REMATCH[1]}/DetSSB.bin 
		fi 
	done
# find coincydence
ls v_trig*.bin >list_148.info
/home/orest/cp_clean/cp_simple -scale CELL list_148.info > 0000_148.txt
/home/orest/cp_clean/cp_simple -ascn -scale CELL list_148.info > 0001_148.txt
/home/orest/cp_clean/cp_simple -decl -scale CELL list_148.info > 0010_148.txt
/home/orest/cp_clean/cp_simple -decl -ascn -scale CELL list_148.info > 0011_148.txt
/home/orest/cp_clean/cp_simple -sdwn -scale CELL list_148.info > 0100_148.txt
/home/orest/cp_clean/cp_simple -sdwn -ascn -scale CELL list_148.info > 0101_148.txt
/home/orest/cp_clean/cp_simple -sdwn -decl -scale CELL list_148.info > 0110_148.txt
/home/orest/cp_clean/cp_simple -sdwn -decl -ascn -scale CELL list_148.info > 0111_148.txt
/home/orest/cp_clean/cp_simple -freq -scale CELL list_148.info > 1000_148.txt
/home/orest/cp_clean/cp_simple -freq -ascn -scale CELL list_148.info > 1001_148.txt
/home/orest/cp_clean/cp_simple -freq -decl -scale CELL list_148.info > 1010_148.txt
/home/orest/cp_clean/cp_simple -freq -decl -ascn -scale CELL list_148.info > 1011_148.txt
/home/orest/cp_clean/cp_simple -freq -sdwn -scale CELL list_148.info > 1100_148.txt
/home/orest/cp_clean/cp_simple -freq -sdwn -ascn -scale CELL list_148.info > 1101_148.txt
/home/orest/cp_clean/cp_simple -freq -sdwn -decl -scale CELL list_148.info > 1110_148.txt
/home/orest/cp_clean/cp_simple -freq -sdwn -decl -ascn -scale CELL list_148.info > 1111_148.txt

# find max coincydence
max=0
for f in $(ls *.txt)
  do
    if [[ $f =~ ([0-1]{4})_([0-9]{3}).+ ]]
    then 
      sh=${BASH_REMATCH[1]}
      st=$(grep 'Scalin' $f)
      if [[ $st =~ .+\ ([0-9]+)$ ]]
      then 
        if [ ${BASH_REMATCH[1]} -gt $max ]
        then
          max=${BASH_REMATCH[1]}
          bdn=$sh
        fi
      fi
    fi
  done

cd ../
rm -fr state* candidates/* *.e* *.o*

# write down summary 
echo $sim_num $sig_pars $snr $max $bdn CELL >> H0_BAND_${sim_num}.sum

exit 0
