#!/bin/bash

#-------------------
# Signal generation
#-------------------

sig_pars=''
until [ -n "$sig_pars" ]; do  
	sig_pars=$(./sigen -amp H0 -band BAND -dt DT -gsize GSIZE -reffr REFFR) 
done 

# Recovering the hemisphere number
hemi=($sig_pars); hemi=${hemi[2]}

# Saving injection details ''just in case''
echo $sig_pars > sig_BAND

num_of_frames=0 

#rm BAND.snr 

#----------------------------------
# Adding signal to data and search  
#----------------------------------

for frame in $(cat LOF|sort -r); do 

	LD_LIBRARY_PATH=LDLP ./search -data DATA -output . -ident ${frame} -band BAND -dt DT -addsig <(echo $sig_pars) --nocheckpoint -threshold THRESH 2>> /dev/null 1>> BAND.out
	let "num_of_frames += 1"
	# ./search exits with exit(137) when the signal goes out of band 
	if [[ $? -eq 137 ]]; then
		last_frame=${frame}
		break
	fi	

done 

## this loop is used when the signal goes out of the initial band
#if [[ $last_frame ]]; then 
#
#	for frame in $(awk '{if($1<='${last_frame}') print $1}' LONF|sort -r); do
#
#	        ./search -d DATA -i ${frame} -b NEXTB --whitenoise -a <(echo $sig_pars)	2>> BAND.snr
#		let "num_of_frames += 1"
#
#	done 
#fi 

# signal-to-noise estimate (mean of the values generated during the search)
# snr=$(grep "SNR:" BAND.snr | awk '{sum +=$2*$2} END {print sqrt(sum/'$num_of_frames')}')

sim_num=${PWD##*/} 

#--------------
# Coincidences
#--------------

mkdir BAND 

# Calculate band frequency from band number
fpo=$(echo "BAND DT"|awk '{print 10 + 0.96875*$1/(2.0*$2)}')

# Searching for coincidences 
for shi in {0..1}{0..1}{0..1}{0..1}; do
	./coincidences -mincoin MINCOIN -snrcutoff SNRCUT -data . -trigname triggers_ -refloc DATA/REFFR -fpo $fpo -shift $shi -scale CELL -refr REFFR -dt DT -output BAND 2>> BAND/summary 1>> BAND.out  
done

# Selecting maximal coincidence for each pulsar among all 16 shifts (col. 5) 
# if many, select the one with highest SNR (col. 10)  
sort -rgk5 -gk10 BAND/summary | head -1 >> summary 
sumvar=$(<summary)

# Cleanup (archiving and deleting the coi files) 
for r in $(ls BAND/*.coi); do 
	LZ4PATH/lz4 -6 ${r} ${r}.lz4 && gzip ${r}.lz4
	if [ -f ${r}".lz4.gz" ]; then 
    	rm $r
    fi 
done 

# cleanup 
rm -fr triggers_* *.e* *.o* BAND

cd ../

# write down the summary with the following info:
# simulation_number injection_parameters coincidences_summary 
echo $sim_num $sig_pars $sumvar >> H0_BAND_${sim_num}.sum

exit 0
