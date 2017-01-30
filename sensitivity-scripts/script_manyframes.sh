#!/bin/bash

# Auxiliary temporary dir for each simulation (locally on the nodes) 
tmpdir="/tmp/${PWD#"${PWD%/*/*}/"}"
parenttmpdir="$(dirname "${tmpdir}")"

rm -fr ${parenttmpdir}/*
mkdir -p ${tmpdir} 

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
echo $sig_pars > ${tmpdir}/sig_BAND

num_of_frames=0 

#rm BAND.snr 

#----------------------------------
# Adding signal to data and search  
#----------------------------------

for frame in $(cat LOF|sort -r); do 

	LD_LIBRARY_PATH=LDLP ./search -data DATA -output . -ident ${frame} -band BAND -dt DT -addsig <(echo $sig_pars) --nocheckpoint -threshold THRESH -usedet USEDET 2>> /dev/null 1>> BAND.out
	let "num_of_frames += 1"
	# ./search exits with exit(137) when the signal goes out of band 
	if [[ $? -eq 137 ]]; then
		last_frame=${frame}
		break
	fi	

done 

sim_num=${PWD##*/} 

#--------------
# Coincidences
#--------------

mkdir ${tmpdir}/BAND 

# Calculate band frequency from band number
fpo=$(echo "BAND DT"|awk '{printf("%.6f", 10 + 0.96875*$1/(2.0*$2))}')

# Searching for coincidences 
for shi in {0..1}{0..1}{0..1}{0..1}; do
	./coincidences -mincoin MINCOIN -snrcutoff SNRCUT -data ${tmpdir} -trigname triggers_ -refloc DATA/REFFR -fpo $fpo -shift $shi -scale CELL -refr REFFR -dt DT -output ${tmpdir}/BAND 2>> ${tmpdir}/BAND/summary 1>> ${tmpdir}/BAND.out  
done

# Selecting maximal coincidence for each pulsar among all 16 shifts (col. 5) 
# if many, select the one with highest SNR (col. 10)  
sort -rgk5 -gk10 ${tmpdir}/BAND/summary | head -1 >> summary 
sumvar=$(<summary)

# Cleanup (archiving and deleting the coi files) 
#for r in $(ls BAND/*.coi); do 
#	LZ4PATH/lz4 -6 ${r} ${r}.lz4 && gzip ${r}.lz4
#	if [ -f ${r}".lz4.gz" ]; then 
#    	rm $r
#    fi 
#done 

# cleanup 
rm -fr ${tmpdir}
rm -fr *.e* *.o* wisdom*

cd ../

# write down the summary with the following info:
# simulation_number injection_parameters coincidences_summary 
echo $sim_num $sig_pars $sumvar >> H0_BAND_${sim_num}.sum

exit 0
