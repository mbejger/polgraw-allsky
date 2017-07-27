#!/bin/bash

# Auxiliary temporary dir for each simulation (locally on the nodes) 
tmpdir="/tmp/${PWD#"${PWD%/*}/"}"
rm -fr ${tmpdir}
mkdir ${tmpdir} 

#-------------------
# Signal generation
#-------------------

sig_pars=''
until [ -n "$sig_pars" ]; do  
	sig_pars=$(SIGEN -amp H0 -band BAND -dt DT -gsize GSIZE -reffr REFFR -nod NOD) 
done 

# Simulation number (passed by job.sub) 
sim_num=$1

#-----------------------------------------------
# Adding signal to data and performing a search  
#-----------------------------------------------

for frame in $(cat LOF|sort -r); do 
	LD_LIBRARY_PATH=LDLP SEARCH -nod NOD -data DATA -output ${tmpdir} -ident ${frame} -band BAND -dt DT -addsig <(echo $sig_pars) --nocheckpoint -threshold THRESH -usedet USEDET 2>> ${tmpdir}/BAND_${sim_num}.out 1>> ${tmpdir}/BAND_${sim_num}.out
done 

#--------------
# Coincidences
#--------------

mkdir ${tmpdir}/coin

# Calculate band frequency from band number
fpo=$(echo "BAND DT"|awk '{printf("%.6f", 10 + 0.96875*$1/(2.0*$2))}')

# Searching for coincidences 
for shi in {0..1}{0..1}{0..1}{0..1}; do
	COINCID -nod NOD -mincoin MINCOIN -snrcutoff SNRCUT -data ${tmpdir} -trigname triggers_ -refloc DATA/REFFR -fpo $fpo -shift $shi -scale CELL -refr REFFR -dt DT -output ${tmpdir}/coin 2>> ${tmpdir}/coin/summary 1>> ${tmpdir}/BAND_${sim_num}.out
done

#-----------------
# Post-processing
#-----------------

# Selecting maximal coincidence for each pulsar among all 16 shifts (col. 5) 
# if many, select the one with highest SNR (col. 10)  
sumvar=$(sort -rgk5 -gk10 ${tmpdir}/coin/summary | head -1)

# cleanup 
cp ${tmpdir}/BAND_${sim_num}.out . 
rm -fr ${tmpdir}/triggers* ${tmpdir}/coin ${tmpdir}/*.out

# Writing down the summary with the following info:
# simulation_number injection_parameters coincidences_summary 
echo $sim_num $sig_pars $sumvar >> H0_BAND.sum

exit 0
