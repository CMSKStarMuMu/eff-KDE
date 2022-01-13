#!/bin/bash

# Same configuration as in sub_composeEff_rooKeys.sh
par=${1}
[ -z "${par}" ] && par=1	# set default
year=${2}
[ -z "${year}" ] && year=2016	# set default

xbin=50
ybin=50
zbin=50

# Create directory for KDE and efficiency histograms
if [ ! -d files ]; then mkdir files; fi

while read -a line; do
    
    bin=${line[0]}

    for seed in {0..9}; do

	wid0=${line[1]}
	wid1=${line[2]}
	wid2=${line[3]}

	for indx in {0..3}; do

	    njobs=5
	    if [ ${indx} -eq 2 ]; then njobs=50; fi

	    # Merge output files from parallel jobs
	    echo 'root -b -q mergeToyEff_parSub.cc('${bin}','${indx}','${par}','${seed}','${wid0}','${wid1}','${wid2}','${xbin}','${ybin}','${zbin}','${njobs}','${year}')'
	    root -b -q 'mergeToyEff_parSub.cc('${bin}','${indx}','${par}','${seed}','${wid0}','${wid1}','${wid2}','${xbin}','${ybin}','${zbin}','${njobs}','${year}')'

	done

	indx=4
	wid0=${line[4]}
	wid1=${line[5]}
	wid2=${line[6]}
	njobs=5

	# Merge output files from parallel jobs
	echo 'root -b -q mergeToyEff_parSub.cc('${bin}','${indx}','${par}','${seed}','${wid0}','${wid1}','${wid2}','${xbin}','${ybin}','${zbin}','${njobs}','${year}')'
	root -b -q 'mergeToyEff_parSub.cc('${bin}','${indx}','${par}','${seed}','${wid0}','${wid1}','${wid2}','${xbin}','${ybin}','${zbin}','${njobs}','${year}')'

    done

done < ../confSF/KDE_SF.list
