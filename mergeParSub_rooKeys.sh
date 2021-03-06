#!/bin/bash

# Same configuration as in sub_composeEff_rooKeys.sh
par=1
xbin=50
ybin=50
zbin=50
year=2016

# Create directory for KDE and efficiency histograms
if [ ! -d files ]; then mkdir files; fi

while read -a line; do
    
    bin=${line[0]}
    wid0=${line[1]}
    wid1=${line[2]}
    wid2=${line[3]}
    
    for indx in {0..3}; do

	njobs=50
	if [ ${indx} -eq 2 ]; then njobs=500; fi

	# Merge output files from parallel jobs
	echo 'root -b -q mergeParSub_rooKeys.cc('${bin}','${indx}','${par}','${wid0}','${wid1}','${wid2}','${xbin}','${ybin}','${zbin}','${njobs}','${year}')'
	root -b -q 'mergeParSub_rooKeys.cc('${bin}','${indx}','${par}','${wid0}','${wid1}','${wid2}','${xbin}','${ybin}','${zbin}','${njobs}','${year}')'
	
    done

    indx=4
    wid0=${line[4]}
    wid1=${line[5]}
    wid2=${line[6]}
    njobs=50
    
    # Merge output files from parallel jobs
    echo 'root -b -q mergeParSub_rooKeys.cc('${bin}','${indx}','${par}','${wid0}','${wid1}','${wid2}','${xbin}','${ybin}','${zbin}','${njobs}','${year}')'
    root -b -q 'mergeParSub_rooKeys.cc('${bin}','${indx}','${par}','${wid0}','${wid1}','${wid2}','${xbin}','${ybin}','${zbin}','${njobs}','${year}')'
    
done < ../confSF/KDE_SF.list
