#!/bin/bash

# Set parity of the dataset to use for efficiency computation
# (1->procedure definition and validation, 0->fit data and systematics evaluation)
par=${1}
[ -z "${par}" ] && par=1	# set default
year=${2}
[ -z "${year}" ] && year=2016	# set default

# Number of steps in the grid used to sample the KDE functions
xbin=50
ybin=50
zbin=50

# Create directory for log files
if [ ! -d logs_test ]; then mkdir logs_test; fi

# Use external configuration file with list of bin numbers (to be coherent with the numbers in effDataset_b*.root files)
# and corresponding scale factor configuration
if make composeEff_rooKeys_parSub_test
then

    while read -a line; do
	
	bin=${line[0]}
	wid0=${line[1]}
	wid1=${line[2]}
	wid2=${line[3]}
	
	# Submit jobs for genDen, genNum, recoDen, and correct-tag recoNum
	# which use the same scale-sactor configuration
	for indx in {0..3}; do
	    
	    # Number of parallel jobs for each KDE description
	    njobs=50
	    if [ ${indx} -eq 2 ]; then njobs=500; fi # recoDen sample has very large statistics
	    
	    ./composeEff_rooKeys_parSub_test ${bin} ${indx} ${par} ${wid0} ${wid1} ${wid2} ${xbin} ${ybin} ${zbin} 0 ${njobs} ${year} \
		&> logs_test/composeEff_rooKeys_${bin}_${indx}_${par}_${wid0}_${wid1}_${wid2}_${xbin}_${ybin}_${zbin}_0_${njobs}_${year}.out

	done

	# Submit jobs for wrong-tag recoNum
	# which uses a specific scale-factor configuration
	indx=4
	wid0=${line[4]}
	wid1=${line[5]}
	wid2=${line[6]}
	njobs=50

	./composeEff_rooKeys_parSub_test ${bin} ${indx} ${par} ${wid0} ${wid1} ${wid2} ${xbin} ${ybin} ${zbin} 0 ${njobs} ${year} \
	    &> logs_test/composeEff_rooKeys_${bin}_${indx}_${par}_${wid0}_${wid1}_${wid2}_${xbin}_${ybin}_${zbin}_0_${njobs}_${year}.out

    done < ../confSF/KDE_SF.list #File containing scale-factor configuration

fi
