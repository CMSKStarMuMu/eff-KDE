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

if make composeEff_rooKeys_parSub; then

    # Use external configuration file with list of bin numbers (to be coherent with the numbers in effDataset_b*.root files)
    # and corresponding scale factor configuration
    while read -a line; do
	
	bin=${line[0]}
	wid0=${line[1]}
	wid1=${line[2]}
	wid2=${line[3]}
	vers=${line[7]}

	# Create directory for log files
	[ -z "${vers}" ] && vers="-1"
	if [ ! -d logs_parSub_v${vers} ]; then mkdir logs_parSub_v${vers}; fi

	if [ ! -d tmpdir_v${vers} ]; then mkdir tmpdir_v${vers}; fi
	cd tmpdir_v${vers}

	# Submit jobs for genDen, genNum, recoDen, and correct-tag recoNum
	# which use the same scale-sactor configuration
	for indx in {0..3}; do
	    
	    # Number of parallel jobs for each KDE description
	    njobs=20
	    if [ ${indx} -eq 0 ]; then njobs=100; fi
	    if [ ${indx} -eq 2 ]; then njobs=200; fi
	    if [ ${indx} -eq 3 ]; then njobs=5; fi
	    
	    sbatch -a 0-$((${njobs}-1)) \
		   --mem 2900 \
		   --export=bin=${bin},indx=${indx},par=${par},wid0=${wid0},wid1=${wid1},wid2=${wid2},xbin=${xbin},ybin=${ybin},zbin=${zbin},totdiv=${njobs},year=${year},vers=${vers} \
		   ../run_composeEff_rooKeys.sh
	    
	done

	# Submit jobs for wrong-tag recoNum
	# which uses a specific scale-factor configuration
	indx=4
	wid0=${line[4]}
	wid1=${line[5]}
	wid2=${line[6]}
	njobs=5
	
	sbatch -a 0-$((${njobs}-1)) \
	       --mem 2900 \
	       --export=bin=${bin},indx=${indx},par=${par},wid0=${wid0},wid1=${wid1},wid2=${wid2},xbin=${xbin},ybin=${ybin},zbin=${zbin},totdiv=${njobs},year=${year},vers=${vers} \
	       ../run_composeEff_rooKeys.sh

	cd ..

    done < ../confSF/KDE_SF.list #File containing scale-factor configuration

fi
