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

vnjobs=(65 4 90 2 1)
# (100 20 200 5 5) previous run

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

	for indx in {0..4}
	do

	    njobs=${vnjobs[$indx]}
	    if [ $indx == 4 ]
	    then
		wid0=${line[4]}
		wid1=${line[5]}
		wid2=${line[6]}
	    fi

	    sbatch -a 0-$((${njobs}-1)) \
		   --mem 3900 \
		   --export=bin=${bin},indx=${indx},par=${par},wid0=${wid0},wid1=${wid1},wid2=${wid2},xbin=${xbin},ybin=${ybin},zbin=${zbin},totdiv=${njobs},year=${year},vers=${vers} \
		   ../run_composeEff_rooKeys_SLURM.sh
		   # -x "wn155" \
	    
	done

	cd ..

    done < ../confSF/KDE_SF.list #File containing scale-factor configuration

fi