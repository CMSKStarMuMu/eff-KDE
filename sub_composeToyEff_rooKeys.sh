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

seed=${3}
nseeds=10

if make composeToyEff_rooKeys_parSub
then cp composeToyEff_rooKeys_parSub /eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-Swave/
else exit 1
fi

# Use external configuration file with list of bin numbers (to be coherent with the numbers in effDataset_b*.root files)
# and corresponding scale factor configuration
while read -a line; do
    
    bin=${line[0]}
    wid0=${line[1]}
    wid1=${line[2]}
    wid2=${line[3]}

    # [ "${bin}" == "0" ] && continue
    
    # Create directory for log files
    if [ ! -d logs_parSub_b${bin}_toy${seed} ]; then mkdir logs_parSub_b${bin}_toy${seed}; fi

    # Submit jobs for genDen, genNum, recoDen, and correct-tag recoNum
    # which use the same scale-sactor configuration
    for indx in {0..3}; do
	
	# Number of parallel jobs for each KDE description
	njobs=5
	if [ ${indx} -eq 2 ]; then njobs=50; fi # recoDen sample has very large statistics
	
	# Creation of the submit HTCondor file
	cat << EOF > temp_sub_composeToyEff_rooKeys_oneBin.sub
Executable  = run_composeToyEff_rooKeys.sh
bin         = ${bin}
indx        = ${indx}
par	    = ${par}
seed        = ${seed} + \$(ProcId) / ${njobs}
wid0        = ${wid0}
wid1        = ${wid1}
wid2        = ${wid2}
xbin        = ${xbin}
ybin        = ${ybin}
zbin        = ${zbin}
ndiv        = \$(ProcId) % ${njobs}
totdiv      = ${njobs}
year        = ${year}
Arguments   = \$INT(bin) \$INT(indx) \$INT(par) \$INT(seed) \$(wid0) \$(wid1) \$(wid2) \$INT(xbin) \$INT(ybin) \$INT(zbin) \$INT(ndiv) \$INT(totdiv) \$INT(year)
Log         = logs_parSub_b${bin}_toy${seed}/sub_\$(ClusterId).log
Output      = logs_parSub_b${bin}_toy${seed}/composeToyEff_rooKeys_\$INT(bin)_\$INT(indx)_\$INT(par)_\$INT(seed)_\$(wid0)_\$(wid1)_\$(wid2)_\$INT(xbin)_\$INT(ybin)_\$INT(zbin)_\$INT(ndiv)_\$INT(totdiv)_\$INT(year).out
Error       = logs_parSub_b${bin}_toy${seed}/composeToyEff_rooKeys_\$INT(bin)_\$INT(indx)_\$INT(par)_\$INT(seed)_\$(wid0)_\$(wid1)_\$(wid2)_\$INT(xbin)_\$INT(ybin)_\$INT(zbin)_\$INT(ndiv)_\$INT(totdiv)_\$INT(year).err
transfer_output_files = ""
+JobFlavour = "testmatch"
EOF
        if [ "${USER}" == "fiorendi" ]; then
            echo '+AccountingGroup = "group_u_CMST3.all"'>>temp_sub_composeToyEff_rooKeys_oneBin.sub
        fi
        echo "Queue $((${njobs}*${nseeds}))">>temp_sub_composeToyEff_rooKeys_oneBin.sub

        # Submission and file removal
        condor_submit temp_sub_composeToyEff_rooKeys_oneBin.sub
	rm temp_sub_composeToyEff_rooKeys_oneBin.sub
	
    done

    # Submit jobs for wrong-tag recoNum
    # which uses a specific scale-factor configuration
    indx=4
    wid0=${line[4]}
    wid1=${line[5]}
    wid2=${line[6]}
    njobs=5
    
    # Creation of the submit HTCondor file
    cat << EOF > temp_sub_composeToyEff_rooKeys_oneBin.sub
Executable  = run_composeToyEff_rooKeys.sh
bin         = ${bin}
indx        = ${indx}
par	    = ${par}
seed        = ${seed} + \$(ProcId) / ${njobs}
wid0        = ${wid0}
wid1        = ${wid1}
wid2        = ${wid2}
xbin        = ${xbin}
ybin        = ${ybin}
zbin        = ${zbin}
ndiv        = \$(ProcId) % ${njobs}
totdiv      = ${njobs}
year        = ${year}
Arguments   = \$INT(bin) \$INT(indx) \$INT(par) \$INT(seed) \$(wid0) \$(wid1) \$(wid2) \$INT(xbin) \$INT(ybin) \$INT(zbin) \$INT(ndiv) \$INT(totdiv) \$INT(year)
Log         = logs_parSub_b${bin}_toy${seed}/sub_\$(ClusterId).log
Output      = logs_parSub_b${bin}_toy${seed}/composeToyEff_rooKeys_\$INT(bin)_\$INT(indx)_\$INT(par)_\$INT(seed)_\$(wid0)_\$(wid1)_\$(wid2)_\$INT(xbin)_\$INT(ybin)_\$INT(zbin)_\$INT(ndiv)_\$INT(totdiv)_\$INT(year).out
Error       = logs_parSub_b${bin}_toy${seed}/composeToyEff_rooKeys_\$INT(bin)_\$INT(indx)_\$INT(par)_\$INT(seed)_\$(wid0)_\$(wid1)_\$(wid2)_\$INT(xbin)_\$INT(ybin)_\$INT(zbin)_\$INT(ndiv)_\$INT(totdiv)_\$INT(year).err
transfer_output_files = ""
+JobFlavour = "testmatch"
EOF
    if [ "${USER}" == "fiorendi" ]; then
        echo '+AccountingGroup = "group_u_CMST3.all"'>>temp_sub_composeToyEff_rooKeys_oneBin.sub
    fi
    echo "Queue $((${njobs}*${nseeds}))">>temp_sub_composeToyEff_rooKeys_oneBin.sub

    # Submission and file removal
    condor_submit temp_sub_composeToyEff_rooKeys_oneBin.sub
    rm temp_sub_composeToyEff_rooKeys_oneBin.sub

done < ../confSF/KDE_SF.list #File containing scale-factor configuration
