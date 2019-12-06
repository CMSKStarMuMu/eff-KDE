#!/bin/bash

# Set parity of the dataset to use for efficiency computation
# (1->procedure definition and validation, 0->fit data and systematics evaluation)
par=1

# Number of steps in the grid used to sample the KDE functions
xbin=50
ybin=50
zbin=50
year=2017

# Create directory for log files
if [ ! -d logs_parSub ]; then mkdir logs_parSub; fi

# Use external configuration file with list of bin numbers (to be coherent with the numbers in effDataset_b*.root files)
# and corresponding scale factor configuration
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
	
	# Creation of the submit HTCondor file
	cat << EOF > temp_sub_composeEff_rooKeys_oneBin.sub
Executable  = run_composeEff_rooKeys.sh
bin         = ${bin}
indx        = ${indx}
par	    = ${par}
wid0        = ${wid0}
wid1        = ${wid1}
wid2        = ${wid2}
xbin        = ${xbin}
ybin        = ${ybin}
zbin        = ${zbin}
ndiv        = \$(ProcId)
totdiv      = ${njobs}
year        = ${year}
Arguments   = \$INT(bin) \$INT(indx) \$INT(par) \$(wid0) \$(wid1) \$(wid2) \$INT(xbin) \$INT(ybin) \$INT(zbin) \$INT(ndiv) \$INT(totdiv) \$INT(year)
Log         = logs_parSub/sub_\$(ClusterId).log
Output      = logs_parSub/composeEff_rooKeys_\$INT(bin)_\$INT(indx)_\$INT(par)_\$(wid0)_\$(wid1)_\$(wid2)_\$INT(xbin)_\$INT(ybin)_\$INT(zbin)_\$INT(ndiv)_\$INT(totdiv)_\$INT(year).out
Error       = logs_parSub/composeEff_rooKeys_\$INT(bin)_\$INT(indx)_\$INT(par)_\$(wid0)_\$(wid1)_\$(wid2)_\$INT(xbin)_\$INT(ybin)_\$INT(zbin)_\$INT(ndiv)_\$INT(totdiv)_\$INT(year).err
+JobFlavour = "testmatch"
EOF
        if [ "${USER}" == "fiorendi" ]; then
            echo '+AccountingGroup = "group_u_CMST3.all"'>>temp_sub_composeEff_rooKeys_oneBin.sub
        fi
        echo 'Queue ${njobs}'>>prova.sub

        # Submission and file removal
        condor_submit temp_sub_composeEff_rooKeys_oneBin.sub
	rm temp_sub_composeEff_rooKeys_oneBin.sub
	
    done

    # Submit jobs for wrong-tag recoNum
    # which uses a specific scale-factor configuration
    indx=4
    wid0=${line[4]}
    wid1=${line[5]}
    wid2=${line[6]}
    njobs=50
    
    # Creation of the submit HTCondor file
    cat << EOF > temp_sub_composeEff_rooKeys_oneBin.sub
Executable  = run_composeEff_rooKeys.sh
bin         = ${bin}
indx        = ${indx}
par	    = ${par}
wid0        = ${wid0}
wid1        = ${wid1}
wid2        = ${wid2}
xbin        = ${xbin}
ybin        = ${ybin}
zbin        = ${zbin}
ndiv        = \$(ProcId)
totdiv      = ${njobs}
year        = ${year}
Arguments   = \$INT(bin) \$INT(indx) \$INT(par) \$(wid0) \$(wid1) \$(wid2) \$INT(xbin) \$INT(ybin) \$INT(zbin) \$INT(ndiv) \$INT(totdiv) \$INT(year)
Log         = logs_parSub/sub_\$(ClusterId).log
Output      = logs_parSub/composeEff_rooKeys_\$INT(bin)_\$INT(indx)_\$INT(par)_\$(wid0)_\$(wid1)_\$(wid2)_\$INT(xbin)_\$INT(ybin)_\$INT(zbin)_\$INT(ndiv)_\$INT(totdiv)_\$INT(year).out
Error       = logs_parSub/composeEff_rooKeys_\$INT(bin)_\$INT(indx)_\$INT(par)_\$(wid0)_\$(wid1)_\$(wid2)_\$INT(xbin)_\$INT(ybin)_\$INT(zbin)_\$INT(ndiv)_\$INT(totdiv)_\$INT(year).err
+JobFlavour = "testmatch"
EOF
    if [ "${USER}" == "fiorendi" ]; then
        echo '+AccountingGroup = "group_u_CMST3.all"'>>temp_sub_composeEff_rooKeys_oneBin.sub
    fi
    echo 'Queue ${njobs}'>>prova.sub

    # Submission and file removal
    condor_submit temp_sub_composeEff_rooKeys_oneBin.sub
    rm temp_sub_composeEff_rooKeys_oneBin.sub

done < ../confSF/KDE_SF.list #File containing scale-factor configuration
