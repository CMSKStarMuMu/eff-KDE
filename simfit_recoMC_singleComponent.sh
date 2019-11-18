#!/bin/bash

par=1
tag=1
year1=2016
year2=2017
year3=0

# Create directories for fit logs, results and plots
if [ ! -d logs_simFit ]; then mkdir logs_simFit; fi
if [ ! -d simFitResults ]; then mkdir simFitResults; fi
if [ ! -d plotSimFit_d ]; then mkdir plotSimFit_d; fi

# Compile dictionary and macro
make AngDict
make simfit_recoMC_singleComponent

while [ ${tag} -ge 0 ]; do
    while read -a line; do
	bin=${line[0]}
	
	./simfit_recoMC_singleComponent ${bin} ${par} ${tag} 1 1 ${year1} ${year2} ${year3}\
	    &> logs_simFit/simfit_recoMC_singleComponent_${bin}_${par}_${tag}_${year1}_${year2}_${year3}.out &

    done < ../confSF/KDE_SF.list
    tag=$(( ${tag} - 1 ))
done
