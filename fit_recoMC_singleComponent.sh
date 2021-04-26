#!/bin/bash

par=1
tag=1
year=2016

# Create directories for fit logs, results and plots
if [ ! -d logs_fit ]; then mkdir logs_fit; fi
if [ ! -d fitResults ]; then mkdir fitResults; fi
if [ ! -d plotFit_d ]; then mkdir plotFit_d; fi

# Compile dictionary and macro
make AngDict
make fit_recoMC_singleComponent

while [ ${tag} -ge 0 ]; do
    while read -a line; do
	bin=${line[0]}
	
	./fit_recoMC_singleComponent ${bin} ${par} ${tag} 1 1 ${year}\
	    &> logs_fit/fit_recoMC_singleComponent_${bin}_${par}_${tag}_${year}.out &

    done < ../confSF/KDE_SF.list
    tag=$(( ${tag} - 1 ))
done
