#!/bin/bash

par=1

# Create directories for fit logs, results and plots
if [ ! -d logs_fit ]; then mkdir logs_fit; fi
if [ ! -d fitResults ]; then mkdir fitResults; fi
if [ ! -d plotFit_d ]; then mkdir plotFit_d; fi

# Compile dictionary and macro
make AngDict
make fit_recoMC_fullAngular

while read -a line; do
    bin=${line[0]}
    
    ./fit_recoMC_fullAngular ${bin} ${par} 0 1 \
	&> logs_fit/fit_recoMC_fullAngular_${bin}_${par}.out &
    
done < ../confSF/KDE_SF.list
