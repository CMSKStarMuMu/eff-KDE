#!/bin/bash

par=1
closure=1
year=2016

# Create directories for logs and plots
if [ ! -d logs_plot ]; then mkdir logs_plot; fi
if [ ! -d plotEff_d ]; then mkdir plotEff_d; fi

# Compile macro
make plotEff

while read -a line; do
    bin=${line[0]}
    
    ./plotEff ${bin} ${par} ${closure} ${year} \
	&> logs_plot/plotEff_${bin}_${par}_${year}.out &
    
done < ../confSF/KDE_SF.list
