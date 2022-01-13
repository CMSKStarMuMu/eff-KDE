#!/bin/bash

par=${1}
[ -z "${par}" ] && par=1	# set default
year=${2}
[ -z "${year}" ] && year=2016	# set default

closure=1

# Create directories for logs and plots
if [ ! -d logs_plot ]; then mkdir logs_plot; fi
if [ ! -d plotToyEff_d ]; then mkdir plotToyEff_d; fi

# Compile macro
if make plotToyEff
then

    while read -a line; do
	bin=${line[0]}
    
	./plotToyEff ${bin} ${par} ${closure} ${year} \
	    &> logs_plot/plotToyEff_${bin}_${par}_${year}.out &
    
    done < ../confSF/KDE_SF.list

fi
