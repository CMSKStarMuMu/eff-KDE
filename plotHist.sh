#!/bin/bash

par=1

# Create directories
if [ ! -d logs_plot ]; then mkdir logs_plot; fi
if [ ! -d plotHist_d ]; then mkdir plotHist_d; fi

while read -a line; do
    
    bin=${line[0]}
    
    for indx in {0..0}; do

	# Plot KDE histograms
	root -b -q 'plotHist.cc('${bin}','${indx}','${par}')' \
	    &> logs_plot/plotHist_${bin}_${indx}_${par}.out &
	
    done

done < ../confSF/KDE_SF.list
