#!/bin/bash

par=${1}
[ -z "${par}" ] && par=1	# set default
year=${2}
[ -z "${year}" ] && year=2016	# set default
vers=${3}
[ -z "${vers}" ] && vers=0	# set default

# Create directories
if [ ! -d logs_plot ]; then mkdir logs_plot; fi
if [ ! -d plotHist_d ]; then mkdir plotHist_d; fi

while read -a line; do
    
    bin=${line[0]}

    [ "${line[7]}" != "${vers}" ] && continue

    for indx in {0..4}; do

	# Plot KDE histograms
	root -b -q 'plotHist.cc('${bin}','${indx}','${par}','${year}','${vers}')' \
	    &> logs_plot/plotHist_${bin}_${indx}_${par}_${year}_${vers}.out &
	
    done

done < ../confSF/KDE_SF.list
