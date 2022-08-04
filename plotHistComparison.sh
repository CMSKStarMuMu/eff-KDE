#!/bin/bash

par=${1}
[ -z "${par}" ] && par=1	# set default
year=${2}
[ -z "${year}" ] && year=2016	# set default
vers=${3}
[ -z "${vers}" ] && vers=0	# set default
qbin=${4}
[ -z "${qbin}" ] && qbin=-1	# set default

# Create directories
if [ ! -d logs_plot ]; then mkdir logs_plot; fi
if [ ! -d plotHist_d ]; then mkdir plotHist_d; fi

if [ $qbin == "-1" ]
then

    while read -a line; do
    
	bin=${line[0]}

	# [ "${line[7]}" != "${vers}" ] && continue
	# [ $bin != 4 ] && continue

	for indx in {3..4}; do

	    # Plot KDE histograms
	    root -b -q 'plotHistComparison.cc('${bin}','${indx}','${par}','${year}','${vers}')' \
		&> logs_plot/plotHistComparison_${bin}_${indx}_${par}_${year}_${vers}.out &
	
	done

    done < ../confSF/KDE_SF.list

else

    bin=$qbin

    for indx in {3..4}; do

	# Plot KDE histograms
	root -b -q 'plotHistComparison.cc('${bin}','${indx}','${par}','${year}','${vers}')' \
	    &> logs_plot/plotHistComparison_${bin}_${indx}_${par}_${year}_${vers}.out &
	
    done

fi
