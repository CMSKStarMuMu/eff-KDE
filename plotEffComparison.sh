#!/bin/bash

par=${1}
[ -z "${par}" ] && par=1	# set default
year=${2}
[ -z "${year}" ] && year=2016	# set default
vers=${3}
[ -z "${vers}" ] && vers=-1	# set default

closure=1

# Create directories for logs and plots
if [ ! -d logs_plot ]; then mkdir logs_plot; fi
if [ ! -d plotEff_d ]; then mkdir plotEff_d; fi

# Compile macro
if make plotEffComparison
then
    
    while read -a line; do
	bin=${line[0]}

	# [ "${line[7]}" != "${vers}" ] && continue
	# [ $bin != 0 ] && continue
    
	./plotEffComparison ${bin} ${par} ${closure} ${year} ${vers} \
	    &> logs_plot/plotEffComparison_${bin}_${par}_${year}_${vers}.out &
    
    done < ../confSF/KDE_SF.list

fi
