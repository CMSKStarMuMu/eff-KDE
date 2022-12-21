#!/bin/bash

par=${1}
[ -z "${par}" ] && par=1	# set default
year=${2}
[ -z "${year}" ] && year=2016	# set default

tag=1

while [ ${tag} -ge 0 ]; do
    while read -a line; do
	
	bin=${line[0]}
	vers=${line[7]}

	if [ -r files/KDEeffVers_${year}_${bin}_${par}_${vers}.log ]; then
	    if [ -r versfiles/sub_preComp_${year}_${bin}_${par}_${vers}.log ]; then
		if [ "$(cat files/KDEeffVers_${year}_${bin}_${par}_${vers}.log)" != "$(cat versfiles/sub_preComp_${year}_${bin}_${par}_${vers}.log)" ]; then
		    echo "Error, version files don't match: files/KDEeffVers_${year}_${bin}_${par}_${vers}.log and versfiles/sub_preComp_${year}_${bin}_${par}_${vers}.log"
		    continue
		fi
	    else
		echo "Warning, version file not found: versfiles/sub_preComp_${year}_${bin}_${par}_${vers}.log"
	    fi
	else
	    echo "Warning, version file not found: files/KDEeffVers_${year}_${bin}_${par}_${vers}.log"
	fi

	root -q -b 'mergeParSub_preComp_Integrals_MC.cc('${bin}','${par}','${tag}',50,1,'${year}','${vers}')'
	
    done < ../confSF/KDE_SF.list
    tag=$(( ${tag} - 1 ))
done
