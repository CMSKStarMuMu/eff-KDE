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

	for toy in {0..9}
	do
	    root -q -b 'mergeParSub_preComp_toyIntegrals_MC.cc('${bin}','${par}','${tag}',5,1,'${year}','${toy}','${vers}')'
	done
	
    done < ../confSF/KDE_SF.list
    tag=$(( ${tag} - 1 ))
done
