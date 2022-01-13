#!/bin/bash

par=${1}
[ -z "${par}" ] && par=1	# set default
year=${2}
[ -z "${year}" ] && year=2016	# set default

tag=1

while [ ${tag} -ge 0 ]; do
    while read -a line; do
	
	bin=${line[0]}
	
	root -q -b 'mergeParSub_preComp_Integrals_MC_Swave.cc('${bin}','${par}','${tag}',500,1,'${year}')'
	
    done < ../confSF/KDE_SF.list
    tag=$(( ${tag} - 1 ))
done
