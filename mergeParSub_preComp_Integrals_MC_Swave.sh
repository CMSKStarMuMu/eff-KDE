#!/bin/bash

par=1
tag=1
year=${1}

while [ ${tag} -ge 0 ]; do
    while read -a line; do
	
	bin=${line[0]}
	
	root -q -b 'mergeParSub_preComp_Integrals_MC_Swave.cc('${bin}','${par}','${tag}',100,1,'${year}')'
	
    done < ../confSF/KDE_SF.list
    tag=$(( ${tag} - 1 ))
done
