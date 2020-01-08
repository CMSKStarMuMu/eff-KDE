#!/bin/bash

par=1
tag=1
year=2016

while [ ${tag} -ge 0 ]; do
    while read -a line; do
	
	bin=${line[0]}
	
	root -q -b 'mergeParSub_preComp_Integrals_MC.cc('${bin}','${par}','${tag}',50,1,'${year}')'
	
    done < ../confSF/KDE_SF.list
    tag=$(( ${tag} - 1 ))
done
