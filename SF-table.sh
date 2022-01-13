#!/bin/bash

par=0

# Number of steps in the grid used to sample the KDE functions
xbin=50
ybin=50
zbin=50

while read -a line; do
	
    bin=${line[0]}
    wid0=${line[1]}
    wid1=${line[2]}
    wid2=${line[3]}

    echo -n "${bin}"

    for iVar in {0..2}
    do
	ls logs_test/composeEff_rooKeys_${bin}_*_${par}_${wid0}_${wid1}_${wid2}_${xbin}_${ybin}_${zbin}_*.out | grep -v "s_${bin}_4" | xargs cat |  grep "j: ${iVar}" | sed 's!.*weight=!!g' | sort | awk '{a[NR]=$0} END{printf(" & %.2f - %.2f", a[1], a[NR])}'
    done

    wid0=${line[4]}
    wid1=${line[5]}
    wid2=${line[6]}

    for iVar in {0..2}
    do
	grep "j: ${iVar}" logs_test/composeEff_rooKeys_${bin}_4_${par}_${wid0}_${wid1}_${wid2}_${xbin}_${ybin}_${zbin}_*.out | sed 's!.*weight=!!g' | sort | awk '{a[NR]=$0} END{printf(" & %.2f - %.2f", a[1], a[NR])}'
    done
    echo "\\\\"
	
done < ../confSF/KDE_SF.list #File containing scale-factor configuration
