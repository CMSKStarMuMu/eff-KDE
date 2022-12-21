#!/bin/bash

# Same configuration as in sub_composeEff_rooKeys.sh
par=${1}
[ -z "${par}" ] && par=1	# set default
year=${2}
[ -z "${year}" ] && year=2016	# set default

xbin=50
ybin=50
zbin=50

while read -a line; do

    bin=${line[0]}
    wid0=${line[1]}
    wid1=${line[2]}
    wid2=${line[3]}
    wid3=${line[4]}
    wid4=${line[5]}
    wid5=${line[6]}
    vers=${line[7]}

    # Compose efficiency descriptions using KDE histograms tagged with scale-factor configuration
    # The first set refers to wrong-tag recoNum, the others to the other terms
    root -b -q 'extractEff.cc('${bin}','${par}','\
${wid3}','${wid4}','${wid5}','\
${wid0}','${wid1}','${wid2}','\
${wid0}','${wid1}','${wid2}','\
${wid0}','${wid1}','${wid2}','\
${wid0}','${wid1}','${wid2}','\
${xbin}','${ybin}','${zbin}','\
${year}','${vers}')'

    if [ -r versfiles/merge_composeEff_${year}_${bin}_${par}_${vers}.log ]; then
	cp  versfiles/merge_composeEff_${year}_${bin}_${par}_${vers}.log files/KDEeffVers_${year}_${bin}_${par}_${vers}.log
    else
	echo "Warning, version file not found: versfiles/merge_composeEff_${year}_${bin}_${par}_${vers}.log"
    fi

done < ../confSF/KDE_SF.list
