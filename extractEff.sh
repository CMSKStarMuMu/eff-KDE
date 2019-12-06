#!/bin/bash

# Same configuration as in sub_composeEff_rooKeys.sh
par=1
xbin=50
ybin=50
zbin=50
year=2016

while read -a line; do

    bin=${line[0]}
    wid0=${line[1]}
    wid1=${line[2]}
    wid2=${line[3]}
    wid3=${line[4]}
    wid4=${line[5]}
    wid5=${line[6]}

    # Compose efficiency descriptions using KDE histograms tagged with scale-factor configuration
    # The first set refers to wrong-tag recoNum, the others to the other terms
    root -b -q 'extractEff.cc('${bin}','${par}','\
${wid3}','${wid4}','${wid5}','\
${wid0}','${wid1}','${wid2}','\
${wid0}','${wid1}','${wid2}','\
${wid0}','${wid1}','${wid2}','\
${wid0}','${wid1}','${wid2}','\
${xbin}','${ybin}','${zbin}','\
${year}')'

done < ../confSF/KDE_SF.list
