#!/bin/bash

rwv=${1}
bopt=${2}
yopt=${3}

if [ -z "$rwv" ]
then
    rwv=0
elif [ ! "$rwv" -ge "0" ]
then
    echo "Illegal reweighting version: ${rwv}"
    echo "Empty or non-negative values allowed"
    return 1
fi

if [ -z "$bopt" ]
then
    blist=("-1" "4" "6")
elif [ "$bopt" == "a" ]
then
    blist=("-1" "4" "6")
elif [ "$bopt" == "n" ]
then
    blist=("-1")
elif [ "$bopt" == "j" ]
then
    blist=("4")
elif [ "$bopt" == "p" ]
then
    blist=("6")
elif [ "$bopt" == "r" ]
then
    blist=("4" "6")
else
    echo "Illegal bin configuration: ${bopt}"
    echo 'Good values: "a" (all bins, default), "n" (non-resonant), "j" (J/psi), "p" (psi2S), "r" (both resonant bins)'
    return 1
fi

if [ -z "$yopt" -o "$yopt" == "a" ]
then
    ylist=("6" "7" "8")
elif [ "$yopt" -ge "6" -a "$yopt" -le "8" ]
then
    ylist=("$yopt")
else
    echo "Illegal year configuration: ${yopt}"
    echo 'Good values: "a" (all years, default), or a number from 6 to 8'
    return 1
fi

p="-1"
theta="false"
plot="false"

cd /eos/user/a/aboletti/BdToKstarMuMu/fileIndex/
if [ -z "$(git status --porcelain | grep "^[ ]*M")" ]
then
    ntu="$(git rev-parse --short HEAD)"
else
    echo "Not-commited changes in file-index repo. Abort!"
    return 1
fi
cd -

for b in "${blist[@]}"
do
    for y in "${ylist[@]}"
    do
	if [ ${b} == "-1" ]
	then
	    for ib in {0..7}
	    do
		if [ ${ib} != "4" -a ${ib} != "6" ]
		then
		    if [ ${p} == "-1" ]
		    then
			echo "$ntu" > "/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-Swave/effDatasetVers_201${y}_${ib}_0_${theta}_${rwv}.log"
			echo "$ntu" > "/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-Swave/effDatasetVers_201${y}_${ib}_1_${theta}_${rwv}.log"
		    else
			echo "$ntu" > "/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-Swave/effDatasetVers_201${y}_${ib}_${p}_${theta}_${rwv}.log"
		    fi
		fi
	    done
	else
	    if [ ${p} == "-1" ]
	    then
		echo "$ntu" > "/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-Swave/effDatasetVers_201${y}_${b}_0_${theta}_${rwv}.log"
		echo "$ntu" > "/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-Swave/effDatasetVers_201${y}_${b}_1_${theta}_${rwv}.log"
	    else
		echo "$ntu" > "/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-Swave/effDatasetVers_201${y}_${b}_${p}_${theta}_${rwv}.log"
	    fi
	fi
	echo "root -b -q -l 'createDataset.cc('${y}','${b}','${p}','${theta}','${rwv}','${plot}')' \
	    &> createDataset_${y}_${b}_${p}_${theta}_${rwv}.log &"
	root -b -q -l 'createDataset.cc('${y}','${b}','${p}','${theta}','${rwv}','${plot}')' \
	    &> createDataset_${y}_${b}_${p}_${theta}_${rwv}.log &
    done
done
