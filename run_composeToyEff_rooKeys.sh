#!/bin/bash

bin=${1}
indx=${2}
par=${3}
seed=${4}
wid0=${5}
wid1=${6}
wid2=${7}
xbin=${8}
ybin=${9}
zbin=${10}
ndiv=${11}
totdiv=${12}
year=${13}

if [ "${USER}" == "fiorendi" ]; then
    export HOME=/afs/cern.ch/work/f/fiorendi/private/effKDE/eff-KDE
    export CMSSWDIR=/afs/cern.ch/work/f/fiorendi/private/effKDE/CMSSW_10_4_0/src
    export SAMPLEDIR=/eos/cms/store/user/fiorendi/p5prime/effKDE/${year}/lmnr/
elif [ "${USER}" == "aboletti" ]; then
    export HOME=/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-Swave
    export CMSSWDIR=/eos/user/a/aboletti/BdToKstarMuMu/CMSSW_10_4_0/src
    export SAMPLEDIR=/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-Swave
else
    echo no user found
    exit 1
fi

echo setting HOME to $HOME 
echo setting CMSSWDIR to $CMSSWDIR

export WORKDIR=$PWD
cd $CMSSWDIR
source  /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`

cd $WORKDIR

echo 'now submitting for bin ' ${bin}

if [ ! -r $SAMPLEDIR/effDataset_b${bin}_${year}.root ]; then
    echo $SAMPLEDIR/effDataset_b${bin}_${year}.root not found
    exit 1
fi
if [ ! -r $HOME/composeEff_rooKeys_parSub.cc ]; then
    echo $HOME/composeEff_rooKeys_parSub.cc not found
    exit 1
fi

cp $SAMPLEDIR/effDataset_b${bin}_${year}.root .
cp $HOME/composeToyEff_rooKeys_parSub .

if [ "${par}" == "0" ]; then
    cp $SAMPLEDIR/files/KDEhist_b${bin}_ev_${year}.root .
else
    cp $SAMPLEDIR/files/KDEhist_b${bin}_od_${year}.root .
fi

echo "./composeToyEff_rooKeys_parSub ${bin} ${indx} ${par} ${seed} ${wid0} ${wid1} ${wid2} ${xbin} ${ybin} ${zbin} ${ndiv} ${totdiv} ${year}"
./composeToyEff_rooKeys_parSub ${bin} ${indx} ${par} ${seed} ${wid0} ${wid1} ${wid2} ${xbin} ${ybin} ${zbin} ${ndiv} ${totdiv} ${year}

if [ ! -d $HOME/tmp_b${bin}_toy${seed} ]; then mkdir $HOME/tmp_b${bin}_toy${seed}; fi
cp KDEhist_*_rooKeys* $HOME/tmp_b${bin}_toy${seed}/

rm composeToyEff_rooKeys_parSub
rm effDataset_b*
rm KDEhist*
