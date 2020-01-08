#!/bin/bash

bin=${1}
indx=${2}
par=${3}
wid0=${4}
wid1=${5}
wid2=${6}
xbin=${7}
ybin=${8}
zbin=${9}
ndiv=${10}
totdiv=${11}
year=${12}

if [ "${USER}" == "fiorendi" ]; then
    export HOME=/afs/cern.ch/work/f/fiorendi/private/effKDE/eff-KDE
    export CMSSWDIR=/afs/cern.ch/work/f/fiorendi/private/effKDE/CMSSW_10_4_0/src
    export SAMPLEDIR=/eos/cms/store/user/fiorendi/p5prime/effKDE/${year}/lmnr/
elif [ "${USER}" == "aboletti" ]; then
    export HOME=/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-parSub
    export CMSSWDIR=/eos/user/a/aboletti/BdToKstarMuMu/CMSSW_10_4_0/src
    export SAMPLEDIR=/eos/user/a/aboletti/BdToKstarMuMu/datasets/PUweight
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
cp $HOME/composeEff_rooKeys_parSub.cc .

echo 'root -l -q -b composeEff_rooKeys_parSub.cc($bin,$indx,$par,$wid0,$wid1,$wid2,$xbin,$ybin,$zbin,$ndiv,$totdiv,$year)'
root -l -q -b 'composeEff_rooKeys_parSub.cc('${bin}','${indx}','${par}','${wid0}','${wid1}','${wid2}','${xbin}','${ybin}','${zbin}','${ndiv}','${totdiv}','${year}')'

if [ ! -d $HOME/tmp ]; then mkdir $HOME/tmp; fi
cp KDEhist* $HOME/tmp/

rm composeEff_rooKeys_parSub.cc
rm effDataset_b*
rm KDEhist*
