#!/bin/bash

export HOME=/eos/user/a/aboletti/BdToKstarMuMu/efficiency
export CMSSWDIR=/eos/user/a/aboletti/BdToKstarMuMu/CMSSW_10_4_0/src
export WORKDIR=$PWD

cd $CMSSWDIR
source  /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`

cd $WORKDIR

bin=${1}
indx=${2}
par=${3}
wid=${4}
xbin=${5}
ybin=${6}
zbin=${7}
ndiv=${8}
totdiv=${9}

cp $HOME/effDataset_b${bin}.root .
cp $HOME/composeEff_rooKeys_parSub.cc .

root -l -q -b 'composeEff_rooKeys_parSub.cc('${bin}','${indx}','${par}','${wid}','${xbin}','${ybin}','${zbin}','${ndiv}','${totdiv}')'

cp KDEhist* $HOME/
