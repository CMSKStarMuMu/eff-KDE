#!/bin/bash

export HOME=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/eff-KDE
export CMSSWDIR=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/CMSSW_10_4_0/src
# export HOME=/eos/user/a/aboletti/BdToKstarMuMu/efficiency
# export CMSSWDIR=/eos/user/a/aboletti/BdToKstarMuMu/CMSSW_10_4_0/src
export WORKDIR=$PWD

cd $CMSSWDIR
source  /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`

cd $WORKDIR

bin=${1}
par=${2}
tag=${3}
cnt_hit=${4}
seed=${5}

if [ ${par} -eq 0 ]; then cp $HOME/files/KDEeff_b${bin}_ev.root $WORKDIR ; fi
if [ ${par} -eq 1 ]; then cp $HOME/files/KDEeff_b${bin}_od.root $WORKDIR ; fi
cp $HOME/preComp_Integrals_MC .

./preComp_Integrals_MC ${bin} ${par} ${tag} ${cnt_hit} ${seed}

cp PreIntMC_* $HOME/
rm PreIntMC_*
rm KDEeff_b*
rm preComp_Integrals_MC
