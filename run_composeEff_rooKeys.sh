#!/bin/bash

cd /afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/CMSSW_10_4_0/src
source  /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`

cd /afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/eff-KDE/
pwd

bin=${1}
tag=${2}
xbin=${3}
ybin=${4}
zbin=${5}
wid=${6}
ndiv=${7}
totdiv=${8}

root -q -b 'composeEff_rooKeys_parSub.cc('${bin}','${tag}','${xbin}','${ybin}','${zbin}','${wid}','${ndiv}','${totdiv}')' \
    &>logs/log_run_composeEff_rooKeys_${bin}_${tag}_${xbin}_${ybin}_${zbin}_${wid}_${ndiv}_${totdiv}.out
