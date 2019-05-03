#!/bin/bash

cd /afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/CMSSW_10_4_0/src
source  /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`

cd /afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/eff-KDE/
pwd

bin=${1}
tag=${2}
kern=${3}
doMi=${4}
xbin=${5}
ybin=${6}
zbin=${7}
ndiv=${8}
divx=${9}
divy=${10}
divz=${11}

root -q -b 'composeEff_customKDE_uniform_parSub.cc('${bin}','${tag}','${kern}','${doMi}','${xbin}','${ybin}','${zbin}','${ndiv}','${divx}','${divy}','${divz}')' \
    &>logs/log_run_composeEff_customKDE_uniform_${bin}_${tag}_${kern}_${doMi}_${xbin}_${ybin}_${zbin}_${ndiv}_${divx}_${divy}_${divz}.out
