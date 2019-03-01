#!/bin/bash

cd /afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/CMSSW_10_4_0/src
source  /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`

cd /afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/eff-KDE/
pwd

bin=${1}
tag=${2}
kern=${3}
xbin=${4}
ybin=${5}
zbin=${6}
ndiv=${7}
divx=${8}
divy=${9}
divz=${10}

root -q -b 'composeEff_customKDE_uniform_parSub.cc('${bin}','${tag}','${kern}','${xbin}','${ybin}','${zbin}','${ndiv}','${divx}','${divy}','${divz}')' \
    &>logs/log_run_composeEff_customKDE_uniform_${bin}_${tag}_${kern}_${xbin}_${ybin}_${zbin}_${ndiv}_${divx}_${divy}_${divz}.out
