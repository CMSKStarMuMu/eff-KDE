#!/bin/bash

if [ "${USER}" == "fiorendi" ]; then
    export HOME=/afs/cern.ch/work/f/fiorendi/private/effKDE/eff-KDE
    export CMSSWDIR=/afs/cern.ch/work/f/fiorendi/private/effKDE/CMSSW_10_4_0/src
elif [ "${USER}" == "aboletti" ]; then
    export HOME=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/eff-KDE-Swave
    export SAMPLEDIR=/eos/cms/store/user/fiorendi/p5prime/effKDE
    export CMSSWDIR=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/CMSSW_10_4_0/src
else
    echo no user found
fi

echo setting HOME to $HOME 
echo setting CMSSWDIR to $CMSSWDIR

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
year=${6}

if [ ${par} -eq 0 ]
then
    if [ ! -r $SAMPLEDIR/${year}/lmnr/newphi/KDEeff_b${bin}_ev_${year}.root ]; then
	echo $SAMPLEDIR/${year}/lmnr/newphi/KDEeff_b${bin}_ev_${year}.root not found
	exit 1
    fi
    cp $SAMPLEDIR/${year}/lmnr/newphi/KDEeff_b${bin}_ev_${year}.root $WORKDIR
else
    if [ ! -r $SAMPLEDIR/${year}/lmnr/newphi/KDEeff_b${bin}_od_${year}.root ]; then
	echo $SAMPLEDIR/${year}/lmnr/newphi/KDEeff_b${bin}_od_${year}.root not found
	exit 1
    fi
    cp $SAMPLEDIR/${year}/lmnr/newphi/KDEeff_b${bin}_od_${year}.root $WORKDIR
fi

if [ ! -r $HOME/preComp_Integrals_MC_Swave ]; then
    echo $HOME/preComp_Integrals_MC_Swave not found
    exit 1
fi
cp $HOME/preComp_Integrals_MC_Swave .

./preComp_Integrals_MC_Swave ${bin} ${par} ${tag} ${cnt_hit} ${seed} ${year}

if [ ! -d $HOME/tmp ]; then mkdir $HOME/tmp; fi
cp PreIntMC_* $HOME/tmp/

rm PreIntMC_*
rm KDEeff_b*
rm preComp_Integrals_MC_Swave