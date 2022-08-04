#!/bin/bash

#SBATCH -p lipq

seed=$((${SLURM_ARRAY_TASK_ID} + 1))

echo "bin=${bin},par=${par},tag=${tag},seed=${seed},cnt_hit=${cnt_hit},year=${year},vers=${vers}"

HOME=/lstore/cms/boletti/Run2-BdToKstarMuMu/eff-KDE-theta
SAMPLEDIR=/lstore/cms/boletti/Run2-BdToKstarMuMu/eff-KDE-theta

module load root/6.12.06

if [ ${par} -eq 0 ]
then
    if [ ! -r $SAMPLEDIR/files/KDEeff_b${bin}_ev_${year}_v${vers}.root ]; then
	echo $SAMPLEDIR/files/KDEeff_b${bin}_ev_${year}_v${vers}.root not found
	exit 1
    fi
else
    if [ ! -r $SAMPLEDIR/files/KDEeff_b${bin}_od_${year}_v${vers}.root ]; then
	echo $SAMPLEDIR/files/KDEeff_b${bin}_od_${year}_v${vers}.root not found
	exit 1
    fi
fi

if [ ! -r $HOME/preComp_Integrals_MC ]; then
    echo $HOME/preComp_Integrals_MC not found
    exit 1
fi

echo "$HOME/preComp_Integrals_MC ${bin} ${par} ${tag} ${cnt_hit} ${seed} ${year} ${vers}" \
    &> $HOME/logs_preComp_v${vers}/preComp_Integrals_MC_${bin}_${par}_${tag}_${cnt_hit}_${seed}_${year}.out

$HOME/preComp_Integrals_MC ${bin} ${par} ${tag} ${cnt_hit} ${seed} ${year} ${vers} \
    &>> $HOME/logs_preComp_v${vers}/preComp_Integrals_MC_${bin}_${par}_${tag}_${cnt_hit}_${seed}_${year}.out
