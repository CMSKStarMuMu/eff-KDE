#!/bin/bash

#SBATCH -p lipq

ndiv=$SLURM_ARRAY_TASK_ID

echo "bin=${bin},indx=${indx},par=${par},wid0=${wid0},wid1=${wid1},wid2=${wid2},xbin=${xbin},ybin=${ybin},zbin=${zbin},ndiv=${ndiv},totdiv=${totdiv},year=${year},vers=${vers}"

HOME=/lstore/cms/boletti/Run2-BdToKstarMuMu/eff-KDE-theta
SAMPLEDIR=/lstore/cms/boletti/Run2-BdToKstarMuMu/eff-KDE-theta

module load root/6.12.06

if [ ! -r $SAMPLEDIR/effDatasetTheta_b${bin}_${year}.root ]; then
    echo $SAMPLEDIR/effDatasetTheta_b${bin}_${year}.root not found
    exit 1
fi
if [ ! -r $HOME/composeToyEff_rooKeys_parSub.cc ]; then
    echo $HOME/composeToyEff_rooKeys_parSub.cc not found
    exit 1
fi

[ -z "${vers}" ] && vers="-1"
if [ "${par}" == "0" ]; then
    parstr="ev"
else
    parstr="od"
fi
if [ ! -r $SAMPLEDIR/files/KDEhist_b${bin}_${parstr}_${year}_v${vers}.root ]; then
    echo $SAMPLEDIR/files/KDEhist_b${bin}_${parstr}_${year}_v${vers}.root not found
    exit 1
fi
if [ ! -d $SAMPLEDIR/tmptoy_v${vers} ]; then mkdir $SAMPLEDIR/tmptoy_v${vers}; fi

echo "$HOME/composeToyEff_rooKeys_parSub $bin $indx $par $wid0 $wid1 $wid2 $xbin $ybin $zbin $ndiv $totdiv $year $vers" \
     &> $HOME/logs_parSub_v${vers}/composeToyEff_rooKeys_parSub_${bin}_${indx}_${par}_${wid0}_${wid1}_${wid2}_${xbin}_${ybin}_${zbin}_${ndiv}_${totdiv}_${year}.out

$HOME/composeToyEff_rooKeys_parSub ${bin} ${indx} ${par} ${wid0} ${wid1} ${wid2} ${xbin} ${ybin} ${zbin} ${ndiv} ${totdiv} ${year} ${vers} \
    &>> $HOME/logs_parSub_v${vers}/composeToyEff_rooKeys_parSub_${bin}_${indx}_${par}_${wid0}_${wid1}_${wid2}_${xbin}_${ybin}_${zbin}_${ndiv}_${totdiv}_${year}.out

echo "run end"
