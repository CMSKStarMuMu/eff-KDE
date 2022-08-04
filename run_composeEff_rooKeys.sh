#!/bin/bash

#SBATCH -p lipq
#SBATCH -t 1200

ndiv=$SLURM_ARRAY_TASK_ID

echo "bin=${bin},indx=${indx},par=${par},wid0=${wid0},wid1=${wid1},wid2=${wid2},xbin=${xbin},ybin=${ybin},zbin=${zbin},ndiv=${ndiv},totdiv=${totdiv},year=${year},vers=${vers}"

HOME=/lstore/cms/boletti/Run2-BdToKstarMuMu/eff-KDE-theta
SAMPLEDIR=/lstore/cms/boletti/Run2-BdToKstarMuMu/eff-KDE-theta

module load root/6.12.06

if [ $vers -ge "10" ]; then
    XGBv="_XGBv$(($vers / 10))"
fi

if [ ! -r $SAMPLEDIR/effDatasetTheta_b${bin}_${year}${XGBv}.root ]; then
    echo $SAMPLEDIR/effDatasetTheta_b${bin}_${year}${XGBv}.root not found
    exit 1
fi
if [ ! -r $HOME/composeEff_rooKeys_parSub.cc ]; then
    echo $HOME/composeEff_rooKeys_parSub.cc not found
    exit 1
fi
[ -z "${vers}" ] && vers="-1"
if [ ! -d $SAMPLEDIR/tmp_v${vers} ]; then mkdir $SAMPLEDIR/tmp_v${vers}; fi

echo "$HOME/composeEff_rooKeys_parSub $bin $indx $par $wid0 $wid1 $wid2 $xbin $ybin $zbin $ndiv $totdiv $year $vers" \
     &> $HOME/logs_parSub_v${vers}/composeEff_rooKeys_parSub_${bin}_${indx}_${par}_${wid0}_${wid1}_${wid2}_${xbin}_${ybin}_${zbin}_${ndiv}_${totdiv}_${year}.out

$HOME/composeEff_rooKeys_parSub ${bin} ${indx} ${par} ${wid0} ${wid1} ${wid2} ${xbin} ${ybin} ${zbin} ${ndiv} ${totdiv} ${year} ${vers} \
    &>> $HOME/logs_parSub_v${vers}/composeEff_rooKeys_parSub_${bin}_${indx}_${par}_${wid0}_${wid1}_${wid2}_${xbin}_${ybin}_${zbin}_${ndiv}_${totdiv}_${year}.out
