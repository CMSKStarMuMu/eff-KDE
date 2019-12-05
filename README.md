# eff-KDE
Set of root macros used to create a description of selection efficiency, by means of Kernel Density Estimators 

# Quick run
Use this list of commands to produce a result as quickly as possible

## Setup working area
Make sure to be in an area with ROOT v6.12 (or later) available. The code could work with previous versions, but it is not guaranteed.
If you are on CERN lxplus machines, you can achieve it by running:
```sh
export SCRAM_ARCH=slc7_amd64_gcc820
cmsrel CMSSW_10_4_0
cd CMSSW_10_4_0/src/ && cmsenv && cd ../..
```
Clone this branch in the working directory:
```sh
git clone -b working-basicMacroImplementation git@github.com:CMSKStarMuMu/eff-KDE.git
cd eff-KDE
```
## Create datasets
If needed, change the [location of the ntuples](createDataset.cc#L42-L55), which need to be produced with the code in the [B0KstMuMuNtuple repository](https://github.com/CMSKStarMuMu/B0KstMuMuNtuple).
Then produce files with all the needed datasets (example for 2017 ntuples):
```sh
root -q -b 'createDataset.cc(7)'
```

## Produce KDE description of numerators and denominators
This code is configured to submit parallel computation of the KDE, using HTCondor (available at CERN and accessible from lxplus machines).
Feel free to adapt the code to run on other kind of infrastructure or, discouraged, to run it locally.
Before running, you need to create a file `../confSF/KDE_SF.list`, containing the q2-bins to process and the corresponding KDE scale factors.
```sh
if [ ! -d ../confSF ]; then mkdir ../confSF; fi
cat << EOF > ../confSF/KDE_SF.list
0       0.50    1.00    1.00    0.50    1.00    1.00
1       0.30    1.00    0.70    0.30    1.00    0.70
2       0.40    0.80    0.50    0.40    0.80    0.50
3       0.60    1.00    0.80    0.60    1.00    0.80
5       0.40    1.00    0.50    0.40    1.00    0.50
7       0.60    1.00    0.60    0.60    1.00    0.60
EOF
```
Adapt the paths to CMSSW area, working directory and directory for output files in [run_composeEff_rooKeys.sh](run_composeEff_rooKeys.sh#L3-L5).
Submit 4200 jobs by running:
```sh
source sub_composeEff_rooKeys.sh
```
when all the jobs have finished (you can check with `condor_q`) you can merge them:
```sh
source mergeParSub_rooKeys.sh
```

## Compose efficiency histograms
Compose the numerators and denominators to create efficiency descriptions:
```sh
source extractEff.sh
```
and find the efficiency histograms in the `files/KDEeff_b*_*.root` files.
