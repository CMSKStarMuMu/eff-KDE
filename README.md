# eff-KDE
Set of root macros used to create a description of selection efficiency, by means of Kernel Density Estimators 
and perform maximum likelihood fits with them

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
git clone -b working-fitWorkflow git@github.com:CMSKStarMuMu/eff-KDE.git
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
Adapt the paths to CMSSW area, working directory and directory for output files in [run_composeEff_rooKeys.sh](run_composeEff_rooKeys.sh#L3-L4), and in [mergeParSub_rooKeys.cc](mergeParSub_rooKeys.cc).
Here 50 jobs are submitted per bin, per parity, per efficiency term (5400 jobs with current q2 binning):
```sh
mkdir logs_parSub
condor_submit sub_composeEff_rooKeys.sub
```
when all the jobs have finished (you can check with `condor_q`) you can merge them:
```sh
mkdir files
root -q -b 'mergeParSub_rooKeys.cc'
```

## Compose efficiency histograms
Compose the numerators and denominators to create efficiency descriptions:
```sh
root -q -b 'extractEff.cc'
```
and find the efficiency histograms in the `KDEeff_b*_*.root` files.

## Run partial-integral numeric computation
This code is configured to submit parallel computation of the KDE, using HTCondor (available at CERN and accessible from lxplus machines).
Feel free to adapt the code to run on other kind of infrastructure or, discouraged, to run it locally.
Adapt the paths to CMSSW area and working directory in [run_preComp_Integrals_MC.sh](run_preComp_Integrals_MC.sh#L3-L6).
```sh
mkdir logs_preComp
make preComp_Integrals_MC
condor_submit sub_preComp_Integrals_MC.sub
```
when all the jobs have finished (you can check with `condor_q`) you can merge them:
```sh
root -b -q 'mergeParSub_preComp_Integrals_MC.cc'
```

## Fit generator-level distributions
Compile and run with:
```sh
mkdir fitResults
mkdir plotFit_d
make AngDict
make fit_genMC
./fit_genMC
```
This produces a root file `fitResults/fitResult_genMC.root` containing the RooFitResult objects, and fit projection plots in `plotFit_d/fitResult_genMC_*.pdf`.

## Fit single flavour-tagged components of post-selection distributions
Compile and run with:
```sh
mkdir fitResults
mkdir plotFit_d
make AngDict
make fit_recoMC_singleComponent
./fit_recoMC_singleComponent
```
This produces a root file `fitResults/fitResult_recoMC_singleComponent.root` containing the RooFitResult objects, and fit projection plots in `plotFit_d/fitResult_recoMC_singleComponent_*.pdf`.

## Fit post-selection distributions
Compile and run with:
```sh
mkdir fitResults
mkdir plotFit_d
make AngDict
make fit_recoMC_fullAngular
./fit_recoMC_fullAngular
```
This produces a root file `fitResults/fitResult_recoMC_fullAngular.root` containing the RooFitResult objects, and fit projection plots in `plotFit_d/fitResult_recoMC_fullAngular_*.pdf`.

## Plot and compare fit results
```sh
root -b -q 'plotFitResults.cc(0)' # for fit with even efficiency on odd dataset
root -b -q 'plotFitResults.cc(1)' # for fit with odd efficiency on even dataset
```
This produces one plot for each parameter in `plotFit_d/fitResult_*.pdf
