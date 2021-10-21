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
git clone -b working-basicMacroImplementation git@github.com:CMSKStarMuMu/eff-KDE.git
cd eff-KDE
```
## Create datasets
If needed, change the [location of the ntuples](createDataset.cc#L49-L60), which need to be produced with the code in the [B0KstMuMuNtuple repository](https://github.com/CMSKStarMuMu/B0KstMuMuNtuple).
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
1       0.30    1.00    0.70    0.30    0.80    0.30
2       0.40    0.80    0.50    0.40    0.80    0.50
3       0.60    1.00    0.80    0.40    1.00    0.50
5       0.40    1.00    0.50    0.40    1.20    1.10
7       0.60    1.00    0.60    0.40    0.60    0.40
EOF
```
The KDE scale factors that have been used for the normalisation channels are
```
4       0.60    1.00    0.80    0.40    1.00    0.50
6       0.60    1.00    0.60    0.40    1.00    1.00
```
Adapt the paths to CMSSW area, working directory and directory for output files in [run_composeEff_rooKeys.sh](run_composeEff_rooKeys.sh#L16-L19).
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

## Plot projections and slices
To plot the individual numerator and denominator terms of the final efficiency, and the derived KDE description, use
```
source plotHist.sh
```
If necessary, change the [input dataset location](plotHist.cc#L60).


To plot the efficiency slices and projections and the derived KDE description and to run the efficiency closure test, use
```
source plotEff.sh
```
If necessary, change the [input dataset location](plotEff.cc#L52).

## Run partial-integral numeric computation
This code is configured to submit parallel computation of the KDE, using HTCondor (available at CERN and accessible from lxplus machines).
Feel free to adapt the code to run on other kind of infrastructure or, discouraged, to run it locally.
Adapt the paths to CMSSW area and working directory in [run_preComp_Integrals_MC.sh](run_preComp_Integrals_MC.sh#L3-L6), then submit the jobs with:
```sh
source sub_preComp_Integrals_MC.sh
```
when all the jobs have finished (you can check with `condor_q`) you can merge them:
```sh
source mergeParSub_preComp_Integrals_MC.sh
```

## Fit generator-level distributions
Compile and run with:
```sh
source fit_genMC.sh
```
This produces a root file `fitResults/fitResult_genMC.root` containing the RooFitResult objects, and fit projection plots in `plotFit_d/fitResult_genMC_*.pdf`.

## Fit single flavour-tagged components of post-selection distributions
Compile and run with:
```sh
source fit_recoMC_singleComponent.sh
```
This produces a root file `fitResults/fitResult_recoMC_singleComponent.root` containing the RooFitResult objects, and fit projection plots in `plotFit_d/fitResult_recoMC_singleComponent_*.pdf`.

## Fit post-selection distributions
Compile and run with:
```sh
source fit_recoMC_fullAngular.sh
```
This produces a root file `fitResults/fitResult_recoMC_fullAngular.root` containing the RooFitResult objects, and fit projection plots in `plotFit_d/fitResult_recoMC_fullAngular_*.pdf`.

## Plot and compare fit results
```sh
root -b -q 'plotFitResults.cc(0)' # for fit with even efficiency on odd dataset
root -b -q 'plotFitResults.cc(1)' # for fit with odd efficiency on even dataset
```
This produces one plot for each parameter in `plotFit_d/fitResult_*.pdf`
If necessary, change the [input dataset location](plotEff.cc#L52).
