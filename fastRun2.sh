#!/bin/bash

# merge output of parralel jobs
root -l -q -b mergeParSub_rooKeys.cc

# plot KDE functions
mkdir plotHist_d
root -l -q -b plotHist.cc

# create efficiency from KDE functions
root -l -q -b extractEff.cc

# plot efficiency functions
make plotEff && ./plotEff
