# MicroBooNE_NeutronStudy_Reco

source setup.sh


## Under Constants.h, Choose to apply(not apply) PMissingCut by setting
bool PMissingCut = true (false);


## Choose Event File to run on:
emacs NeutrinoSelectionFilter.h
### ctrl+s: Open
### Change TFile and TDirectory


## Run the Neutron Selection Filter
root -l script_MicroBooNE_Reco.C

## Run the 1D Plot:
### Choose whether a stacked (unstacked) histogram
### Under PlotRoot.cpp
bool stackedHist = true (false);

root -l PlotRoot.cpp


## Run the 2D Plot:
root -l PlotRoot2D.cpp



# Using High Statistics Version
cd HighStats
### All of the scripts run in the same fashion, but must access Constants.h file through ../Constants.h

