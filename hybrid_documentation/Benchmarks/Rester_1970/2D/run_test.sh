#!/bin/bash
#run_test.sh
# This script reproduces the 2.8 MeV Au photon spectrum from Rester (1970) Fig.
# 17. When plotting, we compare against the experimental results, and find that
# the low energy photon spectrum is wrong (photoelectric effect has been
# neglected, which fails to attenuate low energy photons). The results are also
# compared to Geant4 with the photoelectric effect disabled, showing better
# agreement with EPOCH

inputDir=`pwd`
coreNo=4
../../quick2d.sh $inputDir $coreNo
matlab -nosplash -nodesktop -r 'plot_rester; quit;'


