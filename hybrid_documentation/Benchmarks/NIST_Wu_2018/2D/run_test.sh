#!/bin/bash
#run_test.sh
# This script reproduces Fig 4 in Wu for electrons on an Al target. This shows
# the total stopping power due to both radiative (bremsstrahlung) and
# collisional (ionisation) losses for e- between 10 keV and 1 GeV. This is
# compared to the NIST data, stored in Nist_dat.mat.

inputDir=`pwd`
coreNo=4

# Run Davies first
../../quick2d.sh $inputDir $coreNo
matlab -nosplash -nodesktop -r 'plot_NIST; quit;'
rm *.sdf


