#!/bin/bash
#run_test.sh
# A Geant4 simulation was performed with all the physics libraries switched off
# apart from electron ionisation energy loss. A beam of 100,000 50MeV (KE) e- 
# was injected into a 100 um target, and the momenta of electrons escaping the
# rear was output and saved to G4_100um_Al_50MeV.mat. Note that this contains
# over 100,000 e- due to delta-ray production. This script runs an equivalent 
# simulation in EPOCH, and compares the resulting angular and energy spectra. 

inputDir=`pwd`
coreNo=4
../../quick3d.sh $inputDir $coreNo
matlab -nosplash -nodesktop -r 'plot_G4_delta; quit;'


