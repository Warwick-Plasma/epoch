#!/bin/bash
#run_test.sh
# We have written an input deck which fires two bunches of different energies.
# One bunch has energy over the user-defined TNSA energy cutoff, and the other
# should reflect losing 50% of it's energy. This script runs the simulation
# along with a quick MATLAB analysis to ensure the boundaries are working as
# expected.

inputDir=`pwd`
coreNo=4
../../quick2d.sh $inputDir $coreNo
matlab -nosplash -nodesktop -r 'plot_matlab_test; quit;'


