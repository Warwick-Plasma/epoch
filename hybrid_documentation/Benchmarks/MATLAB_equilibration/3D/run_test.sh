#!/bin/bash
#run_test.sh
# We have written equilibration scripts in both MATLAB and EPOCH. This script 
# tests to see if both codes agree with each other. This is a very basic test, 
# which only tests if I have been consistent with myself - this is not a
# a rigorous benchmark against the physics

inputDir=`pwd`
coreNo=4
../../quick3d.sh $inputDir $coreNo
matlab -nosplash -nodesktop -r 'plot_matlab_test; quit;'


