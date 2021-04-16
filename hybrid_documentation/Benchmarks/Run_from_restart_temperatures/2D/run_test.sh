#!/bin/bash
#run_test.sh
# This script tests the ability to reload temperatures from a restart dump.
#
# input_first.deck and input_second.deck should be run back to back, and should
# produce the same results as input_full.deck being run alone.
#
# The test used is identical to the MATLAB equilibration test.

inputDir=`pwd`
coreNo=4

# Run the first half of the test
cp input_first.deck input.deck
../../quick2d.sh $inputDir $coreNo

# Run the second half
cp input_second.deck input.deck
../../quick2d.sh $inputDir $coreNo

# Inspect data
matlab -nosplash -nodesktop -r 'plot_matlab_test; quit;'
mv Temperature_evolution.fig restart_plot.fig
mv Temperature_evolution.png restart_plot.png

# Run without restart
rm *.sdf
cp input_full.deck input.deck
../../quick2d.sh $inputDir $coreNo
matlab -nosplash -nodesktop -r 'plot_matlab_test; quit;'
mv Temperature_evolution.fig full_plot.fig
mv Temperature_evolution.png full_plot.png

# Compare figures
matlab -nosplash -nodesktop -r 'combine_figs; quit;'

