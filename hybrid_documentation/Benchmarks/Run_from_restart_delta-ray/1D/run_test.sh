#!/bin/bash
#run_test.sh
# This script tests the ability to reload delta-ray optical depths
#
# input_first.deck and input_second.deck should be run back to back, and should
# produce the same results as input_full.deck being run alone.
#
# The test used is identical to the MATLAB equilibration test.

inputDir=`pwd`
coreNo=4

# Run the first half of the test
cp input_first.deck input.deck
../../quick1d.sh $inputDir $coreNo

# Run the second half
cp input_second.deck input.deck
../../quick1d.sh $inputDir $coreNo

# Inspect data
matlab -nosplash -nodesktop -r 'plot_G4_delta; quit;'
mv Energy_spectrum.fig restart_plot1.fig
mv Energy_spectrum.png restart_plot1.png
mv Angular_dist.fig restart_plot2.fig
mv Angular_dist.png restart_plot2.png

# Run without restart
rm *.sdf
cp input_full.deck input.deck
../../quick1d.sh $inputDir $coreNo
matlab -nosplash -nodesktop -r 'plot_G4_delta; quit;'
mv Energy_spectrum.fig full_plot1.fig
mv Energy_spectrum.png full_plot1.png
mv Angular_dist.fig full_plot2.fig
mv Angular_dist.png full_plot2.png

# Compare figures
matlab -nosplash -nodesktop -r 'combine_figs; quit;'

