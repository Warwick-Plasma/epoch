#!/bin/bash
#run_test.sh
# This script reproduces the 18.66 mg/cm² result from Hanson 1951, Figure 3,
# using the Davies algorithm. In this benchmark, we are testing the ability of
# the code to continue from a restart dump. One input deck is run for the first 
# half, and a second input deck runs from 0001.sdf.
#
# We then repeat the run without a restart, and compare the final figures. This
# verifies hybrid particles are correctly restarted

inputDir=`pwd`
coreNo=4

# Run up to 25 fs
cp input_first.deck input.deck
../../quick2d.sh $inputDir $coreNo

# Run remaining
cp input_second.deck input.deck
../../quick2d.sh $inputDir $coreNo

# Inspect data
matlab -nosplash -nodesktop -r 'plot_hanson; quit;'
mv Hanson_plot.fig restart_plot.fig
mv Hanson_plot.png restart_plot.png

# Run without restart
rm *.sdf
cp input_full.deck input.deck
../../quick2d.sh $inputDir $coreNo
matlab -nosplash -nodesktop -r 'plot_hanson; quit;'
mv Hanson_plot.fig full_plot.fig
mv Hanson_plot.png full_plot.png

# Final formatting
matlab -nosplash -nodesktop -r 'combine_figs; quit;'

