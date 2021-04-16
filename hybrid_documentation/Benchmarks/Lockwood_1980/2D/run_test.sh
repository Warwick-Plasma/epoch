#!/bin/bash
#run_test.sh
# This script reproduces the figures from Lockwood 1980 (report pg. 100 & 101,
# PDF pg. 92, 93), which shows a depth/dose curve for 0.5 MeV electrons in Ta

inputDir=`pwd`
coreNo=4

# Run Davies first
cp input_Davies.deck input.deck
../../quick2d.sh $inputDir $coreNo
matlab -nosplash -nodesktop -r 'plot_lockwood; quit;'
mv Lockwood_plot.fig Davies_plot.fig
mv Lockwood_plot.png Davies_plot.png

# Run Urban second
cp input_Urban.deck input.deck
../../quick2d.sh $inputDir $coreNo
matlab -nosplash -nodesktop -r 'plot_lockwood; quit;'
mv Lockwood_plot.fig Urban_plot.fig
mv Lockwood_plot.png Urban_plot.png

# Final formatting
matlab -nosplash -nodesktop -r 'combine_figs; quit;'

