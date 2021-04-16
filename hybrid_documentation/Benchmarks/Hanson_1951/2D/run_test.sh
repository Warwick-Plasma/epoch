#!/bin/bash
#run_test.sh
# This script reproduces the 18.66 mg/cm² result from Hanson 1951, Figure 3.
# Two simulations are performed, one using the Davies elastic scatter
# algorithm, the other using Geant4's Urban algorithm.

inputDir=`pwd`
coreNo=4

# Run Davies first
cp input_Davies.deck input.deck
../../quick2d.sh $inputDir $coreNo
matlab -nosplash -nodesktop -r 'plot_hanson; quit;'
mv Hanson_plot.fig Davies_plot.fig
mv Hanson_plot.png Davies_plot.png

# Run Urban second
cp input_Urban.deck input.deck
../../quick2d.sh $inputDir $coreNo
matlab -nosplash -nodesktop -r 'plot_hanson; quit;'
mv Hanson_plot.fig Urban_plot.fig
mv Hanson_plot.png Urban_plot.png

# Final formatting
matlab -nosplash -nodesktop -r 'combine_figs; quit;'

