#!/bin/bash
#run_test.sh
# This script attempts to reproduce the Milchberg resistivity using the reduced
# Lee-More model. The starting electron temperature varies linearly between 0.1 
# and 200 eV, and the corresponding resistivity is dumped to file

inputDir=`pwd`
coreNo=4

# Run with default fit first
cp input_fit.deck input.deck
../../quick3d.sh $inputDir $coreNo
matlab -nosplash -nodesktop -r 'plot_resistivity; quit;'
mv reduced_Lee_More.fig with_fit.fig
mv reduced_Lee_More.png with_fit.png

# Run without fit second
cp input_nofit.deck input.deck
../../quick3d.sh $inputDir $coreNo
matlab -nosplash -nodesktop -r 'plot_resistivity; quit;'
mv reduced_Lee_More.fig without_fit.fig
mv reduced_Lee_More.png without_fit.png

# Final formatting
matlab -nosplash -nodesktop -r 'combine_figs; quit;'
