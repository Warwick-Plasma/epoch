#!/bin/bash
#run_test.sh
# This script attempts to reproduce Fig 4. of the Evans paper, which simulates
# the temperature increase due to a ~3e20 W/cm² Vulcan shot on a plastic target
# with an aluminium tracer layer. We also plot the temperature heatmap, like
# in Fig 2 (in the central z=0 plane).

inputDir=`pwd`
coreNo=4
../../quick1d.sh $inputDir $coreNo
matlab -nosplash -nodesktop -r 'plot_Evans; quit;'


