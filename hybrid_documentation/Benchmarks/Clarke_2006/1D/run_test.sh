#!/bin/bash
#run_test.sh
# This script attempts to reproduce Fig 13. of the Clarke paper, which simulates
# the X-ray emission due to a ~4e20 W/cm² Vulcan shot on a 3 mm gold target.

inputDir=`pwd`
coreNo=4
../../quick1d.sh $inputDir $coreNo
matlab -nosplash -nodesktop -r 'plot_Clarke; quit;'


