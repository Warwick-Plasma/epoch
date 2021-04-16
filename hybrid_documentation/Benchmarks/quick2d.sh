#!/bin/bash
#quick2d.sh
#This script takes two arguments
#
# 1) The directory of the input deck
# 2) The number of cores to run the job on 

# Path to the EPOCH directory containing the bin dir
epochDir=/path/to/epoch/epoch2d

# Run EPOCH
cd $epochDir
mpirun -np $2 ./bin/epoch2d <<< $1
