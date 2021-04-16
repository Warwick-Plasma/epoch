#!/bin/bash
#cmake format: cmake -DGeant4_DIR=<path to geant4's install folder with "Geant4Config.cmake" in> <path to source code>
cmake -DGeant4_DIR=$HOME/geant4.10.04.p01-install/lib/Geant4-10.4.1 $HOME/bethe-bloch-G4/Au
make -j4

#Run this before executing to make sure the tables are linked correctly
#NOTE for some reason these lines need to be executed line by line, they will not work if ran from this shell script
cd $HOME/geant4.10.04.p01-install/bin
. geant4.sh
cd -
