#! /bin/bash

COMPILER=''

MPIF90_VERS=`mpif90 --version 2>/dev/null`

# Check for gfortran
echo -n $MPIF90_VERS | grep -iwE "GNU|gfortran" &> /dev/null
[ $? == 0 ] && COMPILER='gfortran'

# Check for intel
echo -n $MPIF90_VERS | grep -iE "ifort|intel" &> /dev/null
[ $? == 0 ] && COMPILER='intel'

# Check for PGI
echo -n $MPIF90_VERS | grep -iwE "PGI|pgfortran" &> /dev/null
[ $? == 0 ] && COMPILER='pgi'

# Check for G95
echo -n $MPIF90_VERS | grep -iw "G95" &> /dev/null
[ $? == 0 ] && COMPILER='g95'

# Check for Cray, i.e. hector/archer
ftn -V &> /dev/null
[ $? == 0 ] && COMPILER='archer'

# Check for ibm according to current Makefile
mpixlf90_r &> /dev/null
[ $? == 0 ] && COMPILER='ibm'

echo -n $COMPILER
