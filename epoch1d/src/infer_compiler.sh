#!/bin/bash

COMPILER=''

MPIF90_VERS=`mpif90 --version 2>/dev/null`
#check for gfortran
echo -n $MPIF90_VERS | grep -iE "GNU|gfortran" &> /dev/null
if [ $? ==  0 ] ; then COMPILER='gfortran'; fi
#Check for intel
echo -n $MPIF90_VERS | grep -iE "ifort|intel" &> /dev/null
if [ $? ==  0 ] ; then COMPILER='intel'; fi
#Check for Cray, i.e. hector/archer
ftn -V &> /dev/null
if [ $? ==  0 ] ; then COMPILER='archer'; fi
#Check for ibm according to current Makefile
mpixlf90_r &> /dev/null
if [ $? ==  0 ] ; then COMPILER='ibm'; fi

echo -n $COMPILER
