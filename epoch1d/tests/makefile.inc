
# this makefile should be included into each makefile in the examples folder.
# its a template adding some useful targets. only in exceptions those targets
# can be changed inside the actual makefile for the example. this ensures,
# that all examples work the same, so they can be used for testing.
#
# Stephan Kuchel, 150818


EPOCH=../../bin/epoch1d
# MPRIPROCS can be set as an environmental variable, using
# $export MPRIPROCS=8
MPIPROCS?=2

default: epoch1d.dat

epoch1d.dat: input.deck $(EPOCH)
	echo . | mpirun -n $(MPIPROCS) $(EPOCH)

clean:
	rm -rf *.png *.sdf deck.status epoch{1,2,3}d.dat *.visit
	rm -rf __pycache__
