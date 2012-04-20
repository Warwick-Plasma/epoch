import os, sys
from distutils.sysconfig import get_python_lib
from distutils.core import setup, Extension

def get_numpy_dir():
  for r,d,fl in os.walk(get_python_lib(standard_lib=1)):
    for f in fl:
      if f == 'arrayobject.h':
        return os.path.realpath(os.path.join(r,'..'))
  print 'Unable to build python module. Numpy directory not found.'
  sys.exit(1)


sdffiles = ['sdf_control.c', 'sdf_input.c', 'sdf_input_cartesian.c',
            'sdf_input_point.c']

srcfiles = ['sdf_python.c'] + sdffiles

incdirs = [get_numpy_dir()]

setup(name="sdf", version="1.0",
      ext_modules=[Extension("sdf", srcfiles, include_dirs=incdirs)])
