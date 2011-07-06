import os
from distutils.sysconfig import get_python_lib
from distutils.core import setup, Extension

def getname():
  for r,d,fl in os.walk(get_python_lib()):
    for f in fl:
      if f == 'arrayobject.h':
        return os.path.realpath(os.path.join(r,'..'))
  return NULL

incdir=getname()

setup(name="sdf", version="1.0",
      ext_modules=[Extension("sdf", ["sdf_python.c", "sdf_control.c",
      "sdf_input.c", "sdf_input_cartesian.c", "sdf_input_point.c"],
      include_dirs=[incdir])])
