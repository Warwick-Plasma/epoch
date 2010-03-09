#ifndef CFD_CONSTANTSH
#define CFD_CONSTANTSH

#include <stdio.h>

  #define TYPE_SCRIBBLE -1
  #define TYPE_ADDITIONAL 0
  #define TYPE_MESH 1
  #define TYPE_MESH_VARIABLE 2
  #define TYPE_SNAPSHOT 3
  #define TYPE_STITCHED_VECTOR 4
  #define TYPE_STITCHED_MAGNITUDE 5

  #define MESH_CARTESIAN 0
  #define MESH_PARTICLE 1

  #define PARTICLE_CARTESIAN 0
  #define PARTICLE_POLAR 1
  #define PARTICLE_CYLINDRICAL 2

  #define VAR_CARTESIAN  0
  #define VAR_PARTICLE 1

  #define DIMENSION_IRRELEVANT 0
  #define DIMENSION_1D 1
  #define DIMENSION_2D 2
  #define DIMENSION_3D 3

  #define CFD_VERSION 1
  #define CFD_REVISION 0

#endif
