##
## Programs and options...
##
include /Users/christopherbrady/visit/1.11.2/darwin-i386/include/make-variables
TOPDIR=/Users/christopherbrady/visit/1.11.2/darwin-i386


##
## Libraries and includes...
##
VTK_INCLUDE= \
  -I/Users/christopherbrady/visit/1.11.2/darwin-i386/include/vtk \
  -I/Users/christopherbrady/visit/1.11.2/darwin-i386/include/vtk/Common \
  -I/Users/christopherbrady/visit/1.11.2/darwin-i386/include/vtk/Filtering \
  -I/Users/christopherbrady/visit/1.11.2/darwin-i386/include/vtk/Graphics \
  -I/Users/christopherbrady/visit/1.11.2/darwin-i386/include/vtk/Hybrid \
  -I/Users/christopherbrady/visit/1.11.2/darwin-i386/include/vtk/IO \
  -I/Users/christopherbrady/visit/1.11.2/darwin-i386/include/vtk/Imaging \
  -I/Users/christopherbrady/visit/1.11.2/darwin-i386/include/vtk/Rendering \
  -I/Users/christopherbrady/visit/1.11.2/darwin-i386/include/vtk/Utilities
CXXFLAGS=$(CXXFLAGSORIG)
CPPFLAGS=$(CPPFLAGSORIG) $(VTK_INCLUDE) -I. -I/Users/christopherbrady/visit/1.11.2/darwin-i386/include -I/Users/christopherbrady/visit/1.11.2/darwin-i386/include/mesa -I/Users/christopherbrady/visit/1.11.2/darwin-i386/include/visit $(PY_CXXFLAGS)
LDFLAGS=$(LDFLAGSORIG) $(PLUGIN_LDFLAGS) -L/Users/christopherbrady/.visit/darwin-i386/plugins/databases

##
## Files...
##
PLUGINDIR=databases
PLUGINNAME=CFDDatabase

WIDGETS=
VIEWERWIDGETS=
ISRC=CFDPluginInfo.C
COMMONSRC=CFDPluginInfo.C CFDCommonPluginInfo.C
MSRC=CFDMDServerPluginInfo.C avtCFDFileFormat.C CommonCommands.C Block.C BlockHandler.C BlockReader.C MeshReader.C Mesh_Var_Reader.C Stitched_Vector_Reader.C Stitched_Magnitude_Reader.C
ESRC=CFDEnginePluginInfo.C avtCFDFileFormat.C CommonCommands.C Block.C BlockHandler.C BlockReader.C MeshReader.C Mesh_Var_Reader.C Stitched_Vector_Reader.C Stitched_Magnitude_Reader.C


SRC=$(ISRC) $(COMMONSRC) $(MSRC) $(ESRC)

SHLIB_FORCED=

ELIBS_FOR_MACOSX_PREBINDING=$(ZLIB_LIB) $(GLEW_LIBS) $(MESA_LIBS) $(GL_LIBS)
MLIBS_FOR_MACOSX_PREBINDING=$(GLEW_LIBS) $(ZLIB_LIB) $(MESA_LIBS) $(GL_LIBS)
ILIBS=
MLIBS=-lpipeline_ser -ldatabase_ser -lmir_ser -lavtmath_ser -ldbatts -lavtexceptions -lstate -lmisc -lcomm -lexpr -lparser -lutility -lplugin -lvisit_vtk -llightweight_visit_vtk  $(MLIBS_FOR_MACOSX_PREBINDING) $(VTK_LIBS)
ESERLIBS=-lpipeline_ser -ldatabase_ser -lmir_ser -lavtmath_ser -ldbatts -lavtexceptions -lstate -lcomm -lmisc -lexpr -lparser -lplugin -lutility -lvisit_vtk -llightweight_visit_vtk $(ELIBS_FOR_MACOSX_PREBINDING) $(VTK_LIBS)
EPARLIBS=-lpipeline_par -ldatabase_par -lmir_par -lavtmath_par -ldbatts -lavtexceptions -lstate -lcomm -lmisc -lexpr -lparser -lplugin -lutility -lvisit_vtk -llightweight_visit_vtk $(ELIBS_FOR_MACOSX_PREBINDING) $(VTK_LIBS) $(SHLIB_MPI_LIBS)

IDSO=/Users/christopherbrady/.visit/darwin-i386/plugins/databases/libICFDDatabase.dylib
MDSO=/Users/christopherbrady/.visit/darwin-i386/plugins/databases/libMCFDDatabase.dylib
ESERDSO=/Users/christopherbrady/.visit/darwin-i386/plugins/databases/libECFDDatabase_ser.dylib
EPARDSO=/Users/christopherbrady/.visit/darwin-i386/plugins/databases/libECFDDatabase_par.dylib


DISTRIB=

##
## Derived objects
##
IOBJ_tmp=$(ISRC:.C=.o)
GOBJ_tmp=$(COMMONSRC:.C=.o) $(GSRC:.C=.o)
SOBJ_tmp=$(COMMONSRC:.C=.o) $(SSRC:.C=.o)
VOBJ_tmp=$(COMMONSRC:.C=.o) $(VSRC:.C=.o)
MOBJ_tmp=$(COMMONSRC:.C=.o) $(MSRC:.C=.o) $(MSPECIFICSRC:.C=_mds.o)
ESEROBJ_tmp=$(COMMONSRC:.C=.o) $(ESRC:.C=.o)     $(ESPECIFICSRC:.C=_eng.o)
EPAROBJ_tmp=$(COMMONSRC:.C=.o) $(ESRC:.C=_par.o) $(ESPECIFICSRC:.C=_par_eng.o)

IOBJ=$(IOBJ_tmp:.c=.o)
GOBJ=$(GOBJ_tmp:.c=.o)
SOBJ=$(SOBJ_tmp:.c=.o)
VOBJ=$(VOBJ_tmp:.c=.o)
MOBJ=$(MOBJ_tmp:.c=.o)
ESEROBJ=$(ESEROBJ_tmp:.c=.o)
EPAROBJ=$(EPAROBJ_tmp:.c=.o)

MOCSRC = $(WIDGETS:.h=_moc.C)
MOCOBJ = $(MOCSRC:.C=.o)

VIEWERMOCSRC = $(VIEWERWIDGETS:.h=_moc.C)
VIEWERMOCOBJ = $(VIEWERMOCSRC:.C=.o)

##
## Standard targets...
##
all: message $(IDSO) $(MDSERVERLIB) $(ENGINELIBSER) $(ENGINELIBPAR)

clean:
	$(RM) $(IOBJ) $(COMMONOBJ)
	$(RM) $(GOBJ) $(SOBJ) $(VOBJ) $(MOBJ) $(ESEROBJ) $(EPAROBJ)
	$(RM) $(MOCSRC) $(MOCOBJ) $(VIEWERMOCSRC) $(VIEWERMOCOBJ)
	$(RM) $(IDSO) $(GDSO) $(SDSO) $(VDSO) $(MDSO) $(ESERDSO) $(EPARDSO)
	$(RM) $(JAVAOBJ)

##
## Other targets...
##
message:
	@echo
	@echo "****************************************************************************"
	@echo "*** Building Compact Fluid Datatype Database Plugin"
	@echo "****************************************************************************"

##
## Automatic dependency stuff
##
include /Users/christopherbrady/visit/1.11.2/darwin-i386/include/make-targets
