##
## Programs and options...
##
include /home/space/phsgbd/visit/1.7.1/linux-x86_64/include/make-variables
TOPDIR=/home/space/phsgbd/visit/1.7.1/linux-x86_64


##
## Libraries and includes...
##
VTK_INCLUDE= \
  -I/home/space/phsgbd/visit/1.7.1/linux-x86_64/include/vtk \
  -I/home/space/phsgbd/visit/1.7.1/linux-x86_64/include/vtk/Common \
  -I/home/space/phsgbd/visit/1.7.1/linux-x86_64/include/vtk/Filtering \
  -I/home/space/phsgbd/visit/1.7.1/linux-x86_64/include/vtk/Graphics \
  -I/home/space/phsgbd/visit/1.7.1/linux-x86_64/include/vtk/Hybrid \
  -I/home/space/phsgbd/visit/1.7.1/linux-x86_64/include/vtk/IO \
  -I/home/space/phsgbd/visit/1.7.1/linux-x86_64/include/vtk/Imaging \
  -I/home/space/phsgbd/visit/1.7.1/linux-x86_64/include/vtk/Rendering \
  -I/home/space/phsgbd/visit/1.7.1/linux-x86_64/include/vtk/Utilities
CXXFLAGS=$(CXXFLAGSORIG)
CPPFLAGS=$(CPPFLAGSORIG) $(VTK_INCLUDE) -I. -I/home/space/phsgbd/visit/1.7.1/linux-x86_64/include -I/home/space/phsgbd/visit/1.7.1/linux-x86_64/include/visit $(PY_CXXFLAGS)
LDFLAGS=$(LDFLAGSORIG)  -L/home/space/phsgbd/.visit/linux-x86_64/plugins/databases

##
## Files...
##
PLUGINDIR=databases

WIDGETS=
VIEWERWIDGETS=
ISRC=CFDPluginInfo.C
COMMONSRC=CFDPluginInfo.C CFDCommonPluginInfo.C
MSRC=CFDMDServerPluginInfo.C avtCFDFileFormat.C CommonCommands.C Block.C BlockHandler.C BlockReader.C MeshReader.C Mesh_Var_Reader.C Stitched_Vector_Reader.C Stitched_Magnitude_Reader.C
ESRC=CFDEnginePluginInfo.C avtCFDFileFormat.C CommonCommands.C Block.C BlockHandler.C BlockReader.C MeshReader.C Mesh_Var_Reader.C Stitched_Vector_Reader.C Stitched_Magnitude_Reader.C


SRC=$(ISRC) $(COMMONSRC) $(MSRC) $(ESRC)

SHLIB_FORCED=

ELIBS_FOR_MACOSX_PREBINDING=$(BZIP2_LIBS) $(GLEW_LIBS) $(MESA_LIBS) $(GL_LIBS)
MLIBS_FOR_MACOSX_PREBINDING=$(GLEW_LIBS) $(BZIP2_LIBS) $(MESA_LIBS) $(GL_LIBS)
ILIBS=
MLIBS=-lpipeline_ser -ldatabase_ser -lmir_ser -lavtmath_ser -ldbatts -lavtexceptions -lstate -lmisc -lcomm -lexpr -lparser -lutility -lplugin -lvisit_vtk -llightweight_visit_vtk  $(MLIBS_FOR_MACOSX_PREBINDING) $(VTK_LIBS)
ESERLIBS=-lpipeline_ser -ldatabase_ser -lmir_ser -lavtmath_ser -ldbatts -lavtexceptions -lstate -lcomm -lmisc -lexpr -lparser -lplugin -lutility -lvisit_vtk -llightweight_visit_vtk $(ELIBS_FOR_MACOSX_PREBINDING) $(VTK_LIBS)
EPARLIBS=-lpipeline_par -ldatabase_par -lmir_par -lavtmath_par -ldbatts -lavtexceptions -lstate -lcomm -lmisc -lexpr -lparser -lplugin -lutility -lvisit_vtk -llightweight_visit_vtk $(ELIBS_FOR_MACOSX_PREBINDING) $(VTK_LIBS) $(SHLIB_MPI_LIBS)

IDSO=/home/space/phsgbd/.visit/linux-x86_64/plugins/databases/libICFDDatabase.so
MDSO=/home/space/phsgbd/.visit/linux-x86_64/plugins/databases/libMCFDDatabase.so
ESERDSO=/home/space/phsgbd/.visit/linux-x86_64/plugins/databases/libECFDDatabase_ser.so
EPARDSO=/home/space/phsgbd/.visit/linux-x86_64/plugins/databases/libECFDDatabase_par.so


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
include /home/space/phsgbd/visit/1.7.1/linux-x86_64/include/make-targets
