# Set the compiler flags
#FFLAGS = -fast#-fast #-arch pn4 -tpp7 -tune pn4 -ax
FFLAGS = -O3

# Set some of the build parameters
TARGET = epoch2d

# Set any precompiler options
#NONMPIIO = -DNONMPIIO   # Uncomment for non-MPI file access.

# Uncomment one of the following lines if on a cluster.
#COPSON = -compiler intel
#MHDCLUSTER = -f90=pgf90 -DMHDCLUSTER -fpic



# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------

SRCDIR = src
OBJDIR = obj
BINDIR = bin
MODULEFLAG = -module
MACHINEFLAGS = $(COPSON) $(MHDCLUSTER)
OPFLAGS = $(QMONO)
FC = mpif90 $(MACHINEFLAGS) $(OPFLAGS)
PREPROFLAGS = $(NONMPIIO)

OBJFILES = shared_data.o boundary.o iocommon.o output.o output_arb.o output_cartesian.o output_particle.o input.o inputfunctions.o input_particle.o input_cartesian.o \
 input_arb.o iocontrol.o calc_df.o diagnostics.o setup.o fields.o particles.o  \
 strings.o mpi_routines.o initial_conditions.o deck_io_block.o deck_boundaries_block.o deck_control_block.o deck_species_block.o deck.o welcome.o epoch2d.o
FULLTARGET = $(BINDIR)/$(TARGET)

#vpath %.f90 $(SRCDIR)
#vpath %.o $(OBJDIR)
VPATH = $(SRCDIR):$(OBJDIR)

# Rule to build the fortran files

%.o: %.f90
	@mkdir -p $(BINDIR) $(OBJDIR)
	$(FC) -c $(FFLAGS)  $(MODULEFLAG) $(OBJDIR) -o $(OBJDIR)/$@ $<

%.o: %.F90
	@mkdir -p $(BINDIR) $(OBJDIR) 
	$(FC) -c $(FFLAGS)  $(MODULEFLAG) $(OBJDIR) -o $(OBJDIR)/$@ $(PREPROFLAGS) $<

$(FULLTARGET): $(OBJFILES)
	$(FC) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILES))

.PHONEY: clean
clean:
	@rm -rf *~ $(BINDIR) $(OBJDIR) *.pbs.* *.sh.* $(SRCDIR)/*~ *.log

.PHONEY: tidy
tidy:
	@rm -rf $(OBJDIR) *.pbs.* *.sh.* $(SRCDIR)/*~ *.log

.PHONEY: visit
visit:
	@cd VisIT;xml2makefile -clobber cfd.xml;make

.PHONEY: visitclean
visitclean:
	@cd VisIT;make clean;rm -f .depend *.d
# All the dependencies
shared_data.o:shared_data.f90
mpi_routines.o:mpi_routines.f90 shared_data.o
setup.o:setup.f90 shared_data.o input_arb.o input_cartesian.o input_particle.o strings.o
boundary.o:boundary.f90 shared_data.o 
iocommon.o:iocommon.f90 shared_data.o
inputfunctions.o: inputfunctions.f90 iocommon.o shared_data.o
input_cartesian.o: input_cartesian.f90 inputfunctions.o iocommon.o shared_data.o
input_particle.o: input_particle.f90 inputfunctions.o iocommon.o shared_data.o
input_arb.o: input_arb.f90 inputfunctions.o iocommon.o shared_data.o input.o iocontrol.o inputfunctions.o
input.o:input.f90 iocommon.o shared_data.o inputfunctions.o input_cartesian.o
output.o:output.f90 iocommon.o shared_data.o
output_arb.o:output_arb.f90 iocommon.o shared_data.o
output_cartesian.o:output_cartesian.f90 iocommon.o shared_data.o output.o
output_particle.o:output_particle.f90 iocommon.o shared_data.o output.o
iocontrol.o:iocontrol.f90 iocommon.o input.o output.o shared_data.o
calc_df.o:calc_df.f90 shared_data.o boundary.o
diagnostics.o:diagnostics.f90 shared_data.o boundary.o input.o output.o iocontrol.o iocommon.o output_particle.o output_cartesian.o output_arb.o
fields.o:fields.f90 shared_data.o boundary.o
particles.o: particles.f90 shared_data.o boundary.o
strings.o:strings.f90 shared_data.o
initial_conditions.o:initial_conditions.f90 strings.o shared_data.o
deck_io_block.o:deck_io_block.f90 strings.o shared_data.o
deck_control_block.o:deck_control_block.f90 strings.o shared_data.o
deck_boundaries_block.o:deck_boundaries_block.f90 strings.o shared_data.o
deck_species_block.o:deck_species_block.f90 strings.o shared_data.o
deck.o:deck.f90 shared_data.o initial_conditions.o strings.o deck_control_block.o deck_boundaries_block.o deck_species_block.o
welcome.o: welcome.f90 shared_data.o
epoch2d.o:epoch2d.f90 shared_data.o setup.o boundary.o diagnostics.o mpi_routines.o welcome.o initial_conditions.o strings.o deck.o mpi_routines.o fields.o particles.o
