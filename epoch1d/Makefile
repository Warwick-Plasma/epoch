# Set the compiler flags
#FFLAGS = -fast#-fast #-arch pn4 -tpp7 -tune pn4 -ax
FFLAGS = -O3


# Set some of the build parameters
TARGET = epoch1d

# Set any precompiler options
PERPARTICLEWEIGHT = -DPER_PARTICLE_WEIGHT   # Uncomment for per particle weighting
#PERPARTICLECHARGEMASS = -DPER_PARTICLE_CHARGEMASS   # Uncomment for per particle weighting
INPUTFILE = -DINPUTFILE

# Uncomment one of the following lines if on a cluster.


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
PREPROFLAGS = $(PERPARTICLEWEIGHT) $(PERPARTICLECHARGEMASS) $(INPUTFILE)

OBJFILES = shared_data.o boundary.o iocommon.o output.o output_cartesian.o output_particle.o input.o iocontrol.o calc_df.o diagnostics.o setup.o fields.o particles.o  \
 output_arb.o inputfunctions.o input_arb.o input_cartesian.o input_particle.o partlist.o\
 strings.o balance.o mpi_routines.o helper.o initial_conditions.o deck_io_block.o deck_boundaries_block.o deck_control_block.o deck_species_block.o deck.o welcome.o epoch1d.o
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

# All the dependencies
shared_data.o:shared_data.F90
partlist.o:shared_data.o partlist.F90
balance.o:balance.f90 shared_data.o
mpi_routines.o:mpi_routines.f90 balance.o shared_data.o
setup.o:setup.F90 shared_data.o input_arb.o input_cartesian.o input_particle.o strings.o 
boundary.o:boundary.F90 shared_data.o partlist.o
iocommon.o:iocommon.f90 shared_data.o
inputfunctions.o: inputfunctions.f90 iocommon.o shared_data.o
input_cartesian.o: input_cartesian.f90 inputfunctions.o iocommon.o shared_data.o
input_particle.o: input_particle.f90 inputfunctions.o iocommon.o shared_data.o
input_arb.o: input_arb.f90 inputfunctions.o iocommon.o shared_data.o input.o iocontrol.o inputfunctions.o
input.o:input.f90 iocommon.o shared_data.o inputfunctions.o input_cartesian.o
utput.o:output.f90 iocommon.o shared_data.o
output_cartesian.o:output_cartesian.f90 iocommon.o shared_data.o output.o
output_particle.o:output_particle.f90 iocommon.o shared_data.o output.o
ouput_arb.o:output_arb.f90 iocommon.o shared_data.o output.o
iocontrol.o:iocontrol.f90 iocommon.o input.o output.o shared_data.o
calc_df.o:calc_df.F90 shared_data.o
diagnostics.o:diagnostics.F90 shared_data.o boundary.o input.o output.o iocontrol.o iocommon.o output_particle.o output_cartesian.o balance.o output_arb.o
fields.o:fields.f90 shared_data.o boundary.o
particles.o: particles.F90 shared_data.o boundary.o
strings.o:strings.f90 shared_data.o
helper.o:helper.F90 shared_data.o partlist.o
initial_conditions.o:initial_conditions.F90 strings.o shared_data.o helper.o
deck_io_block.o:deck_io_block.f90 strings.o shared_data.o
deck_control_block.o:deck_control_block.F90 strings.o shared_data.o
deck_boundaries_block.o:deck_boundaries_block.f90 strings.o shared_data.o
deck_species_block.o:deck_species_block.f90 strings.o shared_data.o
deck.o:deck.F90 shared_data.o initial_conditions.o strings.o deck_control_block.o deck_boundaries_block.o deck_species_block.o
welcome.o: welcome.f90 shared_data.o
epoch1d.o:epoch1d.F90 shared_data.o setup.o boundary.o diagnostics.o mpi_routines.o welcome.o initial_conditions.o strings.o deck.o mpi_routines.o fields.o particles.o balance.o
