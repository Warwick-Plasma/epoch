# Set the compiler flags

FFLAGS = -O3 #-heap-arrays 0#-fast

# Set some of the build parameters
TARGET = epoch1d

#uncomment this line to make each particle have an individual weight
DEFINES = -DPER_PARTICLE_WEIGHT -DPART_IONISE -DTRACER_PARTICLES -DPARTICLE_PROBES

#Copy this line onto the DEFINES line to turn on particle cell division
#-DPARTICLE_CELL_DIVISION
#Copy this line ont
#-DPART_DEBUG -DFIELD_DEBUG
# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------

SRCDIR = src
OBJDIR = obj
BINDIR = bin
MODULEFLAG = -module
MACHINEFLAGS = $(COPSON) $(MHDCLUSTER)
OPFLAGS = $(QMONO)
FC = mpif90 $(MACHINEFLAGS) $(OPFLAGS) $(DEFINES)
PREPROFLAGS = $(NONMPIIO)

OBJFILES = shared_data.o boundary.o iocommon.o output.o output_arb.o output_cartesian.o output_particle.o input.o inputfunctions.o input_particle.o input_cartesian.o partlist.o\
 input_arb.o iocontrol.o calc_df.o diagnostics.o setup.o fields.o particles.o ionise.o\
 strings.o mpi_routines.o initial_conditions.o deck_io_block.o deck_boundaries_block.o helper.o\
 split_particles.o strings_advanced.o shunt.o evaluate.o tokenizer_blocks.o evaluator_blocks.o\
 custom_deck.o custom_parser.o stack.o deck_ic_species_block.o deck_constant_block.o deck_deo_block.o\
 deck_window_block.o deck_eio_dist_fn_block.o deck_eio_particle_probe_block.o\
 iterators.o probes.o mpi_subtype_control.o deck_ic_external_block.o simple_io.o\
 deck_ic_fields_block.o laser.o deck_ic_laser_block.o custom_laser.o window.o dist_fn.o control.o\
 particlepointeradvance.o deck_control_block.o deck_species_block.o deck.o welcome.o balance.o \
 shape_functions.o solve_gauss.o FMC_MPI.o epoch1d.o
FULLTARGET = $(BINDIR)/$(TARGET)

#vpath %.f90 $(SRCDIR)
#vpath %.o $(OBJDIR)
VPATH = $(SRCDIR):$(SRCDIR)/deck:$(SRCDIR)/io:$(SRCDIR)/parser:$(SRCDIR)/housekeeping:$(SRCDIR)/user_interaction/:$(SRCDIR)/multigrid:$(SRCDIR)/physics_packages:$(OBJDIR)

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
shared_data.o:shared_data.F90
mpi_routines.o:mpi_routines.f90 shared_data.o partlist.o
balance.o: balance.F90 shared_data.o partlist.o boundary.o mpi_subtype_control.o
partlist.o:partlist.F90 shared_data.o
setup.o:setup.F90 shared_data.o input_arb.o input_cartesian.o input_particle.o strings.o partlist.o balance.o mpi_subtype_control.o
boundary.o:boundary.f90 shared_data.o partlist.o
simple_io.o:simple_io.F90 shared_data.o mpi_subtype_control.o boundary.o
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
calc_df.o:calc_df.F90 shared_data.o boundary.o
particlepointeradvance.o:particlepointeradvance.F90 shared_data.o
dist_fn.o:dist_fn.F90 shared_data.o partlist.o output_cartesian.o output.o iocontrol.o
iterators.o:iterators.F90 shared_data.o particlepointeradvance.o
probes.o:probes.F90 shared_data.o partlist.o output_particle.o output.o iocontrol.o
diagnostics.o:diagnostics.F90 shared_data.o boundary.o input.o output.o iocontrol.o iocommon.o output_particle.o \
	output_cartesian.o output_arb.o balance.o iterators.o dist_fn.o calc_df.o probes.o mpi_subtype_control.o
laser.o:laser.f90 shared_data.o evaluate.o shunt.o partlist.o custom_laser.o
fields.o:fields.f90 shared_data.o boundary.o laser.o
shape_functions.o: shape_functions.F90 shared_data.o
particles.o: particles.F90 shared_data.o boundary.o partlist.o shape_functions.o
strings.o:strings.F90 shared_data.o
strings_advanced.o:strings_advanced.F90 shared_data.o strings.o shunt.o evaluate.o
shunt.o:shunt.F90 shared_data.o strings.o tokenizer_blocks.o
tokenizer_blocks.o:tokenizer_blocks.f90 shared_data.o strings.o
evaluate.o:evaluate.f90 shared_data.o stack.o custom_parser.o evaluator_blocks.o
stack.o:stack.f90 shared_data.o
evaluator_blocks.o:evaluator_blocks.f90 shared_data.o strings.o custom_parser.o stack.o
custom_parser.o:custom_parser.f90 shared_data.o tokenizer_blocks.o stack.o
custom_deck.o:custom_deck.f90 shared_data.o strings.o strings_advanced.o shunt.o evaluate.o
custom_laser.o:custom_laser.f90 shared_data.o
window.o:window.F90 shared_data.o boundary.o partlist.o helper.o
helper.o:helper.F90 shared_data.o boundary.o
initial_conditions.o:initial_conditions.f90 strings.o shared_data.o helper.o shunt.o laser.o dist_fn.o probes.o simple_io.o
deck_io_block.o:deck_io_block.F90 strings.o shared_data.o strings_advanced.o
deck_control_block.o:deck_control_block.f90 strings.o shared_data.o strings_advanced.o
deck_boundaries_block.o:deck_boundaries_block.f90 strings.o shared_data.o strings_advanced.o
deck_species_block.o:deck_species_block.F90 strings.o shared_data.o strings_advanced.o setup.o
deck_ic_species_block.o:deck_ic_species_block.f90 shared_data.o strings_advanced.o strings.o shunt.o evaluate.o
deck_ic_fields_block.o:deck_ic_fields_block.f90 shared_data.o strings_advanced.o strings.o shunt.o evaluate.o
deck_ic_external_block.o:deck_ic_external_block.f90 shared_data.o strings_advanced.o simple_io.o
deck_constant_block.o:deck_constant_block.f90 shared_data.o strings.o strings_advanced.o
deck_deo_block.o:deck_deo_block.f90 shared_data.o strings.o strings_advanced.o shunt.o
deck_ic_laser_block.o:deck_ic_laser_block.f90 shared_data.o laser.o strings.o strings_advanced.o shunt.o evaluate.o
deck_window_block.o:deck_window_block.f90 shared_data.o strings_advanced.o
deck_eio_dist_fn_block.o:deck_eio_dist_fn_block.f90 shared_data.o strings_advanced.o dist_fn.o
deck_eio_particle_probe_block.o:deck_eio_particle_probe_block.F90 shared_data.o probes.o strings_advanced.o
deck.o:deck.F90 shared_data.o strings.o deck_control_block.o deck_boundaries_block.o deck_species_block.o\
	strings_advanced.o custom_deck.o deck_ic_species_block.o deck_constant_block.o deck_deo_block.o deck_ic_laser_block.o deck_window_block.o\
	deck_eio_dist_fn_block.o deck_eio_particle_probe_block.o deck_ic_external_block.o
welcome.o: welcome.F90 shared_data.o
split_particles.o: split_particles.F90 partlist.o shared_data.o helper.o
FMC_MPI.o: FMC_MPI.f90
solve_gauss.o: solve_gauss.F90 shared_data.o FMC_MPI.o
ionise.o: ionise.F90 shared_data.o partlist.o calc_df.o helper.o
control.o: control.F90 shared_data.o setup.o
epoch1d.o:epoch1d.F90 shared_data.o setup.o boundary.o diagnostics.o mpi_routines.o welcome.o initial_conditions.o strings.o deck.o mpi_routines.o fields.o particles.o boundary.o balance.o split_particles.o custom_deck.o FMC_MPI.o solve_gauss.o window.o ionise.o control.o
