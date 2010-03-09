# Set the compiler flags
FFLAGS = -O3 #-heap-arrays 0 -fast

# Set some of the build parameters
TARGET = epoch2d

# Set pre-processor defines
DEFINES = -DPER_PARTICLE_WEIGHT -DTRACER_PARTICLES

# Copy this line onto the DEFINES line to turn on particle cell division
# -DPARTICLE_CELL_DIVISION
# Copy this line onto the DEFINES line to turn on debugging
# -DPART_DEBUG -DFIELD_DEBUG

MODULEFLAG = -module 
# Uncomment for the g95 compiler
#MODULEFLAG = -fmod=
# Uncomment for the gfortran compiler
#MODULEFLAG = -J

# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------

SRCDIR = src
OBJDIR = obj
BINDIR = bin
FC = mpif90
PREPROFLAGS = $(DEFINES)

OBJFILES = balance.o boundary.o calc_df.o control.o current_smooth.o \
  custom_deck.o custom_laser.o custom_parser.o deck_boundaries_block.o \
  deck_constant_block.o deck_control_block.o deck_deo_block.o \
  deck_eio_dist_fn_block.o deck_eio_particle_probe_block.o \
  deck_ic_external_block.o deck_ic_fields_block.o deck_ic_laser_block.o \
  deck_ic_species_block.o deck_io_block.o deck.o deck_species_block.o \
  deck_window_block.o diagnostics.o dist_fn.o epoch2d.o evaluate.o \
  evaluator_blocks.o fields.o helper.o ic_module.o input_cartesian.o \
  input_functions.o input.o input_particle.o iocommon.o iocontrol.o ionise.o \
  iterators.o laser.o mpi_routines.o mpi_subtype_control.o output_cartesian.o \
  output.o output_particle.o particle_pointer_advance.o particles.o partlist.o \
  probes.o setup.o shape_functions.o shared_data.o shunt.o simple_io.o \
  split_particles.o stack.o strings_advanced.o strings.o tokenizer_blocks.o \
  welcome.o window.o
FULLTARGET = $(BINDIR)/$(TARGET)

#vpath %.f90 $(SRCDIR)
#vpath %.o $(OBJDIR)
VPATH = $(SRCDIR):$(SRCDIR)/deck:$(SRCDIR)/housekeeping:$(SRCDIR)/io:$(SRCDIR)/parser:$(SRCDIR)/physics_packages:$(SRCDIR)/user_interaction:$(OBJDIR)

# Rule to build the fortran files

%.o: %.f90
	@mkdir -p $(BINDIR) $(OBJDIR)
	$(FC) -c $(FFLAGS) $(MODULEFLAG)$(OBJDIR) -o $(OBJDIR)/$@ $<

%.o: %.F90
	@mkdir -p $(BINDIR) $(OBJDIR)
	$(FC) -c $(FFLAGS) $(MODULEFLAG)$(OBJDIR) -o $(OBJDIR)/$@ $(PREPROFLAGS) $<

$(FULLTARGET): $(OBJFILES)
	$(FC) $(FFLAGS) $(MODULEFLAG)$(OBJDIR) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILES))

.PHONEY: clean
clean:
	@rm -rf *~ $(BINDIR) $(OBJDIR) *.pbs.* *.sh.* $(SRCDIR)/*~ *.log

.PHONEY: tidy
tidy:
	@rm -rf $(OBJDIR) *.pbs.* *.sh.* $(SRCDIR)/*~ *.log

.PHONEY: visit
visit:
	@cd ../VisIt;xml2makefile -clobber cfd.xml;make

.PHONEY: visitclean
visitclean:
	@cd ../VisIt;make clean;rm -f .depend *.d

# All the dependencies

balance.o: balance.F90 boundary.o mpi_subtype_control.o partlist.o
boundary.o: boundary.f90 partlist.o
calc_df.o: calc_df.F90 boundary.o shape_functions.o shared_data.o
control.o: control.F90 setup.o shared_data.o
current_smooth.o: current_smooth.F90 shape_functions.o
custom_deck.o: custom_deck.f90 shared_data.o
custom_laser.o: custom_laser.f90 shared_data.o
custom_parser.o: custom_parser.f90 shared_data.o
deck_boundaries_block.o: deck_boundaries_block.f90 strings.o
deck_constant_block.o: deck_constant_block.f90 strings_advanced.o
deck_control_block.o: deck_control_block.f90 strings_advanced.o
deck_deo_block.o: deck_deo_block.f90 shunt.o strings.o
deck_eio_dist_fn_block.o: deck_eio_dist_fn_block.f90 dist_fn.o \
  strings_advanced.o
deck_eio_particle_probe_block.o: deck_eio_particle_probe_block.F90 probes.o \
  strings_advanced.o
deck_ic_external_block.o: deck_ic_external_block.f90 simple_io.o \
  strings_advanced.o
deck_ic_fields_block.o: deck_ic_fields_block.f90 strings_advanced.o
deck_ic_laser_block.o: deck_ic_laser_block.f90 laser.o strings_advanced.o
deck_ic_species_block.o: deck_ic_species_block.f90 strings_advanced.o
deck_io_block.o: deck_io_block.F90 strings_advanced.o
deck.o: deck.F90 custom_deck.o deck_boundaries_block.o deck_constant_block.o \
  deck_control_block.o deck_deo_block.o deck_eio_dist_fn_block.o \
  deck_eio_particle_probe_block.o deck_ic_external_block.o \
  deck_ic_fields_block.o deck_ic_laser_block.o deck_ic_species_block.o \
  deck_io_block.o deck_species_block.o deck_window_block.o strings.o
deck_species_block.o: deck_species_block.F90 setup.o strings_advanced.o
deck_window_block.o: deck_window_block.f90 strings_advanced.o
diagnostics.o: diagnostics.F90 calc_df.o dist_fn.o iocontrol.o iterators.o \
  mpi_subtype_control.o output_cartesian.o output_particle.o probes.o
dist_fn.o: dist_fn.F90 output_cartesian.o
epoch2d.o: epoch2d.F90 balance.o boundary.o control.o deck.o diagnostics.o \
  fields.o ic_module.o ionise.o mpi_routines.o particles.o setup.o \
  split_particles.o welcome.o window.o
evaluate.o: evaluate.f90 evaluator_blocks.o stack.o
evaluator_blocks.o: evaluator_blocks.f90 custom_parser.o stack.o
fields.o: fields.F90 boundary.o laser.o
helper.o: helper.F90 boundary.o shape_functions.o strings.o
ic_module.o: ic_module.f90
input_cartesian.o: input_cartesian.f90 input_functions.o
input_functions.o: input_functions.f90 iocommon.o
input.o: input.f90 input_functions.o
input_particle.o: input_particle.f90 input_functions.o shared_data.o
iocommon.o: iocommon.f90 shared_data.o
iocontrol.o: iocontrol.f90 input.o output.o
ionise.o: ionise.F90 calc_df.o helper.o
iterators.o: iterators.F90 particle_pointer_advance.o
laser.o: laser.f90 custom_laser.o evaluate.o
mpi_routines.o: mpi_routines.f90 partlist.o
mpi_subtype_control.o: mpi_subtype_control.F90 shared_data.o
output_cartesian.o: output_cartesian.f90 iocommon.o output.o
output.o: output.f90 iocommon.o
output_particle.o: output_particle.f90 iocommon.o output.o shared_data.o
particle_pointer_advance.o: particle_pointer_advance.F90 shared_data.o
particles.o: particles.F90 boundary.o current_smooth.o shape_functions.o
partlist.o: partlist.F90 shared_data.o
probes.o: probes.F90 mpi_subtype_control.o output_particle.o partlist.o
setup.o: setup.F90 input_cartesian.o input_particle.o iocontrol.o \
  mpi_subtype_control.o partlist.o strings.o
shape_functions.o: shape_functions.F90 shared_data.o
shared_data.o: shared_data.F90
shunt.o: shunt.F90 tokenizer_blocks.o
simple_io.o: simple_io.F90 boundary.o mpi_subtype_control.o
split_particles.o: split_particles.F90
stack.o: stack.f90 shared_data.o
strings_advanced.o: strings_advanced.F90 evaluate.o shunt.o
strings.o: strings.F90 shared_data.o
tokenizer_blocks.o: tokenizer_blocks.f90 strings.o
welcome.o: welcome.F90 strings.o
window.o: window.F90 helper.o
