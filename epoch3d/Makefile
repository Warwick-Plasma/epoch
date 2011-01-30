# Set the compiler flags
FFLAGS = -O3 #-heap-arrays 0 -fast

# The default MODULEFLAG is set for the Intel and PGI compilers
MODULEFLAG = -module 
# Uncomment for the g95 compiler
#MODULEFLAG = -fmod=
# Uncomment for the gfortran compiler
#MODULEFLAG = -I/usr/include -Iobj -J

# Set some of the build parameters
TARGET = epoch3d

# Set pre-processor defines
DEFINES = -DPER_PARTICLE_WEIGHT -DPARTICLE_PROBES -DTRACER_PARTICLES

# The following are a list of pre-processor defines which can be added to
# the above line modifying the code behaviour at compile time.

# Allow each pseudoparticle to represent a different number of real particles.
# -DPER_PARTICLE_WEIGHT

# Specify charge and mass per particle rather than per species.
# -DPER_PARTICLE_CHARGE_MASS

# Turn on particle cell division. Required for collisions.
# -DSPLIT_PARTICLES_AFTER_PUSH

# Turn on debugging.
# -DFIELD_DEBUG -DPARSER_DEBUG -DPARTICLE_DEBUG

# Turn on particle ionisation model. (Beta).
# -DPARTICLE_IONISE

# Keep global particle counts up to date.
# -DPARTICLE_COUNT_UPDATE

# Enable particles which behave like normal particles but do not affect
# the field, current or each other.
# -DTRACER_PARTICLES

# Include diagnostics which gather information on particles which cross a
# plane in the simulation.
# -DPARTICLE_PROBES

# Use second order particle weighting (default is third order).
# -DPARTICLE_SHAPE_TOPHAT

# Use fifth order particle weighting (default is third order).
# -DPARTICLE_SHAPE_BSPLINE3

# If you are having trouble compiling the encoded_source module or
# don't want embedded source code in restart dumps, uncommenting the
# following line.
#ENCODED_SOURCE = dummy_encoded_source.o



# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------

#BUILDOPT += opt_dont_build_sha1
#BUILDOPT += opt_dont_build_embedded_source

-include COMMIT

SRCDIR = src
OBJDIR = obj
BINDIR = bin
FC = mpif90
DATE := $(shell date +%s)
MACHINE := $(shell uname -n)
PREPROFLAGS = $(DEFINES) -D_COMMIT='"$(COMMIT)"' -D_DATE=$(DATE) -D_MACHINE='"$(MACHINE)"'

SRCFILES = balance.F90 boundary.f90 calc_df.F90 custom_deck.f90 \
  custom_laser.f90 custom_parser.f90 deck.F90 deck_boundaries_block.f90 \
  deck_constant_block.f90 deck_control_block.f90 deck_eio_dist_fn_block.f90 \
  deck_eio_particle_probe_block.F90 deck_ic_fields_block.f90 \
  deck_ic_laser_block.f90 deck_io_block.F90 deck_species_block.F90 \
  deck_window_block.f90 diagnostics.F90 dist_fn.F90 epoch3d.F90 evaluate.f90 \
  evaluator_blocks.f90 fields.f90 helper.F90 ic_module.f90 ionise.F90 \
  iterators.F90 laser.f90 mpi_routines.F90 mpi_subtype_control.f90 \
  particle_pointer_advance.f90 particle_temperature.F90 particles.F90 \
  partlist.F90 probes.F90 sdf.f90 sdf_common.f90 sdf_control.f90 sdf_input.f90 \
  sdf_input_cartesian.f90 sdf_input_point.f90 sdf_input_util.f90 \
  sdf_job_info.f90 sdf_output.f90 sdf_output_cartesian.f90 \
  sdf_output_point.f90 sdf_output_util.f90 setup.F90 shape_functions.F90 \
  shared_data.F90 shunt.F90 simple_io.f90 split_particles.F90 stack.f90 \
  strings.f90 strings_advanced.f90 tokenizer_blocks.f90 utilities.f90 \
  version_data.F90 welcome.F90 window.F90

OBJFILES := $(SRCFILES:.f90=.o)
OBJFILES := $(OBJFILES:.F90=.o)

ENCODED_SOURCE ?= encoded_source.o
OBJFILES := $(OBJFILES) $(ENCODED_SOURCE)

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

opt_dont_build_sha1: ;
opt_dont_build_embedded_source: ;
$(OBJDIR)/encoded_source.f90: $(BUILDOPT) $(SRCFILES) $(SRCDIR)/include \
  $(SRCDIR)/gen_src_module $(SRCDIR)/gen_commit_string Makefile COMMIT
	./$(SRCDIR)/gen_src_module $@ $^

encoded_source.o: $(OBJDIR)/encoded_source.f90
	@mkdir -p $(BINDIR) $(OBJDIR)
	$(FC) -c $(FFLAGS)  $(MODULEFLAG)$(OBJDIR) -o $(OBJDIR)/$@ $<

COMMIT: $(SRCFILES)
	./$(SRCDIR)/gen_commit_string

clean: tidy
	@rm -rf *~ $(BINDIR)

tidy:
	@rm -rf $(OBJDIR) *.pbs.* *.sh.* $(SRCDIR)/*~ *.log

visit:
	@cd ../VisIt; ./build sdf

visitclean:
	@cd ../VisIt; make clean; \
	  rm -rf .depend *.d *Info.C *Info.h CMake* cmake* Makefile

.PHONY: clean tidy visit visitclean

# All the dependencies

balance.o: balance.F90 boundary.o mpi_subtype_control.o partlist.o
boundary.o: boundary.f90 deck_io_block.o particle_temperature.o partlist.o
calc_df.o: calc_df.F90 boundary.o shape_functions.o shared_data.o
custom_deck.o: custom_deck.f90 shared_data.o
custom_laser.o: custom_laser.f90 shared_data.o
custom_parser.o: custom_parser.f90 shared_data.o
deck.o: deck.F90 custom_deck.o deck_boundaries_block.o deck_constant_block.o \
  deck_control_block.o deck_eio_dist_fn_block.o \
  deck_eio_particle_probe_block.o deck_ic_fields_block.o deck_ic_laser_block.o \
  deck_io_block.o deck_species_block.o deck_window_block.o sdf.o shared_data.o \
  strings.o version_data.o
deck_boundaries_block.o: deck_boundaries_block.f90 strings.o
deck_constant_block.o: deck_constant_block.f90 shunt.o strings.o
deck_control_block.o: deck_control_block.f90 fields.o strings_advanced.o
deck_eio_dist_fn_block.o: deck_eio_dist_fn_block.f90 dist_fn.o \
  strings_advanced.o
deck_eio_particle_probe_block.o: deck_eio_particle_probe_block.F90 probes.o \
  strings_advanced.o
deck_ic_fields_block.o: deck_ic_fields_block.f90 simple_io.o strings_advanced.o
deck_ic_laser_block.o: deck_ic_laser_block.f90 laser.o strings_advanced.o
deck_io_block.o: deck_io_block.F90 strings_advanced.o
deck_species_block.o: deck_species_block.F90 setup.o simple_io.o \
  strings_advanced.o utilities.o
deck_window_block.o: deck_window_block.f90 strings_advanced.o
diagnostics.o: diagnostics.F90 calc_df.o deck.o dist_fn.o iterators.o \
  mpi_subtype_control.o probes.o sdf.o shared_data.o version_data.o \
  $(ENCODED_SOURCE)
dist_fn.o: dist_fn.F90 mpi_subtype_control.o sdf.o shared_data.o
epoch3d.o: epoch3d.F90 balance.o boundary.o custom_parser.o deck.o \
  diagnostics.o fields.o helper.o ic_module.o ionise.o mpi_routines.o \
  particles.o sdf_job_info.o setup.o split_particles.o welcome.o window.o
evaluate.o: evaluate.f90 evaluator_blocks.o stack.o
evaluator_blocks.o: evaluator_blocks.f90 custom_parser.o stack.o
fields.o: fields.f90 boundary.o laser.o
helper.o: helper.F90 boundary.o strings.o
ic_module.o: ic_module.f90 helper.o shared_data.o
ionise.o: ionise.F90 calc_df.o helper.o
iterators.o: iterators.F90 particle_pointer_advance.o
laser.o: laser.f90 custom_laser.o evaluate.o
mpi_routines.o: mpi_routines.F90 partlist.o
mpi_subtype_control.o: mpi_subtype_control.f90 shared_data.o
particle_pointer_advance.o: particle_pointer_advance.f90 shared_data.o
particle_temperature.o: particle_temperature.F90 shape_functions.o
particles.o: particles.F90 boundary.o shape_functions.o
partlist.o: partlist.F90 shared_data.o
probes.o: probes.F90 mpi_subtype_control.o partlist.o sdf.o
sdf.o: sdf.f90 sdf_common.o sdf_control.o sdf_input.o sdf_input_cartesian.o \
  sdf_input_point.o sdf_input_util.o sdf_output.o sdf_output_cartesian.o \
  sdf_output_point.o sdf_output_util.o
sdf_common.o: sdf_common.f90 sdf_job_info.o
sdf_control.o: sdf_control.f90 sdf_common.o sdf_input.o sdf_output.o \
  sdf_output_util.o
sdf_input.o: sdf_input.f90 sdf_common.o
sdf_input_cartesian.o: sdf_input_cartesian.f90 sdf_common.o sdf_input.o
sdf_input_point.o: sdf_input_point.f90 sdf_common.o sdf_input.o
sdf_input_util.o: sdf_input_util.f90 sdf_input.o sdf_input_cartesian.o \
  sdf_input_point.o
sdf_job_info.o: sdf_job_info.f90
sdf_output.o: sdf_output.f90 sdf_common.o
sdf_output_cartesian.o: sdf_output_cartesian.f90 sdf_common.o sdf_output.o
sdf_output_point.o: sdf_output_point.f90 sdf_common.o sdf_output.o
sdf_output_util.o: sdf_output_util.f90 sdf_output.o sdf_output_cartesian.o \
  sdf_output_point.o
setup.o: setup.F90 fields.o mpi_subtype_control.o partlist.o sdf.o \
  shared_data.o strings.o version_data.o $(ENCODED_SOURCE)
shape_functions.o: shape_functions.F90 shared_data.o
shared_data.o: shared_data.F90 sdf.o sdf_job_info.o
shunt.o: shunt.F90 tokenizer_blocks.o
simple_io.o: simple_io.f90 boundary.o mpi_subtype_control.o
split_particles.o: split_particles.F90 helper.o shared_data.o
stack.o: stack.f90 shared_data.o
strings.o: strings.f90 shared_data.o
strings_advanced.o: strings_advanced.f90 evaluate.o shunt.o
tokenizer_blocks.o: tokenizer_blocks.f90 strings.o
utilities.o: utilities.f90 shared_data.o
version_data.o: version_data.F90
welcome.o: welcome.F90 strings.o version_data.o
window.o: window.F90 boundary.o helper.o partlist.o
