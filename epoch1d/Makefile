# Set the compiler flags
FFLAGS = -O3 #-heap-arrays 0 -fast

# Set some of the build parameters
TARGET = epoch1d

# Set pre-processor defines
DEFINES = -DPER_PARTICLE_WEIGHT -DPARTICLE_PROBES -DTRACER_PARTICLES

# Copy this line onto the DEFINES line to turn on particle cell division
# -DSPLIT_PARTICLES_AFTER_PUSH
# Copy this line onto the DEFINES line to turn on debugging
# -DPARTICLE_DEBUG -DFIELD_DEBUG

MODULEFLAG = -module 
# Uncomment for the g95 compiler
#MODULEFLAG = -fmod=
# Uncomment for the gfortran compiler
#MODULEFLAG = -I/usr/include -Iobj -J

# If you are having trouble compiling the encoded_source module, try
# uncommenting the following line.
#ENCODED_SOURCE = dummy_encoded_source.o

#BUILDOPT += opt_dont_build_sha1
#BUILDOPT += opt_dont_build_embedded_source

# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------

-include COMMIT

SRCDIR = src
OBJDIR = obj
BINDIR = bin
FC = mpif90
DATE := $(shell date +%s)
MACHINE := $(shell uname -n)
PREPROFLAGS = $(DEFINES) -D_COMMIT='"$(COMMIT)"' -D_DATE=$(DATE) -D_MACHINE='"$(MACHINE)"'

SRCFILES = balance.F90 boundary.f90 calc_df.F90 cfd.f90 cfd_common.f90 \
  cfd_control.f90 cfd_input.f90 cfd_input_cartesian.f90 \
  cfd_input_functions.f90 cfd_input_particle.f90 cfd_job_info.f90 \
  cfd_output.f90 cfd_output_cartesian.f90 cfd_output_particle.f90 \
  custom_deck.f90 custom_laser.f90 custom_parser.f90 deck.F90 \
  deck_boundaries_block.f90 deck_constant_block.f90 deck_control_block.f90 \
  deck_deo_block.f90 deck_eio_dist_fn_block.f90 \
  deck_eio_particle_probe_block.F90 deck_ic_fields_block.f90 \
  deck_ic_laser_block.f90 deck_io_block.F90 deck_species_block.F90 \
  deck_window_block.f90 diagnostics.F90 dist_fn.F90 epoch1d.F90 evaluate.f90 \
  evaluator_blocks.f90 fields.f90 helper.F90 ic_module.f90 ionise.F90 \
  iterators.F90 laser.f90 mpi_routines.F90 mpi_subtype_control.f90 \
  particle_pointer_advance.f90 particle_temperature.F90 particles.F90 \
  partlist.F90 probes.F90 setup.F90 shape_functions.F90 shared_data.F90 \
  shunt.F90 simple_io.f90 split_particles.F90 stack.f90 strings.f90 \
  strings_advanced.f90 tokenizer_blocks.f90 utilities.f90 version_data.F90 \
  welcome.F90 window.F90

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
$(OBJDIR)/encoded_source.f90: $(BUILDOPT) $(SRCFILES) Makefile gen_src_module COMMIT
	./gen_src_module $@ $^

encoded_source.o: $(OBJDIR)/encoded_source.f90
	@mkdir -p $(BINDIR) $(OBJDIR)
	$(FC) -c $(FFLAGS)  $(MODULEFLAG)$(OBJDIR) -o $(OBJDIR)/$@ $<

COMMIT: $(SRCFILES)
	./gen_commit_string

clean: tidy
	@rm -rf *~ $(BINDIR)

tidy:
	@rm -rf $(OBJDIR) *.pbs.* *.sh.* $(SRCDIR)/*~ *.log

visit:
	@cd ../VisIt; xml2makefile -clobber cfd.xml; make

visit2:
	@cd ../VisIt; xml2info -clobber cfd2.xml; \
	  xml2cmake -clobber cfd2.xml; cmake .; make

visitclean:
	@cd ../VisIt; make clean; \
	  rm -rf .depend *.d CFD*Info.C CFD*Info.h CMake* cmake* Makefile

.PHONY: clean tidy visit visitclean

# All the dependencies

balance.o: balance.F90 boundary.o mpi_subtype_control.o partlist.o
boundary.o: boundary.f90 deck_io_block.o particle_temperature.o partlist.o
calc_df.o: calc_df.F90 boundary.o shape_functions.o shared_data.o
cfd.o: cfd.f90 cfd_common.o cfd_control.o cfd_input.o cfd_input_cartesian.o \
  cfd_input_functions.o cfd_input_particle.o cfd_output.o \
  cfd_output_cartesian.o cfd_output_particle.o
cfd_common.o: cfd_common.f90 cfd_job_info.o
cfd_control.o: cfd_control.f90 cfd_common.o cfd_input.o cfd_output.o
cfd_input.o: cfd_input.f90 cfd_common.o cfd_input_functions.o shared_data.o
cfd_input_cartesian.o: cfd_input_cartesian.f90 cfd_common.o \
  cfd_input_functions.o shared_data.o
cfd_input_functions.o: cfd_input_functions.f90 cfd_common.o
cfd_input_particle.o: cfd_input_particle.f90 cfd_common.o \
  cfd_input_functions.o shared_data.o
cfd_job_info.o: cfd_job_info.f90
cfd_output.o: cfd_output.f90 cfd_common.o cfd_job_info.o shared_data.o \
  version_data.o
cfd_output_cartesian.o: cfd_output_cartesian.f90 cfd_common.o cfd_output.o \
  shared_data.o
cfd_output_particle.o: cfd_output_particle.f90 cfd_common.o cfd_output.o \
  shared_data.o
custom_deck.o: custom_deck.f90 shared_data.o
custom_laser.o: custom_laser.f90 shared_data.o
custom_parser.o: custom_parser.f90 shared_data.o
deck.o: deck.F90 cfd.o custom_deck.o deck_boundaries_block.o \
  deck_constant_block.o deck_control_block.o deck_deo_block.o \
  deck_eio_dist_fn_block.o deck_eio_particle_probe_block.o \
  deck_ic_fields_block.o deck_ic_laser_block.o deck_io_block.o \
  deck_species_block.o deck_window_block.o shared_data.o strings.o \
  version_data.o
deck_boundaries_block.o: deck_boundaries_block.f90 strings.o
deck_constant_block.o: deck_constant_block.f90 strings_advanced.o
deck_control_block.o: deck_control_block.f90 fields.o strings_advanced.o
deck_deo_block.o: deck_deo_block.f90 shunt.o strings.o
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
diagnostics.o: diagnostics.F90 calc_df.o cfd.o deck.o dist_fn.o iterators.o \
  mpi_subtype_control.o probes.o shared_data.o version_data.o $(ENCODED_SOURCE)
dist_fn.o: dist_fn.F90 cfd.o mpi_subtype_control.o shared_data.o
epoch1d.o: epoch1d.F90 balance.o boundary.o cfd_job_info.o custom_parser.o \
  deck.o diagnostics.o fields.o helper.o ic_module.o ionise.o mpi_routines.o \
  particles.o setup.o split_particles.o welcome.o window.o
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
probes.o: probes.F90 cfd.o mpi_subtype_control.o partlist.o
setup.o: setup.F90 cfd.o fields.o mpi_subtype_control.o partlist.o \
  shared_data.o strings.o version_data.o $(ENCODED_SOURCE)
shape_functions.o: shape_functions.F90 shared_data.o
shared_data.o: shared_data.F90 cfd_job_info.o
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
