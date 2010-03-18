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
PREPROFLAGS = $(DEFINES) -D_COMMIT='"$(COMMIT)"' -D_DATE=$(DATE) -D_MACHINE='"$(MACHINE)"' -D_FLAGS='"$(FFLAGS)"'

SRCFILES = balance.F90 boundary.f90 calc_df.F90 control.F90 current_smooth.F90 \
  custom_deck.f90 custom_laser.f90 custom_parser.f90 deck.F90 \
  deck_boundaries_block.f90 deck_constant_block.f90 deck_control_block.f90 \
  deck_deo_block.f90 deck_eio_dist_fn_block.f90 \
  deck_eio_particle_probe_block.F90 deck_ic_external_block.f90 \
  deck_ic_fields_block.f90 deck_ic_laser_block.f90 deck_ic_species_block.f90 \
  deck_io_block.F90 deck_species_block.F90 deck_window_block.f90 \
  diagnostics.F90 dist_fn.F90 epoch2d.F90 evaluate.f90 evaluator_blocks.f90 \
  fields.f90 helper.F90 ic_module.f90 input.f90 input_cartesian.f90 \
  input_functions.f90 input_particle.f90 iocommon.f90 iocontrol.f90 ionise.F90 \
  iterators.F90 job_info.f90 laser.f90 mpi_routines.f90 \
  mpi_subtype_control.F90 output.f90 output_cartesian.f90 output_particle.f90 \
  particle_pointer_advance.F90 particles.F90 partlist.F90 probes.F90 setup.F90 \
  shape_functions.F90 shared_data.F90 shunt.F90 simple_io.F90 \
  split_particles.F90 stack.f90 strings.F90 strings_advanced.F90 \
  tokenizer_blocks.f90 version_data.F90 welcome.F90 window.F90

OBJFILES := $(SRCFILES:.f90=.o)
OBJFILES := $(OBJFILES:.F90=.o) encoded_source.o

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

.PHONEY: clean
clean: tidy
	@rm -rf *~ $(BINDIR)

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
deck.o: deck.F90 custom_deck.o deck_boundaries_block.o deck_constant_block.o \
  deck_control_block.o deck_deo_block.o deck_eio_dist_fn_block.o \
  deck_eio_particle_probe_block.o deck_ic_external_block.o \
  deck_ic_fields_block.o deck_ic_laser_block.o deck_ic_species_block.o \
  deck_io_block.o deck_species_block.o deck_window_block.o strings.o \
  version_data.o
deck_boundaries_block.o: deck_boundaries_block.f90 strings.o
deck_constant_block.o: deck_constant_block.f90 strings_advanced.o
deck_control_block.o: deck_control_block.f90 fields.o strings_advanced.o
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
deck_species_block.o: deck_species_block.F90 setup.o strings_advanced.o
deck_window_block.o: deck_window_block.f90 strings_advanced.o
diagnostics.o: diagnostics.F90 calc_df.o deck.o dist_fn.o encoded_source.o \
  iocontrol.o iterators.o mpi_subtype_control.o output_cartesian.o \
  output_particle.o probes.o
dist_fn.o: dist_fn.F90 output_cartesian.o
epoch2d.o: epoch2d.F90 balance.o boundary.o control.o deck.o diagnostics.o \
  fields.o ic_module.o ionise.o mpi_routines.o particles.o setup.o \
  split_particles.o welcome.o window.o
evaluate.o: evaluate.f90 evaluator_blocks.o stack.o
evaluator_blocks.o: evaluator_blocks.f90 custom_parser.o stack.o
fields.o: fields.f90 boundary.o laser.o
helper.o: helper.F90 boundary.o shape_functions.o strings.o
ic_module.o: ic_module.f90
input.o: input.f90 input_functions.o
input_cartesian.o: input_cartesian.f90 input_functions.o
input_functions.o: input_functions.f90 iocommon.o
input_particle.o: input_particle.f90 input_functions.o shared_data.o
iocommon.o: iocommon.f90 job_info.o shared_data.o
iocontrol.o: iocontrol.f90 input.o output.o
ionise.o: ionise.F90 calc_df.o helper.o
iterators.o: iterators.F90 particle_pointer_advance.o
job_info.o: job_info.f90
laser.o: laser.f90 custom_laser.o evaluate.o
mpi_routines.o: mpi_routines.f90 partlist.o
mpi_subtype_control.o: mpi_subtype_control.F90 shared_data.o
output.o: output.f90 iocommon.o version_data.o
output_cartesian.o: output_cartesian.f90 iocommon.o output.o
output_particle.o: output_particle.f90 iocommon.o output.o shared_data.o
particle_pointer_advance.o: particle_pointer_advance.F90 shared_data.o
particles.o: particles.F90 boundary.o current_smooth.o shape_functions.o
partlist.o: partlist.F90 shared_data.o
probes.o: probes.F90 mpi_subtype_control.o output_particle.o partlist.o
setup.o: setup.F90 fields.o input_cartesian.o input_particle.o iocontrol.o \
  mpi_subtype_control.o partlist.o strings.o version_data.o
shape_functions.o: shape_functions.F90 shared_data.o
shared_data.o: shared_data.F90 job_info.o
shunt.o: shunt.F90 tokenizer_blocks.o
simple_io.o: simple_io.F90 boundary.o mpi_subtype_control.o
split_particles.o: split_particles.F90
stack.o: stack.f90 shared_data.o
strings.o: strings.F90 shared_data.o
strings_advanced.o: strings_advanced.F90 evaluate.o shunt.o
tokenizer_blocks.o: tokenizer_blocks.f90 strings.o
version_data.o: version_data.F90
welcome.o: welcome.F90 strings.o version_data.o
window.o: window.F90 helper.o
