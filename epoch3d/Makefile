# Specify a particular compiler with "make COMPILER=pgi", etc.
# Specify debugging flags with "make MODE=debug"


# Compiler specific flags

# Use intel compatible flags by default
FFLAGS = -O3
#FFLAGS = -O3 -heap-arrays 64 -ipo -xHost # Optimised (B)
#FFLAGS = -O3 -heap-arrays 64 -ipo -xAVX  # Optimised (W)
#FFLAGS = -O0 -fpe0 -nothreads -traceback -fltconsistency \
#         -C -g -heap-arrays 64 -warn -save-temps -fpic -Wl,-no_pie # Debug
MODULEFLAG = -module $(OBJDIR)

ifeq ($(strip $(MODE)),debug)
  # Autodetect OSX
  SYSTEM := $(shell uname -s)
  FFLAGS = -O0 -fpe0 -nothreads -traceback -fltconsistency \
           -C -g -heap-arrays 64 -warn -fpic
  ifeq ($(strip $(SYSTEM)),Darwin)
    FFLAGS += -Wl,-no_pie
  endif
endif

MPIF90 ?= mpif90
D = -D

# PGI
# ===
ifeq ($(strip $(COMPILER)),pgi)
  FFLAGS = -r8 -fast -fastsse -O3 -Mipa=fast,inline -Minfo # Optimised
  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -Mbounds -g                                     # Debug
  endif
  MODULEFLAG = -module $(OBJDIR)
endif

# Intel
# =====
ifeq ($(strip $(COMPILER)),intel)
  FFLAGS = -O3
  #FFLAGS = -O3 -heap-arrays 64 -ipo -xHost # Optimised (B)
  #FFLAGS = -O3 -heap-arrays 64 -ipo -xAVX  # Optimised (W)
  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -O0 -fpe0 -nothreads -traceback -fltconsistency \
             -C -g -heap-arrays 64 -warn -fpic
    ifeq ($(strip $(SYSTEM)),Darwin)
      FFLAGS += -Wl,-no_pie
    endif
  endif
  MODULEFLAG = -module $(OBJDIR)
endif

# gfortran
# ========
ifeq ($(strip $(COMPILER)),gfortran)
  FFLAGS = -O3
  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -O0 -g -Wall -Wextra -pedantic -fbounds-check \
             -ffpe-trap=invalid,zero,overflow -Wno-unused-parameter \
             -ffpe-trap=underflow,denormal

    GNUVER := $(shell gfortran -dumpversion | head -1 \
        | sed 's/[^0-9\.]*\([0-9\.]\+\).*/\1/')
    GNUMAJOR := $(shell echo $(GNUVER) | cut -f1 -d\.)
    GNUMINOR := $(shell echo $(GNUVER) | cut -f2 -d\.)

    # gfortran-4.3
    GNUGE43 := $(shell expr $(GNUMAJOR) \>= 4 \& $(GNUMINOR) \>= 3)
    ifeq "$(GNUGE43)" "1"
      FFLAGS += -fbacktrace -fdump-core

      # gfortran-4.6
      GNUGE46 := $(shell expr $(GNUMINOR) \>= 6)
      ifeq "$(GNUGE46)" "1"
        FFLAGS += -Wno-unused-dummy-argument

        # gfortran-4.8
        GNUGE48 := $(shell expr $(GNUMINOR) \>= 8)
        ifeq "$(GNUGE48)" "1"
          FFLAGS += -Wno-target-lifetime
        endif
      endif
    endif
  endif
  MODULEFLAG = -I/usr/include -I$(OBJDIR) -J$(OBJDIR)
  INFO_FLAGS = -Wno-conversion -fno-range-check
endif

# g95
# ========
ifeq ($(strip $(COMPILER)),g95)
  FFLAGS = -O3
  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -O0 -g                                        # Debug
  endif
  MODULEFLAG = -fmod=$(OBJDIR)
endif

# IBM Bluegene
# ============
ifeq ($(strip $(COMPILER)),ibm)
  FFLAGS = -O5 -qhot -qipa # Optimised
  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -O0 -C -g -qfullpath -qinfo #-qkeepparm -qflttrap \
          -qnosmp -qxflag=dvz -Q! -qnounwind -qnounroll # Debug
    #FFLAGS = -O0 -qarch=qp -qtune=qp
    #FFLAGS = -qthreaded -qsmp=noauto -qsmp=omp # Hybrid stuff
  endif
  MODULEFLAG = -I$(OBJDIR) -qmoddir=$(OBJDIR)
  MPIF90 = mpixlf90_r

  # IBM compiler needs a -WF to recognise preprocessor directives
  D = -WF,-D
endif

# HECToR
# ========
ifeq ($(strip $(COMPILER)),hector)
  FFLAGS = -O3
  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -O0 -g -ea -ec -eC -eD -eI -en -hfp_trap -Ktrap=fp -m0 -M1438,7413
  endif
  MODULEFLAG = -em -I/usr/include -I$(OBJDIR) -J$(OBJDIR)
  MPIF90 = ftn
endif

# Don't compile encoded source if MODE=debug or ENC=no
# Do compile encoded source if MODE=debug and ENC=yes
ifeq ($(strip $(ENC)),no)
  ENCODED_SOURCE = epoch_source_info_dummy.o
else
  ifneq ($(strip $(ENC)),yes)
    ifeq ($(strip $(MODE)),debug)
      ENCODED_SOURCE = epoch_source_info_dummy.o
    endif
  endif
endif

FFLAGS += -I$(SDF)/include
FFLAGS += $(MODULEFLAG)
LDFLAGS = $(FFLAGS) -L$(SDF)/lib -lsdf

# Set some of the build parameters
TARGET = epoch3d

# Set pre-processor defines
DEFINES := $(DEF)
DEFINES += $(D)PER_PARTICLE_WEIGHT
DEFINES += $(D)PARTICLE_PROBES
DEFINES += $(D)TRACER_PARTICLES

# The following are a list of pre-processor defines which can be added to
# the above line modifying the code behaviour at compile time.

# Allow each pseudoparticle to represent a different number of real particles.
#DEFINES += $(D)PER_PARTICLE_WEIGHT

# Enable particles which behave like normal particles but do not affect
# the field, current or each other.
#DEFINES += $(D)TRACER_PARTICLES

# Include diagnostics which gather information on particles which cross a
# plane in the simulation.
#DEFINES += $(D)PARTICLE_PROBES

# Use second order particle weighting (default is third order).
#DEFINES += $(D)PARTICLE_SHAPE_TOPHAT

# Use fifth order particle weighting (default is third order).
#DEFINES += $(D)PARTICLE_SHAPE_BSPLINE3

# Include a unique global particle ID. The first flag defines the ID using
# an 8-byte integer, the second uses 4-bytes.
#DEFINES += $(D)PARTICLE_ID
#DEFINES += $(D)PARTICLE_ID4

# Keep global particle counts up to date.
#DEFINES += $(D)PARTICLE_COUNT_UPDATE

# Include QED routines
#DEFINES += $(D)PHOTONS

# Use the Trident process for pair production
#DEFINES += $(D)TRIDENT_PHOTONS

# Use Intel-specific 'mm_prefetch' calls to load next particle in the list
# into cache ahead of time. Gives ~10% speed-up with Intel 12 compilers.
#DEFINES += $(D)PREFETCH

# Turn on debugging.
#DEFINES += $(D)PARSER_DEBUG $(D)PARTICLE_DEBUG $(D)MPI_DEBUG $(D)SIMPLIFY_DEBUG

# Don't generate any output at all. Useful for benchmarking.
#DEFINES += $(D)NO_IO

# Bypass the main simulation and only perform collision tests.
# This flag should only ever be used for testing purposes.
#DEFINES += $(D)COLLISIONS_TEST

# Specify charge and mass per particle rather than per species.
# This is a legacy option which will be removed soon.
#DEFINES += $(D)PER_PARTICLE_CHARGE_MASS

# Perform checks on evaluated deck expressions.
# This slows down the code but may be required if floating point exceptions
# are enabled.
#DEFINES += $(D)PARSER_CHECKING


# If you are having trouble compiling the encoded_source module or
# don't want embedded source code in restart dumps, uncomment the
# following line.
#ENCODED_SOURCE = epoch_source_info_dummy.o



# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------

ENCODED_SOURCE ?= epoch_source_info.o

all: epoch

SDF := ../SDF/FORTRAN
SDFMOD = $(SDF)/include/sdf.mod
SRCDIR = src
OBJDIR = obj
BINDIR = bin
INCDIR = $(SRCDIR)/include
FC = $(MPIF90)
DATE := $(shell date +%s)
MACHINE := $(shell uname -n)
PREPROFLAGS = $(DEFINES) $(D)_COMMIT='"$(COMMIT)"' $(D)_DATE=$(DATE) \
  $(D)_MACHINE='"$(MACHINE)"' -I$(INCDIR)

FC_INFO := $(shell ${FC} --version 2>/dev/null \
    || ${FC} -V 2>&1 | grep '[a-zA-Z]' | head -n 1)

SRCFILES = balance.F90 boundary.f90 calc_df.F90 collisions.F90 \
  current_smooth.F90 custom_deck.f90 custom_laser.f90 custom_parser.f90 \
  deck.F90 deck_boundaries_block.f90 deck_collision_block.F90 \
  deck_constant_block.f90 deck_control_block.F90 deck_dist_fn_block.f90 \
  deck_fields_block.f90 deck_io_block.F90 deck_io_global_block.F90 \
  deck_laser_block.f90 deck_particle_probe_block.F90 deck_qed_block.F90 \
  deck_species_block.F90 deck_subset_block.f90 deck_window_block.f90 \
  diagnostics.F90 dist_fn.F90 epoch3d.F90 evaluate.F90 evaluator_blocks.F90 \
  fields.f90 finish.f90 helper.F90 ic_module.f90 ionise.F90 iterators.F90 \
  laser.f90 md5.f90 mpi_routines.F90 mpi_subtype_control.f90 numerics.f90 \
  particle_migration.F90 particle_pointer_advance.f90 particle_temperature.F90 \
  particles.F90 partlist.F90 photons.F90 prefetch.F90 probes.F90 \
  random_generator.f90 redblack_module.f90 setup.F90 shape_functions.F90 \
  shared_data.F90 shunt.F90 simple_io.f90 split_particle.F90 stack.f90 \
  strings.f90 strings_advanced.f90 timer.f90 tokenizer_blocks.f90 \
  utilities.f90 version_data.F90 welcome.F90 window.F90

OBJFILES := $(SRCFILES:.f90=.o)
OBJFILES := $(OBJFILES:.F90=.o)

INCLUDES = $(INCDIR)/particle_head.inc $(INCDIR)/particle_to_grid.inc \
  $(INCDIR)/bspline3/b_part.inc $(INCDIR)/bspline3/e_part.inc \
  $(INCDIR)/bspline3/gx.inc $(INCDIR)/bspline3/gxfac.inc \
  $(INCDIR)/bspline3/hx_dcell.inc \
  $(INCDIR)/tophat/b_part.inc $(INCDIR)/tophat/e_part.inc \
  $(INCDIR)/tophat/gx.inc $(INCDIR)/tophat/gxfac.inc \
  $(INCDIR)/tophat/hx_dcell.inc \
  $(INCDIR)/triangle/b_part.inc $(INCDIR)/triangle/e_part.inc \
  $(INCDIR)/triangle/gx.inc $(INCDIR)/triangle/gxfac.inc \
  $(INCDIR)/triangle/hx_dcell.inc

OBJFILES := $(OBJFILES) $(ENCODED_SOURCE)

FULLTARGET = $(BINDIR)/$(TARGET)

SOURCE_ALL := Makefile $(SRCFILES) $(INCLUDES) pack.py pack.sh \
  $(SRCDIR)/COMMIT epoch_source_info_dummy.f90 gen_commit_string

#vpath %.f90 $(SRCDIR)
#vpath %.o $(OBJDIR)
VPATH = $(SRCDIR):$(SRCDIR)/deck:$(SRCDIR)/housekeeping:$(SRCDIR)/io:\
  $(SRCDIR)/parser:$(SRCDIR)/physics_packages:$(SRCDIR)/user_interaction:\
  $(SDF)/src:$(OBJDIR)

PACK_SDF = $(SRCDIR)/pack.sh

GIT_WORK_TREE = ".."
PACK_PREFIX = epoch
PACK_SOURCE_CODE = 1
PACK_GIT_DIFF = 1
PACK_GIT_DIFF_FROM_ORIGIN = 1
GENERATE_CHECKSUM = 1
F77_OUTPUT = 0
PACK_OPTS = $(GIT_WORK_TREE) $(PACK_PREFIX) $(PACK_SOURCE_CODE) \
    $(PACK_GIT_DIFF) $(PACK_GIT_DIFF_FROM_ORIGIN) $(GENERATE_CHECKSUM) \
    $(F77_OUTPUT)

$(SRCDIR)/COMMIT: FORCE
	@./$(SRCDIR)/gen_commit_string || $(MAKE) $(MAKECMDGOALS)
-include $(SRCDIR)/COMMIT

# Rule to build the fortran files

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $(OBJDIR)/$@ $<

%.o: %.F90
	$(FC) -c $(FFLAGS) -o $(OBJDIR)/$@ $(PREPROFLAGS) $<

epoch: $(FULLTARGET)
$(FULLTARGET): $(OBJFILES)
	@mkdir -p $(BINDIR)
	$(FC) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILES)) $(LDFLAGS)

$(OBJDIR)/epoch_source_info.f90: $(SOURCE_ALL)
	sh $(PACK_SDF) $(PACK_OPTS) $@ "$(FC_INFO)" "$(FFLAGS)" $^
epoch_source_info.o: $(OBJDIR)/epoch_source_info.f90 $(SOURCE_ALL)
	$(FC) -c $(FFLAGS) $(INFO_FLAGS) -o $(OBJDIR)/$@ $<

$(SDFMOD):
	$(MAKE) -C $(SDF)

clean:
	@rm -rf $(BINDIR) $(OBJDIR)

cleanall: tidy

tidy:
	@rm -rf $(OBJDIR) *~ *.pbs.* *.sh.* $(SRCDIR)/*~ *.log
	$(MAKE) -C $(SDF) cleanall

tarball:
	@cd ..; sh make_tarball.sh

visit:
	@cd ../SDF/VisIt; ./build

visitclean:
	@cd ../SDF/VisIt; make clean; ./build -c; \
	  rm -rf .depend *.d *Info.C *Info.h CMake* cmake* Makefile

$(OBJFILES): | $(OBJDIR)

$(OBJDIR):
	@mkdir -p $(OBJDIR)

.PHONY: clean cleanall tidy visit visitclean epoch FORCE

# All the dependencies

balance.o: balance.F90 boundary.o mpi_subtype_control.o redblack_module.o \
  timer.o
boundary.o: boundary.f90 laser.o mpi_subtype_control.o particle_temperature.o \
  partlist.o
calc_df.o: calc_df.F90 boundary.o shared_data.o
collisions.o: collisions.F90 calc_df.o prefetch.o
current_smooth.o: current_smooth.F90 shape_functions.o shared_data.o
custom_deck.o: custom_deck.f90 shared_data.o
custom_laser.o: custom_laser.f90 shared_data.o
custom_parser.o: custom_parser.f90 shared_data.o
deck.o: deck.F90 custom_deck.o deck_boundaries_block.o deck_collision_block.o \
  deck_constant_block.o deck_control_block.o deck_dist_fn_block.o \
  deck_fields_block.o deck_io_block.o deck_io_global_block.o \
  deck_laser_block.o deck_particle_probe_block.o deck_qed_block.o \
  deck_species_block.o deck_subset_block.o deck_window_block.o md5.o photons.o \
  strings.o version_data.o
deck_boundaries_block.o: deck_boundaries_block.f90 strings_advanced.o
deck_collision_block.o: deck_collision_block.F90 collisions.o strings_advanced.o
deck_constant_block.o: deck_constant_block.f90 evaluate.o shunt.o
deck_control_block.o: deck_control_block.F90 fields.o strings_advanced.o timer.o
deck_dist_fn_block.o: deck_dist_fn_block.f90 dist_fn.o strings_advanced.o
deck_fields_block.o: deck_fields_block.f90 simple_io.o strings_advanced.o
deck_io_block.o: deck_io_block.F90 strings_advanced.o
deck_io_global_block.o: deck_io_global_block.F90 strings_advanced.o
deck_laser_block.o: deck_laser_block.f90 laser.o strings_advanced.o
deck_particle_probe_block.o: deck_particle_probe_block.F90 probes.o \
  strings_advanced.o
deck_qed_block.o: deck_qed_block.F90 strings_advanced.o
deck_species_block.o: deck_species_block.F90 partlist.o setup.o simple_io.o \
  strings_advanced.o utilities.o
deck_subset_block.o: deck_subset_block.f90 strings_advanced.o utilities.o
deck_window_block.o: deck_window_block.f90 strings_advanced.o
diagnostics.o: diagnostics.F90 calc_df.o deck.o dist_fn.o iterators.o probes.o \
  random_generator.o setup.o shared_data.o strings.o timer.o version_data.o \
  window.o $(SDFMOD) $(ENCODED_SOURCE)
dist_fn.o: dist_fn.F90 mpi_subtype_control.o
epoch3d.o: epoch3d.F90 balance.o collisions.o deck.o diagnostics.o fields.o \
  finish.o helper.o ic_module.o ionise.o mpi_routines.o particle_migration.o \
  particles.o photons.o setup.o split_particle.o welcome.o window.o
evaluate.o: evaluate.F90 evaluator_blocks.o shunt.o
evaluator_blocks.o: evaluator_blocks.F90 custom_parser.o stack.o
fields.o: fields.f90 boundary.o
finish.o: finish.f90 collisions.o deck.o diagnostics.o laser.o partlist.o \
  setup.o shared_data.o shunt.o window.o
helper.o: helper.F90 boundary.o calc_df.o partlist.o strings.o
ic_module.o: ic_module.f90 helper.o shared_data.o
ionise.o: ionise.F90 numerics.o partlist.o random_generator.o
iterators.o: iterators.F90 particle_pointer_advance.o partlist.o
laser.o: laser.f90 custom_laser.o evaluate.o
md5.o: md5.f90
mpi_routines.o: mpi_routines.F90 helper.o
mpi_subtype_control.o: mpi_subtype_control.f90 shared_data.o
numerics.o: numerics.f90 shared_data.o
particle_migration.o: particle_migration.F90 calc_df.o partlist.o prefetch.o
particle_pointer_advance.o: particle_pointer_advance.f90 shared_data.o
particle_temperature.o: particle_temperature.F90 random_generator.o \
  shared_data.o
particles.o: particles.F90 boundary.o current_smooth.o partlist.o prefetch.o
partlist.o: partlist.F90 random_generator.o shared_data.o
photons.o: photons.F90 partlist.o
prefetch.o: prefetch.F90 shared_data.o
probes.o: probes.F90 partlist.o
random_generator.o: random_generator.f90
redblack_module.o: redblack_module.f90 partlist.o
setup.o: setup.F90 fields.o helper.o laser.o mpi_subtype_control.o \
  shared_data.o shunt.o split_particle.o timer.o version_data.o welcome.o \
  window.o
shape_functions.o: shape_functions.F90 shared_data.o
shared_data.o: shared_data.F90 $(SDFMOD)
shunt.o: shunt.F90 tokenizer_blocks.o
simple_io.o: simple_io.f90 boundary.o mpi_subtype_control.o
split_particle.o: split_particle.F90 boundary.o
stack.o: stack.f90 shared_data.o
strings.o: strings.f90 shared_data.o
strings_advanced.o: strings_advanced.f90 evaluate.o shunt.o
timer.o: timer.f90 shared_data.o
tokenizer_blocks.o: tokenizer_blocks.f90 strings.o
utilities.o: utilities.f90 shared_data.o
version_data.o: version_data.F90 $(SRCDIR)/COMMIT
welcome.o: welcome.F90 strings.o version_data.o
window.o: window.F90 boundary.o evaluate.o partlist.o
