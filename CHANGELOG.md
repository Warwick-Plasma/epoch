## v4.16.0 to v4.17.0

 * Added volume correction sampling Zenitani 2015 DOI 10.1063/1.4919383

 * Added diagnostic for the average particle momentum per cell. This adds
   the flags "average_p{x,y,z}" to the output block.

 * Added low temperature correction to the Nanbu collision algorithm
   (Pérez 2012 DOI 10.1063/1.4742167)

 * Added "deck_warnings_fatal" flag to the control block. If set to "T" then the
   code will abort if any warnings are encountered during the deck parse. The
   default value is "F".

 * Added the DECK_DEBUG compiler flag. If enabled then this forces
   "deck_warnings_fatal" to always be set to "T"

 * Added missing "nproc_[xyz]" and isotropic "temp" parameters. These were
   mentioned in the documentation but had been omitted from the code.

 * Added "drift_p{x,y,z}" and "drift_p{x,y,z}_back" as aliases for
   "drift_{x,y,z}" and "drift_{x,y,z}_back"

 * Added the ability to specify "temp" and "temp_[xyz]" in electronvolts when
   setting up injectors.

 * Added parameter "number_density_max" to the injectors, matching the
   functionality available in the species block.

 * Added non-abbreviated aliases for the following deck keywords:
   "temperature" in place of "temp"
   "background" for "back"
   "nparticles" for "npart"
   "fraction" for "frac"
   "energy" for "en"
   "gamma_minus_one" for "gamma_m1"
   "average_particle_energy" for "ekbar"
   "particle_energy_flux" for "ekflux"
   "poynting_flux" for "poynt_flux"
   "polarisation" for "pol"
   "breit_wheeler" for "bw"

 * Allow subset restrictions to be time-varying functions

 * Added "heat_bath" particle boundary conditions. These use the particle
   injectors to create an improved form of thermal boundary conditions.

 * Add model for simulating bremsstrahlung radiation. This is disabled by
   default and can be enabled by specifying the -DBREMSSTRAHLUNG flag at
   compile time. The model is configured using a new "bremsstrahlung" block
   which takes the following parameters:
     enable
     start_time
     produce_photons
     photon_energy_min
     photon_weight
     photon_dynamics
     use_plasma_screening
     use_radiation_reaction
     table_location
   It also adds the "bremsstrahlung_optical_depth" output block parameter
   and "atomic_number", "background_species" and "brem_photon" species
   parameters. See the manual for further details.

 * Allow total_energy_sum diagnostic to be per-species

 * Removed the restrictions on minimum domain size. The domain can now be as
   small as a single cell.

 * Allow use_offset_grid to be specified on a per-output block basis.

 * Enable "use_exact_restart" by default


## v4.15.0 to v4.16.0

 * Replaced USE_ISATTY compiler flag with NO_USE_ISATTY

 * Add NO_MPI3 compiler flag to disable MPI-3 features such as the
   MPI_TYPE_SIZE_X routine. This allows the code to be compiled against
   older versions of the MPI library.

 * Added collisions block option "use_nanbu". If set to true then the
   scattering algorithm of Nanbu, with relativistic corrections by
   Perez will be used in the collisions module. If false, the previous
   Sentoku-Kemp algorithm will be used. Due to known issues in the Sentoku-Kemp
   method for some test problems, the Nanbu method is now the default.

 * Bugfixes for Sentoku-Kemp collisions. These changes will affect the results
   for some problems. Note that some test problems continue to demonstrate
   unexpected behaviour and users are advised to use the Nanbu method instead.

 * Added "number_density" aliases for "density" in the species and injector
   blocks.
   These include:
     - number_density for density
     - promote_number_density for promote_density
     - demote_number_density for demote_density
     - number_density_min for density_min
     - number_density_max for density_max
     - number_density_back for density_back

 * Added a "zero_current" alias for "tracer" in the species blocks. The use of
   "tracer" has now been deprecated and will be removed in version 5.0.
   At that time, the compiler flag will also be renamed.


## v4.14.0 to v4.15.0 (2019-01-15)

 * Enabled averaging of poynt_flux and ekflux variables

 * Added control block option "use_pre_balance". If set to true then the
   load balancer attempts to calculate the load balance prior to loading
   any particles

 * Added control block option "use_optimal_layout". If set to true then the
   load balancer attempts to find the best nprocx/y/z split of the domain
   prior to the initial load balance

 * Removed static parameter arrays from the input deck block parsers to
   simplify the process of adding new entries

 * Added control block option "use_more_setup_memory". If set to false then
   only one set of arrays will be used for storing temperature, density and
   drift during species loading. This can be a significant memory saving but
   comes at the expense of recalculating grid quantities multiple times.
   Setting the flag to true enables one set of arrays per species.
   The default value is false.

 * Replaced non-standard ISATTY intrinsic function with a call to the
   POSIX C-library function "isatty"

 * Changed FORTRAN standard to 2003

 * Added persistent subsets
   This adds the flags "persist_start_time" and "persist_start_step" to the
   subset block.  If either of these flags is present then once either the time
   specified by "persist_start_time" or the step specified by
   "persist_start_step" has been reached, the particles that have been selected
   using the other criteria in the subset block will be recorded. Each
   subsequent output for this subset will then use the particle list selected
   at the specified start time.

 * Added relativistic Maxwellians (Maxwell-Juttner)
   Adds the "use_maxwell_juttner" logical flag to the species block. If set
   to true then within the particle loader (both initially and for the moving
   window) the particles will be loaded following a Maxwell-Juttner distribution
   rather than a simple Maxwellian. If set to "F" then the standard Box-Muller
   loader for a Maxwellian distribution is used. The default value is "F".

 * Added arbitrary distribution functions in loader
   Adds the "dist_fn", "dist_fn_px_range", "dist_fn_py_range" and
   "dist_fn_pz_range" keys. The latter three keys set then range of momentum
   over which to sample the distribution function. The first key specifies the
   acceptance function. This should be a function having the maximum value of 1
   but the same shape as the true distribution function. It is combined with
   the density to calculate the full distribution function. The function is
   specified using the new deck keys "px", "py" and "pz".  The arbitrary
   distribution function can be combined with the "drift" keys and the
   function is shifted up by the drift momentum specified. The temperature
   key is not compatible and will be silently ignored if specified.

 * Added temperature_{x,y,z} output diagnostics


## v4.12.0 to v4.14.0 (2018-08-13)

 * Added the Higuera-Cary relativistic particle push. This is enabled by
   compiling with the flag -DHC_PUSH

 * Added sanity check for binary file size

 * Fixed the ETA string on restart

 * Added control block option "reset_walltime" to reset the walltime
   to zero when restarting

 * Added the ability to request output dumps at runtime by touching a file
   named "DUMP"

 * Made the x-coordinate a time-dependent input deck variable when the moving
   window is enabled

 * Added "window_stop_time" parameter to the window block

 * Added example SDF reader code

 * Added "atan2" function to the deck parser

 * Omit redundant dimension in dist_fn output

 * Modified load balance calculation so that a perfectly balanced problem has
   a value of one

 * Only redistribute the simulation if it would result in an improvement to
   the load balance

 * Improvements to the load balancing calculations

 * Added dlb_maximum_interval parameter to the input deck

 * Added dlb_force_interval parameter to the input deck

 * Added balance_first parameter to the input deck

 * Added y and z versions of the bc_x_{min,max}_after_move parameters to
   the input deck

 * Added "dump_at_walltimes" flag to the io block

 * Added "walltime_interval", "walltime_start" and "walltime_stop" flags to
   the io block


## v4.11.0 to v4.12.0 (2018-06-16)

 * Fixes for exact restarts

 * Apply reflecting particle BCs when conduct is specified

 * Fix moving window for fractional particles per cell

 * Don't average ppc diagnostic over shape function

 * Added average particle weight diagnostic

 * Replaced PARTICLE_COUNT_UPDATE Makefile DEFINE with
   "use_particle_count_update" control block flag

 * Added "use_accurate_n_zeros" control block flag to enable finding the
   exact number of output that are generated

 * Prevent duplicated final dump

 * Added "use_flux_maxwellian" option to the injector block

 * Added warnings for BCs that are incompatible with moving windows

 * Added y,z versions of the Lehe Maxwell solver

 * Added "stencil" block for specifying a Maxwell solver with custom
   coefficients

 * Added "work_{x,y,z}" and "work_{x,y,z}_total" dumpmasks for tracking the
   work done on each particle by the electric field.

 * Various bugfixes

 * Updated the SDF submodule:
   SDF/FORTRAN
   - Improve Fortran error handling

   SDF/utilities
   - Added support for SDF_BLOCKTYPE_POINT_DERIVED to sdffilter
   - Make sdffilter print all point mesh components
   - Return non-zero code for unsupported blocks in sdffilter
   - Fixed station file column order in sdffilter
   - Return an error for failed python builds

   SDF/VisIt
   - Add an option for disabling obstacle boundaries in VisIt
   - Don't generate boundary meshes when requested in VisIt


## v4.10.0 to v4.11.0 (2018-04-13)

 * Added time dependent moving window

 * Added multiple RNG states

 * Fixed F2003 extension in shared_data

 * Fixed errors in derived variable ouput for photons

 * Fix bug in grow_array in 2D

 * If print_const is used, output constant values into a separate file

 * Allow CPML with moving window

 * Added COMPILER=auto option to automatically detect compiler

 * Added time dependent moving window

 * Fix pack.sh script for non-bash shells

 * Added issue templates

 * Changed the Make rule from hector to archer

 * Fixed lasers so that positions work as expected

 * Fixed #1691, compilation error with older gfortran

 * Updated mediawiki links to the new location

 * Updated SDF submodule
   SDF/FORTRAN
   - Fix pack.sh script for non-bash shells

   SDF/IDL
   - Functions to read key-val pair file into struct
   - Rename main function
   - Added licensing header
   - Uppercase keywords
   - Correct indentation
   - General tidy

   SDF/Matlab
   - Read name-value pairs into structure
   - Removed DOS linefeed characters
   - Added license header
   - Tidied up a little

   SDF/utilities
   - Script to read name-val pairs from file into dict
   - Fix float vs int ordering
   - Made PEP8 compliant

   SDF/VisIt
   - Fix build script for VisIt-2.13 and 2.13.1 on macOS


## v4.9.0 to v4.10.0 (2018-03-13)

 * Added time varying particle injectors

 * Add per species particle boundaries You can now specify bc_x_min and
   bc_x_max to a species block. This overrides the global boundaries for that
   species

 * Fixes #1567 with moving window and offset_grid Bug 1567: UNLOCATABLE
   PARTICLE with restart dump and moving window Shift particles back to true
   position and adjust offset grid

 * Bugfix for assigning particle IDs

 * Various fixes for compiling with DEFINES flags

 * Added particles_per_cell diagnostic

 * Fixed missing ierror argument to MPI_ALLGATHER This fixes issue #1636

 * Fix PER_PARTICLE_CHARGE and fractional particles Previously 'fractional'
   particles were not given mass or charge giving segfault if ppc did not
   evenly divide by number of cells

 * Added CONTRIBUTING.md guide

 * Updated SDF
   SDF/FORTRAN
   - Fix PGI build on OSX
   - Work around a PGI floating point bug.
   - Fixed unix_seconds for Jan and Feb of leap years
   - Fully initialise sdf_file_handle
   - Changed FPE generation in error handler
   - Fix count of leap years since 1970
   - Check for broken OpenMPI 2.1.1 and 2.1.2

   SDF/utilities
   - Added support for string-type namevalue pairs
   - Avoid duplicated dictionary names
   - Made build of sdf2ascii optional (off by default)
   - Removed unhelpful dimension squashing
   - Removed mmap option from python reader

   SDF/VisIt
   - Added test for compatible g++ version


## v4.8.0 to v4.9.0 (2017-07-28)

 * Add alternative field solvers for the Maxwell equations.
   The following stencils are implemented lehe: R. Lehe
   et al., Phys. Rev. Special Topics 16, 021301 (2013) pukhov: A. Pukhov,
   Journal of Plasma Physics 61, 425-433 (1999) cowan: Cowan et al., Phys.
   Rev. ST Accel. Beams 16, 041303 (2013) The solvers can be selected in the
   control section of the input deck via
     maxwell_solver = {yee, lehe, cowan. pukhov}

   The default value is
     maxwell_solver = yee

   In the 1D-case Pukhov and Cowan reduce to Yee. In the 2D-case Pukhov and
   Cowan are identical, so Cowan reduces to Pukhov because Pukhov is older. In
   the 3D-case all solvers are implemented separately.

 * Added ability to read initial particle setup from a file
   This adds a "particles" custom particle loader block block works in
   conjunction with species block and must target a valid species block in
   order to load particles.  Particle data is loaded from RAW binary files,
   similar to the "fields" block, with one file per particle parameter and
   ordering assumed preserved across files.

 * Added a delta-f implementation of EPOCH.

 * Added automated continuous integration testing

 * Improved accuracy of relativistic kinetic energy calculation

 * Drive lasers using the integrated phase history rather than sin(omega*t).
   This is the same for first order in time but differs for higher order
   terms. See the discussion in issue #1468

 * Added the capability to set the frequency of lasers as a time dependent
   parameter in the input deck.

 * Added method for subtracting initial current (issues #1422 and #347)
   This adds the "use_current_correction" flag to the control block which will subtract any initial current from the simulation.

 * Add the input deck aliases for the alternative spellings:
   "field_ionization", "ionization_energies", "ionization_electron_species"
   and "collisional_ionization" (issue #1361)

 * Add |p| (mod_p) as distribution function direction (issue #1356)

 * Allow subsets to be applied to field variables (issue #1523)

 * Allow distribution functions to output spatial subsections if the range
   keyword is set. (issue #1525)

 * Add colour output for the EPOCH banner

 * Allow {y,z}_data variable to be present in 1d and 2d decks

 * Enable build of the VisIt SDF reader on Windows machines

Bugfixes:

 * Fixed a typo in deallocation of yb,zb

 * Fixed regression in restart code (issue #1358)

 * Fixed SEGFAULT when using field ionisation (issue #1407)

 * Use Rayleigh distribution for sampling in the thermal boundary conditions

 * Ensure plasma frequency calculation respects density_max (issue #1475)

 * Allow restarting of code with time varying laser frequency

 * Fix bug that in recording Trident electrons (issue #1512)

 * Disable load-balancer for dlb_threshold <= 0 (issue #1417)

 * Fix long filename path lengths (issue #1452)

 * Fix parallel loading of particles for the TOPHAT shape function (issue #287)

 * Fix current boundary conditions for thermal boundaries (issue #1542)

 * Prevent SEGFAULT when averaging with multiple output blocks by only
   allowing one average per variable. (issue #1529)

 * Change phase to inject laser properly at angle (issue #1538)

 * Fix secondary_list size when using field ionisation

 * Fixed missing c_dir_z in evalutor_blocks

 * Various minor bugfixes

 * Close files in VisIt after use (issue #1357)

 * Fixed the VisIt SDF reader build script for VisIt 2.7.{2,3}, 2.12.x


## v4.7.6 to v4.8.0 (2016-04-17)

 * Corrections for various issues in EPOCH's collisions

 * Updated SDF submodules

 * Added Makefile rule for making SDF utilities

 * Fixed operator precedence and associativity (issue #1295)

 * Fixed ETA timestring to be less than 80 columns

 * Warn about the use of -DPREFETCH (issue #1318)

 * Enforce nproc{x,y,z} if supplied (issues #1321 and #1322)

 * Added c_const_{xb,yb,zb} to the deck parser

 * Stagger grid when evaluating field components (issue #1337)

 * Add missing "restrict_<flag>" properties to dist_fn block, where <flag>
   includes the following: en, gamma_m1, xy_angle, yz_angle, zx_angle


## v4.7.0 to v4.7.6 (2016-02-12)

 * Fixed typo in check_qed_variable() routine

 * Updated README.md to prevent use of "Download zip"

 * Fixed compilation performance witn Intel compiler
   The intel compiler doesn't like inheriting the SDF module.

 * Handle unknown entries during stack sanity check
   Commit aa45aeb introduced sanity checking of the stack but not all block
   types are checked. This meant that the evaluation stack could be depleted
   and lead to out of bounds errors.

 * Fixed overflow when calculating particle count

 * Added GPL licensing headers

 * Moved the Changelog to the top-level directory


## v4.6.0 to v4.7.0 (2015-11-23)

 * Fixed build status URL

 * Updated submodules to point to GitLab server

 * Change MPI_ABORT exit codes to be less than 255

 * Disable dump_first flag when restarting

 * Added "dump_first_after_restart" flag to io block
   In the past, a "dump_first" flag in the output block would cause an output
   dump immediately after restarting. Since this is rarely the desired
   behaviour, the flag is now ignored when restarting. To force a dump to
   occur immediately after restart, set
   "dump_first_after_restart = T" in the output block.

 * Apply current boundary conditions after ionisation
   This fixes issue #1308 in which artefacts are visible at CPU boundaries
   when ionisation is enabled.

 * Fixed makefile to remove warning messages
   We now test for the existence of /usr/local before adding it to the include
   path.

 * Don't use fixed sizes for stack arrays
   This was causing SEGFAULTS when the input deck had a large number of terms
   on the right-hand side. Since we already know how many elements there are,
   we can just dynamically allocate them.

 * Use plain arrays for the eval_stack

 * Changed eval_stack to be dynamically sized

 * Perform sanity check on interpolation parameters

 * Issue a warning when interpolation range exceeded
   Previous versions of the code would exit with an error if the value of a
   point lay outside the range of values provided for interpolation. Now, the
   interpolation is performed and a warning message printed.

 * Issue a warning message when changing n_zeros
   EPOCH automatically adjusts the number of digits used for SDF files to
   accommodate more than 10000 dumps. This is now accompanied by a warning
   message so that the user knows why this is happening.

 * Store value of n_zeros specified from input deck

 * Print warning if requested n_zeros is invalid

 * Print warning if n_zeros is insufficient for t_end / dt_snapshots

 * Increase n_zeros if needed for further output

 * Added warning for when CPU split is too large

 * Give useful warning when domain split is too small

 * Output per-species charge and mass as scalar variables
   It was previously only possible to output these quantities as per-particle
   variables.

 * Print the name of invalid restart file.

 * Fixes for exact restarting from restart dumps

 * Allow forced write on same timestep

 * Updated SDF submodules
  SDF/C
   - Added support for the station_file flag
   - Bumped file and library revision numbers
   - Fixed some compiler warnings
   - Removed uthash submodule
   - Updated uthash.h symbolic link for new location
   - Use nblocks from header info when building summary
   - Added a flag for ignoring the nblocks header value
   - Fixed some coding-style issues
   - Fixed some memory leaks
   - Added routine for printing extension info.
   - Fixed double free SEGFAULT
  SDF/FORTRAN
   - Set station_file flag when a station block is found
   - Add symbolic names for the header offsets
   - Added station_file flag to the SDF header
   - Fixed makefile to remove warning messages
  SDF/IDL
   - Added support for station_flag
   - Fixed station_file typo
   - Bumped SDF file compatible revision
   - Allow for files with other than 4 digits
   - Fixed format specifier for previous patch
   - Specify minimum number of zeros for SDF file format
   - Fixed IDL load of 0th SDF snapshot
  SDF/VisIt
   - Added support for station_flag
   - Fixed station_file typo
   - Fix a concatenated string.
  SDF/utilities
   - fix segfaults when building against python3
   - Eliminated memory leak from previous change
   - Fix compilation on Ubuntu and add exit code
   - Fixed build script for python3
   - Fixed compiler warnings
   - Added support for station_flag
   - Fixed show-stopping typo
   - Added option to ignore header value of nblocks
   - Don't look for dictionary item by name
   - Removed incorrect Py_DECREF for Py_{True,False}
   - Added option for printing the loaded extension module
   - Added printing of file header to sdffilter
   - Add open limits to block ranges
   - Don't print header block when specifying block id
   - Store string properties as unicode for python3


## v4.5.0 to v4.6.0 (2015-09-18)

 * Exit with error if no valid COMPILER specified

 * Added basic gitlab continuous integration script.

 * Added more build targets to gitlab-ci script

 * Fixed missing script entries

 * Added VisIt to the gitlab-ci build script

 * Added test_all.sh script

 * Added separate gfortran and intel builds

 * Added -DPARSER_CHECKING for divide in input.deck parser

 * Show -DPARSER_CHECKING in welcome message

 * Added CI build status to README.md

 * Correctly initialise particles in ionisation routines

 * Properly initialise all newly created particles.

 * Revert "Fixed field setup when restarting."

 * Removed INTENT attribute to POINTER argument

 * Fixed pack.py when the argparse module is not found

 * Removed the "manuals" submodule

 * Removed manuals submodule from make_tarball.sh

 * Only match version tags in gen_commit_string.sh

 * Synced make_tarball script with ODIN.

 * Remove .git files and directories from tarball

 * add keyword 'VERSION' to exit after welcome message.
   In order to check the EPOCH version, the code will now exit with status 0,
   if 'VERSION' is given as data directory.

 * Renamed VERSION flag to VERSION_INFO

 * Correctly exit code in parallel

 * Reorganised setup routines

 * Change default behaviour of "make clean"
   In the old Makefile, "make clean" only cleans up the EPOCH object files and
   leaves the SDF module untouched. There is a separate "make cleanall"
   command which cleans up everything, but most users don't know about it and
   find the default behaviour unintuitive.
   "make clean" now does the same as "make cleanall" and there is a new rule
   "make rm" which performs the minimal clean.

 * Exit with a non-zero status when an error is found

 * Fixed warning messages about string truncation.

 * Update SDF submodules
  SDF/C
   - Add submodule uthash
   - Add members to sdf_block and sdf_file for hashing.
   - Add functions sdf_hash_block and sdf_hash_block_list
   - Use hash lookup table to retrieve blocks by id.
   - Use hash table lookup to find duplicate block ids.
   - Add newly read blocks to the hash table.
   - Add function sdf_delete_hash_block
   - Delete the blocks from the hash table when they are removed.
   - Add derived blocks to the hash table.
   - Delete blocks from the hash table before their ID gets mangled.
   - Define HASH_DEBUG when debugging.
   - Hash blocks by name, as well as id.
   - Added -fPIC flag to the Makefile.
   - Increase the length of the arrays dims and local_dims.
   - Fix some memory leaks.
   - Free extension data when closing a file.
   - Avoid freeing memory twice.
   - Add an option for compiling with valgrind support.
   - Increase the length of the arrays dims and local_dims.
   - Re-order sdf_block members.
   - Changed array dimensions to SDF_MAXDIMS
   - Fixed Makefile flags
   - Updated Doxygen version string
   - Fixed declarations
  SDF/FORTRAN
   - Fixed pack.py when the argparse module is not found
   - Fix the file path 'incfile', used by pack.py.
  SDF/VisIt
   - Added some missing library files to the build process
  SDF/utilities
   - Allow selective building of Python interface
   - Fix compiler flags.
   - Fix compiler option that only applies to gcc.
   - Use a dict object to find mesh blocks by id, instead of iterating through
     a list of items.
   - Add the function 'Block_setdata'
   - Fix some memory leaks.


## v4.4.2 to v4.5.0 (2015-06-19)

 * Fixed COMMIT file generation scripts.

 * Possible fix for recursive make.

 * Fixes for output of field variables.
   The recent addition of subset grid reduction broke normal ouput of field
   variables. This ought to be fixed now.

 * Prevent git error messages in pack.sh

 * Synced gen_commit_string.sh script with other repos

 * Added fix for the pack.sh shell script

 * Added functions for querying the SDF string length
   Also increased the default maximum size.

 * Check for name truncation due to the SDF library.

 * Fixed the MPI subarrays used for field subsets

 * Fixed reduced grid calculation.

 * Fixes for final memory deallocation.

 * Removed unnecessary repetition of subset name.

 * Removed gfortran floating underflow compilation flags.

 * Fixed GIT_DIR when CWD does not equal GIT_WORK_TREE

 * Fixed sdf_write_namevalue when h%string_length changes.
   Care must be taken when writing the array of strings in the
   sdf_write_namevalue routine since the strings are of length
   c_max_string_length which might be longer than the current value of
   h%string_length.

   It is also necessary to truncate the strings being passed into this routine
   in order to prevent warning messages being printed to screen.

 * Used the stack allocator in sdf_derived

 * Fixed building of SDF utility man pages on OSX

 * Added an error when no COMPILER is set for "make"

 * Added README.txt with full build instructions.

 * Disable the -ipo flag for Intel by default.

 * Fixed reading of integer arrays in parallel.

 * Removed unused variables.

 * Attempt to calculate GIT_DIR as an absolute pathname.

 * Changed mmap of data to be on a per-block basis.

 * Fix VisIt build script when PATH contains spaces

 * Added type conversion for integer arrays.

 * Added better load balancing heuristics

 * Removed the "error" spam from the io block parser.

 * Added some missing deallocations during finalise.

 * Use a persistent buffer for packed particle data

 * Always put .mod files into the include/ directory

 * Added control of n_zeros to the input deck

 * Added documentations to the SDF submodule.

 * Bugfixes for mmap'ed file access

 * Fixed memory corruption bug in python reader.
   The python reader called sdf_close() once the reference count of the SDF
   object reached zero. This would then deallocate all of the memory allocated
   by the SDF library, including the arrays containing the data. For now, I
   have removed the sdf_close() call entirely so there is now a massive memory
   leak. The obvious fix is to close the file once the last reference to the
   file's contents has reached zero, but first we need to return a data
   structure which tracks the returned items.

 * Fix memory leak in python reader.
   The memory leak introduced by the previous patch has now been fixed. The
   numpy arrays are now created with the sdf module as a base class. This way,
   the reference count for the sdf module only becomes zero once all the numpy
   arrays have gone out of scope.

 * Fixed lines longer than 80 columns

 * Minor fixes and tidying to the photons module.

 * Moved README.txt to README.md

 * Abort execution when the input deck is invalid.

 * Fixed field setup when restarting.

 * Fixed a compiler warning about formatting.

 * Added a floating-point consistency option.
   If the code uses 80-bit precision arithmetic for the current arrays, the
   answer will be inconsistent when run in parallel since only 64-bits of
   precision is communicated to neighbouring processes. This change adds a
   flag "CONS=1" that can be passed to to "make" which enforces 64-bit
   arithmetic in the field solver. Currently, it is only used for the Intel
   compiler.

 * strip non-ascii before packing source
   non-ascii chars are not supported in fortran. The code doesnt compile, if
   i.e. 'gfortran --version' outputs non-ascii characters, which it does in
   certain languages (i.e. German). This has now been fixed.

 * Added "manuals" repository as a submodule.

 * Avoid duplicate restart dumps.
   Fixes issue #384

 * A few small changes to the README file.
   References to CCPForge have been changed to GitLab now that the hosting
   location has changed

 * Fix a 'raise' statement.
   'Raise' should return an instance, not the class itself.

 * Made Downloads location more explicit in README.md

 * Fixed typos introduced by commit 761e415c
   This fixes issue #1259

 * Fixes for make_tarball.sh
   Added generation of commit_info.h for SDF/C and SDF/utilities submodules.
   Replaced submodule search paths with local versions.

 * Updated manuals submodule to latest master
   * Merge branch 'kuschel/manuals-master'
   - use hyperref to set links in pdf

 * Updated pack.py to match the version in SDF/FORTRAN

 * Call sdf_write_datablock() from all CPUs
   Commit 920465f removed MPI_BCAST from the sdf_write_datablock() routines
   but this had unintended consequences. This has now been reverted and the
   call is a collective operation so must be called by all CPUs in the
   communicator. This provides an alternative fix to issue #1256

 * Fixed gfortran detection for versions > 4.x

 * Apply commit ad6a808b to 1d and 3d

 * Updated SDF submodules
  SDF/C
   - Used the stack allocator in sdf_derived
   - Fixed a few integer sizes.
   - Changed mmap of data to be on a per-block basis.
   - Bugfixes for mmap'ed file access
   - Fixed lines longer than 80 columns
   - Bugfixes for setting up current CPU mesh
   - Added a routine for destroying the SDF stack.
   - Fixed the free'ing of NAMEVALUE and RUN_INFO blocks.
   - Fix copying of info_length.
   - Replaced a couple of numeric constants with symbols.
   - Fix for allocation of derived station blodk
   - Added a derived flag to the block type.
   - Added a handle to the stack allocator.
   - Fixed a memory leak for derived variables.
   - Fixed a few more memory leaks
   - Changed doxygen version to match library version.
   - Add a mechanism for including a layer of ghost cells for the VisIt reader.
   - Tidy-up of ghost-cell changes.
   - Added a few fields to the blocktype.
   - Use the field "offset" for station file offsets.
   - Fixed generation of CPU meshes with internal boundaries.
   - Fix a compiler option that only applies to gcc.
   - Fix compilation using clang
   - Don't mark data as read in contiguous stitched blocks.
   - Ensure that blocklist contains no duplicate IDs.
   - Always populate required variables
   - Don't modify local_dims for station variables.
   - Ensure that newly created block IDs are unique.
   - Added commit info to the SDF C library
   - Added generation of commit_info.h to CMakeLists.txt
   - Changed C library name to libsdfc.a
  SDF/FORTRAN
   - Fixed COMMIT file generation scripts.
   - Prevent git error messages in pack.sh
   - Synced gen_commit_string.sh script with other repos
   - Added fix for the pack.sh shell script
   - Added functions for querying the SDF string length
   - Fixed GIT_DIR when CWD does not equal GIT_WORK_TREE
   - Fixed sdf_write_namevalue when h%string_length changes.
   - Disable the -ipo flag for Intel by default.
   - Attempt to calculate GIT_DIR as an absolute pathname.
   - Force GIT_DIR to be an absoulte path.
   - Fix a spelling mistake.
   - Guarantee a clean GIT_DIR1, when finding the absolute GIT_DIR.
   - Added type conversion for integer arrays.
   - Always put .mod files into the include/ directory
   - Fixed lines longer than 80 columns
   - Fix the file view when reading in 2D and 3D Cartesian meshes.
   - Added SDF minor revision number.
   - Fix compiler options.
   - strip non-ascii before packing source
   - Fix a 'raise' statement.
   - Ensure pack.py complies with PEP8 guidelines.
   - Removed MPI_BCAST calls from serial routines.
   - Revert "Removed MPI_BCAST calls from serial routines."
   - Fixed gfortran detection for versions > 4.x
  SDF/Matlab
   - Added support for Lagrangian meshes to Matlab
   - Added a missing semicolon to the Matlab reader
   - Added missing support for INTEGER datatypes
  SDF/VisIt
   - Fixed COMMIT file generation scripts.
   - Fix VisIt build script when PATH contains spaces
   - Fixed lines longer than 80 columns
   - Added a handle to the stack allocator.
   - Remove unnecessary memory allocation.
   - Move a #include statement.
   - Include a layer of ghost cells at mesh boundaries when vieweing in parallel.
   - Allow freeing the stack when there is no filehandle.
   - Stop build after specified build type
   - Fixes to build script for VisIt 2.9.0
   - Fill in ghost cells for lagrangian meshes.
   - Tidy-up of ghost-cell changes.
   - Fixed generation of CPU meshes with internal boundaries.
   - Add an option to configure the installation directory.
   - Changed the names of commit info defines
   - Updated VisIt patches for the latest version.
  SDF/utilities
   - Fixed building of SDF utility man pages on OSX
   - Fixed memory corruption bug in python reader.
   - Fix memory leak in python reader.
   - Added NPY_1_7_API_VERSION to python reader.
   - Added backwards compatibility for older versions of numpy
   - Fixed indentation levels.
   - Added "mmap" flag to SDF_open
   - Added Array type for tracking numpy arrays
   - Free block data when numpy arrays are destroyed.
   - Added a python datatype for returning SDF blocks.
   - Return block structures instead of arrays.
   - Reorganised sdf_python file.
   - Added a set of python types to represent SDF block types.
   - Added more metadata fields to the block type.
   - Added BlockMesh members
   - Added __version__ atribute to python module
   - Removed redundant file handle from Array type
   - Replaced file handle with pointer to sdf in block type.
   - Export Block type
   - Removed a bunch of typecasts
   - Renamed SDF_type
   - Changed the python SDF read call to be a class method.
   - Added better docstring for the SDF read method
   - Removed some of the temporary "sub" objects
   - Made setup_* calls more consistent with each other.
   - Added a reference counter for the Block structure.
   - Changed mesh block to be closer to that of SDF
   - Added a getter routine for BlockArray variables.
   - Added a getter routine for the BlockMesh variables.
   - Add links to the meshes for mesh-based variables.
   - Return blocklist as member names instead of dictionary.
   - Add derived variables to the python reader.
   - Added more metadata to the Block entries.
   - Added support for namevalue blocks
   - Changed namevalue entries to be members.
   - Fix data population of array blocks
   - Fixed calculation of median mesh.
   - Removed broken stpncpy call
   - Added a routine for destroying the SDF stack.
   - Fixed a few reference counting errors in the python module
   - Free up all memory in sdffilter
   - Added an sdf_legacy module for backwards compatibility.
   - Add sdf_legacy to cmake build.
   - Perform a type check before trying to recast 'PyObject *value' as a 'Block *'.
   - Remove an incorrect Py_DECREF(value).
   - Removed incorrect incref to the blocklist dictionary.
   - Removed unnecessary python keyword substitutions.
   - Only modify the dictionary if an entry has changed.
   - Removed stack_destroy() call from SDF_dealloc.
   - Removed unnecessary DECREFs when error handling.
   - Added stitched material as a block type.
   - Added material block to sdf_legacy.py
   - Added grid_id block member.
   - Added {PLAIN,POINT}_DERIVED blocktype support
   - Added a handle to the stack allocator.
   - Fixed a memory leak for derived variables.
   - Fix some variable names.
   - Fix Python setup when building.
   - Removed unused and confusing return values.
   - Removed unnecessary NULL tests before Py_XDECREF
   - Use the field "offset" for station file offsets.
   - Fix compilation using clang
   - Added BlockList type to the module
   - Don't add arrays which contain no data
   - Speed up location of numpy header file
   - Increment the reference count when adding to list.
   - Added an option for purging duplicate block IDs
   - Fixed creation of python install prefix.
   - A few fixes for building the reader with python3.
   - Fix build script for python module
   - Build using libsdf.a
   - Added git commit info to python reader
   - Added git commit info to sdffilter and sdf2ascii
   - Added commit info from SDF library.
   - Added generation of commit_info.h to CMakeLists.txt
   - Changed C library name to libsdfc.a
   - Added error message to build script.


## v4.3.0 to v4.4.2 (2015-01-10)

 * Added data reduction for grid-based variables.
   This change adds the integer parameters "skip_{x,y,z}" and "skip" to the
   subset blocks. If not present the default value is 0. If a parameter is set
   to a positive integer then all grid-based variables using the subset
   restriction will be reduced when being written to file. This is achieved by
   skipping by the specified number of cells in each of the specified
   directions. The "skip" parameter provides a quick method for setting the
   same number of cells to skip in all directions.
   This currently only applies to grid-based variables and is ignored for data
   averages.

 * Added "print_eta_string" for writing time to complete.
   If "print_eta_string" is set to "T" in the control block of an input deck,
   then the current estimated time to completion will be appended to the
   status updates.

 * Added "ndims" constant to input.deck
   The input.deck now accepts an additional constant "ndims" returning the
   total number of dimensions in the simulation. So this is 1, if the
   input.deck was started using epoch1d; 2, if started using epoch2d and 3, if
   started using epoch3d.

 * added variable 'r_xyz' to input.deck parser.
   'r_xyz' will always retrun the distance to coordinate origin.

 * Added "allow_missing_restart" flag.
   This adds the flag "allow_missing_restart" to the control block of input
   decks. The default value is "F". When "restart_snapshot" is specified then
   the simulation first checks that the specified restart dump is valid. If
   the restart dump exists and is valid then it is used to provide initial
   conditions for the simulation. However, if the restart dump does not exist
   or is not usable for some reason then by default the simulation will abort.
   If "allow_missing_restart" is set to "T" then the simulation will not abort
   but will continue to run and use the initial conditions contained in the
   input deck to initialise the simulation.

 * Added "cc" constant for converting volume to cm^3.

 * Changed compiler flags so that none are specified by default.
   The following flags have been renamed in the Makefile:
     PER_PARTICLE_WEIGHT is replaced with PER_SPECIES_WEIGHT
     PARTICLE_PROBES is replaced with NO_PARTICLE_PROBES
     TRACER_PARTICLES is replaced with NO_TRACER_PARTICLES

 * Moved SDF libraries and utilities into a separately maintained repository.

 * Build and link SDF as a standalone library

 * EPOCH now gets all version information from git.

 * Simplified COMMIT string and moved into the src directory.

 * Changed the embedded source and input decks to a newer, more general method

 * Deallocate all memory at the end of the simulation.

 * Added IDL helper scripts.

 * Added a number of enhancements to the sdffilter command-line utility,
   including support for more block types, formatting options for output
   and methods for filtering the data to be output such as 2d slices, etc.

 * Added --derived option to sdffilter.
   If this option is used then the sdffilter routine attempts to load the
   external library and display derived data.

 * Added --exclude flag to the sdffilter utility.
   This flag allows all variables to be printed except for the ones specified.

 * Added asciidoc documentation for the sdffilter and sdf2ascii -line utilities.

 * Added Doxygen comments to the SDF C-library public interface.

 * Changed particle weight calculation to be per cell.
   The particle weights are now calculated on a per cell basis rather than
   trying work out a smoothed version based on the particle shapes. This makes
   the loading slightly more noisy but avoids discrepencies at the boundaries
   when performing particle injection.

 * Set thermal bc arrays when restarting.

 * Added PARSER_CHECKING to the deck parser.

 * Added support for DATABLOCK ant NAMEVALUE blocks to SDF.

 * Added support for calculating derived data from stations.

 * Added routines to the SDF C-library for adding/removing materials in-place,
   adding/removing data blocks in-place and modifying arrays in-place.

 * Added preliminary C-library write support.

 * Added repository and compile-time information to the SDF library.
   This adds a packed byte array of repository differences along with various
   other useful compile-time information. The information is compiled into the
   library file and can be included in output dumps generated by any program
   which links with the library.

 * Added name/value pair, generic, checksummed arbitrary data and
   station_derived block types to SDF.

 * Read/write particle IDs to SDF as integers

 * Added input/output of integer point variables.

 * Added sanity checks for particle species.

 * Added current accumulation optimisations.

 * Added routine that checks for block name truncation.

Bugfixes:

 * Fixed linking stage on some platforms.

 * Added fixes for compiling on ARCHER

 * Added fix for broken Platform-MPI implementation.

 * Fixed compilation when using the Intel compiler.

 * Removed unrequired and troublesome recursive make.

 * Fixed incorrect implementation of Box-Muller.
   If the two random numbers fail to meet the criteria then both should be
   thrown away and regenerated. If one of the previous random numbers is
   re-used then this will skew the distribution. Thanks to Mikhail Garasyov
   who explained the issue on the CCPForge open-discussion list, 2014-12-04.

 * Fixed total energy diagnostic.

 * Fixed laser absorption diagnostics.

 * Bugfix for calculating plasma frequency in 3D.

 * Fixed particle loading when using TOPHAT shape.

 * Fixed bug for grid not being written to file.

 * Add fallback for the python argparse module.

 * Deallocate the prefix_first_call variable.

 * Fixed array dimensions for point mesh.

 * Prevent double deallocations of memory.

 * Fixed typo in 3d current calculation.

 * Fixed a bug when no output blocks are defined.

 * Bugfix for "dump_first" when using file prefixes.

 * Fixed maximum string length in MD5 calculation.

 * Fixed commit string generation.

 * Fixed bug in initialising dist_fn counter.

 * Fixes for parsing unary operators in the input deck.

 * Fixed calculation of initial time.

 * Fixed restarting from dumps with "use_exact_restart"

 * Fixed bug when using "r_xy", "r_yz" or "r_xz".
   These variables worked fine in epoch3d, but using them in a lower
   dimensional version "r_xz" for example became simplified to "x**2" instead
   of "x". This is now been fixed.

 * fixed warning message
   'r_xy' does no longer yield a warning message if used in epoch2d.

 * Bugfix for the simplify stack parser.
   In the previous version, the stack parser would omit any time or space
   varying parameter and substitute it with an arbitrary number. This had some
   unintended consequences, so now the expression is evaluated in full and
   discarded later on.

 * Fixed typos in dist_fn block parser.

 * Fixed typos in basic_evaluate_standard.

 * Fixed typo in 3D laser absorption diagnostics.

 * Fixed typo in 1D probe block parser.

 * Fix time-varying variables when evaluating deck.
   The recent additions which simplify the output from the deck parser
   inadvertantly stomp on the "is_time_varying" flag. This change fixes
   things.

 * Fixed bug with unallocated dump_point_grid.
   dump_point_grid was being assigned to before it was allocated.

 * Fixed a typo in creating id_max restricted subsets.
   This patch was contributed by Matthew Tooley on CCPForge.

 * Fixed sanity checking in the I/O block parser.

 * Fixed the '-DNO_IO' compile-time flag.
   When this flag was enabled, the deck parsing routines were not being called
   properly. This has now been fixed.

 * Fixed compilation when not using PER_PARTICLE_WEIGHT

 * Various small bugfixes and corrections.


## v4.2.0 to v4.3.0 (2014-01-21)

 * Added 'restart' dumpmask constant.

 * Ionization species inherit "dumpmask".
   The species which are automatically created for each ionization
   level are now given the same "dumpmask" or "dump" parameter as
   the parent species. This behaviour can be overridden by explicitly
   adding a species block of the same name with a differing dumpmask.
   eg.

   begin:species
     name = Helium
     ionisation_energies = (24.6*ev,54.4*ev)
     dump = F
   end:species

   begin:species
     name = Helium1
     dump = T
   end:species

 * Rationalised the usage of multiple "output" blocks.
   This set of changes introduces a clear distinction between the
   old style of output block and the new style.

   With the old method, there is a single unnamed output block. It
   has three levels of dump frequency: "normal", "full" and
   "restart". There is a single dump interval specified by either
   "dt_snapshot" or "nstep_snapshot" and the frequency between full
   and restart dumps is specified in relation to this interval using
   "full_dump_every" and "restart_dump_every".

   With the new method, there can be multiple output blocks and each
   one has a "name" field to distinguish it. It is the existence or
   absence of this field which determines whether it is a new or old
   style output block. Each output block has a single dump frequency.
   In order to have multiple dump levels you must create a separate
   output block for each level. Since there is no longer the notion
   of multiple dump levels, the "full_dump_every", etc. flags are
   not supported.
   There are a few output parameters which only make sense globally.
   When using the new style of output block, these are now specified
   by placing them in their own "output_global" block.
   The list of fields which are now specified globally are:
     force_first_to_be_restartable
     force_final_to_be_restartable
     use_offset_grid

 * Made "dump_source_code" a per-output block parameter.
   The flags "dump_source_code" and "dump_input_decks" parameters
   can now be specified individually for each named output block.
   If the flags are not present, they default to "T" for restartable
   output blocks and "F" for non-restartable ones.

 * Added start and stop times to the "output" blocks.
   "output" blocks in the input deck now accept the parameters
   "time_start", "time_stop", "nstep_start" and "nstep_stop".
   These determine the time or step number at which to start
   considering output for a block and similarly the time or step
   number at which to stop considering output for a block.
   For example:
   begin:output
     nstep_snapshot = 4
     nstep_start = 5
     nstep_stop = 12
   end:output

   With this block, the first output will occur at step number 8
   since this is the first multiple of 4 to occur after 5. The last
   dump will occur at step number 12 regardless of how many steps
   the simulation continues to run for.
   Note that if "dump_first" or "dump_last" are set to true for this
   block then dumps will occur first or last regardless of the values
   of "time_start", etc.
   Note also that the flags are applied on a per-output block basis
   in the case of multiple output blocks. A subsequent patch will
   add global applied variants.

 * Added global versions of "time_start", etc. flags.
   This changeset adds "time_start", "time_stop", "nstep_start"
   and "nstep_stop" parameters to the "output_global" input deck
   block.
   These have the same semantics as the per-output block versions
   added in the previous changeset, except that they apply to all
   output blocks in the input deck.
   The per-output block parameters are still honoured and the most
   restrictive version applies.
   For example:
   begin:output_global
     nstep_start = 5
     nstep_stop = 8
   end:output_global

   begin:output
     name = block1
     nstep_snapshot = 1
     nstep_start = 4
     nstep_stop = 7
   end:output

   begin:output
     name = block2
     nstep_snapshot = 1
     nstep_start = 6
     nstep_stop = 9
   end:output

   In this example, the first "block1" output will occur at step 5
   due to the global restriction and the first "block2" output will
   occur at step 6 due to the per-block restriction.
   Meanwhile, the last "block1" output will occur at step 7 due to
   the per-block restriction and the last "block2" output will occur
   at step 8 due to the global restriction.

 * Added global variants of "dump_first", "dump_last"
   This adds "dump_first" and "dump_last" parameters to the
   "output_global" block of the input deck.
   If specified, they override the values of any per-output block
   or default values of these parameters.

 * Added "t_end" and "nsteps" to the deck parser.
   This changeset adds the "t_end" and "nsteps" parameters from the
   "control" block to the input deck parser.
   This enables them to be used within maths expressions elsewhere
   in the input deck. For example, you can now specify that there
   should be 10 output dumps at equal intervals using:

   begin:output
     dt_snapshot = t_end / 10
     nstep_snapshot = nstep / 10
   end:output

   The change also exposes the "nprocx/y/z" parameters although since
   these may be redefined when the MPI grid is setup, they may not
   always give the correct values.

 * Added dump_at_times flag to the "output" block.
   This change adds both the "dump_at_times" and "dump_at_nsteps"
   parameters to "output" blocks in the input deck.
   These specify a list of times/nsteps at which to write the
   current output block. They act in addition to other output
   frequency parameters.
   eg.
   begin:output
     name = test
     nstep_snapshot = 10
     dump_at_nsteps = 3,9,10,18,22
     nstep_start = 7
     ex = always
   end:output

   This will output the "test" block at nsteps 9,10,18,20,22,30,40...

 * Added "dump_cycle" parameter to the "output" block.
   This adds the "dump_cycle" integer parameter to the "output"
   block in the input deck.
   If set to a positive integer then the output file number will
   be reset to zero after the specified cycle number is reached.
   The default is "0", so dump cycling never occurs.
   eg.
   block:output
     nstep_snapshot = 10
     dump_cycle = 2
     ex = always
   end:output

   With this block, output will be written at step 10 to the output
   file "0000.sdf", step 20 to "0001.sdf", step 30 to "0002.sdf",
   step 40 to "0000.sdf", step 50 to "0001.sdf", etc.

 * Allow "restart_snapshot" to accept filenames.
   The "restart_snapshot" parameter in the "control" block of the
   input deck will now accept a filename instead of a cycle
   number.
   If you want to restart from "0012.sdf" then it can either be
   specified using:
     restart_snapshot = 12

   as in previous versions, or alternatively it can be specified
   using:
     restart_snapshot = 0012.sdf

   This syntax will be useful once output file prefixes have been
   implemented.

 * Added the "file_prefix" parameter to output blocks.
   It is now possible to separate output files generated by
   different output blocks by specifying a prefix.
   This is done using the "file_prefix" parameter in an "output"
   block of the input deck.
   If "file_prefix = aa" is specified then the files generated by
   the output block will be named "aa0000.sdf", etc. instead of
   just "0000.sdf". The prefix is not applied to output files
   generated by other output blocks.

   The reason for introducing this change, is that you can now
   write different variables to separate files which was not
   previously possible.

   For example, here are two output blocks which do not use file
   prefixes:
   begin:output
     name = o1
     nstep_snapshot = 1
     charge_density = always
   end:output

   begin:output
     name = o2
     dump_at_nsteps = 10
     restartable = T
   end:output

   With this input deck, we want to have the "charge_density"
   derived variable at every snapshot and then periodically write a
   restart dump. The problem is that the dump file "0010.sdf"
   contains both the restart information and the "charge_density"
   variable. At the end of the run we can't just delete the large
   restart dumps without losing the smaller variables at that
   time step.

   With the new version we would add a prefix to one or both
   blocks:
   begin:output
     name = o1
     file_prefix = small
     nstep_snapshot = 1
     charge_density = always
   end:output

   begin:output
     name = o2
     nstep_snapshot = 10
     restartable = T
   end:output

   Now the "charge_density" will be written to "small0000.sdf", etc.
   At step 10, two files will be written: "small0010.sdf" containing
   just the charge_density and "0000.sdf" containing all the restart
   variables.

   Note that some care must be taken, since if the same variable is
   in the output block for multiple file prefixes then multiple
   copies will be written to file. This obviously uses more disk
   space and is more time consuming than necessary.

   It should also be noted that if multiple output blocks use the
   same file stem then their output will be combined. eg:
   begin:output
     name = o1
     file_prefix = a
     dump_at_nsteps = 2,4
     ex = always
   end:output

   begin:output
     name = o2
     file_prefix = a
     dump_at_nsteps = 3,4
     ey = always
   end:output

   begin:output
     name = o3
     file_prefix = b
     dump_at_nsteps = 4
     ez = always
   end:output

   In this example, at step 2 a0000.sdf contains ex, step 3 a0001.sdf
   contains ey, step 4 a0002.sdf contains ex,ey and b0000.sdf
   contains ez.

 * Added "rolling_restart" flag to the output block.
   If the "rolling_restart" flag is set to "T" in a named output
   block then it automatically sets the output block to have the
   properties required for performing rolling restarts.
   It is currently just a shorthand for setting the following flags:
     dump_cycle = 1
     restartable = T
     file_prefix = roll

   In the future, there may be some additional properties attached
   to rolling restart dumps and these will also be enabled using
   this flag.

 * Added "dump_cycle_first_index" flag to output block.
   This parameter can be used to specify the starting number used
   when output numbering is reset due to the "dump_cycle" parameter.

 * Added command-line options for Makefile.

 * Allow the "phase" parameter to be time varying.

 * Allow the "profile" parameter to be time varying.

 * Added "sdffilter" utility for manipulating SDF files.
   The "sdf2ascii" utility has been developed primarily as a debugging
   tool. Its purpose is to spit out the exact contents of a file.
   Rather than extend this to handle more complex post-processing
   tasks, I have instead added the "sdffilter" utility.
   At the moment it is fairly basic and just handles station data.
   In the future it will handle tasks such as rewriting SDF files,
   performing lineouts, generating derived variables, etc.

 * Added support for constants to the Matlab reader.

 * Added a makefile flag for disabling encoded source.

 * Reduce the number of CPUs if necessary.
   The previous patches enforce a minimum domain size. In some
   situations it may not be possible to divide the simulation
   amongst all the processors and still meet this restriction.
   In such circumstances, we keep reducing the number of processes
   until we find a suitable split.
   On a typical large unattended run this behaviour is better than
   exiting and losing allocation slot.

   To have the simulation exit under such circumstances, add the
   flag "allow_cpu_reduce = F" to the control block.

 * Fixed compilation on IBM Bluegene machines.
   To compile on a Bluegene machine, use the command
   "make COMPILER=ibm"

 * Added default values for ignorable directions.
   This change allows submitting 3D or 2D input decks to a 1D
   version of EPOCH and 3D input decks to a 2D version of EPOCH.
   Any references to y/z will be set equal to zero unless overridden
   by a deck constant. Other y/z values also assume sensible
   defaults, eg. 1 grid cell, 1 metre thick, etc.

 * Added automatic byteswapping to the SDF library.
   The library now checks the endianness of the SDF file and
   byteswaps the data if required.

 * Fixed unpacking script on big-endian machines.

 * Added "sdf_buffer_size" parameter to output_global.
   This adds an input deck parameter to the "output_global" block,
   allowing users to control the size of the buffer used when
   writing point data to SDF files.

 * Added "filesystem" option to output blocks.
   Some filesystems can be unreliable when performing parallel I/O.
   Often this is fixable by prefixing the filename with 'ufs' or
   'nfs'. This change adds an option to both the "output" block
   and "output_global" block enabling the filesystem prefix to be
   used. eg:

   begin:output_global
     filesystem = ufs
   end:output_global

 * Added Makefile flag for compiling on HECToR

 * Added particle gamma and relativistic mass output.
   It is now possible to output per-particle gamma and per-particle
   relativistic mass. These are specified in the same way as other
   particle variables in the "output" block using the parameters
   "gamma" and "relativistic_mass". eg.

   begin:output
     dt_snapshot = 1e-15
     particles = always
     gamma = always
     relativistic_mass = always
   end:output

   There is also now a "rest_mass" parameter which is just an alias
   for "mass", controlling the output of per-particle rest mass.

 * Allow "qed" blocks when -DPHOTONS is not used.
   Previous versions of the code would halt execution if a "qed"
   block was found in the input deck for code which was not compiled
   using the -DPHOTONS preprocessor flag. This can make it difficult
   to use the same input deck with and without QED enabled.
   The code now only halts if "use_qed=T" inside the "qed" block.

 * Added option for QED radiation reaction.
   The "qed" block now accepts the logical flag
   "use_radiation_reaction". If set to "T", the radiation reaction
   force is calculated. If set to "F", it is not calculated.
   The default value is "T".

 * Output boundary thickness to SDF dumps.

 * Added constants to the SDF python reader.

 * Added "disabled" flag for output blocks.
   If "disabled=T" is set in an output block then the block is
   ignored and never generates any output.

 * Add reading of the data directory name from a file.
   EPOCH now first checks for the existence of a file named
   "USE_DATA_DIRECTORY" in the current working directory. If such
   a file exists, it reads it to obtain the name of the data
   directory to use. If no such file exists, it prompts for a
   data directory name as before.
   This is useful for cluster setups in which it is difficult or
   impossible to pipe in the directory name using a job script.

 * Check for the existence of a "STOP" file.
   EPOCH will now check for the existence of a file named either
   "STOP" or "STOP_NODUMP" in the simulation output directory.
   The check is performed at regular intervals and if such a file
   is found then the code exits immediately. If "STOP" is found
   then a restart dump is written before exiting. If "STOP_NODUMP"
   is found then no I/O is performed.

   The interval between checks is controlled by the integer parameter
   "check_stop_frequency" which can be specified in the "control"
   block of the input deck. If it is less than or equal to zero then
   the check is never performed.
   The default value is 10.

 * Stop execution after a given elapsed walltime.
   This adds a two parameters to the "control" block of the input
   deck which allow the code to halt execution and create a
   restart dump after a specified walltime has elapsed.

   Two parameters have been added:

   "stop_at_walltime" is the walltime at which to halt the simulation.
   This is specified as the number of seconds since the start of the
   simulation.

   "stop_at_walltime_file" is the filename from which to read the
   value for "stop_at_walltime". Since the walltime will often be
   found by querying the queuing system in a job script, it may be
   more convenient to pipe this value into a text file rather than
   modifying the input deck.

 * Added "immobile" particle species property.
   "species" blocks now accept the logical flag "immobile". If set
   to "T", the species will be ignored during the particle push.
   The default value is "F".

 * Changed the default value of "dump_first" to "T".

 * Added collisional ionisation

 * Added an option for outputting deck constants.
   This adds a logical flag "print_constants" to the "control"
   block of an input deck. If set to "T", deck constants are printed
   to the "deck.status" file as they are parsed.

 * Automatically write the field grid if needed.
   If any grid variables are written to file then the grid will
   automatically written as well. This can be disabled by adding
   a "grid = never" entry in the output block.

 * Automatically write the particle grid if needed.
   As with field variables, if any particle variable is written to
   file then the corresponding particle grid will automatically
   written as well. This can be disabled by adding a
   "particles = never" entry in the output block.

 * Added "total_energy_sum" diagnostic.

Bugfixes:

 * Bug fix for reading in minor revision numbers.

 * Added SDF compatibility interface for old codes.
   The addition of the minor revision field to sdf_write_run_info()
   meant that codes could not switch between the two versions.
   This patch adds an interface so that old codes do not have to be
   changed.

 * Fix compiler warnings related to the TRANSFER() function.

 * Updated the ASCII EPOCH logo.

 * Reduce module USE's to make ifort-13.x happy.
   Recent versions of the Intel Fortran compiler take hours to
   compile the code if there are recursively included modules.
   Although this is quite clearly a compiler bug, it is not too
   painful to work around.
   It may turn out that this causes issues with other compilers and
   may need to be reverted.

 * Fixed the gen_src_module on OSX 10.8
   The new version of "base64" on OSX-10.8 machines segfaults
   when using the command-line option for folding long lines.
   This change works around the issue by piping the output through
   the "fold" command instead.

 * Fixes for building the COMMIT file.

 * Fix building the reader with VisIt-2.6.1

 * Added commit info to SDF reader plugin.

 * Changed multiplication factor in the Bohm-Gross relation.

 * Made date commands in build scripts more portable.

 * Ignore non-dumped species for derived variables.
   If a particle species has its dumpmask set to "never" or dump=F
   then it should be ignored when writing per-species derived
   variables.

 * Added a couple of aliases to the output blocks.

 * Dump all the information required for QED restart.

 * Various fixes to make restart dumps work correctly.
   Dumps are now performed at the half-timestep to ensure that we
   have all the information required for applying the laser boundary
   conditions.
   We also dump some extra state information for moving windows.

 * Fixed bug in which restart variables were ommitted.
   The list of variables which need to be written into restart files
   was not always being set. This would occasionally lead to restart
   dumps which did not contain all the variables required for a
   restart to occur.

 * Write a unique list to the "*.visit" files.
   Each time an output dump is written, the name of the output
   dump gets written to a "<name>.visit" file, where <name> is
   the name of the corresponding output block.
   When dump cycling is enabled, this could lead to multiple
   repetitions of the same output files in the list.
   The diagnostic routine now keeps a list of the filenames already
   written and checks for duplicates before writing to the *.visit
   file.

 * Don't allow multiple "output" blocks with same name.
   It is an error for the input deck to contain multiple output
   blocks which share the same name. The current changeset properly
   enforces this restriction.

 * Added a more efficient particle shuffling algorithm.

 * Made date commands in build scripts more portable.

 * Ignore non-dumped species for derived variables.
   If a particle species has its dumpmask set to "never" or dump=F
   then it should be ignored when writing per-species derived
   variables.

 * Only write to the status file from rank zero.

 * Added a couple of aliases to the output blocks.

 * Possible fix for generating correct COMMIT string.

 * Ensure that dist_fn gets the right species count.

 * Check for valid species names during deck parse.

 * Prevent possible infinite looping in dist_fn.

 * Fixed bug setting up initial density in parallel.

 * Made laser driven boundaries slightly faster.

 * Fixed per-species field averaging.

 * Stop the output dump times from drifting.

 * Fix per-output block averaging start times.

 * Only read SDF serial input using one process.
   Many of the SDF input routines read in the same data on all
   processes. However, it is much more efficient to read the data
   using one process and then broadcast it to all the others.

 * Added support for station files to the SDF library and VisIt reader.

 * Fix re-init of CPML boundaries during load balance.

 * Fixed the per-particle kinetic energy diagnostic.

 * IOR subset dumpmasks with non-subset ones.
   This change allows 'var = subset_name + species' to work in an
   output block. However, there may be unintended consequences.

 * Fixed compilation issues with MPICH.

 * Check for a minimum size when splitting the domain.
   Communication of the ghost cells is incorrect if the simulation
   domain is smaller than the number of ghost cells. This is fixable
   but it is easier just to ensure that the domain sizes are never
   that small.

 * Fixed missing factor of 2 from BSI ionization.

 * Fixed compilation with some compilers.
   It would seem that it is against the F90 standard to declare
   variables with both POINTER and INTENT attributes.

 * In collisions, skip particles with same momentum.

 * Use normalised momentum for collision tests.

 * Honour the minimum domain size when load balancing.

 * Added integer point variable types to the IDL reader.

 * Exit on unparsable input decks.

 * Fix the units of momenta in distribution functions.

 * Output simulation time in a more user-friendly way.

 * Use Buneman's form of current smoothing.

 * Bugfixes for the collision tests.

 * Avoid issues due to interpolating particle weights.

 * Only convert point data when necessary.
   There was a bug when the SDF library was instructed to convert
   point arrays from double to single precision before writing to
   file. The entire buffer array was converted even though it is
   only necessary to convert the data being written. This led to
   output dumps being much slower than they should be.

 * Add extra debugging flags on OSX.

 * Bug fixes for the lagrangian mesh SDF blocktype.

 * Added lagrangian mesh support to the IDL reader.

 * Fixed the use of mmap by command-line tools.

 * Fix "restrict_z" parameter in EPOCH3D.

 * Fixed typo when sanity checking species block.
   The sanity check for the "identify:" parameter in the species
   block contained a typo which meant that an error was never raised.
   This has now been fixed.

 * Added an additional warning message for QED.
   The code now prints a warning message if QED is turned on and
   both electron and positron species are empty.

 * Only issue a warning on an unrecognised "identify"
   In previous versions, if an unrecognised "identify:" parameter
   was found in the species block, the code would issue an error
   message and halt execution. This can be a bit annoying when
   turning on/off flags such as -DPHOTONS. The code now just issues
   a warning and continues execution.

 * Fix reference counting in the SDF python reader

 * Fixes for angular distribution functions.

 * Added better sanity checks to dist_fn parser.

 * Autodetect gfortran compiler version in Makefile.

 * Allow the VisIt reader to cope with missing grids.

 * Improved performance of CPML boundaries.

 * Improved performance of the field solver.
   The generic form for the arbitrary-order field solver was a fair
   bit slower than explicitly having separate loops for each.

 * Fix {time/nstep}_start for io blocks.
   This change fixes an issue where {time,nstep}_start were not
   both being tested for when calculating output dump times.
   Similarly for {time,nstep}_stop

 * Fixed occasional spikes in autoloaded particles.
   The calculation of per-particle weights was being carried out
   incorrectly at the corners where parallel interfaces intersected
   with the domain boundary. This has now been fixed.

 * Fix derived output variables at domain boundaries.

 * Fixed a few serious regressions in v4.2.12
   Some of the tests for equality which had been changed in commit
   3cd93e3 (SVN 1080) were the wrong way around. This completely broke
   a few things.

 * Added "USE mpi" to ionise.F90, needed by some compilers.

 * Added "species_id" field to SDF point blocks.
   This adds read/write support to the fortran SDF library for
   adding a "species_id" field to the metadata of all point blocks.
   It also adds read support for this field to the C-library.

 * Fixed EPOCH to work with new "species_id" fields.
   The SDF library now requires the use of a "species_id" field
   on all point blocks.

 * Added stack simplifier routines.

 * Added "c_io_never" dumpmask constant.
   Previously, the dumpmask flag "never" was assigned to zero.
   The problem with this is it makes it impossible to tell whether
   the "never" flag was explicitly used in a dumpmask.
   This change might lead to some breakage.

 * Added generic timing routines.

 * Hack to work around issues on filesystems without flock().
   If a filesystem is mounted without support for the flock()
   system call then non-interleaved parallel writes will fail.
   This patch enables a ROMIO hint which is known to work around
   the issue on some systems.
   It is not a true fix because the underlying problem is with the
   filesystem and not EPOCH. It is a bit of a cludge and is not
   guaranteed to work on all MPI implementations.

 * Speeded up encoded source for certain compilers.

 * Fix laser absorption diagnostic.

 * Updated MatLab and IDL readers for new file revision.


## v4.1.0 to v4.2.0 (2013-02-23)

 * Allow the data directory to contain a '/' character.

 * Added unary plus to the deck parser.

 * Added per-particle kinetic energy diagnostic.
   This can be requested using "ek = <dumpmask>".

 * Added an optional error handler for MPI debugging.

 * Allow dump frequency to be reset during restart.
   The code now outputs the time an step number of the last dump for each
   output block, rather than the next one. This allows the calculation of
   the next dump time to be reset immediately after restarting rather than
   only being reset once the next dump occurs.

Bugfixes:

 * Fixed preventing EPOCH3D from running in parallel.
   This fixes a typo which was introduced by the patch applied on Mon Oct
   22 05:41:01 2012.

 * Don't calculate dist_fn for empty species.

 * Fix dumping of particle species subsets.

 * Trim whitespace from the end of ID strings.

 * Added a "long_id" field to the SDF block type.
   If a block ID is longer than "c_id_length" then it will be truncated.
   Once this occurs, it is possible that the shortened ID string does not
   contain enough information to find it again. With this change, the
   original string is stored in the block structure so that up until the
   SDF file is closed the block can still be retrieved.

 * Use an untruncated ID string to find point meshes.
   When writing point variables, we also pass through the ID string of the
   associated point mesh. If this was too long and has been truncated, the
   shortened version may no longer uniquely identify the correct mesh
   block. To avoid this problem, we use the untruncated version of the
   string to request the shortened version from the SDF library


## v4.0.0 to v4.1.0 (2013-01-24)

 * Changed the default for "dump_first" back to "T".

 * Made sample_dist_function significantly faster.

 * Added an option randomise particles in VisIt.

 * Revert to previous default behaviour for IDL.
   The new GUI routines for IDL make use of a more complex hierarchical
   data structure. This was the new default behaviour and broke many
   peoples existing scripts. With this patch, the old behaviour is now the
   default. Two new functions have been introduced:
   "getstruct()" and "explore_struct()" are the same as "getdata()" and
   "explore_data()" except that they return the new style data structure.

 * Allow sdf2ascii to accept multiple variable ids.
   When filtering the output of sdf2ascii using the "-v" flag, it is now
   possible to use the flag multiple times to enable multiple variables to
   be output. Also added a few bugfixes and some tidying.

 * Enable reading of constants in the IDL SDF reader.

 * Added node-centred grid to the python interface.

 * Added drift to the moving window species.

 * Added current CPU split to VisIt reader

 * Added the "use_exact_restart" control block parameter.
   This logical flag may be set to "T" or "F". By default it is "F". If
   this flag is set in conjunction with "restart_snapshot", EPOCH will
   attempt to restart the code using the exact configuration it was using
   when the restart dump was produced. The domain split amongst processors
   will be identical, as will the particle lists. The code will also reset
   the status of the random number generator on each process. Note that the
   flag will be ignored if the number of processors does not match that
   used in the original run.

Bugfixes:

 * Fix final output dump to be non-empty.
   The recent multiple output block support broke the final output dump
   such that it was written but contained no data.

 * Fixed floating overflow bug in collisions routine.

 * Allow restart to work with missing particle species.

 * Rebalance ghost cells as well as computational domain.
   This fixes several problems with boundary conditions, including the
   ability to have an imposed electromagnetic field and also correctly
   calculating laser driven boundaries.

 * Fix current BCs for BSPLINE3 particle weighting.

 * Use correct grid stagger for reflecting BCs.
   When particles are set to reflect from the boundary, the current arrays
   are flipped in the direction perpendicular to that boundary. The
   routines which do this were not taking account of grid staggering
   correctly which led to electron depletion at the boundary. This has now
   been fixed.

 * Fixed VisIt reader build script for older versions.
   Also removed a few compiler warnings and added a fallback to build the
   serial version of the reader if all else fails.

 * Bugfix for rebalancing ghost cells.

 * Fix the "dump_first" flag in "output" blocks.

 * Bugfixes for per-output block data averaging.

 * Fix load balancing in only one direction.

 * Improved logging messages in VisIt reader.

 * Added "cpu_split" block to the output dumps.

 * Allow ionisation to compile without PER_PARTICLE_WEIGHT

 * Sanity checks for multiphoton ionisation.
   The code exits with an error when multiphoton ionisation is requested
   and there are either no lasers or multiple lasers with different
   frequencies.

 * Fix typo in 3d rebalancing of ghost cells.
   This fixes the bug which was found and diagnosed by Anupam Karmakar on
   the CCPForge help forum, 2012-09-11.

 * Fix MPI_Send/Recv so that they use an integer count.
   The send/receive count in the particle transfer calls was using an
   8-byte integer, whereas the MPI standard states that it is a plain
   integer. This broke things on Juelich's BlueGene machine.

 * Ensure correct output dump times for restarted runs.
   Information about when the next output dump should occur must be written
   to the restart dump in order for restarts to behave correctly.

 * Read domain extents from dump file when restarting.
   Restarting from simulations with a moving window broke because the
   domain extents differ from that found in the input deck. This change
   fixes it.

 * Fix broken diagnostics output.
   Laser absorption and poynting flux arrays were not using the iomask
   array which has replaced dumpmask. This patch addresses the issue.

 * Fix for per-species averaging.

 * Fixes for restarting from a previous snapshot.
   When restarting, a load balance call is not immediately followed by the
   particle push so we have to explicitly rebalance the current arrays.

 * Fix for overlapping regions in MPI_Sendrecv.
   When the simulation domain is smaller than the number of ghostcells, the
   memory region being sent from the simulation domain overlaps with the
   memory region occupied by the receiving ghost cells. This was causing
   random errors during communications for some problems. (eg. Bug #782 on
   CCPForge)

 * Removed a few MPI-related memory leaks.

 * Species block "temp" now sets temperature in 3D.

 * Removed "dof" (degrees of freedom) variable.
   The calculation of temperature depends on the degrees of freedom of the
   code. The presence of collisions, initial temperature conditions, etc.
   mean that the degrees of freedom is not always the same as the
   dimensions of the code. To avoid confusion, the parameter has now been
   removed entirely.

 * Fixed the IDL quick_view for multiple timesteps.

 * Fixed the "full_dump_every" flag.

 * Added script for generating tarballs from git.

 * Read initial conditions when restarting.

 * Fix for restarting with moving windows.

 * Added a more helpful error message to the autoloader.

 * Allow species with no particles specified.
   If no particles are specified for a species, the code now prints a
   message indicating that no particles have been loaded and continues to
   run. The previous behaviour was to exit with an error message.

 * Fix automatically calculated dist_fn ranges.

 * Various other minor bugfixes and tweaks


## v3.1.0 to v4.0.0 (2012-07-24)

Additions and changes:

 * Added Sentoku and Kemp collision routine.
   This adds a new output block named "collisions" which accepts the
   following three parameters.

     use_collisions - This is a logical flag which determines whether or not
                      to call the collision routine.
                      If omitted, the default is "true" if any of the
                      frequency factors are non-zero (see below) and "false"
                      otherwise.

     coulomb_log - This may either be set to a real value, specifying the
                   Coulomb logarithm to use when scattering the particles
                   or to the special value "auto". If "auto" is used then
                   the routine will calculate a value based on the properties
                   of the two particles being scattered.
                   If omitted, the default value is "auto".

     collide - This sets up a symmetric square matrix of size nspecies*nspecies
               containing the collision frequency factors to use between
               particle species. The element (s1,s2) gives the frequency
               factor used when colliding species s1 with species s2.
               If the factor is less than zero, no collisions are performed.
               If it is equal to one, collisions are performed normally.
               For any value between zero and one, the collisions are performed
               using a frequency multiplied by the given factor.
               If "collide" has a value of "all" then all elements of the
               matrix are set to one. If it has a value of "none" then all
               elements are set to minus one.
               If the syntax "species1 species2 <value>" is used, then the
               (species1,species2) element of the matrix is set to
               the factor "<value>". This may either be a real number, or
               the special value "on" or "off". The "collide" parameter may
               be used multiple times.
               The default value is "all" (ie. all elements of the matrix
               are set to one).

     For example:

       begin:collisions
         use_collisions = T
         coulomb_log = auto
         collide = all
         collide = spec1 spec2 off
         collide = spec2 spec3 0.1
       end:collisions

     With this block, collisions are turned on and the Coulomb logarithm is
     automatically calculated. All values of the frequency array are set
     to one except (spec1,spec2) is set to minus one (and also (spec2,spec1))
     and (spec2,spec3) is set to 0.1

 * Added QED pair production as described in Duclous et al. PPCF 53 (2011)
   It is enabled using the compiler flag "-DPHOTONS". Additionally,
   the Trident process is enabled using "-DTRIDENT_PHOTONS".

   A new input deck block named "qed" has been added which accepts the
   following parameters:

     use_qed - Logical flag which turns QED on or off. The default is "F".

     qed_start_time - Floating point value specifying the time after which QED
                      effects should be turned on. The default is 0.

     produce_photons - Logical flag which specifies whether to track the
                       photons generated by synchrotron emission. If this
                       is "F" then the radiation reaction force is calculated
                       but the properties of the emitted photons are not
                       tracked. The default is "F".

     photon_energy_min - Minimum energy of produced photons. Radiation
                         reaction is calculated for photons of all energies,
                         but photons with energy below this cutoff are not
                         tracked.  The default is 0.

     photon_dynamics - Logical flag which specifies whether to push photons.
                       If "F" then photons are generated, but their motion
                       through the domain is not simulated and they stay where
                       they were generated. The default is "F".

     produce_pairs - Logical flag which determines whether or not to simulate
                     the process of pair generation from gamma ray photons.
                     Both produce_photons and photon_dynamics must be "T" for
                     this to work. The default is "F".

     qed_table_location - EPOCH's QED routines use lookup tables to calculate
                          gamma ray emission and pair production. If you want
                          to use tables in a different location from the
                          default, specify the new location using this
                          parameter.
                          The default is "src/physics_packages/TABLES".

   QED also requires that the code now know which species are electrons,
   positrons and photons. The species type is specified using a single
   "identify" tag in a species block. To specify an electron the block in the
   deck would look like

     begin:species
       name = electron
       frac = 0.5
       rho = 7.7e29
       identify:electron
     end:species

   Once the identity of a species is set then the code automatically assigns
   mass and charge states for the species.
   Possible identities are:

     electron - A normal electron species. All species of electrons in the
                simulation must be identified in this way or they will not
                generate photons.

     positron - A normal positron species. All species of positron in the
                simulation must be identified in this way or they will not
                generate photons.

     photon   - A normal photon species. One species of this type is needed
                for photon production to work. If multiple species are present
                then generated photons will appear in the first species of
                this type.

     bw_electron - The electron species for pair production by the
                   Breit-Wheeler process. If a species of this type exists
                   then electrons from the pair production module will be
                   created in this species. If no species of this type is
                   specified then pair electrons will be generated in the
                   first electron species.

     bw_positron - As above but for positrons.

     trident_electron - The electron species for pair production by the
                        Trident process. If a species of this type exists then
                        electrons from the pair production module will be
                        created in this species. If no species of this type is
                        specified then pair electrons will be generated in the
                        first electron species.

     trident_positron - As above but for positrons.

     proton - A normal proton species. This is for convenience only and is not
              required by the pair production routines.

   A species should be identified only once, so a "bw_electron" species does
   not need to also be identified as an "electron" species. If the code is
   running with "produce_photons=T" then a photon species must be created by
   user and identified. If the code is running with "produce_pairs=T" then the
   code must specify at least one electron (or bw_electron) species and one
   positron (or bw_positron) species. The code will fail to run if the needed
   species are not specified.

 * Added Alistair Lawrence-Douglas's ionisation routines.
   This adds the following two parameters to the input deck control block:

     use_multiphoton - Logical flag which turns on modelling ionisation by
                       multiple photon absorption. This should be set to "F"
                       if there is no laser attached to a boundary as it
                       relies on laser frequency. The default is "T".

     use_bsi - Logical flag which turns on barrier suppression ionisation
               correction to the tunneling ionisation model for high intensity
               lasers. The default is "T".

   Also two additional parameters to the species block:

     ionisation_energies - When specified, this turns on ionisation modeling.
                           This is an array of ionisation energies starting
                           from the outermost shell. This expects to be given
                           all energies down to the fully ionised ion; if the
                           user wishes to exclude some inner shell ionisation
                           for some reason they need to give this a very large
                           number. Note the ionisation model assumes that the
                           outermost electron ionises first always, and that
                           the orbitals are filled assuming ground state.

     electron_species - Name of the electron species. This can be specified as
                        an array in the event that the user wishes some levels
                        to have a different electron species which can be
                        handy for monitoring ionisation at specific levels.

 * Added Holger's CPML boundary conditions.
   This implementation closely follows that outlined in the book
   "Computational Electrodynamics. The Finite-Difference Time-Domain
   Method.", A. Taflove and S. Hagness (2005).
   See also J. A. Roden and S. D. Gedney, "Convolutional PML (CPML):
   an efficient FDTD implementation of the CFS-PML for arbitrary
   media", Microwave Optical Technol. Lett., vol. 27, 2000.

   CPML boundaries are specified in the input deck by specifying
   either "cpml_outflow" or "cpml_laser" in the boundaries block.
   There are also four configurable parameters.
   "cpml_thickness" is the thickness of the CPML boundary.
   "cpml_kappa_max", "cpml_a_max" and "cpml_sigma_max" are tunable
   parameters which affect the behaviour of the absorbing media.
   See the above references for details.

   eg.

   begin:boundaries
     cpml_thickness = 16
     cpml_kappa_max = 20
     cpml_a_max = 0.2
     cpml_sigma_max = 0.7
     bc_x_min = cpml_laser
     bc_x_max = cpml_outflow
     bc_y_min = cpml_outflow
     bc_y_max = cpml_outflow
   end:boundaries

 * Added continuation lines to the input deck.
   If the input deck contains a '\' character then the rest of
   the line is ignored and the next line becomes a continuation
   of the current one.

 * Added "dump_first" logical flag to the output block.
   Set "dump_first=F" to prevent EPOCH from generating an
   output dump immediately after initialising the simulation.
   The default value is "F".

 * Added "dump_last" logical flag to the output block.
   Set "dump_last=F" to prevent EPOCH from generating an
   output dump just before exiting. The default value is "T"
   if an output block exists in the input deck and "F" otherwise.

 * Added "force_first_to_be_restartable" flag.
   This is a logical flag which can be added to an "output"
   block to determine whether or not the first output dump
   should be a restart dump. The default value is "F".

 * Implemented the "ejected_particles" output parameter.
   If the "ejected_particles" variable is requested in the
   output block then all the particles which have left the
   simulation domain are written to file.
   The parameter is assigned a dumpmask in the same way as
   other variables in the output block.
   Once the data has been written, the ejected particle lists
   are reset and will accumulate particles until the next
   requested output dump.

 * Added global particle IDs.
   Particle IDs are useful if you want to track the progress
   of each particle throughout the simulation. Since they increase
   the size of each particle data structure, they are disabled by
   default and must be enabled using a compiler flag.
   The PARTICLE_ID flag will use an 8-byte integer to represent the
   ID and PARTICLE_ID4 uses a 4-byte integer.
   They are written to file using the "id" flag.
   Eg.
     particles = always
     id = always

   in the "output" block will write particle IDs to file.

   Note: In the current implementation, the particle IDs are
   passed between processors and written to file using REAL numbers.
   This means that in double precision the maximum particle ID
   is 2^53 ~ 10^16. This should be ample for the forseeable future.
   However, if the code is compiled for single precision then the
   maximum ID is 2^24 = 16777216. Probably not big enough.

 * Added a 'subset' block for specifying a subset of particles.
   This adds a new block which can be used to pick a subset of
   particles based on given criteria.
   The subset name can then be used in place of (or in addition to)
   the dumpmask in an 'output' block.

   For example:

   begin:subset
     name = background
     random_fraction = 0.1
     include_species:electron
     include_species:proton
   end:subset

   begin:subset
     name = high_gamma
     gamma_min = 1.3
     include_species:electron
   end:subset

   begin:output
     particles = background + high_gamma + always
     px = background + high_gamma
     py = background
     pz = always
   end:output

   In this example, three 'px' blocks will be written:
   'Particles/background/electron/Px', 'Particles/background/proton/Px'
   and 'Particles/high_gamma/electron/Px'.
   The 'background' blocks will contain 10% of the each species,
   randomly selected.
   The 'high_gamma' block will contain all the electrons with a gamma
   greater than 1.3.
   There will also be 'Particles/background/electron/Py' and
   'Particles/background/proton/Py' block containing y-momentum for the same
   10% random subset of particles.
   Finally, the 'Particles/All/electron/Pz' and 'Particles/All/proton/Pz' will
   contain the z-momentum for all particles.

   The following restrictions can be specified in a 'subset' block:
   random_fraction, {px,py,pz,weight,charge,mass,gamma}_{min,max},
   x_{min,max}, y_{min,max} (2d), z_{min,max} (2d & 3d), dumpmask,
   include_species.

 * Added subset restrictions for particle IDs.
   A range of particle IDs to dump can be specified in the "subset" block
   using "id_min" and "id_max"

 * Added single-precision output for EPOCH.
   This adds a new 'single' field to the dumpmask, allowing
   output variables to be converted to single precision when being
   dumped to file. This is specified on a per-variable basis.

   eg.
   begin:output
     dt_snapshot = 8 * femto

     grid = always
     ex = always
     ey = always + single
   end:output

   In this example, the grid variable 'ex' will be written as a
   double precision array and 'ey' will be converted to single
   precision.

 * Allow averages to be stored in single precision.
   If the 'average_single' dumpmask is supplied to a variable
   then the averaged data will be stored in a single precision
   array throughout the simulation.

 * Enable multiple "output" blocks in the input deck.
   It is now possible to have multiple output blocks, each
   with their own "dt_snapshot" or "nstep_snapshot" and their
   own set of output variables.
   The syntax remains the same as before with the addition of
   "name" and "restartable" fields.
   The "name" field specifies the file name to use for the output
   list. Each time EPOCH generates an output dump, it writes an
   entry into the file "<name>.visit". This can be used to find
   all the output dumps of a specific output block.
   The "restartable" field specifies that the output block should
   generate ouput dumps containing all the information necessary
   to restart a simulation.

 * Added "npart_per_cell" to species blocks.
   This patch adds the ability to specify the number of particles
   per cell to use for the initial particle loading.
   In the current implementation, "npart_per_cell" is a simple
   integer constant. At a later stage this will be extended to
   allow "npart_per_cell" to be a spatially varying function.

   If per-species weighting is used then the value of
   "npart_per_cell" will be the average number of particles per
   cell.  If "npart" or "frac" have also been specified for a species,
   then they will be ignored.

   To avoid confusion, there is no globally used "npart_per_species".
   If you want to have a single value to change in the input deck
   then this can be achieved using a constant block.

 * Added thermal boundaries

 * Ramped density in moving windows.
   When the moving window is active, particles are injected into the
   right-hand side of the domain as the domain moves to the left.
   The density of the particle species being injected may now be a function
   of space and time. The particles will adhere to this functional form.

 * Improve MPI performance for particle redistribution.
   Sending and receiving is now performed using red-black processor
   ordering which is significantly faster.

 * Don't write to disk in the load balancing routine.
   The previous load balancer would write the field variables
   to a file so that each processor can read back its new
   domain. With this patch, the same functionality is obtained
   using MPI calls. This will significantly increase performance.

 * Added per-species and per-subset current diagnostic.

 * Added distribution functions for the angle of particle momentum.
   This adds the flags "dir_xy_angle", "dir_yz_angle" and "dir_zx_angle"
   to the "dist_fn" block. These calculate the distribution of particle
   momentum directions in the X-Y, Y-Z and Z-X planes.

 * Added laser absorption diagnostic.
   This adds a variable named "absorption" to the list of variables that
   can be dumped to file. It accepts a dumpmask in the same manner as other
   output variables. When selected, two numbers will be calculated and
   written to file:
     'Absorption/Laser_enTotal' - The total amount of energy injected into
                                  the simulation by laser boundaries.
     'Absorption/Abs_frac' - The fraction of this laser energy being absorbed
                             by the open boundaries.

 * Added particle prefetching
   Use Intel-specific 'mm_prefetch' calls to load next particle in the list
   into cache ahead of time.
   Enabled using the "-DPREFETCH" compiler flag.

 * Added user-defined delay on particle push.
   This adds the parameter "particle_tstart" to the "control" block which
   specifies the time at which to start pushing particles.

 * Added particle migration between species based on energy.
   This adds two new parameters to the "control" block:

     use_migration - Logical flag which determines whether or not to use
                     particle migration. The default is "F".

     migration_interval - The number of timesteps between each migration event.
                          The default is 1 (migrate at every timestep).

   The following parameters are added to the "species" block:

     migrate - Logical flag which determines whether or not to consider this
               species for migration. The default is "F".

     promote_to - The name of the species to promote particles to.

     demote_to - The name of the species to demote particles to.

     promote_multiplier - The particle is promoted when its energy is greater
                          than "promote_multiplier" times the local average.
                          The default value is 1.

     demote_multiplier  - The particle is demoted when its energy is less
                          than "demote_multiplier" times the local average.
                          The default value is 1.

     promote_density - The particle is only considered for promotion when
                       the local density is less than "promote_density".
                       The default value is HUGE.

     demote_density  - The particle is only considered for demotion when
                       the local density is greater than "demote_density".
                       The default value is 0.

 * Added default compiler flags.

 * Added "milli" to "atto" constants to the maths parser.

 * Added NO_IO pre-processor directive.
   This compiles the code so that no disk I/O is performed.
   Useful for benchmarking.

 * Variable listings in IDL have been moved to a separate routine.
   Instead of the command "data=getdata(0,wkdir='Data',/variables)" you
   should now use:
     list_variables,0,'Data'

 * The global variable "wkdir_global" has been removed from the IDL interface.
   To set the default value for "wkdir" you must now use the command
   "set_wkdir, wkdir". To view the currently set default value, use the
   command "print,get_wkdir()".

 * Added IDL widgets for interactively exploring SDF files.
   Data can now be loaded via a GUI interface using the command:
     data=explore_data('Data',snapshot=5)
   First parameter is optional data directory

   Simple plots can be viewed using the command:
     quick_view,'Data',snapshot=5
   The first parameter is optional data directory

 * Multiple improvements to SDF C-library routines.

 * Added basic python reader for SDF files.

 * Split the grid in a more optimal way.

 * Removed SPLIT_PARTICLES_AFTER_PUSH directive

 * Added command-line arguments to sdf2ascii

 * Added error handling to SDF files.

Bugfixes:

 * Fix for particle load and derived data calculation

 * Fix dump times

 * Corrected "intensity_w_cm2" conversion factor.

 * Fixed bugs in the load balancer.

 * Fix parsing of input decks with no newline at the end.

 * Changed "gen_src_module" to print less errors.

 * Only regenerate embedded code when source code has changed.

 * Fixed some issues building the embedded source module.

 * Fixes for restarting from snapshots.

 * Fixed enforcement of the "density_max" restriction.

 * Enforce positive mass for particle species.

 * Various fixes for the VisIt reader plugin.

 * Write the correct units in output dumps.

 * Fixes for VisIt build script.

 * Improved portability of VisIt build script.

 * Modify CFL condition for high order field solvers.

 * Output Poynting flux when there are no particles.

 * Don't write particle variables when there are no particles.

 * Fixed calculation of dist_fns for parallel runs.
   The spatial ranges for distribution functions was incorrectly
   calculated, leaving gaps between grids when run in parallel.

 * Fix timestep when restarting.
   The value for "dt_plasma_frequency" is calculated using the
   initial conditions which are never used when restarting.
   To fix this, restart dumps now contain the value of
   "dt_plasma_frequency".

 * Fixed zero-gradient boundary conditions.

 * Output dump files using both nstep_snapshot and dt_snapshot.

 * Only require "dt_snapshot" or "nstep_snapshot"
   The deck parsing routines were mistakenly requiring both
   "dt_snapshot" and "nstep_snapshot" to be specified. Only
   one of the two is required.

 * Only require "dt_average" or "nstep_average"
   The deck parsing routines were mistakenly requiring both
   "dt_average" and "nstep_average" to be specified. Only
   one of the two is required.

 * Adjust averaging for "nsteps" less than "dt_snapshot".
   Averaging failed to work if "nsteps" was specified
   such that the job completed before the first proper output
   dump. This code now adjusts the time at which averaging begins
   accumulating data.

 * Fix boundary conditions on moving windows.
   Fixed the bug reported by Jurgen Boker in the CCPForge
   help forum on 31/01/2012.

 * Added equilibrium fields to bc calculations.
   This fixes outflow and laser conditions when there is a background field.

 * Allow missing particle species in restart dumps.

 * Added missing nint,floor,ceil input deck functions.

 * Always allocate the number of particles specified.
   Before this patch, the code would round to the nearest number
   of particles on each process. In cases where there was less than
   one particle per process, no particles would be assigned.
   In the case of tracer particles and a large number of processors,
   it is not unreasonable to want less than one particle per
   processor.

 * Improved loading of particles for moving window.
   The code now uses a real value for the species' npart_per_cell.
   This allows the injected particle population to match the
   one used in the initial conditions when there is not a fixed
   number of particles per cell.

 * Check for input deck strings which are too long.
   The code now checks to see if an input deck string is too
   long for the parser to handle and issues an error message
   which explains how to fix it.

 * Truncate strings which are too long for SDF.

 * Various fixes for the python SDF reader.
   Dictionary keys now use the name field instead of the ID,
   making them behave more like the VisIt reader.
   Useful fields from the header are now added to the dictionary.
   File is closed and memory deallocated on object destruction.

 * Avoid divide-by-zero in particle loading routine.

 * Corrected the degrees of freedom for temperature calculation.

 * Allow dumping of tracer particles for derived variables.
   With this change, tracer particles are ignored when summing
   over all particle species but not ignored when doing per-species
   output. If the user wants to ignore tracer particles for
   per-species output, this can be achieved using the dumpmask for
   that species.

 * Use more consistent naming for variables in dumps.

 * Don't prepend "derived_" tags derived data.

 * Fixes for outflow boundary conditions.

 * Ensure that CPU splits are the same for all variables.
   This fixes the "striping" bug when running VisIt in parallel.

 * Updated constants with current values from NIST.

 * Fixed potential issues with IDL keyword inheritance.

 * Use consistent tag naming in SDF IDL reader.

 * Added condition to catch potential divide by zero in data averaging.

 * Fix current deposition bug for high order splines.
   Since the particle push advances particles one and a half
   timesteps whilst calculating current, it is sometimes necessary
   to use a larger current array during this calculation.

 * Fix 1D VisIt scatter plots.

 * Fixed an issue with the deck parser.
   When converting from an ASCII string to a numeric value, some
   compilers will translate any occurrence of "t" or "f" into
   true or false. This patch explicitly checks that the string is
   a number before trying to convert it.

 * Perform sanity checks on particle boundary conditions.

 * Various other tidying and bugfixing


## v3.0.0 to v3.1.0 (2011-06-23)

Syntax changes:

 * Added synonyms for averaging block parameters.
   Two new names have been added to improve consistency with
   the rest of the output block.
   "dt_average" can be used in place of "averaging_period".
   "nstep_average" can be used in place of "min_cycles_per_average".

 * Added "intensity" synonym for "irradiance".
   Also added missing *_w_cm2 specifiers for both.

 * Use "dumpmask" for probes and particle species.
   Changed the input deck syntax to be consistent with other
   output blocks. "dumpmask" has the same meaning for probes and
   species as it does for distribution functions and output blocks.

 * Implemented "no_sum" dumpmask parameter.
   By default derived particle variables are summed over all
   species. If the "species" dumpmask parameter is given then
   output is also generated on a per-species basis.
   Adding "no_sum" to the dumpmask prevents the output
   of the variable summed over all species.

Additions and changes:

 * Added current smoothing routine to 1d and 3d versions.

 * Issue a warning, not an error on dist_fn range.
   In dist_fn blocks the range for spatial coordinates is never
   used. Therefore the deck parser may as well allow the entry to
   have invalid values.

 * Tidied up example_decks and added new ones.

 * Issue a warning when deck constants conflict with built-in ones.

 * Changes to the behaviour of data averaging.
   If "dt_average" is longer than "dt_snapshot" then dt_average
   will be set equal to dt_snapshot rather than vice-versa.
   Also, if the timestep is too large to satisfy "nstep_average"
   then a warning will be printed rather than adjusting the timestep.

 * Removed the normalised grids from dist_fn output.

 * Added labels and units to dist_fn output.
   Also added support for units in the VisIt plugin.

 * Updated IDL reader to handle labels correctly.

 * Detect file endianness in the IDL reader.

 * Updated the load balancing routine.
   This change takes into account the workload of the field solver
   as well as the particle push. It also does a slightly better job
   of dividing the work amongst processors.
   Will hopefully prove a little more robust.

 * Print a helpful message for deprecated deck strings.

 * Simplified deck parsing process.
   The input deck is now only parsed twice. The first time
   it is called with deck_state set to "c_ds_first" and this pass
   deals with anything that does not require pre-allocated memory.
   On the second pass, deck_state is set to "c_ds_last" and this
   deals with entries that could not be parsed in the first sweep.

 * More particle push optimisations.

 * Added particle heat flux and Poynting flux diagnostics

 * Output a list of each type of output dump for VisIt.
   The VisIt visualisation tool will accept a list of files
   which belong together and make these available for plotting.
   Since "normal", "full" and "restart" dumps each contain different
   data, it is useful to have a list for each type.
   This may also be useful for the other plotting libraries in future.

Bugfixes:

 * Fixed typo in cone.deck examples.

 * Fixed serious bug in particle loading routines.
   A typo in the routine for testing valid cells meant that no
   particles would get loaded.

 * Fixed bug in maths parser for minus signs.
   Deck constants are now identified using a new block type
   constant "c_pt_deck_constant". This needs to be checked for
   when deciding whether a minus sign is a unary or binary minus.

 * Fixed a typo in the random number generator.

 * Don't deallocate initial conditions until after manual_load.

 * Use thermal particles in the "ramp.deck" example.

 * Fix for calculating gamma-1 in dist_fn routines.

 * Bugfixes for the MatLab reader.
   There was a bug in setting up grid variables in MatLab which
   has now been fixed.
   It has also been changed to reduce the amount of output it
   produces.

 * Fixed typo in the restart dump reading routines.

 * Use negative values to test for output.
   The previous version uses HUGE() values for dt_snapshot and
   nstep_snapshot by default. This can lead to overflow.

 * Ignore cells of zero density in the autoloader.

 * Fix long integer definitions in the maths parser.

 * Start accumulating averaged variables at the correct time.
   Since an output file can now be triggered by either "dt_snapshot"
   or "nstep_snapshot", we need to test which one is due to occur
   next and start accumulating an average in time for this.

 * Make particle probes directional.
   This change fixes particle probes so that they only track
   particles which cross the plane in the direction corresponding
   to the normal direction.

 * Change default kinetic energy limits for probes.
   Particle probes only track particles whose kinetic energy
   falls within the specified range ek_min to ek_max. This fix
   changes the default such that particles of any energy are tracked.

 * Various other tidying and bugfixing


## v2.3.0 to v3.0.0 (2011-02-28)

Syntax changes:

 * Added "eps0", "epsilon0" and "mu0" constants. These can now be used in
   place of "epsilonnaught", etc.

 * Renamed "freq" to "omega". Added "lambda".
   Renamed both the input deck parameter and the name used in the
   deck. Print a warning if "freq" is used in the deck.
   Added "frequency" for specifying non-angular frequency.
   Added "lambda" for specifying the laser wavelength in a vacuum.

 * Made particle probe parameters consistent.
   Previous versions had differing ways of specifying probes in
   the input deck for each of epoch{1,2,3}d. They also had different
   implementations in the code for testing if a particle had crossed
   a probe.
   They are now consistent. They are all specified using
   a point in the plane and the normal vector to the plane. Particles
   are recorded if they cross the plane in the direction of the normal
   vector.
   The points and normals are given using the following syntax:
     point = (1,2.5,3)
     normal = (0.5,1,0)

   In 1d, both the above and scalar form are accepted:
     point = 2.5
     normal = 1

 * Renamed "rho", "minrho", etc. to "density", "density_min", etc.

 * There is now no distinction between the "constant" block and the "deo"
   block. The "deo" block is now deprecated and will be removed at some
   point.

 * Removed the (currently unused) neutral_background option.

Syntax additions:

 * Added "supergauss" function to the maths parser.
   This is identical to the existing "gauss" function except that
   it accepts a fourth parameter which is the power to raise the
   argument to.

 * Added "profile" to epoch1d lasers in the deck for consistency.

 * Added temp_{x,y,z}_ev for specifying temperature in elecronvolts.

 * Added "micron" constant. This is to help readability of units given in
   microns. eg. "lambda = 1.06 * micron"

 * Added "nstep_snapshot" entry to the output block of input.decks
   This allows a user to specify the number of timesteps between
   output dumps as well as the simulation time between dumps.
   Both can be specified and they will both be tested for.
   Also made "dt_snapshot" an optional parameter.
   If "dt_snapshot" is not specified then it is not tested for.
   If "nstep_snapshot" is not specified then it is not tested for.

 * Added "dump_source_code" and "dump_input_decks" options.
   These are logical flags which control whether or not source code
   and input decks are written to restart dumps.

Additions and changes:

 * Improved initial particle loading algorithm.
   There are difficulties in smoothly assigning particle weights when
   using per-particle weights and the density profile is discontinuous.
   To fix this isssue the weight assignments in vacuum cells are now
   reflected back into the non-vacuum area.

 * Write per-species particle data

 * Dump output to the new SDF format.

 * Added SDF VisIt, IDL and MatLab readers.

 * Use a less memory-hungry stack implementation.

 * Adjusted polarization angle and defaults for lasers.
   The polarization angle for a laser is now measured in a systematic
   way with respect to the right-hand triad of propagation direction,
   electric and magnetic fields. The previous version was somewhat random.
   If the laser is on x_min then the default E field is in the y-direction
   and the B field is the z-direction. The polarization angle is measured
   clockwise about the x-axis with zero along the y-axis.
   Similarly, for propagation directions:
    y_min: angle about y-axis, zero along z-axis
    z_min: angle about z-axis, zero along x-axis
    x_max: angle anti-clockwise about x-axis, zero along y-axis
    y_max: angle anti-clockwise about y-axis, zero along z-axis
    z_max: angle anti-clockwise about z-axis, zero along x-axis

   Also set new default values for lasers. The default end time is t_end
   and the default profile is 1.0.

 * Moved unnecessary files out of epochXd directories.
   Also added the missing "include" directory to the source code
   which gets packed by gen_src_module.

 * Added comments to the Makefile

 * Check that laser frequency and amplitude have been specified.

 * Allow the code to still run if npart equals zero.

 * Improved consistency of warning messages.

 * Added constants for status file unit numbers.

 * Only write deck parser diagnostics to one file.

 * Only report missing options if they have been requested.

 * Added a globally shared seed number for random number generation.

 * Replaced random number generator with a better algorithm.

 * Added default values for "nsteps", "t_end" and "dt_multiplier".

 * Normalise particle momentum before doing the push.
   This increases the accuracy of calculations and allows the code to
   run in single precision.

 * Various bugfixes


## v2.2.0 to v2.3.0 (2010-09-28)

 * Added optional "use_random_seed" flag.
   If this boolean flag is set to true in the control block of the
   input deck then the pseudorandom number generation used for
   particle placement will be seed using the system time. Otherwise
   a predefined seed will be used, ensuring that results are exactly
   reproducible.

 * Allow the code to run with no particle species.

 * Updated VisIt reader for VisIt 2.x

 * Define x_min, x_max, etc. at grid boundaries.
   The old version defined these values at cell centres which
   meant that the domain ran from x_min-dx/2 to x_max+dx/2.
   For most users this is unexpected behaviour.
   The current fix also allows the code to run when nx=1.

 * Made more input deck elements optional.

 * Fix bug in creation of 1D MPI subtypes.
   This is actually a bug in the OpenMPI romio implementation. It
   was fixed in OpenMPI v1.4.2. This change is a work around for
   the older versions.

 * Fixed moving window bugs.

 * Changed the random number generation algorithm.

 * Fixed temperature diagnostics

 * Fixed particle probe output.
   This was broken as a result of the file locking patch.

 * Added perfectly conducting boundaries.
   These are specified using "bc_x_min=conduct" in the input deck.
   They set the E-parallel and B-perpendicular to zero and the
   gradient of B-parallel to zero in the perpendicular direction.
   Particles reflect off the boundary.

 * Changed the CFD routines to make them easier to use in other codebases.

 * Added CFD reader for Matlab.


## v2.1.0 to v2.2.0 (2010-08-13)

 * The fields_external and species_external blocks have been removed.
   This functionality is now obtained by passing the external filename
   in single or double quotes. eg:

     rho='rhodata.dat'

   will load rho from the external file named 'rhodata.dat'.

 * The input deck has been changed to handle one species per block.
   The old species, speciesn and species_externaln blocks have been removed.

   Introduced new species block. One species block for each species
   containing both basic properties and initial conditions. Example:

   begin:species

     charge=0.0
     mass=1836.2 * 4
     frac=1.0
     name=helium
     dump=T

     rho=if(x gt 0.5e-5 and x lt 2.5e-5,\
	       den_min*exp((x-0.5e-5)/scale),rho(helium))
     rho=if(x gt 2.5e-5, den_max, rho(helium))
     rho=if(x lt 0.5e-5, 0.0, rho(helium))
     rho_min=den_min
     temp=0

   end:species

   Note that previously defined species properties such as "rho" can now be
   referred to using the species name.

   Syntax for new species block is combined mixture of syntax for old species
   block (without the affixed species numbers) and syntax for old speciesn
   blocks.  If you need to go back and modify a species after creation
   (initial conditions changes etc) then you can specify a new species block
   with the same name.
   For example:

   begin:species

     name=helium
     rho=rho(helium) * 10.0

   end:species

   This will modify the "helium" species, leaving all species unchanged except
   for "rho" which is made ten times bigger.

   species_externaln blocks have been eliminated. To specify a property to be
   initialised by the values in an external file you now just give the filename
   in quotes. For example:

   begin:species

     name=helium
     rho='data.dat'

   end:species

   Species used for calculating distribution function output diagnostics are
   now specified by name rather than number. For example:

   begin:dist_fn

     include_species:helium

   end:dist_fn

   Changed the deck parser to allow names with trailing numbers.

 * Allow separate particle and field boundary conditions.
   These are specified using "bc_x_min_field", "bc_x_min_particle", etc.

 * Added "reflect" and "open" boundary conditions to input deck

 * Added "gamma" to the distribution functions.

 * Allow the writing of timestep info to stdout.
   This is enabled by specifying "stdout_frequency" in the input deck
   control block. The info is written after every "stdout_frequency" number
   of timesteps.

 * Added log functions to the maths parser.
   "loge" and "log10" take one argument and give the natural log and base10
   log, respectively. "log_base" takes two arguments, the second being the
   base to use for the logarithm.

 * Added nprocx/y/z parameters to the input deck. These specify the number
   of processors to use in the x/y/z directions.

 * The example input.deck has been moved to the example_decks directory.

 * The timestep data is now written to the CFD file header rather than
   its own separate block. If you have a custom reader then you will need
   to modify it to account for the change.
   All users will need to recompile the VisIt reader plugin.

 * Removed 'nfs:' from filenames. This allows EPOCH to work on a wider range
   of filesystems.

 * Removed partial writes in output_particle routines.
   This avoids file locking which in turn fixes I/O on
   most Lustre filesystems.

 * Fixed coefficients for 6th order field derivatives.

 * Fixed particle weight coefficients.

 * Fix initial temperature distribution calculation.

 * Fixed laser boundary conditions.

 * Clamp currents to zero whenever E field is clamped.

 * Don't load particles into ghost cells.

 * Fixes to correctly restart from restart dumps.

 * Fix compilation with older versions of gfortran.

 * Fixed bug when writing arrays larger than 2GB.

 * Added current to the restart dumps.

 * Write all particle species into restart dumps and not only those which
   have been asked for.

 * Only write warning messages on rank 0

 * Fixed some non-standard F90

 * Various minor bugfixes

 * More small optimisations to particle push.

 * Made 1d,2d and 3d versions more consistent.

 * Added top-hat particle weighting. Not fully tested.

 * Allow grids which cannot be exactly divided across processes.


## v2.0.0 to v2.1.0 (2010-04-16)

 * Various minor bug fixes.

 * VisIt reader changes.
   Read endianness and other new CFD header items.
   Improve startup speed.
   Added more debug info.

 * Fix bugs in the particle push discovered with the two-stream instability
   problem.

 * Fix compilation when PER_PARTICLE_WEIGHT is not used.

 * Added pol_angle for specifying polarisation angle in radians.

 * Use the names "x_min", "x_max", etc. instead of "left", "right"

 * Made the building of encoded_source slightly more robust.

 * Fixed up the encoded_source routines to comply with Fortran90's limit of
   39 continuation lines.

 * Dump files at end of timestep.
