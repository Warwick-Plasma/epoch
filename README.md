# *** PLEASE READ THIS NOTE ***

If you are obtaining this code from the github server *DO NOT* use the
"Download zip" buttons. They will only download the EPOCH code and not the SDF
subdirectory.

The latest version of the code can be found
[here](https://github.com/Warwick-Plasma/epoch/releases)

For similar reasons, if cloning the repository using git you *MUST* add the
"--recursive" flag.

## DOWNLOADING AND BUILDING THE CODE

The easiest method of downloading the code is to grab a copy of the latest
release package, located in the "Releases" section of the GitHub site.
You can get there by following
[this link](https://github.com/Warwick-Plasma/epoch/releases)

All EPOCH development is carried out using the [git](https://git-scm.com)
revision control system. The master repository is hosted on a web-based
collaborative development site
[here](https://github.com/Warwick-Plasma/epoch)
Further details on this are provided below.

## Compiling the code

The "Releases" section of the GitHub site contains files in the form
"epoch-4.4.5.tar.gz". These are tarred and gzipped packages of the code
which can be unpacked with the command `tar xzf epoch-4.4.5.tar.gz`.
This will create a directory called "epoch-4.4.5".

Within this directory there are various epochXd subdirectories, each of which
has a "Makefile" containing the instructions for compiling the code.

Many people will be used to editing Makefiles by hand in order to set them up
for their own working environment. However, this is not the recommended way
of working with the EPOCH codebase. In theory, all the changes necessary for
compiling EPOCH on any given environment should be possible using command-line
variables.

For most setups, it should only be necessary to set the COMPILER variable to
correspond to the Fortran compiler to be used. This can either be set as an
option to the "make" command or as an environment variable.

For example, to compile the code using Intel's "ifort" Fortran compiler, you
can either type the following:
```
  $> make COMPILER=intel
```

or alternatively:
```
  $> export COMPILER=intel
  $> make
```

In these two examples `$>` represents the terminal's command-line prompt.
After compilation is complete, the binary file will be created in "bin/epochXd",
where X is 1, 2 or 3.

Since most people will always be using the same compiler on a specific machine,
it is often easiest to just add the line `export COMPILER=intel` to your shell
script initialisation file (ie. "$HOME/.bashrc" on most modern UNIX machines).

The actual compiler command used is the MPI fortran wrapper. On nearly all
machines, this is called "mpif90" and so this is what is used by default.
Occasionally, some machines will call it something different. For example,
Intel's MPI wrapper script is sometimes called "mpiifort". If your machine has
such a setup, you can set the name of the wrapper script using the variable
MPIF90. For example:
```
  $> make COMPILER=intel MPIF90=mpiifort
```

Again, it is often easiest to add the line `export MPIF90=mpiifort` to your
$HOME/.bashrc file.

Finally, there are two more variables which can be used to control the options
used when building the code.

Setting "MODE=debug" will build the code with optimisation disabled and
debugging flags turned on. If this variable is not set then the default is to
build a fully optimised version of the code.

There are several pre-processor flags which can be passed at compile-time to
change the behaviour of the code. These flags are described in the Makefile
with the lines beginning "#DEFINES += " and they are commented out by default.
Rather than uncommenting these lines, it is possible to set them on the
command-line using the "DEFINE" variable. For example, to compile a
single-precision version with global particle IDs you would type:
```
  $> make DEFINE="-DPARTICLE_ID -DSINGLE"
```


## COMPILING SDF AND THE VISIT READER

The EPOCH codes use a self-describing file format called SDF. The routines
used in reading and writing such files, along with reader plugins for Matlab,
IDL, python and VisIt are contained in the SDF directory.

The library used by EPOCH for reading and writing the files is automatically
built when you first build EPOCH. However, it is important to note that
whenever you rebuild EPOCH, the SDF library is NOT rebuilt by default. It is
also not removed when you type "make clean". Most of the time, this is what
you want since rebuilding the library adds a significant amount of time to
the compilation of the code. However, occasionally you might want to force the
library to be rebuilt, such as when you switch compilers. To accomplish this
you must first type "make cleanall" which will remove the existing library and
it will then get rebuilt next time you type "make".

In order to visualise data using the VisIt program, you must first build the
SDF VisIt reader plugin. As a pre-requisite, you must have the VisIt binary
in your shell's search path. You can check this by typing:
```
  $> visit -version
```
which should return with a message such as "The current version of VisIt is .."
If instead you get "visit: command not found" then you may need to edit your
PATH environment variable appropriately. Your system administrator should be
able to help.
Next you will need to ensure that you have a C++ compiler (preferably GNU g++)
and CMake. Again, you can check these using `g++ --version` and
`cmake -version`. Note that the appropriate version of these utilities may
depend on the version of VisIt that you have installed.

Once these pre-requisites have been met, you should be able to build the
reader plugin by typing `make visit`. You do not need to do this again unless
you upgrade your version of the VisIt utility. It is rare that any changes to
EPOCH will require an update of the VisIt reader, but if you do start to
encounter errors when reading SDF files then you can try rebuilding the reader
using the commands `make visitclean` followed by `make visit`.

Note that older versions of EPOCH used the CFD format. This is now obsolete
and current versions of the code no longer contain any reader plugin for this
format. However, existing installations of the VisIt CFD plugin will happily
co-exist with the SDF plugin and issuing `make visitclean` will not remove
such plugins.


## WORKING WITH THE GIT REPOSITORY

For more advanced users, the code is also hosted on a git repository. There is
quite a steep learning curve for using git, so using this repository is only
recommended for more advanced users who are comfortable that they can deal with
a "git conflict".

One other added complication, is that the EPOCH repository also uses git
submodules for tracking the SDF file format. This adds an extra source of
possible issues. However, once a workflow is established it can all be quite
straightforward to work with.

To perform an initial checkout of the code using git, you should issue the
following command:

```
  git clone --recursive https://github.com/Warwick-Plasma/epoch.git
```

The "--recursive" flag ensures that not only the "epoch"
repository is checked out, but also the "SDF" submodules.

It is recommended that after checking out a copy of the git repository, users
immediately create a new working branch and leave the default "master" branch
untouched. A new branch can be created and switched to with the command
`git checkout -b work`.

When you wish to update to the latest version, do the following sequence of
actions. First, commit or stash any changes you have made in your "work"
branch. Next, switch to the "master" branch with
`git checkout master`. Now pull the changes with `git pull`,
followed by `git submodule update --recursive`.
At this stage your "master" branch should be fully up to date.

Merging the new version in with your "work" branch is prone to error, so it
is recommended that you create a temporary copy of this branch just in case
everything goes wrong. The command "git branch workold work" will
create a branch named "workold" which is just a copy of "work". This branch
can be deleted once the merge is completed successfully. If everything goes
wrong in the "work" branch, you can reset it back to the original using the
command `git reset --hard workold`.

In order to update your work branch, switch back to it with
`git checkout work` and merge in the changes with `git merge master`.
After issuing this last command, there is a fair chance that you will encounter
conflicts. You must now resolve those conflicts and commit the changes.
After successfully merging in the changes, you can now delete the temporary
copy of your work branch with `git branch -D workold`.

