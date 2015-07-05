#! /bin/sh

failed=""
makeflags="-j8"
build_epoch=0
build_c=0
build_utilities=0
build_visit=0
build_any=0
run_tests=0

while getopts euvth opt
do
   case $opt in
      e) build_epoch=1 ; build_any=1 ;;
      u) build_c=1 ; build_utilities=1 ; build_any=1 ;;
      v) build_visit=1 ; build_any=1 ;;
      t) run_tests=1 ; build_any=1 ;;
      h) cat <<EOF
build_all script options:
  -e: Build EPOCH components
  -u: Build SDF utilities
  -v: Build VisIt reader
EOF
         exit ;;
   esac
done

if [ $build_any -eq 0 ]; then
   build_epoch=1
   build_c=1
   build_utilities=1
   build_visit=1
   run_tests=1
fi

if [ $build_epoch -ne 0 ]; then
   dir=epoch1d
   cd $dir
   echo Building ${dir}...
   make cleanall
   make $makeflags || failed="$failed $dir"
   if [ $run_tests -ne 0 ]; then
      cp example_decks/two_stream.deck Data/input.deck
      echo Data | mpirun -np 4 ./bin/$dir || failed="$failed $dir/tests"
      rm -f Data/*
   fi
   cd -

   dir=epoch2d
   cd $dir
   echo Building ${dir}...
   make clean
   make $makeflags || failed="$failed $dir"
   cd -

   dir=epoch3d
   cd $dir
   echo Building ${dir}...
   make clean
   make $makeflags || failed="$failed $dir"
   cd -
fi

if [ $build_c -ne 0 ]; then
  dir=SDF/C
  cd $dir
  echo Building ${dir}...
  make clean
  make $makeflags || failed="$failed $dir"
  cd -
fi

if [ $build_utilities -ne 0 ]; then
  dir=SDF/utilities
  cd $dir
  echo Building ${dir}...
  sh -x build -r || failed="$failed $dir"
  cd -
fi

if [ $build_visit -ne 0 ]; then
  dir=SDF/VisIt
  cd $dir
  echo Building ${dir}...
  sh -x build -r || failed="$failed $dir"
  cd -
fi


if [ "$failed"x = x ]; then
  echo All components built successfully
  exit 0
else
  echo FAILURE: The following components failed to build - $failed
  exit 1
fi
