#! /bin/sh

# Copyright (C) 2009-2019 University of Warwick
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

GIT_WORK_TREE=$1
GIT_DIR=$(git rev-parse --git-dir 2>/dev/null)

# The following complicated section simulates the use of "readlink -f" to find
# the source directory on platforms which don't support the "-f" flag (eg. OS X)
curdir=$(pwd -P)
target="$GIT_DIR"
cd "$target"
target=$(basename "$target")
# Iterate down a (possible) chain of symlinks
n=0
while [ -L "$target" ]; do
  n=$((n+1))
  if [ $n -gt 1000 ]; then
    echo "ERROR finding source directory."
    exit 1
  fi
  target=$(readlink $target)
  cd $(dirname "$target")
  target=$(basename "$target")
done
GIT_DIR=$(pwd -P)
cd "$curdir"

cd $GIT_WORK_TREE
cd $OLDPWD
export GIT_WORK_TREE GIT_DIR
shift

prefix=$1
pack_source_code=$2
pack_git_diff=$3
pack_git_diff_from_origin=$4
generate_checksum=$5
f77_output=$6

# Use python script by default and then fall back to shell script if that
# fails
BASEDIR=$(dirname $0)
python $BASEDIR/pack.py -- "$@"
if [ $? -eq 0 ]; then
  exit
fi
shift 6

echo WARNING: pack.py script failed. Falling back to shell script.

archive="source_info_archive.tgz"
hexdump="source_info_hexdump.txt"
gitdiff="source_info_gitdiff.txt"
varname="${prefix}_bytes"
diffname="${prefix}_diff_bytes"
module_name="${prefix}_source_info"
outfile="$1"

commitfile="$GIT_WORK_TREE/src/COMMIT"

git_version=$(git describe --always --long --dirty 2>/dev/null)
if [ $? -ne 0 ]; then
  git_version=' '
  if [ -f $commitfile ]; then
    . $commitfile
    git_version=$COMMIT
  fi
  pack_git_diff=0
fi

compile_date=$(date "+%s")
compile_date_string=$(date "+%Y-%m-%d-%H:%M:%S")
compile_machine_info=\
"$(uname -n) $(uname -s) $(uname -r) $(uname -m) $(uname -p)"
compiler_info="$2"
compiler_flags="$3"
if [ "$compiler_info"x = x ]; then
  compiler_info=' '
fi
if [ "$compiler_flags"x = x ]; then
  compiler_flags=' '
fi

ncont=39 # Maximum continuation lines allowed in F95
continuation_lines=39
ncolumns=130
nl=6
nbytes=8
nelements=0
compile_machine_info=$(expr substr "$compile_machine_info" 1 $ncolumns)
compiler_info=$(expr substr "$compiler_info" 1 $ncolumns)
compiler_flags=$(expr substr "$compiler_flags" 1 $ncolumns)

write_data_bytes () {
  filename=$1
  varname=$2

  filesize=$(cat $filename | wc -c)

  if [ $nbytes -eq 4 ]; then
    hexdump -e $((nl-1))'/4 "z'\''%08x'\''," 1/4 "z'\''%08x'\''\n"' \
      $filename | sed "s/,z' .*//" > $hexdump
  else
    # Sadly, some versions of hexdump don't support 8-byte integers
    xxd -c$((8*nl)) -ps $filename | \
      sed "s/\([0-9a-f]\{2\}\)\([0-9a-f]\{2\}\)/\2\1/g;
           s/\([0-9a-f]\{4\}\)\([0-9a-f]\{4\}\)/\2\1/g;
           s/\([0-9a-f]\{8\}\)\([0-9a-f]\{8\}\)/z'\2\1',/g;
           s/,$//; s/,\([^z].*\)/,z'\1'/; s/\(^[^z].*\)/z'\1'/" > $hexdump
  fi

  nlines=$(cat $hexdump | wc -l)
  if [ $nlines -eq 0 ]; then
     nelements=0
     padding=0
  else
    nfull_segments=$(((nlines-1)/ncont))
    nfull_segments=$((nfull_segments*ncont))
    nlast_segment=$((nlines-nfull_segments))
    nlast=$(tail -n 1 $hexdump | tr 'z' '\n' | wc -l)
    nlast=$((nlast-1))
    nelements=$((nl*(nlines-1)+nlast))
    padding=$((nelements*nbytes-filesize))
  fi

  rm -f $filename

cat >> $outfile <<EOF
  CHARACTER(LEN=*), PARAMETER :: ${vname}_mimetype = '$mimetype'
  INTEGER($nbytes) :: $varname($nelements)
  INTEGER, PARAMETER :: ${varname}_padding = $padding
  INTEGER, PARAMETER :: ${varname}_len = $nelements
EOF

  IFS="
"
  i=1
  i0=1
  i1=$((i0+continuation_lines-1))
  n0=1
  n1=$((n0+nl*continuation_lines-1))
  if [ $nfull_segments -gt 0 ]; then
    for ln in `head -n $nfull_segments $hexdump`; do
      if [ $i -eq $i0 ]; then
        echo "DATA($varname(i),i=$n0,$n1)/&" >> $outfile
        echo "$ln,&" >> $outfile
        n0=$((n1+1))
        n1=$((n0+nl*continuation_lines-1))
      elif [ $i -ge $i1 ]; then
        echo "$ln/" >> $outfile
        i0=$((i1+1))
        i1=$((i0+continuation_lines-1))
        [ $i1 -gt $nlines ] && i1=$nlines
      else
        echo "$ln,&" >> $outfile
      fi
      i=$((i+1))
    done
  fi

  nfull=$((nlast_segment-1))
  n1=$((n0+nl*nfull+nlast-1))
  i1=$((i0+nfull))
  for ln in `tail -n $nlast_segment $hexdump`; do
    if [ $i -eq $i0 ]; then
      echo "DATA($varname(i),i=$n0,$n1)/&" >> $outfile
      if [ $i -eq $i1 ]; then
        echo "$ln/" >> $outfile
      else
        echo "$ln,&" >> $outfile
      fi
      n0=$((n1+1))
      n1=$((n0+nl*continuation_lines-1))
    elif [ $i -ge $i1 ]; then
      echo "$ln/" >> $outfile
      i0=$((i1+1))
      i1=$((i0+continuation_lines-1))
      [ $i1 -gt $nlines ] && i1=$nlines
    else
      echo "$ln,&" >> $outfile
    fi
    i=$((i+1))
  done

  rm -f $hexdump
}


shift 3
filelist="$*"
if [ "$filelist"x = x ]; then
  pack_source_code=0
fi


get_bytes_checksum () {
  checksum_type=' '
  checksum=' '
  if [ $generate_checksum -ne 0 ]; then
    files="$*"
    checksum_type='sha256'
    checksum=$(cat $files | shasum -a 256 2>/dev/null | cut -f1 -d' ')
    if [ "$checksum"x = x ]; then
      #echo shasum failed. Trying openssl.
      checksum=$(cat $files | openssl dgst -sha256 2>/dev/null | sed 's/.* //')
      if [ "$checksum"x = x ]; then
        #echo openssl failed. Trying sha256sum.
        checksum=$(cat $files | sha256sum 2>/dev/null | cut -f1 -d' ')
        if [ "$checksum"x = x ]; then
          checksum=' '
        fi
      fi
    fi
    if [ "$checksum"x = x ]; then
      checksum_type=' '
      echo WARNING: unable to create valid SHA-256 checksum
    fi
  fi
}


cat > $outfile <<EOF
MODULE $module_name

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: ${varname}_git_version = '$git_version'
  INTEGER, PARAMETER :: ${varname}_compile_date = $compile_date
  CHARACTER(LEN=*), PARAMETER :: ${varname}_compile_date_string = &
'$compile_date_string'
  CHARACTER(LEN=*), PARAMETER :: ${varname}_compile_machine_info = &
'$compile_machine_info'
  CHARACTER(LEN=*), PARAMETER :: ${varname}_compiler_info = &
'$compiler_info'
  CHARACTER(LEN=*), PARAMETER :: ${varname}_compiler_flags = &
'$compiler_flags'

EOF

if [ $pack_source_code -ne 0 -o $pack_git_diff -ne 0 ]; then
cat >> $outfile <<EOF
! F95 maximum line width of 132 characters, maximum of 39 continuation lines
  INTEGER, PRIVATE :: i
EOF
fi


vname=$varname
checksum_type=' '
checksum=' '
if [ "$filelist"x != x ]; then
  get_bytes_checksum $filelist
fi
cat >> $outfile <<EOF
  CHARACTER(LEN=*), PARAMETER :: ${vname}_checksum_type = '$checksum_type'
  CHARACTER(LEN=*), PARAMETER :: ${vname}_checksum = '$checksum'
EOF
if [ $pack_source_code -eq 0 ]; then
  mimetype=' '
cat >> $outfile <<EOF
  CHARACTER(LEN=*), PARAMETER :: ${vname}_mimetype = '$mimetype'
  INTEGER($nbytes) :: $vname(0)
  INTEGER, PARAMETER :: ${vname}_padding = 0
EOF
else
  tar czf $archive $filelist
  mimetype='application/x-tar-gz'

  write_data_bytes $archive $vname
fi


vname=$diffname
checksum_type=' '
checksum=' '
cat >> $outfile <<EOF
EOF
if [ $pack_git_diff -eq 0 ]; then
  mimetype=' '
cat >> $outfile <<EOF
  CHARACTER(LEN=*), PARAMETER :: ${vname}_checksum_type = '$checksum_type'
  CHARACTER(LEN=*), PARAMETER :: ${vname}_checksum = '$checksum'
  CHARACTER(LEN=*), PARAMETER :: ${vname}_mimetype = '$mimetype'
  INTEGER($nbytes) :: $vname(0)
  INTEGER, PARAMETER :: ${vname}_padding = 0
EOF
else
  if [ $pack_git_diff_from_origin -ne 0 ]; then
    git diff --exit-code origin/main > $gitdiff
  else
    git diff --exit-code > $gitdiff
  fi

  if [ $? -ne 0 ]; then
    cp $gitdiff $gitdiff.tmp
    filterdiff $gitdiff.tmp > $gitdiff
    if [ $? -ne 0 ]; then
      mv $gitdiff.tmp $gitdiff
    fi
    get_bytes_checksum $gitdiff
    gzip -c $gitdiff > ${gitdiff}.gz
    rm $gitdiff
    gitdiff=${gitdiff}.gz
  fi
  mimetype='application/x-gz'
cat >> $outfile <<EOF
  CHARACTER(LEN=*), PARAMETER :: ${vname}_checksum_type = '$checksum_type'
  CHARACTER(LEN=*), PARAMETER :: ${vname}_checksum = '$checksum'
EOF

  write_data_bytes ${gitdiff} $vname
fi

echo "END MODULE $module_name" >> $outfile
