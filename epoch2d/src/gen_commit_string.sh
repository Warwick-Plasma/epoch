#!/bin/sh

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

# Generate a COMMIT file containing a string which describes the
# commit status of the current source tree.

# The commit string has five parts. The first three come from git-describe
# And correspond to tag, commits since tag and commit hash.
# The fourth part is the subversion revision number and the fifth
# part indicates wether there are uncommitted changes

# If any of these fields is unavailable then the tag 'unknown' will be
# used. If none of the tags are known then the COMMIT file is left untouched

COMMIT_FILE=src/COMMIT

LF='
'

state='clean'
commit_string=""

git show-ref > /dev/null 2>&1
if [ $? -eq 0 ]; then
# in a git repo
  gitdescribe=$(git describe --match "v[0-9]*" --long HEAD 2>/dev/null)

  if [ $? -ne 0 ]; then
    always=$(git describe --match "v[0-9]*" --always --long HEAD 2>/dev/null)
    gitdescribe=unknown-unknown-g$always
  fi

  git update-index -q --refresh
  test -z "$(git diff-index --name-only HEAD --)" || state='dirty'

  commit_string=$gitdescribe-$state
else
# not in a git repo
  grep "COMMIT=" $COMMIT_FILE > /dev/null 2>&1
  [ $? -eq 0 ] && exit
  commit_string=unknown-unknown-unknown-unknown
fi

[ -z $commit_string ] && exit

grep "$commit_string" $COMMIT_FILE > /dev/null 2>&1
if [ $? -eq 0 ]; then
  exit
else
  echo "COMMIT=$commit_string" > $COMMIT_FILE
  exit
fi
