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

repo=$(git rev-parse --show-toplevel)
repo=$(basename $repo)
cur=`pwd`
dir=$(mktemp -d -t $repo.XXXXX)

git init -q $dir/$repo
git push -q --tags $dir/$repo HEAD:tmp
cd $dir/$repo
git checkout -q tmp

# Update two levels of submodules
git submodule init
for sub1 in $(git config --local --get-regexp submodule | cut -f2 -d\.); do
  git config --local --replace-all submodule.$sub1.url $cur/$sub1
done
git submodule update

for sub1 in $(git config --local --get-regexp submodule | cut -f2 -d\.); do
  git submodule foreach \
    "git submodule init; \
    for sub2 in \$(git config --local --get-regexp submodule|cut -f2 -d\.); do \
      git config --local --replace-all submodule.\$sub2.url $cur/$sub1/\$sub2; \
    done; \
    git submodule update"
done

cstring=$(git describe --abbrev=0 --match v[0-9]* HEAD | cut -c2-)
fullstring=$(git describe --match v[0-9]* HEAD | cut -c2-)

# Append date if this is not a tagged release
if [ "$cstring"x != "$fullstring"x ]; then
  # Only append date if there are actually changes
  git diff --quiet "$(git describe --abbrev=0 --match v[0-9]* HEAD)"
  if [ $? -ne 0 ]; then
    dt=$(git log --pretty=format:%ci -1 HEAD | cut -f1 -d' ')
    cstring="${cstring}-$dt"
  fi
fi

(cd SDF/VisIt
/bin/sh gen_commit_string.sh)
(cd SDF/FORTRAN
/bin/sh src/gen_commit_string.sh)
(cd SDF/C/src
/bin/sh gen_commit_string.sh)
(cd SDF/utilities
/bin/sh gen_commit_string.sh)
(cd epoch1d
/bin/sh src/gen_commit_string.sh)
cp epoch1d/src/COMMIT epoch2d/src/
cp epoch1d/src/COMMIT epoch3d/src/
rm -rf .git*
find * -name ".git" -o -name ".gitignore" -exec rm -rf {} \;
cd $dir
mv $repo $repo-$cstring
tar -cf - $repo-$cstring | gzip -c > $repo-$cstring.tar.gz
cd $cur
mv $dir/$repo-$cstring.tar.gz .
rm -rf $dir
