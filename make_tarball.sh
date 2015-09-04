#! /bin/sh

repo=epoch
cur=`pwd`
dir=$(mktemp -d -t$repo)

toplevel=$(git rev-parse --show-toplevel)
cd $toplevel/SDF
sdfdir=$(git rev-parse --git-dir)
cd C
subdir=$(git rev-parse --git-dir)
subdir=$(dirname $subdir)

cd $cur
git init -q $dir/$repo
git push -q --tags $dir/$repo HEAD:tmp
cd $dir/$repo
git checkout -q tmp
git submodule init
git config --replace-all submodule.SDF.url $sdfdir
git submodule update
cd SDF
git submodule init
for d in $(git config --get-regexp 'submodule\.*' | cut -f2 -d\.); do
  git config --replace-all submodule.${d}.url $subdir/$d
done
cd ..

git submodule update --init --recursive
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
rm -rf .git
cd $dir
mv $repo $repo-$cstring
tar -cf - $repo-$cstring | gzip -c > $repo-$cstring.tar.gz
cd $cur
mv $dir/$repo-$cstring.tar.gz .
rm -rf $dir
