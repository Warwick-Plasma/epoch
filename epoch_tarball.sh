#! /bin/sh

repo=epoch
cur=`pwd`
dir=$(mktemp -d -tepoch)

git init -q $dir/$repo
git push -q --tags $dir/$repo HEAD:tmp
cd $dir/$repo
git checkout -q tmp
cstring=$(git describe --abbrev=0 --match v[0-9]* HEAD | cut -c2-)
fullstring=$(git describe --match v[0-9]* HEAD | cut -c2-)

# Append date if this is not a tagged release
if [ "$cstring"x != "$fullstring"x ]; then
  # Only append date if there are actually changes
  git diff --quiet "$(git describe --abbrev=0 --match v[0-9]* HEAD)"
  if [ $? -ne 0 ]; then
    sec=$(git log --pretty=format:%ct -1 HEAD)
    dt=$(date -r $sec +"%d_%m_%Y")
    cstring="${cstring}-$dt"
  fi
fi

/bin/sh epoch1d/src/gen_commit_string
cp COMMIT epoch1d/
cp COMMIT epoch2d/
cp COMMIT epoch3d/
rm -rf COMMIT .git
cd $dir
mv $repo $repo-$cstring
tar -cf - $repo-$cstring | gzip -c > $repo-$cstring.tar.gz
cd $cur
mv $dir/$repo-$cstring.tar.gz .
rm -rf $dir
