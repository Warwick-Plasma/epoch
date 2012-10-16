#! /bin/sh

repo=epoch
cur=`pwd`
dir=$(mktemp -d -tepoch)

git init -q $dir/$repo
git push -q --tags $dir/$repo HEAD:tmp
cd $dir/$repo
git checkout -q tmp
sec=$(git log --pretty=format:%ct -1 HEAD)
dt=$(date -r $sec +"%d_%m_%Y")
cstring=$(git describe --abbrev=0 HEAD)-$dt
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
