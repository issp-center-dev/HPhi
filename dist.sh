#!/bin/sh
#
# #############  README  ##############
# This makes an archive file including static copy of submodules (StdFace)
# and PDF formatted documents docs/userguide_HPhi_(ja|en).pdf
# The output filename is HPhi-${vid}.tar.gz,
# where ${vid} is the version number such as 3.5.0 .
# Before using this, install the following python packages:
#   sphinx
#   sphinx_numfig
#   sphinxcontib_spelling
#   git-archive-all
# #####################################

git-archive-all --help > /dev/null 2>&1
if [ $? != 0 ]; then
  echo 'ERROR: git-archive-all is not available'
  echo 'HINT: python3 -m pip install git-archive-all'
  exit 1
fi

#
# Version ID
#
major=`awk '$2=="HPHI_VERSION_MAJOR"{print $3}' src/include/version.h`
minor=`awk '$2=="HPHI_VERSION_MINOR"{print $3}' src/include/version.h`
patch=`awk '$2=="HPHI_VERSION_PATCH"{print $3}' src/include/version.h`
vid=`echo ${major}.${minor}.${patch}`

ROOTDIR=`pwd`

# Build docments
cd $ROOTDIR/doc/ja
make latexpdf
cp ./build/latex/userguide_HPhi_ja.pdf $ROOTDIR/doc
cd $ROOTDIR/doc/en
make latexpdf
cp ./build/latex/userguide_HPhi_en.pdf $ROOTDIR/doc
cd $ROOTDIR/doc/tutorial/en
make latexpdf
cp ./build/latex/tutorial_HPhi_en.pdf $ROOTDIR/doc

# Make a tarball
cd $ROOTDIR
git-archive-all \
  --extra=doc/userguide_HPhi_ja.pdf \
  --extra=doc/userguide_HPhi_en.pdf \
  --extra=doc/tutorial_HPhi_en.pdf \
  --prefix=HPhi-${vid} \
  HPhi-${vid}.tar.gz

# Write the hash of the commit into cmake/git_archive.txt of the tarball,
# unless it is filled in already. "HPhi -v" built from the tarball prints it.
hash=`git rev-parse HEAD`
tmpdir=`mktemp -d`
tar xzf HPhi-${vid}.tar.gz -C ${tmpdir}
if grep -q Format ${tmpdir}/HPhi-${vid}/cmake/git_archive.txt; then
  echo ${hash} > ${tmpdir}/HPhi-${vid}/cmake/git_archive.txt
  COPYFILE_DISABLE=1 tar czf HPhi-${vid}.tar.gz -C ${tmpdir} HPhi-${vid}
fi
rm -rf ${tmpdir}
if [ -n "`git status --porcelain --untracked-files=no`" ]; then
  echo 'WARNING: the source has changes which are not committed.'
  echo '         The hash written into the tarball is that of HEAD.'
fi
