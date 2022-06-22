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
major=`cat src/include/version_major.h`
minor=`cat src/include/version_minor.h`
patch=`cat src/include/version_patch.h`
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
