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
  HPhi-${vid}.tar.gz || exit $?

# Write the hash of the commit into cmake/git_archive.txt of the tarball,
# unless it is filled in already. "HPhi -v" built from the tarball prints it.
#
# A failure of any step below stops the script with the status of the step.
# The temporary files are removed by the trap, also in that case.
hash=`git rev-parse HEAD` || exit $?
tmpdir=`mktemp -d` || exit $?
trap 'rm -rf "${tmpdir}" "HPhi-${vid}.tar.gz.tmp"' EXIT
tar xzf HPhi-${vid}.tar.gz -C ${tmpdir} || exit $?
archive_txt=${tmpdir}/HPhi-${vid}/cmake/git_archive.txt
if [ ! -f ${archive_txt} ]; then
  echo "ERROR: cmake/git_archive.txt is not found in HPhi-${vid}.tar.gz"
  exit 1
fi
if grep -q Format ${archive_txt}; then
  echo ${hash} > ${archive_txt} || exit $?
  # The tarball is replaced only after the new one is written completely
  COPYFILE_DISABLE=1 tar czf HPhi-${vid}.tar.gz.tmp -C ${tmpdir} HPhi-${vid} || exit $?
  mv HPhi-${vid}.tar.gz.tmp HPhi-${vid}.tar.gz || exit $?
fi
if [ -n "`git status --porcelain --untracked-files=no`" ]; then
  echo 'WARNING: the source has changes which are not committed.'
  echo '         The hash written into the tarball is that of HEAD.'
fi
