#!/bin/sh
#
# #############  Note  ################################
# Before packing, you should clean the GIT directory as
# $ git clean -f -d -x
# #####################################################
#
# Version ID
#
major=`cat src/include/version_major.h`
minor=`cat src/include/version_minor.h`
patch=`cat src/include/version_patch.h`
vid=`echo ${major}.${minor}.${patch}`
#
mkdir HPhi-${vid}
#
#cp -rf * HPhi-${vid}
rsync --exclude HPhi-${vid} -a * HPhi-${vid}
#
# Build docments
#
cd HPhi-${vid}/doc/ja
make latexpdf
cp ./build/latex/userguide_HPhi_ja.pdf ../../
cd ../en
make latexpdf
cp ./build/latex/userguide_HPhi_en.pdf ../../
#
# Remove some files
#
find ./ -name ".git*" -delete
#
# Pack
#
cd ../../../
tar czvf HPhi-${vid}.tar.gz HPhi-${vid}
rm -rf HPhi-${vid}
