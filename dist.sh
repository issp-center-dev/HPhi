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
minor=`cat src/include/version_miner.h`
patch=`cat src/include/version_patch.h`
vid=`echo ${major}.${minor}.${patch}`
#
mkdir HPhi-${vid}
#
cp -rf * HPhi-${vid}
#
# Build docments
#
cd HPhi-${vid}/doc/jp
make -f makefile_doc_jp
cd ../en
make -f makefile_doc_en
cd ../fourier/ja
sed -i -e "s/mathjax/pngmath/g" conf.py
make latexpdfja
make html
cd ../en
sed -i -e "s/mathjax/pngmath/g" conf.py
make latexpdfja
make html
cd ../../../
#
# Remove some files
#
find ./ -name ".git*" -delete
rm dist.sh
rm -rf HPhi-${vid}
#
# Pack
#
cd ../
tar czvf HPhi-${vid}.tar.gz HPhi-${vid}
rm -rf HPhi-${vid}
