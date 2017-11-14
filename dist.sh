#!/bin/sh -e
#
echo
echo "===  HPhi distribution tool  ==="
echo
echo "  Did you clean the GIT directory? e.g."
echo "  $ git clean -f -d -x"
echo ""
echo -n "  Start ? yes/no [no] : "
read yesno
if [ -z ${yesno} ] || [ ! ${yesno} = "yes" ]; then
    exit
fi
#
mkdir HPhi-${1}
#
cp -rf * HPhi-${1}
#
# Build docments
#
cd HPhi-${1}/doc/jp
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
#
# Pack
#
cd ../
tar czvf HPhi-${1}.tar.gz HPhi-${1}
rm -rf HPhi-${1}
