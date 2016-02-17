cp input.txt ./MakeDefForHubbard
cd ./MakeDefForHubbard
 sh def_VMC.sh > ../log
 cp *def ../org/
cd ..

cp ./misc/*    ./org
cp ./misc/*pl  ./
cp ./HPhi ./org

echo "Canonical Lanczos begins !"
cp -r org CLanczos
cd ./CLanczos
 perl -w CalMod.pl 0 0  
 ./HPhi EDnamelist.def
cd ..
echo "Canonical Lanczos ends !"

