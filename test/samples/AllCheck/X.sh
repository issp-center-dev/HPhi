cp input.txt ./MakeDefForHubbard
cd ./MakeDefForHubbard
 sh def_VMC.sh > ../log
 cp *def ../org/
cd ..

cp ./misc/* ./org
cp ./HPhi ./org

echo "Canonical Lanczos begins !"
cp -r org CLanczos
cd ./CLanczos
 perl -w CalMod.pl 0 0  
 ./HPhi EDnamelist.def
 #cat output 
cd ..
echo "Canonical Lanczos ends !"

echo "Canonical TPQ begins !"
cp -r org CTPQ
cd ./CTPQ
 perl -w CalMod.pl 1 0  
 ./HPhi EDnamelist.def
 #cat output 
cd ..
echo "Canonical TPQ ends !"

echo "Canonical FullDiag begins !"
cp -r org CFullDiag
cd ./CFullDiag
 perl -w CalMod.pl 2 0  
 ./HPhi EDnamelist.def
 #cat output 
cd ..
echo "Canonical FullDiag ends !"

echo "Grand Canonical Lanczos begins !"
cp -r org GCLanczos
cd ./GCLanczos
 perl -w CalMod.pl 0 3 
 ./HPhi EDnamelist.def
 #cat output 
cd ..
echo "Grand Canonical Lanczos ends !"

echo "Grand Canonical TPQ begins !"
cp -r org GCTPQ
cd ./GCTPQ
 perl -w CalMod.pl 1 3  
 ./HPhi EDnamelist.def
 #cat output 
cd ..
echo "Grand Canonical TPQ ends !"

echo "Grand Canonical FullDiag begins !"
cp -r org GCFullDiag
cd ./GCFullDiag
 perl -w CalMod.pl 2 3  
 ./HPhi EDnamelist.def
 #cat output 
cd ..
echo "Grand Canonical FullDiag ends !"
