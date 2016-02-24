cp input.txt ./MakeDefForHubbard
cd ./MakeDefForHubbard
 sh def_VMC.sh > ../log
 cp *def ../org/
cd ..

cp ./misc/*    ./org
cp ./misc/*pl  ./
cp ./HPhi   ./org

echo "Canonical Lanczos begins !"
cp -r org CLanczos
cd ./CLanczos
 perl -w CalMod.pl 0 0  
 ./HPhi -e  EDnamelist.def
cd ..
echo "Canonical Lanczos ends !"

echo "Canonical TPQ begins !"
cp -r org CTPQ
cd ./CTPQ
 perl -w CalMod.pl 1 0  
 ./HPhi -e  EDnamelist.def
 cd ./output/
   cp ../Ave.pl .
   perl -w Ave.pl
 cd ..
cd ..
echo "Canonical TPQ ends !"

echo "Canonical FullDiag begins !"
cp -r org CFullDiag
cd ./CFullDiag
 perl -w CalMod.pl 2 0  
 ./HPhi -e  EDnamelist.def
 cd ./output/
   cp ../FT.pl .
   perl -w FT.pl
 cd ..
 #cat output 
cd ..
echo "Canonical FullDiag ends !"

echo ""
echo "Energy"
cat CLanczos/output/zvo_CG_energy.dat > tmp_Lanczos.dat
tail -n 1 CTPQ/output/Ave_SS.dat      > tmp_TPQ.dat
head -n 1 CFullDiag/output/Eigenvalue.dat > tmp_FullDiag.dat
perl -w Energy.pl

echo ""
perl -w CheckHop.pl CLanczos/output/zvo_CG_cisajs.dat > CLanczos_Hop.dat
perl -w CheckCor.pl CLanczos/output/zvo_CG_cisajscktalt.dat > CLanczos_Cor.dat

perl -w CheckHop.pl  CTPQ/output/zvo_cisajs_set0step1990.dat  > CTPQ_Hop.dat
perl -w CheckCor.pl  CTPQ/output/zvo_cisajscktalt_set0step1990.dat  > CTPQ_Cor.dat

perl -w CheckHop.pl  CFullDiag/output/zvo_cisajs_eigen0.dat   > CFullDaig_Hop.dat
perl -w CheckCor.pl  CFullDiag/output/zvo_cisajscktalt_eigen0.dat   > CFullDaig_Cor.dat
echo ""

paste CLanczos_Hop.dat CTPQ_Hop.dat CFullDaig_Hop.dat
echo ""
paste CLanczos_Cor.dat CTPQ_Cor.dat CFullDaig_Cor.dat
echo ""
