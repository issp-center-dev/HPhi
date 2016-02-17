cp input.txt ./MakeDefForHubbard
cd ./MakeDefForHubbard
 sh def_VMC.sh > ../log
 cp *def ../org/
cd ..

cp ./misc/*    ./org
cp ./misc/*pl  ./
cp ./HPhi   ./org

echo "Grand Canonical Lanczos begins !"
cp -r org GCLanczos
cd ./GCLanczos
 perl -w CalMod.pl 0 3  
 ./HPhi -e  GCEDnamelist.def
cd ..
echo "Grand Canonical Lanczos ends !"

echo "Grand Canonical TPQ begins !"
cp -r org GCTPQ
cd ./GCTPQ
 perl -w CalMod.pl 1 3  
 ./HPhi -e  GCEDnamelist.def
 cd ./output/
   cp ../Ave.pl .
   perl -w Ave.pl
 cd ..
cd ..
echo "Grand Canonical TPQ ends !"

echo "Grand Canonical FullDiag begins !"
cp -r org GCFullDiag
cd ./GCFullDiag
 perl -w CalMod.pl 2 3  
 ./HPhi -e  GCEDnamelist.def
 cd ./output/
   cp ../FT.pl .
   perl -w FT.pl
 cd ..
 #cat output 
cd ..
echo "Grand Canonical FullDiag ends !"

echo ""
echo "Energy"
cat GCLanczos/output/zvo_CG_energy.dat > tmp_Lanczos.dat
tail -n 1 GCTPQ/output/Ave_SS.dat      > tmp_TPQ.dat
head -n 1 GCFullDiag/output/Eigenvalue.dat > tmp_FullDiag.dat
perl -w Energy.pl

echo ""
perl -w CheckHop.pl  GCLanczos/output/zvo_CG_cisajs.dat > GCLanczos_Hop.dat
perl -w CheckCor.pl  GCLanczos/output/zvo_CG_cisajscktalt.dat > GCLanczos_Cor.dat

perl -w CheckHop.pl  GCTPQ/output/zvo_cisajs_set0step1990.dat  > GCTPQ_Hop.dat
perl -w CheckCor.pl  GCTPQ/output/zvo_cisajscktalt_set0step1990.dat  > GCTPQ_Cor.dat

perl -w CheckHop.pl  GCFullDiag/output/zvo_cisajs_eigen0.dat   > GCFullDaig_Hop.dat
perl -w CheckCor.pl  GCFullDiag/output/zvo_cisajscktalt_eigen0.dat   > GCFullDaig_Cor.dat
echo ""

paste GCLanczos_Hop.dat GCTPQ_Hop.dat GCFullDaig_Hop.dat
echo ""
paste GCLanczos_Cor.dat GCTPQ_Cor.dat GCFullDaig_Cor.dat
echo ""
