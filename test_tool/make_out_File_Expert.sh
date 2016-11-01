#!/bin/sh
source ~/.bash_profile
#Lanczos
echo "Expert calc"
echo ${PWD##*/}
echo "Start: Lanczos"
~/program/qlmpack/build/src/HPhi -e namelist.def > Lanczos.out
rm ./output/*cis*
rm ./output/*Time*
rm -rf ./output_Lanczos
mv ./output ./output_Lancozs
cp ./calcmod.def ./calcmod_.def
#TPQ
echo "Start: TPQ"
sed -i '' s/'CalcType   0'/'CalcType   1'/g calcmod.def  
~/program/qlmpack/build/src/HPhi -e namelist.def > TPQ.out
rm ./output/*cis*
rm ./output/*Time*
rm ./output/*TPQ*
rm -rf ./output_TPQ
mv ./output ./output_TPQ

#FullDiag
#echo "Start: FullDiag"
#sed -i '' s/'CalcType   1'/'CalcType   2'/g calcmod.def  
#~/program/qlmpack/build/src/HPhi -e namelist.def > FullDiag.out
#rm ./output/zvo_cisajs*eigen0*.dat
#rm ./output/zvo_cisajs*eigen1*.dat
#rm ./output/zvo_cisajs*eigen2*.dat
#rm ./output/zvo_cisajs*eigen3*.dat
#rm ./output/zvo_cisajs*eigen4*.dat
#rm ./output/zvo_cisajs*eigen5*.dat
#rm ./output/zvo_cisajs*eigen6*.dat
#rm ./output/zvo_cisajs*eigen7*.dat
#rm ./output/zvo_cisajs*eigen8*.dat
#rm ./output/zvo_cisajs*eigen9*.dat
#rm ./output/*Time*
#rm ./output/*energy*
#rm -rf ./output_FullDiag
#mv ./output ./output_FullDiag
mv ./calcmod_.def ./calcmod.def
