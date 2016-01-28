#!/bin/sh
source ~/.bash_profile
cd ../triangular
#Lanczos
~/program/qlmpack/build/src/HPhi -s ./StdFace.def > Lanczos.out
rm ./output/*cis*
rm ./output/*Time*
rm -rf ./output_Lanczos
mv ./output ./output_Lancozs

#TPQ
sed s/'method = "Lanczos"'/'method = "TPQ"'/g StdFace.def  > StdFace_TPQ.def
~/program/qlmpack/build/src/HPhi -s ./StdFace_TPQ.def > TPQ.out
rm ./output/*cis*
rm ./output/*Time*
rm ./output/*energy*
rm ./output/*TPQ*
rm -rf ./output_TPQ
mv ./output ./output_TPQ
rm ./StdFace_TPQ.def

#FullDiag
sed s/'method = "Lanczos"'/'method = "FullDiag"'/g StdFace.def  > StdFace_FullDiag.def
~/program/qlmpack/build/src/HPhi -s ./StdFace_FullDiag.def > FullDiag.out
rm ./output/zvo_cisajs*eigen0*.dat
rm ./output/zvo_cisajs*eigen1*.dat
rm ./output/zvo_cisajs*eigen2*.dat
rm ./output/zvo_cisajs*eigen3*.dat
rm ./output/zvo_cisajs*eigen4*.dat
rm ./output/zvo_cisajs*eigen5*.dat
rm ./output/zvo_cisajs*eigen6*.dat
rm ./output/zvo_cisajs*eigen7*.dat
rm ./output/zvo_cisajs*eigen8*.dat
rm ./output/zvo_cisajs*eigen9*.dat
rm ./output/*Time*
rm ./output/*energy*
rm -rf ./output_FullDiag
mv ./output ./output_FullDiag
rm ./StdFace_FullDiag.def

cp ./StdFace.def ./StdFace.def_
rm ./*.def
mv ./StdFace.def_ ./StdFace.def

