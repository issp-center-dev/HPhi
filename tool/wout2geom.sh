#!/bin/bash

grep -A 3 " Lattice Vectors " boron.wout | awk 'NR>1{print $2, $3, $4}'
declare -a b0=`grep " b_1 " ${1} | awk '{print "("$2, $3, $4")"}'`
declare -a b1=`grep " b_2 " ${1} | awk '{print "("$2, $3, $4")"}'`
declare -a b2=`grep " b_3 " ${1} | awk '{print "("$2, $3, $4")"}'`

nwan=`grep " |  Number of Wannier Functions               : " ${1} | awk '{print $7}'`

echo ${nwan}
tpi=6.283185307179586476925286766559

for iwan in `seq 1 ${nwan}`
do
    iwan1=`expr ${iwan} + 1`
    declare -a wan=`grep -A ${nwan} " Final State" ${1} | awk 'NR=='${iwan1}'{print "("$7, $8, $9")"}'|sed -e "s/,//g"`
    printf "%d %f %f %f\n" `expr ${iwan} - 1` \
         `echo "scale=10;(${wan[0]}*${b0[0]}+${wan[1]}*${b0[1]}+${wan[2]}*${b0[2]})/$tpi"|bc` \
         `echo "scale=10;(${wan[0]}*${b1[0]}+${wan[1]}*${b1[1]}+${wan[2]}*${b1[2]})/$tpi"|bc` \
         `echo "scale=10;(${wan[0]}*${b2[0]}+${wan[1]}*${b2[1]}+${wan[2]}*${b2[2]})/$tpi"|bc`
done

