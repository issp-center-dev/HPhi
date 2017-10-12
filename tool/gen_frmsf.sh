#!/bin/sh

if test -z ${4}
then
    echo "Usage:"
    echo "$ gen_frmsf.sh filename nk1 nk2 nk3"
    echo "In Standard mode, nk1=W, nk2=L, nk3=Height."
    return
fi

icor[1]=4
icor[2]=6
icor[3]=8
icor[4]=14

cor[1]="up"
cor[2]="down"
cor[3]="density"
cor[4]="s2"

for i in 1 2 3 4
do
    fname="${1}.${cor[i]}.frmsf"
    echo -n "Writing ${fname} ... "
    echo ${4} ${3} ${2} > ${fname}
    
    echo "1" >> ${fname}
    echo "1" >> ${fname}

    awk 'NR==4+'${2}'*'${3}'+1{print '${4}'*$1, '${4}'*$2, '${4}'*$3}' ${1} >> ${fname}
    awk 'NR==4+'${2}'+1       {print '${3}'*$1, '${3}'*$2, '${3}'*$3}' ${1} >> ${fname}
    awk 'NR==4+1+1            {print '${2}'*$1, '${2}'*$2, '${2}'*$3}' ${1} >> ${fname}

    awk 'NR>4{print $'${icor[i]}'}' ${1} >> ${fname}
    awk 'NR>4{print $'${icor[i]}'}' ${1} >> ${fname}

    #for i4 in `seq 1 ${4}`
    #do
    #    ii4=`expr \( ${i4} - 1 + ${4} / 2 \) % ${4} - ${4} / 2`
    #    k4=`echo ${ii4} / ${4} | bc -l`
    #    for i3 in `seq 1 ${3}`
    #    do
    #        ii3=`expr \( $i3 - 1 + $3 / 2 \) % $3 - $3 / 2`
    #        k3=`echo ${ii3} / ${3} | bc -l`
    #        for i2 in `seq 1 ${2}`
    #        do
    #            ii2=`expr \( $i2 - 1 + $2 / 2 \) % $2 - $2 / 2`
    #            k2=`echo ${ii2} / ${2} | bc -l`
    #            echo `echo ${k2}^2 + ${k3}^2 + ${k4}^2 | bc -l` >> ${fname}
    #        done
    #    done
    #done
    echo "done"
done
