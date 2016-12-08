#!/bin/bash

echo ""
if [ -e ${1} ];then
    echo "Input from " ${1}
else
    echo "ERROR ! " ${1} " is not exist."
    return -1
fi

LargeValue=`grep -i largevalue ${1} | awk '{print $2}'`
LanczosMax=`grep -i lanczos_max ${1} | awk '{print $2}'`
NumAve=`find ./ -name "SS_rand*.dat" | wc -l`

echo "LargeValue : " ${LargeValue}
echo "LanczosMax : " ${LanczosMax}
echo "NumAve : " ${NumAve}
LanczosMax=`expr ${LanczosMax} - 2`

FilesSS_rand=`find ./ -name "SS_rand*.dat" | sort`
FilesNorm_rand=`find ./ -name "Norm_rand*.dat" | sort`
echo "SS_rand Files :"
echo ${FilesSS_rand}
echo "Norm_rand Files :"
echo ${FilesNorm_rand}

echo "Read SS_rand"
ii=0
for SS_rand in ${FilesSS_rand}
do

    echo ${SS_rand} " "
    beta0=(`awk 'NR>1{printf "%.16lf ", $1}' ${SS_rand}`)
    energy0=(`awk 'NR>1{printf "%.16lf ", $2}' ${SS_rand}`)
    sqene0=(`awk 'NR>1{printf "%.16lf ", $3}' ${SS_rand}`)

    if [ ${ii} -ne 0 ]; then
        if [ `echo "${beta0[${LanczosMax}]} < ${beta[${LanczosMax}]}" | bc` -eq 1 ]; then
            echo "New beta"
            for iLan in `seq 0 ${LanczosMax}`
            do
                beta[${iLan}]=${beta0[${iLan}]}
            done
        fi
    else
        for iLan in `seq 0 ${LanczosMax}`
        do
            beta[${iLan}]=${beta0[${iLan}]}
        done
    fi

    for iLan in `seq 0 ${LanczosMax}`
    do
        energy[${ii}]=${energy0[${iLan}]}
        sqene[${ii}]=${sqene0[${iLan}]}
        let ++ii
    done
done
echo ""

echo "Read Norm_rand"
ii=0
for Norm_rand in ${FilesNorm_rand}
do

    echo ${Norm_rand} " "
    beta0=(`awk 'NR>1{printf "%.16lf ", $2}' ${Norm_rand}`)
    for iLan in `seq 0 ${LanczosMax}`
    do
        norm[${ii}]=${beta0[${iLan}]}
        let ++ii
    done
done
echo ""

for iAve in `seq 1 ${NumAve}`
do

    for iLan in `seq 0 ${LanczosMax}`
    do
    
        beta1=${beta[${iLan}]}
        coeff=1.0
        bb=0.0
        bHb=0.0
        bHHb=0.0
     
        for jLan in `seq ${LanczosMax} -1 0`
        do

            norm1=${norm[`expr $jLan + \( $iAve - 1 \) \* \( $LanczosMax + 1 \)`]}
            energy1=${energy[${jLan}]}
            sqene1=${sqene[${jLan}]}

            bb=$beta1 \* $beta1 \* $norm1 * &
        &  ($bb + $coeff * ((1 + $beta * ${Nsite} * ${LargeValue} / (2 * ${jLan} + 1)) &
        &                - $beta1 * $energy1 / dble(2*[${jLan}] + 1) &
        &                ) &
        &  )
        bHb = $beta1 * $beta1 * $norm1 * &
        &  ($bHb + $coeff * ((1 + $beta1 * ${Nsite} * ${LargeValue} / (2 * ${jLan} + 1)) * $energy1 &
        &                 - $beta1 \* $sqene1 / (2 * ${jLan} + 1) &
        &                 ) &
        &  )
        bHHb = $beta1 * $beta1 * $norm1 * &
        &  ($bHHb + $coeff * ((1 - $beta1 * ${Nsite} * ${LargeValue} / (2 * ${jLan} + 1)) * $sqene1 &
        &                  + ($beta1 * $Nsite * $Nsite * ${LargeValue} * ${LargeValue} / (2 * ${jLan} + 1) &
        &                    - (2 * ${jLan}) / $beta1 &
        &                    ) * $energy1 &
        &                  ) &
        &  )
        !
        if(bb > 1.0d+10) then
           bb = bb * (1.0d-10)
           bHb = bHb * (1.0d-10)
           bHHb = bHHb * (1.0d-10)
           coeff = coeff * (1.0d-10)
        endif
        !
     enddo
     !
     write(*,*)beta,bHb/bb,bHHb/bb
     !
  enddo

done