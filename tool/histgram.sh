#!/bin/bash

if [ -z ${1} ];then
    echo
    echo "Usage:"
    echo "$ histgram.sh {delta_E}"
    echo
    exit
fi

awk '
NR==1{
    dE='${1}'
    E0=$2
    N=0
}
E0<=$2&&$2<E0+dE{N=N+1}
E0+dE<=$2{
    while(E0+dE<=$2){
        print E0+0.5*dE, N
        N=0
        E0=E0+dE
    }
    N = N+1
}
END{print E0+0.5*dE, N}' output/Eigenvalue.dat

