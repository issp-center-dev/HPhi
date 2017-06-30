#!/bin/bash

awk '
NR==1{
    dE='${2}'
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
END{print E0+0.5*dE, N}' ${1}

