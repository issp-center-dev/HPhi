#!/bin/bash
if [ -z ${1} ] || [ ${1} = "help" ]; then
    echo ""
    echo "Usage:"
    echo "./configure.sh system_name"
    echo "system_name should be choosen from below:"
    echo "  sekirei : ISSP system-B"
    echo "     make : ISSP system-C"
    echo "    intel : Intel compiler + Linux PC"
    echo "      gcc : GCC + Linux"
    echo "  gcc-mac : GCC + Mac"
    echo "   manual : Manual configuration. See below."
    echo ""
    echo "In manual configuretion, it is used as "
    echo "./configure.sh CC=icc LAPACK_FLAGS=\"-Dlapack -mkl=parallel\" \\"
    echo "   FLAGS=\"-qopenmp  -O3 -xCORE-AVX2 -mcmodel=large -shared-intel\""
    echo " where"
    echo "  CC : Compilation command for C"
    echo "  LAPACK_FLAGS : Compile option for LAPACK"
    echo "  FLAGS : Other Compilation options"
    echo ""
else
    if [ ${1} = "sekirei" ]; then
        cat > src/make.sys <<EOF
CC = icc
LAPACK_FLAGS = -Dlapack -mkl=parallel 
FLAGS = -qopenmp  -O3 -xCORE-AVX2 -mcmodel=large -shared-intel
MTFLAGS = -DDSFMT_MEXP=19937 \$(FLAGS)
INCLUDE_DIR=./include
EOF
    elif [ ${1} = "maki" ]; then
        cat > src/make.sys <<EOF
CC = fccpx
LAPACK_FLAGS = -Dlapack -SSL2BLAMP
FLAGS = -Kfast,openmp,SPARC64IXfx,parallel -Kmemalias,alias_const
MTFLAGS = -DDSFMT_MEXP=19937 \$(FLAGS)
INCLUDE_DIR=./include
EOF
    elif [ ${1} = "intel" ]; then
        cat > src/make.sys <<EOF
CC = icc
LAPACK_FLAGS = -Dlapack -mkl=parallel 
FLAGS = -openmp -openmp-report2 -O3 -DHAVE_SSE2
MTFLAGS = -DDSFMT_MEXP=19937 \$(FLAGS)
INCLUDE_DIR=./include
EOF
    elif [ ${1} = "gcc-mac" ]; then
        cat > src/make.sys <<EOF
CC = gcc
LAPACK_FLAGS = -framework Accelerate 
FLAGS = -fopenmp 
MTFLAGS = -DDSFMT_MEXP=19937 \$(FLAGS)
INCLUDE_DIR=./include
EOF
    elif [ ${1} = "gcc" ]; then
        cat > src/make.sys <<EOF
CC = gcc
LAPACK_FLAGS = -Dlapack -llapack -lblas
FLAGS = -fopenmp  -lm
MTFLAGS = -DDSFMT_MEXP=19937 \$(FLAGS)
INCLUDE_DIR=./include
EOF
    else 
        cat > src/make.sys <<EOF
CC = ${CC}
LAPACK_FLAGS = ${LAPACK_FLAGS}
FLAGS = ${FLAGS}
MTFLAGS = -DDSFMT_MEXP=19937 \$(FLAGS)
INCLUDE_DIR=./include
EOF
    fi

    echo "cat src/make.sys"
    cat src/make.sys

    echo
    echo "configure DONE"
    echo

    cat > makefile <<EOF
HPhi:
	cd src;make -f makefile_src

doc:
	cd doc/jp/;make -f makefile_doc_jp
	cd doc/en/;make -f makefile_doc_en

clean:
	cd src; make -f makefile_src clean
	cd doc/jp; make -f makefile_doc_jp clean
	cd doc/en; make -f makefile_doc_en clean

EOF
fi
