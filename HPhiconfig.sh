#!/bin/bash
if [ -z ${1} ] || [ ${1} = "help" ]; then
    echo ""
    echo "Usage:"
    echo "./HPhiconfig.sh system_name"
    echo " system_name should be chosen from below:"
    echo "     sekirei : ISSP system-B"
    echo "        maki : ISSP system-C"
    echo "          sr : SR16000"
    echo "       intel : Intel compiler + Linux PC"
    echo " mpicc-intel : Intel compiler + Linux PC + mpicc"
    echo "         gcc : GCC + Linux"
    echo "   mpicc-gcc : GCC + Linux + mpicc"
    echo "     gcc-mac : GCC + Mac"
    echo "      manual : Manual configuration. See below."
    echo ""
    echo "To configure manually HPhi, please type, for example,  "
    echo "./HPhiconfig.sh CC=icc LAPACK_FLAGS=\"-Dlapack -mkl=parallel\" \\"
    echo "   FLAGS=\"-qopenmp  -O3 -xCORE-AVX2 -mcmodel=large -shared-intel\""
    echo " where"
    echo "            CC : Compilation command for C"
    echo "  LAPACK_FLAGS : Compile option for LAPACK"
    echo "         FLAGS : Other Compilation options"
    echo "           F90 : Complilation command for fortran"
    echo "        FFLAGS : Compilation options for fortran"
    echo ""
else
    if [ ${1} = "sekirei" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
LAPACK_FLAGS = -Dlapack -mkl=parallel 
FLAGS = -qopenmp -O3 -ipo -xHOST -mcmodel=large -shared-intel -D MPI -g -traceback -lifcore
MTFLAGS = -DDSFMT_MEXP=19937 \$(FLAGS)
INCLUDE_DIR=./include
CP = cp -f -v
AR = ar rv
F90 = mpif90
FFLAGS = -fopenmp -fpp -g -traceback -mkl -xHost -O3 -g -traceback -D MPI
EOF
    elif [ ${1} = "maki" ]; then
        cat > src/make.sys <<EOF
CC = mpifccpx
LAPACK_FLAGS = -Dlapack -SSL2BLAMP
FLAGS = -Kfast,openmp,SPARC64IXfx,parallel -Kmemalias,alias_const -D MPI -g -lmpi_f90 -lmpi_f77
MTFLAGS = -DDSFMT_MEXP=19937 \$(FLAGS)
INCLUDE_DIR=./include
CP = cp -f -v
AR = ar rv
F90 = mpifrtpx
FFLAGS = -g -D MPI -Cpp -Kfast,openmp,SPARC64IXfx,parallel
EOF
    elif [ ${1} = "sr" ]; then
        cat > src/make.sys <<EOF
CC = mpcc_r
LAPACK_FLAGS = -L /usr/lib -lessl
FLAGS = -D MPI -D SR -q64 -O3 -lm -qsmp=omp -lxlf90
MTFLAGS = -DDSFMT_MEXP=19937 \$(FLAGS)
INCLUDE_DIR=./include
CP = cp -f
AR = ar -X64 rv 
F90 = mpxlf2003_r
FFLAGS = -WF,-DMPI -qsmp=omp -q64 -O3 -qsuffix=cpp=f90
EOF
    elif [ ${1} = "intel" ]; then
        cat > src/make.sys <<EOF
CC = icc
LAPACK_FLAGS = -Dlapack -mkl=parallel 
FLAGS = -fopenmp -O3 -DHAVE_SSE2 -g -traceback -xHOST -lifcore
MTFLAGS = -DDSFMT_MEXP=19937 \$(FLAGS)
INCLUDE_DIR=./include
CP = cp -f -v
AR = ar rv
F90 = ifort
FFLAGS = -fpp -mkl -xHost -fopenmp -O3 -g -traceback
EOF
    elif [ ${1} = "mpicc-intel" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
LAPACK_FLAGS = -Dlapack -mkl=parallel 
FLAGS = -fopenmp -O3 -DHAVE_SSE2 -D MPI -g -traceback -xHOST -lifcore -lmpifort
MTFLAGS = -DDSFMT_MEXP=19937 \$(FLAGS)
INCLUDE_DIR=./include
CP = cp -f -v
AR = ar rv
F90 = mpif90
FFLAGS = -fopenmp -fpp -g -traceback -mkl -xHost -O3 -g -traceback -D MPI
EOF
    elif [ ${1} = "mpicc-gcc" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
LAPACK_FLAGS = -Dlapack -llapack -lblas 
FLAGS = -fopenmp -O3 -DHAVE_SSE2 -D MPI -g -lm -lgfortran -lmpi_f90 -lmpi_f77
MTFLAGS = -DDSFMT_MEXP=19937 \$(FLAGS)
INCLUDE_DIR=./include
CP = cp -f -v
AR = ar rv
F90 = mpif90
FFLAGS = -fopenmp -cpp -g -lblas -llapack -D MPI
EOF
    elif [ ${1} = "gcc-mac" ]; then
        cat > src/make.sys <<EOF
CC = gcc
LAPACK_FLAGS = -Dlapack -framework Accelerate 
FLAGS = -fopenmp -lm -O3  -D_OSX -lgfortran
MTFLAGS = -DDSFMT_MEXP=19937 \$(FLAGS)
INCLUDE_DIR=./include
CP = cp -f -v
AR = ar rv
F90 = gfortran
FFLAGS = -fopenmp -cpp -g -lblas -llapack
EOF
    elif [ ${1} = "gcc" ]; then
        cat > src/make.sys <<EOF
CC = gcc
LAPACK_FLAGS = -Dlapack -llapack -lblas
FLAGS = -fopenmp -lm -O3 -lgfortran
MTFLAGS = -DDSFMT_MEXP=19937 \$(FLAGS)
INCLUDE_DIR=./include
CP = cp -f -v
AR = ar rv
F90 = gfortran
FFLAGS = -fopenmp -cpp -g -lblas -llapack
EOF
    elif [ ${1} == "manual" ]; then
echo " C compiler ?"
read CC
echo " LAPACK option ?"
read LAPACK_FLAGS
echo " Other compilation flags ?"
read FLAGS
        cat > src/make.sys <<EOF
CC = ${CC}
LAPACK_FLAGS = ${LAPACK_FLAGS}
FLAGS = ${FLAGS}
MTFLAGS = -DDSFMT_MEXP=19937 \$(FLAGS)
INCLUDE_DIR=./include
CP = cp -f -v
AR = ar rv
F90 = ${F90}
FFLAGS = ${FFLAGS}
EOF
    else
        echo ""
        echo "Unsupported system. Please type"
        echo "./HPhiconfig.sh help"
        echo ""
        exit
    fi

    echo "cat src/make.sys"
    cat src/make.sys

    echo
    echo "HPhiconfig DONE"
    echo

    cat > makefile <<EOF
help:
	@echo ""
	@echo "Usage :"
	@echo "make <entry>"
	@echo ""
	@echo "<entry> is chosen from below"
	@echo "      HPhi : Build simulator HPhi in src/ and tool/"
	@echo " userguide : Generate userguid_jp.pdf & userguide_en.pdf in doc/"
	@echo "     clean : Remove all generated files excepting makefile and doc/"
	@echo " veryclean : Remove all generated files including makefile and doc/"
	@echo ""

HPhi:
	cd src;make -f makefile_src
	cd tool;make -f makefile_tool

userguide:
	cd doc/jp/;make -f makefile_doc_jp;mv userguide_jp.pdf ../
	cd doc/en/;make -f makefile_doc_en;mv userguide_en.pdf ../

clean:
	cd src; make -f makefile_src clean
	cd tool; make -f makefile_tool clean

veryclean:
	make clean
	cd doc/jp; make -f makefile_doc_jp clean
	cd doc/en; make -f makefile_doc_en clean
	rm -f doc/userguide_??.pdf
	rm -f src/make.sys makefile
EOF
fi
