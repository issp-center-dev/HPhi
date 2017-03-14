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
    echo ""
    echo "To configure manually HPhi, please type, for example,  "
    echo "./HPhiconfig.sh CC=icc LIBS=\"-Dlapack -mkl=parallel\" \\"
    echo "   FLAGS=\"-qopenmp  -O3 -xCORE-AVX2 -mcmodel=large -shared-intel\""
    echo " where"
    echo "            CC : Compilation command for C"
    echo "           F90 : Complilation command for fortran"
    echo "        CFLAGS : Other Compilation options"
    echo "        FFLAGS : Compilation options for fortran"
    echo "          LIBS : Compile option for LAPACK"
    echo ""
else
    if [ ${1} = "sekirei" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
F90 = mpif90
CFLAGS = -fopenmp -O3 -ipo -xHost -mcmodel=large -shared-intel -g -traceback -D MPI -mkl
FFLAGS = -fopenmp -O3 -ipo -xHost -mcmodel=large -shared-intel -g -traceback -D MPI -mkl -fpp
LIBS = -fopenmp -mkl -lifcore
AR = ar rv
EOF
    elif [ ${1} = "maki" ]; then
        cat > src/make.sys <<EOF
CC = mpifccpx
F90 = mpifrtpx
CFLAGS = -Kfast,openmp,SPARC64IXfx,parallel -g -D MPI -Kmemalias,alias_const
FFLAGS = -Kfast,openmp,SPARC64IXfx,parallel -g -D MPI -Kmemalias,alias_const -Cpp
LIBS = -SSL2BLAMP -lmpi_f90 -lmpi_f77
AR = ar rv
EOF
    elif [ ${1} = "sr" ]; then
        cat > src/make.sys <<EOF
CC = mpcc_rls
F90 = mpxlf2003_r
CFLAGS = -O3 -qsmp=omp -q64 -D SR -D MPI 
FFLAGS = -O3 -qsmp=omp -q64 -qsuffix=cpp=f90 -WF,-DMPI
LIBS = -L /usr/lib -lm -lessl -lxlf90
AR = ar -X64 rv 
EOF
    elif [ ${1} = "intel" ]; then
        cat > src/make.sys <<EOF
CC = icc
F90 = ifort
CFLAGS = -fopenmp -O3 -g -traceback -mkl -xHost -D HAVE_SSE2
FFLAGS = -fopenmp -O3 -g -traceback -mkl -xHost -fpp
LIBS = -fopenmp -mkl -lifcore
AR = ar rv
EOF
    elif [ ${1} = "mpicc-intel" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
F90 = mpif90
CFLAGS = -fopenmp -O3 -g -traceback -D MPI -xHost -mkl -D HAVE_SSE2
FFLAGS = -fopenmp -O3 -g -traceback -D MPI -xHost -mkl -fpp
LIBS = -fopenmp -mkl -lifcore -lmpifort
AR = ar rv
EOF
    elif [ ${1} = "mpicc-gcc" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
F90 = mpif90
CFLAGS = -fopenmp -O3 -D MPI -g -D HAVE_SSE2
FFLAGS = -fopenmp -O3 -D MPI -g -cpp
LIBS = -fopenmp -lm -lgfortran -lmpi_f90 -lmpi_f77 -llapack -lblas 
AR = ar rv
EOF
    elif [ ${1} = "gcc-mac" ]; then
        cat > src/make.sys <<EOF
CC = gcc
F90 = gfortran
CFLAGS = -fopenmp -O3 -g -D_OSX
FFLAGS = -fopenmp -O3 -g -cpp
LIBS = -lm -framework Accelerate -lgfortran
AR = ar rv
EOF
    elif [ ${1} = "gcc" ]; then
        cat > src/make.sys <<EOF
CC = gcc
F90 = gfortran
CFLAGS = -fopenmp -g -O3
FFLAGS = -fopenmp -g -O3 -cpp
LIBS = -fopenmp -lm -llapack -lblas -lgfortran
AR = ar rv
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
