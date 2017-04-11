#!/bin/bash
if [ -z ${1} ] || [ ${1} = "help" ]; then
    echo ""
    echo "Usage:"
    echo "./HPhiconfig.sh system_name"
    echo " system_name should be chosen from below:"
    echo "         sekirei : ISSP system-B (Intel + SGIMPT)"
    echo "         fujitsu : ISSP system-C (FX10)"
    echo "              sr : SR16000"
    echo "           intel : Intel compiler + Linux PC"
    echo "   intel-openmpi : Intel compiler + OpenMPI"
    echo "     intel-mpich : Intel compiler + MPICH2"
    echo "  intel-intelmpi : Intel compiler + IntelMPI"
    echo "             gcc : GCC"
    echo "     gcc-openmpi : GCC + OpenMPI"
    echo "       gcc-mpich : GCC + MPICH2"
    echo "         gcc-mac : GCC + Mac"
    echo ""
    echo "Then src/make.sys is generated."
    echo "  Variables in src/make.sys"
    echo "              CC : C complier"
    echo "             F90 : fortran compiler"
    echo "          CFLAGS : C compiler options"
    echo "          FFLAGS : fortran compiler options"
    echo "            LIBS : Linker option"
    echo ""
else
    if [ ${1} = "sekirei" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
F90 = mpif90
CFLAGS = -fopenmp -O3 -g -traceback -xHost -ipo -mcmodel=large -shared-intel -D MPI -D HAVE_SSE2
FFLAGS = -fopenmp -O3 -g -traceback -xHost -ipo -mcmodel=large -shared-intel -D MPI -fpp
LIBS = -mkl -lifcore
EOF
    elif [ ${1} = "intel-mpich" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
F90 = mpif90
CFLAGS = -fopenmp -O3 -g -traceback -xHost -D MPI -D HAVE_SSE2
FFLAGS = -fopenmp -O3 -g -traceback -xHost -D MPI -fpp
LIBS = -mkl -lifcore -lmpifort
EOF
    elif [ ${1} = "intel-intelmpi" ]; then
        cat > src/make.sys <<EOF
CC = mpiicc
F90 = mpiifort
CFLAGS = -fopenmp -O3 -g -traceback -xHost -D MPI -D HAVE_SSE2
FFLAGS = -fopenmp -O3 -g -traceback -xHost -D MPI -fpp
LIBS = -mkl -lifcore -lmpifort
EOF
    elif [ ${1} = "intel-openmpi" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
F90 = mpif90
CFLAGS = -fopenmp -O3 -g -traceback -xHost -D MPI -D HAVE_SSE2
FFLAGS = -fopenmp -O3 -g -traceback -xHost -D MPI -fpp
LIBS = -mkl -lifcore -lmpi_f90 -lmpi_f77
EOF
    elif [ ${1} = "intel" ]; then
        cat > src/make.sys <<EOF
CC = icc
F90 = ifort
CFLAGS = -fopenmp -O3 -g -traceback -xHost -D HAVE_SSE2
FFLAGS = -fopenmp -O3 -g -traceback -xHost -fpp
LIBS = -mkl -lifcore
EOF
    elif [ ${1} = "gcc-openmpi" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
F90 = mpif90
CFLAGS = -fopenmp -O3 -g -D MPI -D HAVE_SSE2
FFLAGS = -fopenmp -O3 -g -D MPI -cpp
LIBS = -fopenmp -lm -lgfortran -llapack -lblas -lmpi_f90 -lmpi_f77
EOF
    elif [ ${1} = "gcc-mpich" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
F90 = mpif90
CFLAGS = -fopenmp -O3 -g -D MPI -D HAVE_SSE2
FFLAGS = -fopenmp -O3 -g -D MPI -cpp
LIBS = -fopenmp -lm -lgfortran -llapack -lblas -lmpifort
EOF
    elif [ ${1} = "gcc" ]; then
        cat > src/make.sys <<EOF
CC = gcc
F90 = gfortran
CFLAGS = -fopenmp -O3 -g -D HAVE_SSE2
FFLAGS = -fopenmp -O3 -g -cpp
LIBS = -fopenmp -lm -lgfortran -llapack -lblas
EOF
    elif [ ${1} = "gcc-mac" ]; then
        cat > src/make.sys <<EOF
CC = gcc
F90 = gfortran
CFLAGS = -fopenmp -O3 -g -D_OSX -D HAVE_SSE2
FFLAGS = -fopenmp -O3 -g -cpp -D NO_ZDOTC
LIBS = -fopenmp -lm -framework Accelerate -lgfortran
EOF
    elif [ ${1} = "fujitsu" ]; then
        cat > src/make.sys <<EOF
CC = mpifccpx
F90 = mpifrtpx
CFLAGS = -Kfast,openmp,SPARC64IXfx,parallel -g -D MPI -Kmemalias,alias_const -D HAVE_SSE2
FFLAGS = -Kfast,openmp,SPARC64IXfx,parallel -g -D MPI -Cpp -D FUJITSU
LIBS = -Kfast,openmp,parallel -SSL2BLAMP -lmpi_f90 -lmpi_f77
EOF
    elif [ ${1} = "sr" ]; then
        cat > src/make.sys <<EOF
CC = mpcc_r
F90 = mpxlf2003_r
CFLAGS = -O3 -qsmp=omp -q64 -D SR -D MPI 
FFLAGS = -O3 -qsmp=omp -q64 -qsuffix=cpp=f90 -WF,-DMPI -WF,-DSR
LIBS = -qsmp=omp -L /usr/lib -lm -lessl -lxlf90
AROPT = -X64
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
	cd doc/jp/;make -f makefile_doc_jp;
	cd doc/en/;make -f makefile_doc_en;
	cd doc/fourier/ja; make html latexpdfja
	cd doc/fourier/en; make html latexpdfja

clean:
	cd src; make -f makefile_src clean
	cd tool; make -f makefile_tool clean

veryclean:
	make clean
	cd doc/jp; make -f makefile_doc_jp clean
	cd doc/en; make -f makefile_doc_en clean
	cd doc/fourier/ja; make clean
	cd doc/fourier/en; make clean
	rm -f src/make.sys makefile
EOF
fi
