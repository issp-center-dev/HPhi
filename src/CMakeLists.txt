# include guard
cmake_minimum_required(VERSION 2.8)
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of HPhi.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

include_directories(include)
add_definitions(-DDSFMT_MEXP=19937)

set(SOURCES PowerLanczos.c CG_EigenVector.c CalcByFullDiag.c CalcByLanczos.c CalcByTPQ.c FileIO.c FirstMultiply.c HPhiMain.c HPhiTrans.c Lanczos_EigenValue.c Lanczos_EigenVector.c Multiply.c bisec.c bitcalc.c check.c CheckMPI.c dSFMT.c diagonalcalc.c expec_cisajs.c expec_cisajscktaltdc.c expec_energy.c expec_totalspin.c global.c lapack_diag.c log.c makeHam.c matrixlapack.c mltply.c mltplyMPI.c output.c output_list.c phys.c readdef.c sgn.c sz.c vec12.c xsetmem.c ErrorMessage.c LogMessage.c ProgressMessage.c wrapperMPI.c splash.c)

set(SOURCES_STDFACE StdFace/ChainLattice.c StdFace/HoneycombLattice.c StdFace/SquareLattice.c StdFace/StdFace_main.c StdFace/StdFace_ModelUtil.c StdFace/TriangularLattice.c StdFace/Ladder.c)

add_executable(HPhi ${SOURCES} ${SOURCES_STDFACE})
target_link_libraries(HPhi ${LAPACK_LIBRARIES})
if(MPI_FOUND)
  target_link_libraries(HPhi ${MPI_C_LIBRARIES})
endif(MPI_FOUND)
