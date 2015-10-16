include_directories(include)
add_definitions(-DDSFMT_MEXP=19937)

set(SOURCES CG_EigenVector.c CalcByFullDiag.c CalcByLanczos.c CalcByTPQ.c FileIO.c FirstMultiply.c HPhiMain.c HPhiTrans.c Lanczos_EigenValue.c Lanczos_EigenVector.c Multiply.c bisec.c bitcalc.c check.c dSFMT.c diagonalcalc.c expec_cisajs.c expec_cisajscktaltdc.c expec_energy.c expec_totalspin.c global.c lapack_diag.c log.c makeHam.c matrixlapack.c mltply.c output.c output_list.c phys.c readdef.c sgn.c sz.c vec12.c xsetmem.c ErrorMessage.c LogMessage.c ProgressMessage.c wrapperMPI.c)

set(SOURCES_STDFACE StdFace/ChainLattice.c StdFace/HoneycombLattice.c StdFace/SquareLattice.c StdFace/StdFace_main.c StdFace/StdFace_ModelUtil.c StdFace/TriangularLattice.c)

add_executable(HPhi ${SOURCES} ${SOURCES_STDFACE})
target_link_libraries(HPhi ${LAPACK_LIBRARIES} ${MPI_C_LIBRARIES})
