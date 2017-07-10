# include guard
cmake_minimum_required(VERSION 2.8)
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of HPhi.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

add_subdirectory(komega)

include_directories(include)
add_definitions(-D_HPhi)
add_definitions(-DDSFMT_MEXP=19937)
add_definitions(-D_HPhi)

set(SOURCES PowerLanczos.c CG_EigenVector.c CalcByFullDiag.c CalcByLOBPCG.c CalcByLanczos.c CalcByTPQ.c FileIO.c FirstMultiply.c HPhiMain.c HPhiTrans.c Lanczos_EigenValue.c Lanczos_EigenVector.c CalcSpectrum.c CalcSpectrumByBiCG.c CalcSpectrumByLanczos.c CalcSpectrumByTPQ.c CalcSpectrumByFullDiag.c SingleEx.c SingleExHubbard.c PairEx.c PairExHubbard.c PairExSpin.c Multiply.c bisec.c bitcalc.c check.c CheckMPI.c dSFMT.c diagonalcalc.c expec_cisajs.c expec_cisajscktaltdc.c expec_totalspin.c global.c lapack_diag.c log.c makeHam.c matrixlapack.c mltply.c mltplySpin.c mltplyHubbard.c mltplyMPIHubbard.c mltplyMPISpin.c input.c output.c output_list.c phys.c readdef.c sz.c vec12.c xsetmem.c ErrorMessage.c LogMessage.c ProgressMessage.c wrapperMPI.c mltplyMPIBoost.c splash.c expec_energy_flct.c time.c mltplyHubbardCore.c mltplySpinCore.c mltplyMPIHubbardCore.c mltplyMPISpinCore.c )
set(SOURCES PowerLanczos.c CG_EigenVector.c CalcByFullDiag.c CalcByLOBPCG.c CalcByLanczos.c CalcByTPQ.c FileIO.c FirstMultiply.c HPhiMain.c HPhiTrans.c Lanczos_EigenValue.c Lanczos_EigenVector.c CalcSpectrum.c CalcSpectrumByBiCG.c CalcSpectrumByLanczos.c CalcSpectrumByTPQ.c CalcSpectrumByFullDiag.c SingleEx.c SingleExHubbard.c PairEx.c PairExHubbard.c PairExSpin.c Multiply.c bisec.c bitcalc.c check.c CheckMPI.c dSFMT.c diagonalcalc.c expec_cisajs.c expec_cisajscktaltdc.c expec_totalspin.c global.c lapack_diag.c log.c makeHam.c matrixlapack.c mltply.c mltplySpin.c mltplyHubbard.c mltplyMPIHubbard.c mltplyMPISpin.c input.c output.c output_list.c phys.c readdef.c sz.c vec12.c xsetmem.c ErrorMessage.c LogMessage.c ProgressMessage.c wrapperMPI.c mltplyMPIBoost.c splash.c expec_energy_flct.c time.c mltplyHubbardCore.c mltplySpinCore.c mltplyMPIHubbardCore.c mltplyMPISpinCore.c CalcByTEM.c)

set(SOURCES_STDFACE StdFace/ChainLattice.c StdFace/FCOrtho.c StdFace/HoneycombLattice.c StdFace/Kagome.c StdFace/Ladder.c StdFace/Orthorhombic.c StdFace/Pyrochlore.c StdFace/SquareLattice.c StdFace/StdFace_main.c StdFace/StdFace_ModelUtil.c StdFace/TriangularLattice.c StdFace/Wannier90.c)

add_executable(HPhi ${SOURCES} ${SOURCES_STDFACE})
target_link_libraries(HPhi komega ${LAPACK_LIBRARIES} m)
if(MPI_FOUND)
  target_link_libraries(HPhi ${MPI_C_LIBRARIES})
endif(MPI_FOUND)
install(TARGETS HPhi RUNTIME DESTINATION bin)
