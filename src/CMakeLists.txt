# include guard
cmake_minimum_required(VERSION 2.8)
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of HPhi.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

add_subdirectory(komega)

include_directories(include)
add_definitions(-D_HPhi)
add_definitions(-DDSFMT_MEXP=19937)

set(SOURCES FileIO.c HPhiMain.c HPhiTrans.c bitcalc.c check.c CheckMPI.c dSFMT.c diagonalcalc.c global.c log.c input.c output.c output_list.c readdef.c sz.c vec12.c xsetmem.c ErrorMessage.c LogMessage.c ProgressMessage.c wrapperMPI.c splash.c time.c eigenIO.c)

set(SOURCES_MLTPLY Multiply.c mltply.c mltplySpin.c mltplyHubbard.c mltplyMPIHubbard.c mltplyMPISpin.c mltplyMPIBoost.c mltplyHubbardCore.c mltplySpinCore.c mltplyMPIHubbardCore.c mltplyMPISpinCore.c)
set(SOURCES_LANCZOS CalcByLanczos.c CG_EigenVector.c Lanczos_EigenValue.c Lanczos_EigenVector.c)
set(SOURCES_CG CalcByLOBPCG.c )
set(SOURCES_TPQ CalcByTPQ.c FirstMultiply.c)
set(SOURCES_DIAG CalcByFullDiag.c lapack_diag.c makeHam.c matrixlapack.c matrixlapack_magma.c matrixscalapack.c )
set(SOURCES_SPECTRUM CalcSpectrum.c CalcSpectrumByBiCG.c CalcSpectrumByLanczos.c CalcSpectrumByTPQ.c CalcSpectrumByFullDiag.c SingleEx.c SingleExHubbard.c PairEx.c PairExHubbard.c PairExSpin.c )
set(SOURCES_TEM CalcByTEM.c)
set(SOURCES_PHYS expec_cisajs.c expec_cisajscktaltdc.c expec_totalspin.c phys.c expec_energy_flct.c)

set(SOURCES_COMMON common/setmemory.c)

set(SOURCES_STDFACE StdFace/src/ChainLattice.c StdFace/src/FCOrtho.c StdFace/src/HoneycombLattice.c StdFace/src/Kagome.c StdFace/src/Ladder.c StdFace/src/Orthorhombic.c StdFace/src/Pyrochlore.c StdFace/src/SquareLattice.c StdFace/src/StdFace_main.c StdFace/src/StdFace_ModelUtil.c StdFace/src/TriangularLattice.c StdFace/src/Wannier90.c)

add_executable(HPhi ${SOURCES} ${SOURCES_LANCZOS} ${SOURCES_CG} ${SOURCES_TPQ} ${SOURCES_SPECTRUM} ${SOURCES_TEM} ${SOURCES_MLTPLY} ${SOURCES_DIAG} ${SOURCES_PHYS} ${SOURCES_STDFACE} ${SOURCES_COMMON})
target_link_libraries(HPhi komega ${LAPACK_LIBRARIES} m)
if(MPI_FOUND)
  target_link_libraries(HPhi ${MPI_C_LIBRARIES})
endif(MPI_FOUND)

if(MAGMA_FOUND)
  target_link_libraries(HPhi ${MAGMA_C_LIBRARIES})
  add_definitions(-D_MAGMA)
endif(MAGMA_FOUND)
if(CUDA_FOUND)
  target_link_libraries(HPhi ${CUDA_LIBRARIES})
endif(CUDA_FOUND)

if(USE_SCALAPACK)
  string(REGEX REPLACE "-L[ ]+" "-L" sc_libs "${SCALAPACK_LIBRARIES}")
  string(REGEX REPLACE "[ ]+" ";" sc_libs "${sc_libs}")
  foreach(sc_lib IN LISTS sc_libs)
    target_link_libraries(HPhi ${sc_lib})
  endforeach(sc_lib)
  message(STATUS "SCALAPACK_LIBRARIES: ${SCALAPACK_LIBRARIES}")
endif(USE_SCALAPACK)


install(TARGETS HPhi RUNTIME DESTINATION bin)
