# include guard
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should not be executed directly in subdirectory.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

add_subdirectory(komega)

include_directories(include)
include_directories(${STDFACE_DIR})
add_definitions(-D_HPhi)
add_definitions(-DDSFMT_MEXP=19937)

# Create executable
add_executable(HPhi)

# Core sources
target_sources(HPhi PRIVATE
  FileIO.c HPhiMain.c HPhiTrans.c bitcalc.c check.c CheckMPI.c
  dSFMT.c diagonalcalc.c global.c log.c input.c output.c output_list.c
  readdef.c sz.c vec12.c xsetmem.c ErrorMessage.c LogMessage.c
  ProgressMessage.c wrapperMPI.c splash.c time.c eigenIO.c struct.c
  nbody_interall.c anomalous_pair.c green_output.c symmetry_basis.c symmetry_basis_io.c
  symmetry_diagonal.c symmetry_matvec_plan.c symmetry_state_enumerator.c
  symmetry_vector_halo.c
)

# Multiply operations
target_sources(HPhi PRIVATE
  Multiply.c mltply.c mltplySpin.c mltplyHubbard.c
  mltplyMPIHubbard.c mltplyMPISpin.c mltplyMPIBoost.c
  mltplyHubbardCore.c mltplySpinCore.c
  mltplyMPIHubbardCore.c mltplyMPISpinCore.c
  mltplySpinless.c mltplyMPISpinlessFermion.c mltplyMPIBatched.c
  mltplySpinSym.c
)

# Lanczos solver
target_sources(HPhi PRIVATE
  CalcByLanczos.c CG_EigenVector.c
  Lanczos_EigenValue.c Lanczos_EigenVector.c
)

# LOBCG solver
target_sources(HPhi PRIVATE
  CalcByLOBPCG.c
)

# TPQ calculations
target_sources(HPhi PRIVATE
  CalcByTPQ.c FirstMultiply.c CalcByCanonicalTPQ.c MakeIniVec.c
)

# Full diagonalization
target_sources(HPhi PRIVATE
  CalcByFullDiag.c lapack_diag.c makeHam.c
  matrixlapack.c matrixlapack_magma.c matrixscalapack.c
)

# Spectrum calculations
target_sources(HPhi PRIVATE
  CalcSpectrum.c CalcSpectrumByBiCG.c CalcSpectrumByLanczos.c
  CalcSpectrumByTPQ.c CalcSpectrumByFullDiag.c
  SingleEx.c SingleExHubbard.c
  PairEx.c PairExHubbard.c PairExSpin.c PairExSpinless.c
  ExcitationSectorShift.c
)

# Time evolution
target_sources(HPhi PRIVATE
  CalcByTEM.c
)

# Physical quantities
target_sources(HPhi PRIVATE
  expec_cisajs.c expec_cisajscktaltdc.c expec_totalspin.c
  nbody_correlation.c phys.c expec_energy_flct.c
)

# Common utilities
target_sources(HPhi PRIVATE
  common/setmemory.c
)

# StdFace library
target_sources(HPhi PRIVATE
  ${STDFACE_DIR}/src/ChainLattice.c
  ${STDFACE_DIR}/src/FCOrtho.c
  ${STDFACE_DIR}/src/HoneycombLattice.c
  ${STDFACE_DIR}/src/Kagome.c
  ${STDFACE_DIR}/src/Ladder.c
  ${STDFACE_DIR}/src/Orthorhombic.c
  ${STDFACE_DIR}/src/Pyrochlore.c
  ${STDFACE_DIR}/src/SquareLattice.c
  ${STDFACE_DIR}/src/StdFace_main.c
  ${STDFACE_DIR}/src/StdFace_ModelUtil.c
  ${STDFACE_DIR}/src/TriangularLattice.c
  ${STDFACE_DIR}/src/Wannier90.c
)
target_link_libraries(HPhi komega ${LAPACK_LIBRARIES} m)
if(TARGET OpenMP::OpenMP_C)
  target_link_libraries(HPhi OpenMP::OpenMP_C)
endif()
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
