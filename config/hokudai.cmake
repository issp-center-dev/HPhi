# for Intel oneAPI (IntelLLVM) Compiler
set(CMAKE_C_COMPILER "icx" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -DHAVE_SSE2" CACHE STRING "" FORCE)

set(CMAKE_Fortran_COMPILER "ifx" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -DHAVE_SSE2" CACHE STRING "" FORCE)

if(USE_SCALAPACK)
  if(SCALAPACK_LIBRARIES MATCHES "")
    # Use Intel oneAPI MKL with ScaLAPACK (cluster version)
    set(SCALAPACK_LIBRARIES "-qmkl=cluster")
    # set(SCALAPACK_LIBRARIES "-L$ENV{MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64")
  endif()

  message(STATUS "SCALAPACK_LIBRARIES is ${SCALAPACK_LIBRARIES}")
endif()

# for Intel MKL (BLAS/LAPACK detection)
set(BLA_VENDOR "Intel10_64lp" CACHE STRING "" FORCE)
