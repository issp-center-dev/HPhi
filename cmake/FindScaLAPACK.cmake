# FindScaLAPACK.cmake -- locate and validate a ScaLAPACK implementation.
#
# A caller may provide SCALAPACK_LIBRARIES explicitly (including compiler
# driver flags such as -mkl=cluster). Otherwise pkg-config modules commonly
# shipped by Linux distributions are tried. The result is accepted only when
# a link probe for descinit_ succeeds with MPI and LAPACK.
#
# Result variables:
#   ScaLAPACK_FOUND / SCALAPACK_FOUND
#   SCALAPACK_LIBRARIES
#   SCALAPACK_INCLUDE_DIRS

if(NOT SCALAPACK_LIBRARIES)
  find_package(PkgConfig QUIET)
  if(PKG_CONFIG_FOUND)
    foreach(_scalapack_module scalapack scalapack-openmpi scalapack-mpich)
      if(NOT PC_SCALAPACK_FOUND)
        pkg_check_modules(PC_SCALAPACK QUIET ${_scalapack_module})
      endif()
    endforeach()
    if(PC_SCALAPACK_FOUND)
      # PC_*_LDFLAGS is available in all CMake versions supported by HPhi;
      # PC_*_LINK_LIBRARIES was only added in CMake 3.12.
      set(SCALAPACK_LIBRARIES ${PC_SCALAPACK_LDFLAGS})
      set(SCALAPACK_INCLUDE_DIRS ${PC_SCALAPACK_INCLUDE_DIRS})
    endif()
  endif()
endif()

if(NOT SCALAPACK_LIBRARIES)
  # Some installations ship a broken or absent .pc file (for example, a
  # ScaLAPACK pkg-config file whose MPI dependency has a vendor-specific
  # module name). Fall back to the conventional library names; the descinit_
  # probe below still rejects an incompatible implementation.
  find_library(SCALAPACK_LIBRARY_AUTO
    NAMES scalapack scalapack-openmpi scalapack-mpich)
  if(SCALAPACK_LIBRARY_AUTO)
    set(SCALAPACK_LIBRARIES ${SCALAPACK_LIBRARY_AUTO})
  endif()
endif()

if(SCALAPACK_LIBRARIES)
  # Historical config/*.cmake files use a single space-separated flag
  # string. Normalize it once so configure probes and all targets receive
  # separate linker arguments.
  string(REGEX REPLACE "-L[ ]+" "-L" _scalapack_link_args "${SCALAPACK_LIBRARIES}")
  string(REGEX REPLACE "[ \t]+" ";" _scalapack_link_args "${_scalapack_link_args}")
  set(SCALAPACK_LIBRARIES ${_scalapack_link_args})

  include(CheckFunctionExists)
  set(_scalapack_saved_required_libraries ${CMAKE_REQUIRED_LIBRARIES})
  set(_scalapack_saved_required_flags "${CMAKE_REQUIRED_FLAGS}")
  set(CMAKE_REQUIRED_LIBRARIES
      ${SCALAPACK_LIBRARIES} ${LAPACK_LIBRARIES} ${MPI_C_LIBRARIES})
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${MPI_C_LINK_FLAGS}")
  unset(HPHI_SCALAPACK_LINKS CACHE)
  check_function_exists(descinit_ HPHI_SCALAPACK_LINKS)
  set(CMAKE_REQUIRED_LIBRARIES ${_scalapack_saved_required_libraries})
  set(CMAKE_REQUIRED_FLAGS "${_scalapack_saved_required_flags}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ScaLAPACK
  REQUIRED_VARS SCALAPACK_LIBRARIES HPHI_SCALAPACK_LINKS)
set(SCALAPACK_FOUND ${ScaLAPACK_FOUND})

mark_as_advanced(HPHI_SCALAPACK_LINKS SCALAPACK_LIBRARY_AUTO)
