# FindELPA.cmake — locate the (non-threaded) ELPA library.
#
# Search order:
#   1. Explicit cache variables ELPA_INCLUDE_DIR / ELPA_LIBRARY
#   2. pkg-config: module "elpa" or versioned "elpa-<version>"
#   3. ELPA_ROOT hint with versioned include dirs (include/elpa-*/elpa/elpa.h)
#
# The threaded variant (elpa_openmp) is NOT supported: HPhi initializes MPI
# with MPI_Init (MPI_THREAD_SINGLE), which is insufficient for it.
#
# Result variables: ELPA_FOUND, ELPA_INCLUDE_DIRS, ELPA_LIBRARIES

set(ELPA_FOUND FALSE)

if(ELPA_INCLUDE_DIR AND ELPA_LIBRARY)
  if(ELPA_LIBRARY MATCHES "elpa_openmp")
    message(FATAL_ERROR
      "ELPA_LIBRARY points to the threaded ELPA variant (${ELPA_LIBRARY}). "
      "HPhi requires the non-threaded ELPA (MPI_THREAD_SINGLE); "
      "build ELPA without --enable-openmp.")
  endif()
  set(ELPA_INCLUDE_DIRS ${ELPA_INCLUDE_DIR})
  set(ELPA_LIBRARIES ${ELPA_LIBRARY})
  set(ELPA_FOUND TRUE)
endif()

if(NOT ELPA_FOUND)
  find_package(PkgConfig QUIET)
  if(PKG_CONFIG_FOUND)
    pkg_check_modules(PC_ELPA QUIET elpa)
    if(NOT PC_ELPA_FOUND)
      # Versioned .pc files (elpa-YYYY.MM.SSS.pc). Reject elpa_openmp-*.
      execute_process(
        COMMAND ${PKG_CONFIG_EXECUTABLE} --list-all
        OUTPUT_VARIABLE _elpa_pkg_list ERROR_QUIET)
      string(REGEX MATCH "elpa-[0-9][0-9.]*" _elpa_pkg_name "${_elpa_pkg_list}")
      if(_elpa_pkg_name)
        pkg_check_modules(PC_ELPA QUIET ${_elpa_pkg_name})
      endif()
      string(REGEX MATCH "elpa_openmp-[0-9][0-9.]*" _elpa_omp_name "${_elpa_pkg_list}")
      if(NOT PC_ELPA_FOUND AND _elpa_omp_name)
        message(FATAL_ERROR
          "Only the threaded ELPA variant (${_elpa_omp_name}) was found. "
          "HPhi requires the non-threaded ELPA (MPI_THREAD_SINGLE); "
          "build ELPA without --enable-openmp.")
      endif()
    endif()
    if(PC_ELPA_FOUND)
      set(ELPA_INCLUDE_DIRS ${PC_ELPA_INCLUDE_DIRS})
      set(ELPA_LIBRARIES ${PC_ELPA_LINK_LIBRARIES})
      set(ELPA_FOUND TRUE)
    endif()
  endif()
endif()

if(NOT ELPA_FOUND AND ELPA_ROOT)
  file(GLOB _elpa_inc_candidates "${ELPA_ROOT}/include/elpa-*" "${ELPA_ROOT}/include")
  find_path(ELPA_INCLUDE_DIR_AUTO elpa/elpa.h PATHS ${_elpa_inc_candidates} NO_DEFAULT_PATH)
  find_library(ELPA_LIBRARY_AUTO NAMES elpa PATHS "${ELPA_ROOT}/lib" "${ELPA_ROOT}/lib64" NO_DEFAULT_PATH)
  if(ELPA_INCLUDE_DIR_AUTO AND ELPA_LIBRARY_AUTO)
    set(ELPA_INCLUDE_DIRS ${ELPA_INCLUDE_DIR_AUTO})
    set(ELPA_LIBRARIES ${ELPA_LIBRARY_AUTO})
    set(ELPA_FOUND TRUE)
  endif()
endif()

if(ELPA_FOUND)
  message(STATUS "ELPA include: ${ELPA_INCLUDE_DIRS}")
  message(STATUS "ELPA library: ${ELPA_LIBRARIES}")
  # GPU API (ELPA >= 2023.11.001): detect elpa_setup_gpu
  include(CheckSymbolExists)
  set(CMAKE_REQUIRED_INCLUDES ${ELPA_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_LIBRARIES ${ELPA_LIBRARIES})
  check_symbol_exists(elpa_setup_gpu "elpa/elpa.h" ELPA_HAVE_SETUP_GPU)
  unset(CMAKE_REQUIRED_INCLUDES)
  unset(CMAKE_REQUIRED_LIBRARIES)
else()
  message(FATAL_ERROR
    "USE_ELPA=ON but ELPA was not found. Set ELPA_ROOT=<prefix>, or set "
    "ELPA_INCLUDE_DIR and ELPA_LIBRARY explicitly.")
endif()
