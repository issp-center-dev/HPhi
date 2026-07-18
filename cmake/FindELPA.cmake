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
# Result variables: ELPA_FOUND, ELPA_INCLUDE_DIRS, ELPA_LIBRARIES,
#                   ELPA_COMPILE_OPTIONS

set(ELPA_FOUND FALSE)

if(ELPA_INCLUDE_DIR OR ELPA_LIBRARY)
  if(NOT ELPA_INCLUDE_DIR OR NOT ELPA_LIBRARY)
    message(FATAL_ERROR
      "ELPA_INCLUDE_DIR and ELPA_LIBRARY must be specified together.")
  endif()
  if(ELPA_LIBRARY MATCHES "elpa_openmp")
    message(FATAL_ERROR
      "ELPA_LIBRARY points to the threaded ELPA variant (${ELPA_LIBRARY}). "
      "HPhi requires the non-threaded ELPA (MPI_THREAD_SINGLE); "
      "build ELPA without --enable-openmp.")
  endif()
  if(NOT EXISTS "${ELPA_INCLUDE_DIR}/elpa/elpa.h")
    message(FATAL_ERROR
      "ELPA_INCLUDE_DIR does not contain elpa/elpa.h: ${ELPA_INCLUDE_DIR}")
  endif()
  if(NOT EXISTS "${ELPA_LIBRARY}")
    message(FATAL_ERROR "ELPA_LIBRARY does not exist: ${ELPA_LIBRARY}")
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
      # PC_*_LDFLAGS works with CMake 3.0; PC_*_LINK_LIBRARIES requires
      # CMake 3.12 and would break HPhi's declared compatibility range.
      set(ELPA_LIBRARIES ${PC_ELPA_LDFLAGS})
      set(ELPA_COMPILE_OPTIONS ${PC_ELPA_CFLAGS_OTHER})
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
  # Validate the base API before accepting the package. This catches stale
  # cache paths and incomplete transitive dependencies during configure,
  # rather than much later while linking HPhi.
  include(CheckSymbolExists)
  set(_elpa_saved_required_includes ${CMAKE_REQUIRED_INCLUDES})
  set(_elpa_saved_required_libraries ${CMAKE_REQUIRED_LIBRARIES})
  set(_elpa_saved_required_flags "${CMAKE_REQUIRED_FLAGS}")
  string(REPLACE ";" " " _elpa_compile_flags "${ELPA_COMPILE_OPTIONS}")
  set(CMAKE_REQUIRED_INCLUDES ${ELPA_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_LIBRARIES ${ELPA_LIBRARIES} ${SCALAPACK_LIBRARIES}
      ${LAPACK_LIBRARIES} ${MPI_C_LIBRARIES})
  set(CMAKE_REQUIRED_FLAGS
      "${CMAKE_REQUIRED_FLAGS} ${_elpa_compile_flags} ${MPI_C_LINK_FLAGS}")
  unset(ELPA_HAVE_INIT CACHE)
  check_symbol_exists(elpa_init "elpa/elpa.h" ELPA_HAVE_INIT)
  if(NOT ELPA_HAVE_INIT)
    message(FATAL_ERROR
      "ELPA was found, but a test program using elpa_init could not be linked. "
      "Check ELPA, ScaLAPACK, LAPACK, and MPI library compatibility.")
  endif()

  # GPU API (ELPA >= 2023.11.001): detect elpa_setup_gpu.
  # MPI and LAPACK libraries are appended so the check's test link does not
  # fail on unresolved dependency symbols when ELPA is a static library,
  # which would falsely report the GPU API as missing.
  unset(ELPA_HAVE_SETUP_GPU CACHE)
  check_symbol_exists(elpa_setup_gpu "elpa/elpa.h" ELPA_HAVE_SETUP_GPU)
  set(CMAKE_REQUIRED_INCLUDES ${_elpa_saved_required_includes})
  set(CMAKE_REQUIRED_LIBRARIES ${_elpa_saved_required_libraries})
  set(CMAKE_REQUIRED_FLAGS "${_elpa_saved_required_flags}")
else()
  message(FATAL_ERROR
    "USE_ELPA=ON but ELPA was not found. Set ELPA_ROOT=<prefix>, or set "
    "ELPA_INCLUDE_DIR and ELPA_LIBRARY explicitly.")
endif()
