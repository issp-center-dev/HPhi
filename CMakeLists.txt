cmake_minimum_required(VERSION 3.0...3.10)
enable_testing()

project(HPhi NONE)
# Capture whether the user explicitly requested USE_SCALAPACK=OFF (as opposed to
# the default OFF). -D cache vars are defined before option() runs, so a plain
# `if(DEFINED USE_SCALAPACK AND NOT USE_SCALAPACK)` would look identical for
# (a) a fresh `-DUSE_SCALAPACK=OFF` and (b) a *cached* OFF left over from a
# previous configure in the same build dir (e.g. the default from option()
# below, written to CMakeCache.txt on the prior run). Re-running cmake with
# just `-DUSE_ELPA=ON` in an existing build dir would then hit case (b) and
# misreport it as an explicit conflicting USE_SCALAPACK=OFF.
#
# Distinguish them via the cache entry's HELPSTRING: -D command-line cache
# vars get CMake's fixed helpstring "No help, variable specified on the
# command line."; the option() default below carries our own helpstring
# ("Use Scalapack"). Only the command-line form counts as "explicitly set
# OFF" for the USE_ELPA conflict check further down.
#
# Remember across re-configures that the user explicitly disabled ScaLAPACK
# (option() rewrites USE_SCALAPACK's HELPSTRING after the first run, so the
# command-line origin must be recorded in a separate sticky cache entry:
# without it, run 1 `-DUSE_SCALAPACK=OFF` followed by run 2 `-DUSE_ELPA=ON`
# in the same dir would silently auto-enable ScaLAPACK over the user's
# explicit OFF instead of erroring). Also keep a non-cache flag for the
# current configure: CONFIG=elpa may turn USE_SCALAPACK back ON later in this
# same run, but it must not erase an explicit command-line OFF before the
# conflict check has had a chance to report it.
set(_hphi_use_scalapack_explicit_off_this_run FALSE)
if(DEFINED USE_SCALAPACK AND NOT USE_SCALAPACK)
  # A toolchain may define a normal (non-cache) variable. Querying CACHE
  # properties for such a variable is an error, so first use get_property's
  # old-CMake-compatible SET form to confirm that a cache entry exists.
  get_property(_us_is_cache CACHE USE_SCALAPACK PROPERTY TYPE SET)
  if(_us_is_cache)
    get_property(_us_help CACHE USE_SCALAPACK PROPERTY HELPSTRING)
    if(_us_help STREQUAL "No help, variable specified on the command line.")
      set(HPHI_USE_SCALAPACK_EXPLICIT_OFF TRUE CACHE INTERNAL
          "User explicitly passed -DUSE_SCALAPACK=OFF")
      set(_hphi_use_scalapack_explicit_off_this_run TRUE)
    endif()
  endif()
endif()
option(USE_SCALAPACK "Use Scalapack" OFF)
option(USE_ELPA "Use ELPA for FullDiag" OFF)

set(CMAKE_C_FLAGS "-Wall")

if(CONFIG)
  message(STATUS "Loading configration: " ${PROJECT_SOURCE_DIR}/config/${CONFIG}.cmake)
  include(${PROJECT_SOURCE_DIR}/config/${CONFIG}.cmake)
endif(CONFIG)
if(USE_SCALAPACK AND NOT _hphi_use_scalapack_explicit_off_this_run)
  # The command line or a configuration file (re-)enabled ScaLAPACK.
  # Clear a marker left by an earlier configure in this build directory.
  unset(HPHI_USE_SCALAPACK_EXPLICIT_OFF CACHE)
endif()
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(ENABLE_MPI "Enable MPI Parallelization" ON)
option(ENABLE_SCALAPACK "Use Scalapack" ON)
option(GIT_SUBMODULE_UPDATE "Execute `git submodule update --init` automatically" ON)
enable_language(C Fortran)

# HPhi is a C executable, but ScaLAPACK and ELPA are commonly implemented in
# Fortran. Shared libraries carry their runtime dependencies themselves;
# static archives do not. Preserve the Fortran compiler's implicit search
# directories and runtime libraries so configure probes and final/test links
# work for both forms. MPI_Fortran_LIBRARIES is added at each use site after
# find_package(MPI), because a static ScaLAPACK may also call the Fortran MPI
# interface.
set(HPHI_FORTRAN_RUNTIME_LIBRARIES)
foreach(_hphi_fortran_runtime_dir ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES})
  list(APPEND HPHI_FORTRAN_RUNTIME_LIBRARIES
       "-L${_hphi_fortran_runtime_dir}")
endforeach()
list(APPEND HPHI_FORTRAN_RUNTIME_LIBRARIES
     ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})

# Set C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

#set(CMAKE_C_FLAGS "-lscalapack -lgfortran -lm")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_MACOSX_RPATH 1)

if(NOT STDFACE_DIR)
  set(STDFACE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/StdFace")
  # git submodule update
  # ref: https://cliutils.gitlab.io/modern-cmake/chapters/projects/submodule.html
  find_package(Git QUIET)
  if(GIT_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
    # Update submodules as needed
    if(GIT_SUBMODULE_UPDATE)
      message(STATUS "Submodule update")
      execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        RESULT_VARIABLE GIT_SUBMOD_RESULT)
      if(NOT GIT_SUBMOD_RESULT EQUAL "0")
        message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
      endif()
    endif()
  endif()
  if(NOT EXISTS "${STDFACE_DIR}/CMakeLists.txt")
    message(FATAL_ERROR "StdFace is not found. Please `git clone` the repository (https://github.com/issp-center-dev/HPhi) or download an archive file (HPhi-X.Y.Z.tar.gz) from the release page (https://github.com/issp-center-dev/HPhi/releases)")
  endif()
else()
  if(NOT EXISTS "${STDFACE_DIR}/CMakeLists.txt")
    message(FATAL_ERROR "STDFACE_DIR is manually set to ${STDFACE_DIR}, but ${STDFACE_DIR}/CMakeLists.txt does not exists")
  endif()
endif()

# OpenMP detection
find_package(OpenMP)
# OpenMP_C_FOUND is defined by FindOpenMP only since CMake 3.9; older CMake
# (>= our 3.0 minimum) sets the uppercase OPENMP_FOUND instead. Check both so
# OpenMP is not silently disabled for non-Intel compilers on CMake 3.0-3.8.
if(OpenMP_C_FOUND OR OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
elseif(CMAKE_C_COMPILER_ID STREQUAL "Intel")
  # Fallback for older Intel compilers
  if("${CMAKE_C_COMPILER_VERSION}" VERSION_LESS "15.0.0.20140528")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -openmp")
  else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qopenmp")
  endif()
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

if(ENABLE_MPI)
  find_package(MPI)
  if(MPI_FOUND)
    include_directories(${MPI_C_INCLUDE_PATH} ${MPI_C_INCLUDE_DIRS})
    add_definitions(-DMPI)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LINK_FLAGS}")
  endif(MPI_FOUND)
endif(ENABLE_MPI)

# ELPA requires a working MPI implementation and ScaLAPACK. Resolve
# USE_SCALAPACK before dependency discovery, so the block below emits exactly
# one of -D_lapack / -D_SCALAPACK (never both).
if(USE_ELPA)
  if(NOT ENABLE_MPI)
    message(FATAL_ERROR "USE_ELPA=ON requires MPI: configure with -DENABLE_MPI=ON (default).")
  endif()
  if(NOT MPI_FOUND)
    message(FATAL_ERROR
      "USE_ELPA=ON requires a working MPI C implementation, but CMake could not find one.")
  endif()
  if(HPHI_USE_SCALAPACK_EXPLICIT_OFF)
    message(FATAL_ERROR "USE_ELPA=ON requires ScaLAPACK, but USE_SCALAPACK=OFF was set explicitly.")
  endif()
  if(NOT USE_SCALAPACK)
    message(STATUS "USE_ELPA=ON requires ScaLAPACK: enabling USE_SCALAPACK.")
    set(USE_SCALAPACK ON CACHE BOOL "Use Scalapack" FORCE)
  endif()
endif(USE_ELPA)

find_package(LAPACK)
if(USE_SCALAPACK)
  find_package(ScaLAPACK REQUIRED)
  add_definitions(-D_SCALAPACK)
else()
  if(LAPACK_FOUND)
    add_definitions(-D_lapack)
  endif(LAPACK_FOUND)
endif()

if(USE_ELPA)
  find_package(ELPA REQUIRED)
endif(USE_ELPA)

add_subdirectory(src)
add_subdirectory(tool)
add_subdirectory(test)

if (Document)
  add_subdirectory(doc)
endif(Document)
