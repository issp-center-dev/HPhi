cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
project(HPhi NONE)

if(CONFIG)
  message(STATUS "Loading configration: " ${PROJECT_SOURCE_DIR}/config/${CONFIG}.cmake)
  include(${PROJECT_SOURCE_DIR}/config/${CONFIG}.cmake)
endif(CONFIG)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
enable_language(C Fortran)

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif(OPENMP_FOUND)

find_package(MPI)
if(MPI_FOUND)
  include_directories(${MPI_C_INCLUDE_PATH})
  add_definitions(-DMPI)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LINK_FLAGS}")
endif(MPI_FOUND)

find_package(LAPACK)
if(LAPACK_FOUND)
  add_definitions(-Dlapack)
endif(LAPACK_FOUND)

add_subdirectory(src)
