cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
project(HPhi)

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif(OPENMP_FOUND)
find_package(LAPACK)
add_subdirectory(src)

