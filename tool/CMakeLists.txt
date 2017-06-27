# include guard
cmake_minimum_required(VERSION 2.8.0 )

add_definitions(-D_HPhi)
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of HPhi.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

add_library(key2lower key2lower.c)
add_executable(fourier fourier.F90)
add_executable(corplot corplot.F90)
add_executable(cTPQ cTPQ.F90)
target_link_libraries(fourier key2lower ${LAPACK_LIBRARIES})

install(TARGETS fourier RUNTIME DESTINATION bin)
install(TARGETS corplot RUNTIME DESTINATION bin)

