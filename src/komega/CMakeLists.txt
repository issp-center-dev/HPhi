# include guard
cmake_minimum_required(VERSION 2.8)
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of HPhi.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

set(SOURCES_KOMEGA komega_bicg.F90 komega_math.F90 komega_vals.F90)
add_definitions(${MPI_Fortran_COMPILE_FLAGS})
include_directories(${MPI_Fortran_INCLUDE_PATH})
link_directories(${MPI_Fortran_LIBRARIES})
add_library(komega ${SOURCES_KOMEGA})
target_link_libraries(komega ${MPI_Fortran_LIBRARIES} ${LAPACK_LIBRARIES})
install(TARGETS komega ARCHIVE DESTINATION lib LIBRARY DESTINATION lib RUNTIME DESTINATION bin)
