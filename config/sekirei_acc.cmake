# for Intel Compiler
set(CMAKE_C_COMPILER "icc" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS "-O3 -DNDEBUG -xCORE-AVX2 -mcmodel=large -shared-intel" CACHE STRING "" FORCE)
set(CMAKE_Fortran_COMPILER "ifort" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS "-O3 -DNDEBUG -xCORE-AVX2 -mcmodel=large -shared-intel" CACHE STRING "" FORCE)

# for Intel MKL
set(BLA_VENDOR "Intel10_64lp" CACHE STRING "" FORCE)

# for MAGMA
set (MAGMA_FOUND ON)
set (MAGMA_C_LIBRARIES magma)
# set magma installed directory
set (MAGMA /home/issp/materiapps/tool/magma/magma-2.2.0)
include_directories(${MAGMA}/include)
link_directories(${MAGMA}/lib)

# for CUDA
find_package(CUDA)
if(CUDA_FOUND)
  include_directories(${CUDA_INCLUDE_DIRS})
endif(CUDA_FOUND)
