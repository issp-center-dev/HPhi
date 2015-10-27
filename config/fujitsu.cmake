# for Fujitsu Compiler
set(CMAKE_C_COMPILER "mpifccpx" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELEASE "-DNDEBUG -Kfast,parallel -Kmemalias,alias_const" CACHE STRING "" FORCE)
set(OpenMP_C_FLAGS "-Kopenmp" CACHE STRING "" FORCE)

# for SSL2
set(BLAS_LIBRARIES "-SSL2 --linkfortran" CACHE STRING "" FORCE)
