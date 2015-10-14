# for Fujitsu Compiler
set(CMAKE_C_COMPILER "fccpx" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELEASE "-Kfast,SPARC64IXfx,parallel -Kmemalias,alias_const" CACHE STRING "" FORCE)
set(OpenMP_C_FLAGS "-Kopenmp" CACHE STRING "" FORCE)

# for SSL2
set(LAPACK_FOUND TRUE CACHE BOOL "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS "-SSL2 --linkfortran" CACHE STRING "" FORCE)
