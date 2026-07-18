# Example toolchain fragment for building HPhi with ELPA.
#   cmake -DCONFIG=elpa -DELPA_ROOT=/path/to/elpa/prefix ..
# For GPU execution, ELPA >= 2023.11.001 built with CUDA
# (--enable-nvidia-gpu-kernels) is required; HPhi detects elpa_setup_gpu
# automatically and enables the GPU path (_ELPA_GPU).
set(USE_ELPA ON CACHE BOOL "" FORCE)
set(USE_SCALAPACK ON CACHE BOOL "" FORCE)
# ScaLAPACK is discovered with pkg-config or conventional library names. For
# compiler-suite installations requiring special link flags, set
# SCALAPACK_LIBRARIES explicitly.
# MKL-provided ScaLAPACK (adjust BLACS layer to your MPI):
# set(SCALAPACK_LIBRARIES "-L$ENV{MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64")
