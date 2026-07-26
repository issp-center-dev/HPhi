.. highlight:: none

.. _Ch:Prerequisite:

Prerequisite
============

:math:`{\mathcal H}\Phi` requires the following packages:

 * C/fortran compiler (Intel, Fujitsu, GNU, etc. )
 * BLAS/LAPACK library (Intel MKL, Fujitsu, ATLAS, etc.)
 * MPI library (if you do not use MPI, this is not required).
 * ScaLAPACK library (if you do not use it for full diagonalization, this is not required).
 * MAGMA library (if you do not use it for full diagonalization, this is not required).
 * ELPA library (if you do not use it for full diagonalization, this is not required. The non-threaded variant is required: version 2021.11 or later for CPU execution, and a CUDA-enabled build of version 2023.11.001 or later for GPU execution).

.. tip::

 | **E.g. /Settings of Intel compiler**
 | When you use the Intel compiler, you can easily use the scripts attached to the compiler.
 | In the case of the bash in the 64-bit OS, write the following in your ``~/.bashrc``\:
 
 | ``source /opt/intel/bin/compilervars.sh intel64``
 | or
 | ``source /opt/intel/bin/iccvars.sh intel64``
 | ``source /opt/intel/mkl/bin/mklvars.sh``
 
 Please read the manuals of your compiler/library for more information.

