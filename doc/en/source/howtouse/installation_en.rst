Installation
============

You can download :math:`{\mathcal H}\Phi` at the following location.

https://github.com/issp-center-dev/HPhi/releases

You can obtain the :math:`{\mathcal H}\Phi` directory by typing

``$ tar xzvf HPhi-xxx.tar.gz``

:math:`{\mathcal H}\Phi` can be installed by using cmake.
We can compile :math:`{\mathcal H}\Phi` as::

 cd $HOME/build/hphi
 cmake -DCONFIG=gcc $PathTohphi
 make

To use ScaLAPACK library for full diagonalization, the cmake option ``-DUSE_SCALAPACK=ON`` is needed.
To use ELPA library for full diagonalization, the cmake option ``-DUSE_ELPA=ON`` is needed
(ELPA requires ScaLAPACK and MPI; ``USE_SCALAPACK`` is enabled automatically).
If ELPA is installed in a non-standard path, specify its location via ``-DELPA_ROOT=/path/to/elpa``
(individual overrides ``ELPA_INCLUDE_DIR``/``ELPA_LIBRARY`` are also supported).

.. note::

   Platform notes for ELPA builds (verified on an AMD EPYC + Mellanox
   InfiniBand cluster with Intel oneAPI 2022 / Intel MPI 2021.7; both
   issues were actually encountered and diagnosed on the ISSP
   supercomputer kugui):

   * On CPUs **without AVX-512** (e.g. AMD EPYC Rome/Milan), configure
     ELPA itself with ``--disable-avx512 --disable-avx512-kernels``.
     ELPA compiles AVX-512 kernels regardless of the build host and its
     runtime kernel selection may pick one, crashing with an illegal
     instruction (SIGILL) inside the 2-stage solver
     (``single_hh_trafo_*_AVX512_*``) for larger matrices, while small
     matrices (1-stage path) appear to work.
   * With **Intel MPI** on Mellanox (mlx/UCX provider), multi-rank
     FullDiag runs with ``ExpecMode`` 1 or 2 may deadlock in the
     observable-collection step (``MPI_Gatherv``) with the default
     tuned algorithm. Set ``I_MPI_ADJUST_GATHERV=3`` (works multi-node)
     or ``I_MPI_FABRICS=shm`` (single node only) as a workaround.
Here, we set a path to :math:`{\mathcal H}\Phi` as ``$PathTohphi``
and to a build directory as ``$HOME/build/hphi``.
After compilation, ``src`` folder is constructed below a ``$HOME/build/hphi``
folder and we obtain an executable ``HPhi`` in ``src/`` directory.
When no MPI library exists in the system, an executable ``HPhi``
is automatically compiled without an MPI library.

In the above example,
we compile :math:`{\mathcal H}\Phi` by using a gcc compiler.
We can select a compiler by using the following options:

 * ``sekirei`` : ISSP system-B \"sekirei\"
 * ``sekirei_acc``: ISSP system-B \"sekirei\" (for using MAGMA library)
 * ``fujitsu`` : Fujitsu compiler
 * ``intel`` : Intel compiler + Linux PC
 * ``gcc`` : GCC compiler + Linux PC.

An example of compiling :math:`{\mathcal H}\Phi` by using the Intel compiler is shown as follows::

 mkdir ./build
 cd ./build
 cmake -DCONFIG=intel ../
 make

After compilation,
``src`` folder is created below the ``build`` folder and
an execute :math:`{\mathcal H}\Phi` in the ``src`` folder.
Please note that we must delete the ``build`` folder and
repeat the above operations when we change the compiler.


.. tip::

 | CMake automatically trys to download StdFace package (https://github.com/issp-center-dev/StdFace, parser for standard mode),
 | but this may fail due to network problem (e.g., IP unreachable)
 | In such a case, please run ``sh src/StdFace/download.sh`` to download manually.
 |
 | If you want to use StdFace downloaded into another directory,
 | pass ``-DSTDFACE_DIR=<path_to_stdface>`` to cmake command.
