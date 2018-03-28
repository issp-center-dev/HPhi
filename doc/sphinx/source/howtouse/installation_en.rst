Installation
============

You can download :math:`{\mathcal H}\Phi` at the following location.

https://github.com/QLMS/HPhi/releases

You can obtain the :math:`{\mathcal H}\Phi` directory by typing

``$ tar xzvf HPhi-xxx.tar.gz``

There are two types of procedure for installing :math:`{\mathcal H}\Phi`.

Using ``HPhiconfig.sh``
-----------------------

Please run ``HPhiconfig.sh`` script in the :math:`{\mathcal H}\Phi` directory as follows
(for ISSP system-B \"sekirei\"):

``$ bash HPhiconfig.sh sekirei``

Then, the environmental configuration file ``make.sys`` is generated in 
the ``src/`` directory.
The command-line argument of ``HPhiconfig.sh`` is as follows:

 * ``sekirei`` : ISSP system-B \"sekirei\"
 * ``fujitsu`` : ISSP system-C \"maki\"
 * ``sr`` : HITACHI SR16000
 * ``intel`` : Intel compiler
 * ``intel-openmpi`` : Intel compiler + OpenMPI
 * ``intel-mpich`` : Intel compiler + MPICH2
 * ``intel-intelmpi`` : Intel compiler + IntelMPI
 * ``gcc`` : GCC
 * ``gcc-openmpi`` : GCC + OpenMPI
 * ``gcc-mpich`` : GCC + MPICH2
 * ``gcc-mac`` : GCC + Mac.

``make.sys`` is as follows (for ISSP-system-B \"sekirei\")::

 CC = mpicc
 F90 = mpif90
 CFLAGS = -fopenmp -O3 -g -traceback -xHost -ipo -mcmodel=large \
          -shared-intel -D MPI
 FFLAGS = -fopenmp -O3 -g -traceback -xHost -ipo -mcmodel=large \
          -shared-intel -D MPI -fpp
 LIBS = -mkl -lifcore

We explain the macros of this file as: 

 * ``CC`` : C compiler (``icc``, ``gcc``, ``fccpx``).
 * ``F90`` : fortran compiler (``ifort``, ``gfortran``, ``frtpx``)
 * ``CFLAGS`` : C compile options. OpenMP utilization option (``-openmp``, ``-fopenmp``, ``-qopenmp``, etc.) must be specified. When you use MPI, please set ``-D MPI``.
 * ``FFLAGS`` : fortran compile options. Similar to ``CFLAGS``. 
 * ``LIBS`` : Compilation options for LAPACK. ``-Dlapack`` can not be removed.

Then, you are ready to compile HPhi.
Please type

 ``$ make HPhi``

and obtain an executable ``HPhi`` in the ``src/`` directory;
you should add this directory to the ``$PATH``.

If SSE2 is available in your system, please add ``-DHAVE_SSE2`` as an option of CMake.

.. tip::

 | You can make a PATH to :math:`{\mathcal H}\Phi` as follows:
 | ``$ export PATH=${PATH}:HPhi_top_directory/src/``
 | If you retain this PATH, you should write above in ``~/.bashrc`` (for ``bash`` as a login shell)

Using ``cmake``
---------------

.. tip::

 | Before using cmake for sekirei, you must type
 | ``source /home/issp/materiapps/tool/env.sh``
 | while for maki, you must type
 | ``source /global/app/materiapps/tool/env.sh``

We can compile :math:`{\mathcal H}\Phi` as::

 cd $HOME/build/hphi
 cmake -DCONFIG=gcc $PathTohphi
 make

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
 * ``fujitsu`` : Fujitsu compiler (ISSP system-C \"maki\")
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
