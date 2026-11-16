.. highlight:: none

TMcomponents.dat
----------------

This text file stores the tridiagonal-matrix components and the initial
norm used by Lanczos restart and recalculation.

File name
~~~~~~~~~

*  ``##_TMComponents.dat``

``##`` is ``CDataFileHead`` in the ModPara file.

File format
~~~~~~~~~~~

*  Line 1: ``step``
*  Line 2: ``dnorm``
*  The next ``step`` lines: ``alpha beta``
*  Ground-state restart only, final line:
   ``# nproc P nthreads T``

``alpha`` and ``beta`` are the Lanczos tridiagonal components
:math:`\alpha_i` and :math:`\beta_i`. The final comment records the MPI
process and OpenMP thread counts used to write the corresponding
``recalcvec`` files. A ground-state restart rejects a different process
or thread count. Older files without this comment are accepted with a
warning and therefore must be restarted with the same environment.

Spectrum ``CalcSpec`` files retain their existing text format and do not
contain the environment comment.

.. raw:: latex

   \newpage
