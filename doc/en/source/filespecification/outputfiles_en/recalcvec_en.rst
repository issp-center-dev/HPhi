.. highlight:: none

recalcvec.dat
-------------

This per-rank binary file stores the vectors required to continue a
ground-state Lanczos calculation or a shifted-BiCG spectrum calculation.

File name
~~~~~~~~~

*  ``##_recalcvec_rank_$$.dat``

``##`` is ``CDataFileHead`` in the ModPara file, and ``$$`` is the MPI
rank.

File format
~~~~~~~~~~~

Ground-state Lanczos writes:

::

    unsigned long step;
    unsigned long local_size;
    complex double v0[local_size + 1];
    complex double v1[local_size + 1];

The two header fields are written separately, followed by the two
Lanczos vectors. Index zero of each vector is unused.

Shifted BiCG writes a different legacy payload:

::

    int step_or_status;
    unsigned long local_size;
    complex double residual[local_size + 1];
    complex double old_residual[local_size + 1];
    complex double shadow_residual[local_size + 1];
    complex double old_shadow_residual[local_size + 1];

The binary representation is platform-dependent. A restart must use
the same MPI process count and local dimensions as the run that wrote
the files.

.. raw:: latex

   \newpage
