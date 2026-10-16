.. highlight:: none

.. _Subsec:restart:

tmpvec files
------------

The TPQ/cTPQ and LOBPCG (CG) solvers use these binary files for the
legacy ``ReStart`` mechanism. Ground-state Lanczos does not use this
format; it uses the ``recalcvec`` and ``TMComponents`` files instead.

File name
~~~~~~~~~

*  ``tmpvec_set&&_rank_$$.dat``

``&&`` is the zero-based TPQ sample or LOBPCG eigenstate number, and
``$$`` is the MPI rank. The name does **not** contain ``CDataFileHead``.

File format
~~~~~~~~~~~

On the supported 64-bit platforms, the file starts with a 12-byte
header followed by one complex vector:

::

    fwrite(&step, sizeof(int), 1, fp);
    fwrite(&local_size, sizeof(unsigned long int), 1, fp);
    fwrite(&vector[0], sizeof(complex double), local_size + 1, fp);

``step`` is the saved iteration and ``local_size`` is the local vector
dimension. Each file stores ``local_size + 1`` complex values because
index zero is retained for compatibility; ``vector[0]`` is not used in
the calculation. The binary representation is platform-dependent.

.. raw:: latex

   \newpage
