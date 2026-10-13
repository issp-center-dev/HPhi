Restart Mechanism
=================

This section describes the legacy restart and eigenvector-input
mechanisms in :math:`{\mathcal H}\Phi`.

``ReStart`` modes
-----------------

In the ``CalcMod`` input file, ``ReStart`` selects solver continuation:

* ``0``: no restart
* ``1``: write restart data
* ``2``: read restart data and write updated restart data
* ``3``: read restart data

The legacy files depend on the solver:

* Ground-state Lanczos writes ``##_recalcvec_rank_$$.dat`` and
  ``##_TMComponents.dat``. The latter records the MPI process and OpenMP
  thread counts. A mismatch is rejected; an older file without the
  environment line is accepted with a warning.
* LOBPCG (CG), TPQ, and cTPQ write
  ``tmpvec_set&&_rank_$$.dat``. These names do not contain
  ``CDataFileHead``.
* Real-time evolution writes ``##_eigenvec_&&_rank_$$.dat`` when
  ``OutputEigenVec=1``. The header records the next step to execute.

See the output-file reference for the binary layouts. These formats are
platform-dependent.

``InputEigenVec`` modes
-----------------------

``InputEigenVec`` loads a saved eigenvector; it is distinct from
``ReStart``:

* ``0``: do not load an eigenvector.
* ``1``: load the eigenvector and use it directly.
* ``2``: load the eigenvector, apply :math:`H` once, and normalize the
  resulting vector before evaluating observables (legacy Lanczos input
  path).

Loading an eigenvector does not resume the eigensolver iteration.
``OutputEigenVec=1`` writes the resulting eigenvector files.

Iteration limits
----------------

For legacy Lanczos and LOBPCG restart, ``Lanczos_max`` is the number of
additional iterations after the saved point. For TPQ/cTPQ and real-time
evolution it is the total step count, so it must exceed the saved step.

MPI requirements
----------------

Each MPI process reads its own rank file. The number of MPI processes
must match the run that wrote the legacy restart or eigenvector files;
changing the process count is not supported by these formats. Use the
same model, Hilbert-space decomposition, OpenMP thread count where
required, and compatible executable architecture.

Legacy output limitations
-------------------------

The current real-time evolution continuation reopens ``SS.dat`` and
``Norm.dat`` in write mode, so earlier time-series rows are replaced.
The legacy TPQ/cTPQ ``tmpvec`` file does not store
``global_1st_norm``; continued ``Norm_rand`` rows therefore use zero in
that column.

Before a long continuation, keep backup copies of the restart files and
the time-series outputs, and verify the conserved quantities after the
restart.
