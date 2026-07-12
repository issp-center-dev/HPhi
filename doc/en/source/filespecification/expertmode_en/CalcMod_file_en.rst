.. highlight:: none

.. _Subsec:calcmod:

CalcMod file
------------

This file determines the parameters for the calculation method, model, and output mode. The file format is as follows.

::

    CalcType   0
    CalcModel   2
    CalcEigenVec 0

.. _file_format_1:

File format
~~~~~~~~~~~

[string01] [int01]

.. _parameters_1:

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String

   **Description :** Select a word from keywords.

*  [int01]

   **Type :** Int

   | **Description :** A parameter that is correlated with a keyword.

.. _use_rules_1:

Use rules
~~~~~~~~~

*  After setting the keywords at [string 01], a half-width blank is
   needed for setting a parameter.

*  Keywords can be set in random order.

*  If the keywords or filenames are incorrect, the program is
   terminated.

*  The keywords “CalcType" and “CalcModel" are essential.

*  When a head of line is \"#", the line is skipped.

 

Keywords and parameters
~~~~~~~~~~~~~~~~~~~~~~~

The parameters correlated with the keywords are as follows.

*  ``CalcType``

   **Type :** Int

   | **Description :** Select the method for calculation from the
     following list:
   | 0: Lanczos method
   | 1: mTPQ method
   | 2: Full diagonalization method
   | 3: LOBCG for the ground state
   | 4: Time-evolution
   | 5: cTPQ method

*  ``CalcModel``

   **Type :** Int

   | **Description :** Select the model from the following list:
   | 0: Fermion Hubbard model (canonical ensemble: conservation of
     particles or conservation of particles and the component of
     :math:`S_z`)
   | 1: Spin model (canonical ensemble: conservation of the component of
     :math:`S_z`)
   | 2: Kondo lattice model (canonical ensemble: conservation of
     particles, the component of :math:`S_z`)
   | 3: Fermion Hubbard model (grand canonical ensemble)
   | 4: Spin model (grand canonical ensemble)
   | 5: Kondo lattice model (grand canonical ensemble).
   | 7: Spinless fermion model (canonical ensemble: conservation of particles)
   | 8: Spinless fermion model (grand canonical ensemble).
   | 9: :math:`t`-:math:`J` model (canonical ensemble: conservation of
     particles, or conservation of particles and the component of
     :math:`S_z`)
   | 10: :math:`t`-:math:`J` model (grand canonical ensemble)

   For the fermion Hubbard model, you can select the model under the
   conservation of the particles by setting ``NCond`` in the ModPara
   file. When you want to select the model under the conservation of
   particles and the component of :math:`S_z`, set both ``NCond`` and
   ``2Sz`` in the ModPara file.

   The :math:`t`-:math:`J` models (9, 10) follow the same ``NCond`` /
   ``2Sz`` selection as the fermion Hubbard model: for the canonical model
   (9), setting only ``NCond`` conserves the total number of electrons,
   while setting both ``NCond`` and ``2Sz`` also conserves :math:`S_z`.
   Doubly-occupied sites are excluded from the Hilbert space (the local
   dimension per site is 3: empty, up, or down), and the models are
   available only in the expert mode. Note that for MPI runs the number of
   processes must nevertheless be a power of four, the same as for the
   Fermion Hubbard model (**not** a power of three), because the internal
   representation keeps four states per site.

   For the spinless fermion model, only Trans (hopping) and CoulombInter
   (inter-site interaction) terms are valid. CoulombIntra, Hund, Exchange,
   and PairHop cannot be used since there are no spin degrees of freedom.

*  ``CalcEigenVec``

   **Type :** Int (default value: 0)

   | **Description :** Select the method to calculate the eigenvectors:
   | 0: Lanczos+CG methods (when the convergence of eigenvectors is not
     sufficient for using the Lanczos method, the CG method is applied
     to calculate eigenvectors).
   | 1: Lanczos method.

*  ``InitialVecType``

   **Type :** Int (default value: 0)

   | **Description :** Select the type of an initial vector (:math:`v0`):
   | -1: Real part (:math:`{\rm Re}[v0]]`) and imaginary part  (:math:`{\rm Re}[v0]]`) of the initial 
    vector are give as the normally distributed random numbers. Thus, the normalized initial vectors are uniformly distributed
    on the :math:`N_{\rm H}` dimensional super sphere (:math:`N_{\rm H}` is the dimension of the Hilbert space). 
   | 0: Complex type (:math:`{\rm Re}[v0]\in[-1:1]`, :math:`{\rm Im}[v0]\in[-1:1]` ).
   | 1: Real type (:math:`{\rm Re}[v0]\in[-1:1]`, :math:`{\rm Im}[v0]=0`).

*  ``OutputEigenVec``

   **Type :** Int (default value: 0)

   | **Description :** Select the mode of outputting an eigenvector:
   | 0: Not output an eigenvector
   | 1: Output an eigenvector.

*  ``InputEigenVec``

   **Type :** Int (default value: 0)

   | **Description :** Select the mode of inputting an eigenvector:
   | 0: Not input an eigenvector
   | 1: Input an eigenvector.

*  ``ReStart``

   **Type :** Int (default value: 0)

   | **Description :** Select the mode of inputting a restart vector:
   | 0: Not restart calculation
   | 1: Output a restart vector
   | 2: Input a restart vector and output a new restart vector
   | 3: Input a restart vector.

*  ``CalcSpec``

   **Type :** Int (default value: 0)

   | **Description :** Select the mode of calculating dynamical Green’s functions:
   | 0: Not calculate dynamical Green’s functions
   | 1: (not restart) Input an initial vector and files for generating single excited or pair excited states
   | 2: Input components of triangular diagonal matrix
   | 3: Output both components of triangular diagonal matrix and a restart vector
   | 4: Input both components of triangular diagonal matrix and a restart vector
   | 5: Input and output both components of triangular diagonal matrix and a restart vector.

*  ``OutputHam``

   **Type :** Int (default value: 0)

   | **Description :** Full Diag)Select the mode of outputting Hamiltonian:
   | 0: not output Hamiltonian.
   | 1: output Hamiltonian.

*  ``InputHam``

   **Type :** Int (default value: 0)

   | **Description :** (Full Diag)Select the mode of inputting Hamiltonian:
   | 0: not input Hamiltonian.
   | 1: input Hamiltonian.

   | Note: with ``Solver 3`` (ELPA) and more than one MPI process, the
     Hamiltonian is generated in distributed form, so ``OutputHam``/``InputHam``
     are rejected at startup. To output or input the Hamiltonian, run with
     1 MPI process or use a different ``Solver``.

*  ``OutputExcitedVec``

   **Type :** Int (default value: 0)

   | **Description :** Select the mode of outputting an excited vector:
   | 0: Not output an eigenvector
   | 1: Output an eigenvector.
   
*  ``OutputDataHead``

   **Type :** Int (default value: 0)

   | **Description :** Select whether to prefix TPQ/TE physical quantity output filenames (``SS``, ``Norm``, ``Flct``) with the header string defined by ``CDataFileHead`` in the ModPara file:
   | 0: Do not add a prefix (e.g., ``SS_rand0.dat``).
   | 1: Add the ``CDataFileHead`` prefix (e.g., ``zvo_SS_rand0.dat``).
   | When ``OutputGreenFormat=1`` is used for TPQ/cTPQ, the aggregate physical quantity filenames are ``SS_tpq.dat``, ``Norm_tpq.dat``, and ``Flct_tpq.dat``; ``OutputDataHead=1`` prefixes these names in the same way.

*  ``OutputGreenFormat``

   **Type :** Int (default value: 0)

   | **Description :** Select the output format for Green function files and TPQ/cTPQ physical quantity files:
   | 0: Existing split files.
   | 1: Aggregate indexed files for TPQ/cTPQ, real-time evolution, Full diagonalization, and LOBCG.
   | In aggregate mode, TPQ/cTPQ Green function rows and TPQ/cTPQ physical quantity rows start with ``set`` and ``step``, real-time evolution rows start with ``step``, and Full diagonalization/LOBCG rows start with ``eigen``.
   | TPQ/cTPQ physical quantities are written to ``SS_tpq.dat``, ``Norm_tpq.dat``, and ``Flct_tpq.dat`` instead of ``SS_rand*.dat``, ``Norm_rand*.dat``, and ``Flct_rand*.dat``.
   | ``AnomalousG`` in LOBCG keeps the existing non-aggregate output because the existing output is not split by eigen index.

*  ``Solver``

   **Type :** int (default: resolved from legacy keywords)

   | **Description :** (FullDiag)
     Diagonalization backend for the full diagonalization method:
   | 0 (LAPACK, serial), 1 (ScaLAPACK), 2 (MAGMA, single-node multi-GPU),
     3 (ELPA, multi-node CPU/GPU; requires a build with ``USE_ELPA=ON``).
   | When omitted, the backend is resolved from the legacy ``Scalapack``
     and ``NGPU`` keywords so that existing inputs behave as before.

*  ``Scalapack``

   **Type :** Int (default value: 0)

   | **Description :** (Full Diag)Select to use ScaLAPACK library for full diagonalization:
   | 0: not to use ScaLAPACK.
   | 1: use ScaLAPACK.
   | (Deprecated) This keyword is retained for backward compatibility.
     Please use ``Solver 1`` for new inputs.


*  ``NGPU``

   **Type :** Int (default value: see below)

   | **Description :** (FullDiag)
     Number of GPU devices per node for full diagonalization.
   | Default value: when ``Solver`` is not given explicitly, 2 on
     MAGMA-enabled builds and 0 otherwise (legacy behavior, unchanged).
     When ``Solver`` is given explicitly, the default is 2 only for
     ``Solver 2`` (MAGMA); it is 0 for ``Solver 0/1/3`` (ELPA defaults to
     CPU execution; GPU use must be opted into explicitly).
   | For ``Solver 2`` (MAGMA), this specifies the number of GPUs used by a single process.
   | For ``Solver 3`` (ELPA), 0 runs on CPU, and a value >= 1 runs on GPU with automatic GPU assignment (1 process per GPU).
   | It is recommended to match the number of MPI processes per node to the number of GPUs.
   | Note that ``NGPU`` does not physically limit the number of GPU devices.
     To strictly limit GPU count, use the job scheduler (``CUDA_VISIBLE_DEVICES``, etc.).
   | GPU execution with ELPA requires ELPA version >= 2023.11.001.
   | For ``Solver 3``, the Hilbert-space dimension (matrix size) must be at least
     as large as the largest dimension of the MPI process grid (max of nprow, npcol);
     otherwise HPhi exits at startup with an error asking to reduce the number of MPI ranks.
   | With ``Solver 3`` and more than one MPI process, the Hamiltonian is generated
     and stored in distributed (block-cyclic) form, so peak memory per rank scales
     roughly as O(N²/P) rather than O(N²) — larger Hilbert-space dimensions become
     feasible by increasing the process count.

*  ``ExpecMode``

   **Type :** Int (default value: 0)

   | **Description :** (FullDiag)
     Selects the evaluation kernel for full-diagonalization observables
     (energy, N, Sz, S2, doublon, and Green functions):
   | 0: Conventional evaluation (existing behavior).
   | 1: State-task-parallel evaluation. Each MPI rank evaluates the
     observables for its own contiguous block of eigenstates
     independently (no communication during the per-state evaluation
     loop). The aggregate Green-function output is written as
     rank-local partial files and merged into the final aggregate files
     by rank 0 once every rank's manifest reports success.
   | 2: Reserved for a future trace-kernel evaluation mode. Not yet
     implemented; currently runs as ``ExpecMode 1`` and prints
   | ``INFO: ExpecMode 2 kernels are not available in this build; running as ExpecMode 1.``
   | Eligibility: a nonzero ``ExpecMode`` requires ``CalcType`` = 2 (full
     diagonalization) together with ``Solver`` 1 (ScaLAPACK) or 3 (ELPA);
     any other combination (wrong ``CalcType`` or ``Solver``) is rejected
     at startup with an error. With exactly one MPI process, ``ExpecMode``
     is automatically reverted to 0 (results are identical for a single
     process either way), printing
   | ``INFO: ExpecMode reverts to 0 for a single process (results are identical).``
   | Guarantee: ``ExpecMode`` changes only evaluation speed, never the
     physics -- ``ExpecMode`` 0, 1, and 2 produce identical results up to
     floating-point rounding (summation order differs between kernels, so
     agreement is not bit-identical).
   | Memory: ``ExpecMode 1`` additionally holds a state panel roughly the
     same size as the distributed eigenvector storage used by ``Solver 3``
     (O(N2/P) per rank). During the one-time redistribution step both the
     original storage and the new panel coexist, giving a temporary peak
     of roughly 2xO(N2/P) per rank before the original storage is freed;
     steady-state usage afterward is O(N2/P), unchanged from ``ExpecMode 0``.
   | Guidance: prefer ``ExpecMode 1`` for large multi-node ``Solver 3``
     (ELPA) or ``Solver 1`` (ScaLAPACK) runs with many eigenstates and/or
     many Green-function observables, where re-evaluating every eigenstate
     redundantly on every rank (``ExpecMode 0``) becomes the bottleneck.
     Keep the default ``ExpecMode 0`` otherwise, including for small
     systems and single-process runs.
   | Note (behavior fix): as of this phase, distributed FullDiag runs
     (``Solver`` 1 or 3, more than one MPI process) compute S2 and Sz on
     rank 0 for ``ExpecMode 0`` as well -- previously these distributed runs
     zero-filled S2 and Sz and printed a shortened stdout progress line
     without the S2 column. The stdout progress line now always matches
     the single-process (serial) format, with the S2 column included, for
     every ``Solver``/``ExpecMode`` combination. This is an intentional
     correctness fix, independent of ``ExpecMode``'s value.

.. raw:: latex

   \newpage
