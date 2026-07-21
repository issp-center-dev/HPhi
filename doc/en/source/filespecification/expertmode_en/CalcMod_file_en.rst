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
   | When omitted for FullDiag, the backend is resolved from the legacy
     ``Scalapack`` and ``NGPU`` keywords using the existing precedence.
   | **Usage example (expert mode):** write ``Solver`` (and optionally
     ``ExpecMode``) in ``calcmod.def``,

   ::

       CalcType     2
       CalcModel    0
       Solver       3
       ExpecMode    2

   | and run HPhi in expert mode with the desired number of MPI processes:

   ::

       mpiexec -np 4 HPhi -e namelist.def

   | **Usage from standard mode:** the standard-mode input file does not
     accept the ``Solver``/``ExpecMode`` keywords. Generate the expert-mode
     files once with ``HPhi -sdry stan.in``, append the keywords to the
     generated ``calcmod.def``, and **run in expert mode (**\ ``-e``\ **)**
     as above (running in standard mode ``-s`` regenerates the def files,
     so the edits would not take effect). For the same reason, re-running
     ``HPhi -sdry`` after changing ``stan.in`` overwrites ``calcmod.def``
     and discards the appended keywords, so they must be appended again.
     If ``stan.in`` contained the legacy ``Scalapack``/``NGPU`` keywords,
     we recommend removing those lines from the generated ``calcmod.def``
     and specifying ``Solver`` instead.
     See :ref:`Sec:ParallelFullDiag` for the parallel algorithms and a
     performance comparison of the backends.

*  ``Scalapack``

   **Type :** Int (default value: 0)

   | **Description :** (Full Diag)Select to use ScaLAPACK library for full diagonalization:
   | 0: not to use ScaLAPACK.
   | 1: use ScaLAPACK.
   | (Deprecated) This keyword is retained for backward compatibility.
     Please use ``Solver 1`` for new inputs.
   | This legacy keyword applies only to ``CalcType=2`` (FullDiag). For
     another calculation type it is ignored with a warning, preserving the
     normal MPI site decomposition. By contrast, an explicit ``Solver`` value
     1, 2, or 3 outside FullDiag is rejected as an input error.


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
   | 2: Trace-kernel evaluation. Three quantity families use this
     kernel: the one-body (``expec_cisajs``-equivalent) and two-body
     (``expec_cisajscktaltdc``-equivalent) Green functions, and, as of
     this version, the energy/fluctuation family (including the ``var``
     column). For the Green functions, HPhi precomputes, for each
     operator, the basis-state mapping (destination state and amplitude)
     once, then streams every owned eigenstate through that mapping in a
     dense loop, amortizing the per-state operator overhead that
     ``ExpecMode 1`` still pays for these two quantities. For the
     energy/fluctuation family, HPhi instead precomputes the Hamiltonian
     once per rank into a compact CSR (compressed sparse row) matrix and
     streams every owned eigenstate through a CSR sparse
     matrix-vector product, amortizing the per-state
     ``mltply``-style traversal overhead in the same way. ``S2``,
     ``NBodyG``, and ``AnomalousG`` are always evaluated on the
     ``ExpecMode 1`` path in this version -- trace-kernel evaluation for
     these quantities is a candidate for a future phase.
   | Note on terminology: "trace kernel" names this precomputed-mapping,
     streaming evaluation technique -- it does not refer to the matrix
     trace :math:`{\rm Tr}(\cdot)`.
   | Whether the one-body/two-body trace kernel is actually used is
     decided once per run from a per-(model, quantity) capability table
     plus three runtime checks, and reported at the start of the run by a
     rank-0 ``INFO`` line per quantity (indented in the actual log
     output; ``%s`` below stands for ``one-body`` or ``two-body``):
   | ``INFO: ExpecMode 2: %s Green functions use the trace kernel.``
   | ``INFO: ExpecMode 2: %s Green functions use the ExpecMode-1 fallback (unsupported model).``
   | ``INFO: ExpecMode 2: %s Green functions use the ExpecMode-1 fallback (no operators of this kind are defined).``
   | ``INFO: ExpecMode 2: %s Green functions use the ExpecMode-1 fallback (result buffer would exceed HPHI_TRACE_BUF_MAX_MB).``
   | ``INFO: ExpecMode 2: two-body Green functions use the ExpecMode-1 fallback (they share their evaluator with three-/four-/six-body Green functions).``
   | The energy/fluctuation family's outcome is decided independently, by
     three possible reasons (it has no shared-evaluator case -- see below),
     reported by its own rank-0 ``INFO`` line:
   | ``INFO: ExpecMode 2: the energy/fluctuation family uses the trace kernel.``
   | ``INFO: ExpecMode 2: the energy/fluctuation family uses the ExpecMode-1 fallback (the Hamiltonian buffer would exceed HPHI_TRACE_BUF_MAX_MB).``
   | ``INFO: ExpecMode 2: the energy/fluctuation family uses the ExpecMode-1 fallback (the Hamiltonian was read from InputHam).``
   | ``INFO: ExpecMode 2: the energy/fluctuation family uses the ExpecMode-1 fallback (unsupported model).``
   | ``S2``, ``NBodyG``, and ``AnomalousG`` remain unconditionally on the
     ``ExpecMode 1`` path in this version, reported by one fixed line:
   | ``INFO: ExpecMode 2: S2, NBodyG, and AnomalousG always use the ExpecMode-1 path in this version.``
   | Supported models for the one-body/two-body Green-function trace
     kernels (as of this phase): ``Hubbard``/``HubbardGC`` and spin-1/2
     ``Spin``/``SpinGC`` (general-spin models (including
     :math:`S \geq 1`), ``tJ``/
     ``tJGC``, ``Kondo``/``KondoGC``, and ``SpinlessFermion``/
     ``SpinlessFermionGC`` are not yet covered and always print the
     "unsupported model" line above for both quantities).
   | The energy/fluctuation family's trace kernel has broader model
     coverage than the Green-function kernels above: it supports
     ``Hubbard``/``HubbardGC``, ``tJ``/``tJGC``, ``Kondo``/``KondoGC``,
     and ``Spin``/``SpinGC`` (including general spin, :math:`S \geq 1`) --
     not just the four rows supported for the Green functions. It does
     NOT support ``SpinlessFermion``/``SpinlessFermionGC`` (their
     particle-number/spin fluctuation semantics differ and are not
     implemented by the kernel); those models, and any unknown model,
     fall back to the ``ExpecMode 1`` path.
   | For the Green functions, on a supported model each quantity's
     outcome is still decided by up to three further runtime reasons,
     evaluated in a fixed order so that exactly one reason applies (they
     are mutually exclusive by construction, not independently checked):
     (a) a shared-evaluator rule
     for the two-body Green function only, checked first -- it shares its
     evaluator with the ThreeBodyG/FourBodyG/SixBodyG (N-body) Green
     functions, so whenever any of those is requested, the two-body
     quantity falls back together with them (this keeps each output
     file's writer unique: otherwise the fallback loop would silently
     drop the N-body output, or the trace kernel would double-write the
     two-body files); the one-body quantity is unaffected by rule (a).
     (b) a no-operators check, checked next -- if this run defines zero
     operators of that kind (e.g. no ``TwoBodyG``/``CisAjtCkuAlvDC``
     entries), the quantity has nothing to stream and falls back; this
     rule fires BEFORE the memory gate below, so an empty operator table
     is never misreported as exceeding the memory cap. (c) a memory gate,
     checked last (and only reached if neither (a) nor (b) already
     applied) -- if the quantity's result buffer would exceed
     ``HPHI_TRACE_BUF_MAX_MB`` (an environment variable giving the cap in
     MiB, an integer in [1, 1048576], default 1024; the cap applies per
     MPI rank and per quantity to the result buffer only; it is parsed
     from the environment on rank 0, only for ``ExpecMode`` 2 runs, and
     broadcast to every rank, so it only needs to be set in the launch
     environment), that quantity falls back to ``ExpecMode 1``.
   | For the energy/fluctuation family, the outcome is decided by exactly
     three possible reasons, checked in this order: (a) an ``InputHam``
     check, checked first -- if this run's Hamiltonian was read from
     ``InputHam`` (``InputHam 1``), the trace kernel cannot rebuild a
     matching Hamiltonian by re-enumerating the model (doing so would
     re-derive the matrix from the model definition instead of reading
     back the one that was actually diagonalized), so the energy family
     unconditionally falls back to ``ExpecMode 1``. (b) an
     unsupported-model check, checked next -- if the model is not one of
     the supported models listed above (i.e. it is
     ``SpinlessFermion``/``SpinlessFermionGC`` or any unknown model), the
     energy family falls back. (c) the same ``HPHI_TRACE_BUF_MAX_MB``
     memory gate as above, checked last (only reached if neither (a) nor
     (b) applies) -- this cap now also bounds the
     energy family's own per-rank Hamiltonian buffer, a compact CSR
     matrix of approximately :math:`24 \times \mathrm{nnz}` bytes (here
     nnz is the number of matrix entries makeHam emits on this rank,
     which the buffer capacity and the gate are sized for; the merged
     entry count ``rowptr[N]`` after summing duplicates can be smaller);
     if the projected CSR size would exceed the cap, the energy family
     falls back to ``ExpecMode 1``. Unlike the Green-function
     kernels, the energy family has no shared-evaluator case (it does not
     share its evaluator with any other always-fallback quantity).
   | If a run reports the memory-gate fallback but you want the trace
     kernel, raise ``HPHI_TRACE_BUF_MAX_MB`` up to the available per-rank
     memory. Note that, unlike the distributed eigenvector storage, the
     CSR Hamiltonian is replicated in full on every rank, so adding MPI
     ranks does NOT shrink it (see :ref:`Sec:ParallelFullDiag` for why);
     raising the cap, or accepting the reported ``ExpecMode 1`` fallback
     (which produces identical results), are the two options.
   | Eligibility: a nonzero ``ExpecMode`` requires ``CalcType`` = 2 (full
     diagonalization) together with ``Solver`` 1 (ScaLAPACK) or 3 (ELPA);
     any other combination (wrong ``CalcType`` or ``Solver``) is rejected
     at startup with an error. With exactly one MPI process, ``ExpecMode``
     is automatically reverted to 0 (results are identical for a single
     process either way), printing (indented in the actual log output)
   | ``INFO: ExpecMode reverts to 0 for a single process (results are identical).``
   | Guarantee: ``ExpecMode`` changes only evaluation speed, never the
     physics -- ``ExpecMode`` 0, 1, and 2 produce identical results up to
     floating-point rounding (summation order differs between kernels, so
     agreement is not bit-identical). This includes the ``var`` column:
     ``var`` is part of the energy/fluctuation family, so it now uses the
     trace kernel together with the rest of that family whenever the
     kernel is active (see above), and the ``ExpecMode 1`` path
     otherwise; either way ``var`` continues to store
     :math:`\langle H^2 \rangle` (downstream code subtracts
     :math:`\langle H \rangle^2` to obtain the variance), an unchanged
     field contract.
   | Memory: ``ExpecMode 1``/``2`` additionally hold a state panel roughly
     the same size as the distributed eigenvector storage used by
     ``Solver 3`` (O(N²/P) per rank); ``ExpecMode 2`` additionally holds,
     per quantity that uses the trace kernel, a buffer capped by
     ``HPHI_TRACE_BUF_MAX_MB`` above -- a result buffer for the one-body/
     two-body Green functions, or the per-rank CSR Hamiltonian buffer
     described above (approximately :math:`24 \times \mathrm{nnz}`
     bytes) for the energy/fluctuation family. ``HPHI_TRACE_BUF_MAX_MB``
     is a per-quantity threshold: it is checked independently against EACH
     quantity's buffer, not as a budget on the sum of the concurrently-live
     trace buffers -- the one-body result buffer, the two-body result
     buffer, and the energy CSR can each be as large as the cap, so their
     live total can exceed it. During the one-time redistribution
     step both the original storage and the new panel coexist, giving a
     temporary peak of roughly 2xO(N²/P) per rank before the original
     storage is freed; steady-state usage afterward is O(N²/P) plus the
     trace buffers, otherwise unchanged from ``ExpecMode 0``.
   | Guidance: prefer ``ExpecMode 1`` or ``2`` for large multi-node
     ``Solver 3`` (ELPA) or ``Solver 1`` (ScaLAPACK) runs with many
     eigenstates and/or many Green-function observables, where
     re-evaluating every eigenstate redundantly on every rank
     (``ExpecMode 0``) becomes the bottleneck. ``ExpecMode 2`` is
     designed to further amortize per-state overhead for
     one-body/two-body Green-function-heavy workloads on the supported
     models above, and for energy/fluctuation-heavy workloads on any
     model (the energy family's trace kernel is not limited to those
     four rows), and is expected to be at least as fast as
     ``ExpecMode 1`` for such workloads. Keep
     the default ``ExpecMode 0`` otherwise, including for small systems
     and single-process runs.
   | Note (behavior fix): as of phase 3a, distributed FullDiag runs
     (``Solver`` 1 or 3, more than one MPI process) compute S2 and Sz on
     rank 0 for ``ExpecMode 0`` as well -- previously these distributed runs
     zero-filled S2 and Sz and printed a shortened stdout progress line
     without the S2 column. The stdout progress line now always matches
     the single-process (serial) format, with the S2 column included, for
     every ``Solver``/``ExpecMode`` combination. This is an intentional
     correctness fix, independent of ``ExpecMode``'s value.

.. raw:: latex

   \newpage
