MPI Parallelization
===================

This section describes the MPI parallelization strategy used in :math:`{\mathcal H}\Phi`.

Overview
--------

:math:`{\mathcal H}\Phi` uses MPI to parallelize calculations across multiple processes.
The Hilbert space basis states are distributed among MPI processes,
and Hamiltonian matrix-vector multiplication :math:`\hat{H}|\psi\rangle` is performed in parallel.

Site Classification
-------------------

When using :math:`N_{\rm proc}` MPI processes for an :math:`N`-site system,
sites are classified into two categories:

**Local sites** (:math:`N_{\rm local}` sites)
  Sites whose quantum states are stored within each process.
  For these sites, basis state manipulations can be performed locally.

**Inter-process sites** (:math:`N_{\rm inter}` sites)
  Sites whose quantum states span multiple processes.
  Operations on these sites require MPI communication.

The relationship is:

.. math::
   N_{\rm local} = N - \log_2(N_{\rm proc})

.. math::
   N_{\rm inter} = \log_2(N_{\rm proc})

For example, with 4 MPI processes (:math:`N_{\rm proc}=4`) on an 8-site system,
we have :math:`N_{\rm local}=6` local sites and :math:`N_{\rm inter}=2` inter-process sites.

MPI Communication Patterns
--------------------------

Operations in the Hamiltonian are classified by the MPI communication they require:

**Local operations**
  Both sites involved in the operation are local sites.
  No MPI communication is needed.

**MPIsingle operations**
  One site is local, one site is inter-process.
  Requires MPI communication with one partner process.
  The communication partner (origin) is determined by:

  .. math::
     \text{origin} = \text{myrank} \oplus T_{\rm pow}[\text{site}]

  where :math:`\oplus` denotes XOR operation and :math:`T_{\rm pow}` is the power of 2 for each site.

**MPIdouble operations**
  Both sites are inter-process sites.
  Requires MPI communication, but the received data can be processed uniformly.
  The communication partner is:

  .. math::
     \text{origin} = \text{myrank} \oplus (T_{\rm pow}[\text{site}_1] + T_{\rm pow}[\text{site}_2])

Batched MPI Communication
-------------------------

To reduce MPI communication overhead, :math:`{\mathcal H}\Phi` groups multiple operations
that share the same communication partner (origin) into batches:

**Before optimization:**
  Each transfer/interaction term with inter-process sites calls ``MPI_Sendrecv`` individually.
  If :math:`N` terms share the same origin, :math:`N` separate communications occur.

**After optimization:**
  Terms with the same origin are grouped together.
  A single ``MPI_Sendrecv`` retrieves the required data,
  then all :math:`N` terms are processed locally.

This optimization is particularly effective for:

* 2D/3D lattice models where many hopping terms share the same communication partner
* Systems with many transfer terms between neighboring sites

The batched communication is implemented for:

* SpinlessFermion models (MPIsingle)
* Hubbard models (MPIsingle, MPIdouble; InterAll for HubbardGC only)
* Spin models (Exchange MPIsingle; SpinGC also batches PairLift)

In MPI builds, batched communication is enabled by default. No input-file
option is required. The optimization changes the communication/apply path for
supported inter-process terms, but it is designed to reproduce the same
numerical result as the conventional per-term MPI path.

.. note::

   **Limitation in Time-Evolution (TimeEvolution) mode**

   When using ``TETwoBody`` or step-dependent ``TEOneBody``
   (expert mode with ``NTEInterAllMax > 0`` or ``NTETransferMax > 0``),
   the MPI batched communication is automatically disabled and falls
   back to conventional per-term MPI communication.

   This is because the batched group structure is initialized once at
   program startup and cannot track interaction terms that are
   added or modified by ``MakeTEDTransfer``/``MakeTEDInterAll`` at each
   time-evolution step.

   Peierls substitution (AC Laser mode, ``PumpType = "AC Laser"``)
   only updates coefficients of existing transfer entries and can
   therefore benefit from batched communication.

Disabling batched communication
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Set the environment variable ``HPHI_MPI_NOBATCH=1`` to disable MPI
communication batching globally and fall back to the original per-term MPI
exchanges. This is intended for debugging and for verifying that the batched
optimization reproduces the unbatched result. It can also be used to compare
performance with the conventional path on a particular machine and input::

   export HPHI_MPI_NOBATCH=1
   mpirun -np 4 ./HPhi -e namelist.def

The active setting is reported at startup as ``MPI batching : ON`` or
``MPI batching : OFF``. The value is read on rank 0 and broadcast to all ranks,
so all processes always agree.

Distributed symmetry temporary-memory policy
--------------------------------------------

The distributed symmetry-basis path estimates the per-rank temporary-memory
peak for basis distribution, representative-directory batches, and
matrix-vector plan construction. By default, an estimate above 1 GiB per rank
prints a warning but does not stop the calculation.

The following environment variables accept unsigned decimal byte values:

* ``HPHI_SYMMETRY_MEMORY_WARN_BYTES`` sets the warning threshold. Its default
  is ``1073741824`` (1 GiB); ``0`` disables the warning.
* ``HPHI_SYMMETRY_MEMORY_LIMIT_BYTES`` sets an optional hard limit. Its default
  is ``0`` (unlimited). A nonzero value makes the calculation fail
  collectively before a tracked temporary-memory peak is allowed to exceed
  the limit.

For example, the following keeps the default warning threshold and sets a
4 GiB per-rank hard limit::

   export HPHI_SYMMETRY_MEMORY_LIMIT_BYTES=4294967296
   mpirun -np 4 ./HPhi -e namelist.def

Both values are read on rank 0 and broadcast to all ranks. They apply
consistently to the distribution, directory, and matrix-vector plan phases.
Allocation failures, integer-overflow checks, and MPI failures remain fatal
regardless of these settings.

SpinlessFermion Off-diagonal Two-body Green's Function
------------------------------------------------------

For SpinlessFermion,
:math:`\langle c^\dagger_i c_j c^\dagger_k c_l \rangle`
with inter-process sites is supported in MPI paths for both
grand-canonical (GC) and canonical calculations.

This path applies the operator sequence on a global-bit representation and
keeps fermion-sign evaluation consistent across rank communication.

Behavior of MPI-required Tests (ctest)
--------------------------------------

The following tests require an MPI runtime environment:

* ``mpi_consistency_*`` (typically ``-np >= 2``)
* ``lanczos_*_mpidouble`` (typically ``-np = 16``)

If ``MPIRUN`` is not set, or the MPI rank requirement is not met,
these tests are reported as ``Skipped`` by ``ctest``.

Example:

.. code-block:: bash

   MPIRUN='mpirun -np 4 --oversubscribe' ctest -L consistency
   MPIRUN='mpirun -np 16 --oversubscribe' ctest -L batching

Fermion Sign
------------

For fermionic systems (Hubbard, SpinlessFermion), the fermion sign arising from
operator anticommutation must be carefully handled during MPI communication.

The sign is computed using the ``SgnBit`` function, which counts the number of
occupied states between the creation and annihilation operators.
This sign is pre-computed during initialization and stored with each transfer group.

Process Number Requirements
---------------------------

The number of MPI processes must satisfy specific constraints:

**Hubbard/Kondo models:**
  Process number = :math:`4^n` (due to 4 states per site: empty, up, down, double)

**Spin-1/2 models:**
  Process number = :math:`2^n` (due to 2 states per site: up, down)

**General spin models:**
  Process number = product of :math:`(2S_i+1)` for inter-process sites

See :ref:`Subsec:CreatingExpert` for detailed instructions on setting process numbers.
