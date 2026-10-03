.. highlight:: none

.. _Subsec:TransSym:

TransSym file
-------------

This file restricts the calculation to one symmetry sector. It defines a
group :math:`G` of site permutations and a one-dimensional character
:math:`\chi(g)` (a complex number of unit modulus) for every element
:math:`g`. HPhi then works in the basis of symmetry-adapted states

.. math::

   |r;\chi\rangle \propto \sum_{g\in G}\chi(g)^{*}\,T_g|r\rangle ,

where :math:`T_g` moves the content of each site :math:`i` to the site
:math:`g(i)` and :math:`|r\rangle` runs over the representative
configurations. Every state of the sector satisfies
:math:`T_g|\psi\rangle=\chi(g)|\psi\rangle`, and the dimension of the
Hilbert space is reduced by roughly the order of the group. The sector
dimension is printed in the log as ``Symmetry basis: raw_dim=... sector_dim=...``.

The group is not restricted to translations. Reflections, rotations, and
their products with translations can be given in the same way, as long as
the Hamiltonian is invariant under every operation and the representation
is one-dimensional (abelian groups, or a one-dimensional irreducible
representation of a larger group). For a two-dimensional irreducible
representation, restrict the group to an abelian subgroup.

Standard mode writes this file as ``qptransidx.def`` when ``MomentumIndex``
is given (see :doc:`the MomentumIndex parameter <../standardmode_en/Parameters_for_conserved_quantities_en>`).
For a chain of length :math:`L` it contains the :math:`L` translations with
:math:`\chi(T^{g})=\exp(-2\pi i m g/L)`, which the ``MomentumIndex`` section
labels as the momentum :math:`k=2\pi m/L`.

An example of the file format (translations of a 6-site ring with
``MomentumIndex = 1``) is as follows.

::

    # MomentumIndex 1
    =============================================
    NQPTrans          6
    =============================================
    ======== TrIdx_TrWeight_and_TrIdx_i_xi ======
    =============================================
    0  1.000000000000000  0.000000000000000
    1  0.500000000000000 -0.866025403784439
    2 -0.500000000000000 -0.866025403784439
    3 -1.000000000000000  0.000000000000000
    4 -0.500000000000000  0.866025403784439
    5  0.500000000000000  0.866025403784439
    0 0 0 1
    0 1 1 1
    ...
    1 0 1 1
    1 1 2 1
    ...
    5 5 4 1

The following example is the reflection :math:`i\to 5-i` of the same ring,
with the character :math:`-1` (odd parity).

::

    =============================================
    NQPTrans          2
    =============================================
    ======== TrIdx_TrWeight_and_TrIdx_i_xi ======
    =============================================
    0  1.0
    1 -1.0
    0 0 0 1
    0 1 1 1
    0 2 2 1
    0 3 3 1
    0 4 4 1
    0 5 5 1
    1 0 5 1
    1 1 4 1
    1 2 3 1
    1 3 2 1
    1 4 1 1
    1 5 0 1

File format
~~~~~~~~~~~

*  Line 1: Header

*  Line 2: [string01] [int01]

*  Lines 3-5: Header

*  Next [int01] lines: [int02] [double01] [double02]

*  Next [int01] :math:`\times` ``Nsite`` lines: [int02] [int03] [int04] [int05]

Lines that start with ``#`` and empty lines are skipped wherever they
appear and do not count in the line numbers above.

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String (a blank parameter is not allowed)

   **Description :** A keyword for the number of symmetry operations.
   Specify ``NQPTrans`` (case-insensitive).

*  [int01]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** The number of symmetry operations, i.e., the order of
   the group.

*  [int02]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the index of a symmetry operation
   (:math:`0<=` [int02] :math:`<` [int01]).

*  [double01], [double02]

   **Type :** Double ([double02] can be omitted)

   **Description :** The real and imaginary parts of the character
   :math:`\chi(g)` of the operation [int02]. When [double02] is omitted,
   the imaginary part is zero.

*  [int03]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** A site index (:math:`0<=` [int03] :math:`<` ``Nsite``).

*  [int04]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** The site to which the operation [int02] moves the
   site [int03] (:math:`0<=` [int04] :math:`<` ``Nsite``).

*  [int05]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** Reserved for anti-periodic boundary conditions.
   In this version it must be ``1``.

Metadata
~~~~~~~~

A comment line of the form ``# MomentumIndex`` [int06] records the
``MomentumIndex`` that generated the file. [int06] must be a non-negative
integer representable by the C ``int`` type (at most ``INT_MAX``).
Standard mode writes it as the first line, and HPhi reports it in
the log as ``TransSym metadata: MomentumIndex=``\ [int06]. The value is not
used by the calculation in this version. Other comment lines are ignored.

Use rules
~~~~~~~~~

*  Headers cannot be omitted.

*  Every operation must be a bijection of the sites, and every pair of
   [int02] and [int03] must be given exactly once.

*  The operations must form a group: they must contain the identity,
   and the composition of any two operations must be one of them.
   The characters must have unit modulus and be multiplicative,
   :math:`\chi(gh)=\chi(g)\chi(h)`, with :math:`\chi(e)=1` for the identity.

*  The Hamiltonian must be invariant under every operation. Each term is
   mapped by the operation and compared with the original terms; a
   mismatch terminates the program with
   ``TransSym Hamiltonian invariance failed``.

*  For ``SpinlessFermion``, ``Hubbard``, and ``tJ``, the sign of the fermion
   permutation is taken into account automatically. For example, a
   reflection that exchanges occupied orbitals contributes a factor
   :math:`-1`, so the dimensions of the even and odd sectors differ from
   the counting for spins.

*  The ``Lanczos``, ``CG``, and ``TPQ`` (microcanonical TPQ) methods support ``Spin`` with
   :math:`S=1/2` and fixed ``2Sz``, ``SpinlessFermion`` with fixed ``Ncond``,
   and ``Hubbard`` / ``tJ`` with fixed ``Nup`` and ``Ndown``.
   Expert-mode Hamiltonian terms are:

   * ``Spin``: longitudinal ``Trans`` (local diagonal fields), ``Exchange``,
     ``Ising``, ``CoulombInter``, ``Hund``, and fixed-Sz ``InterAll``.
   * ``SpinlessFermion``: ``Trans`` (including on-site potentials),
     ``CoulombInter``, and ``InterAll``.
   * ``Hubbard``: spin-conserving ``Trans``, ``CoulombIntra``, ``CoulombInter``,
     ``Hund``, ``Ising``, ``Exchange``, ``PairHop``, and fixed-spin ``InterAll``.
   * ``tJ``: the same terms as ``Hubbard``, projected onto configurations
     without double occupancy. ``CoulombIntra`` and ``PairHop`` then vanish.
     The raw dimension is
     :math:`\binom{N_{\rm site}}{N_\uparrow}\binom{N_{\rm site}-N_\uparrow}{N_\downarrow}`.
     Both replicated and distributed basis layouts are supported, including
     MPI sizes that cannot be used with raw site decomposition.

   Extended terms are combined after fermionic normal ordering or local Spin
   matrix-unit reduction. This accounts for permutation signs, contractions,
   duplicate terms, and cancellations between families before checking
   invariance and conserved quantum numbers. The coefficient tolerance is
   :math:`10^{-10}`. ``PairLift``, ``NBodyInterAll``, and anomalous terms remain
   unsupported. The raw spinless solver still rejects off-diagonal ``InterAll``;
   this extension applies to ``TransSym``. Standard-mode generation is unchanged.

   Correlation functions, spectrum calculations, restart, and the input and
   output of Hamiltonians and eigenvectors are not supported together with
   this file. Unsupported combinations terminate with an error.

Sector TPQ
~~~~~~~~~~

In expert mode, ``CalcType=1`` with ``TransSym`` evolves a random vector within
the selected symmetry sector using :math:`l-H_q/N_{\rm site}`. The default
basis layout is distributed; ``HPHI_SYMMETRY_BASIS_LAYOUT=replicated`` selects
the reference layout. Empty MPI ranks participate in global normalization.
``Lanczos_max``, ``NumAve``, ``LargeValue``, ``initial_iv``, and ``InitialVecType``
have their usual TPQ meanings. ``exct`` does not restrict the TPQ sector.
Standard mode still limits ``MomentumIndex`` generation to Lanczos and CG.

All existing SS/Norm/Flct columns are supported, including doublon and its
second moment for Hubbard. The SS ``phys_var`` column retains its existing
meaning :math:`\langle H^2\rangle`, not the subtracted variance.
``OutputGreenFormat=1`` selects the usual aggregate SS/Norm/Flct files even
though correlation functions remain unsupported. Restart, vector I/O,
spectrum, and cTPQ remain unavailable with ``TransSym``.

These results estimate the trace within a **single symmetry sector**; they
are not a thermal average over the entire fixed-quantum-number space.
The manifest records ``ensemble=single_symmetry_sector``, ``num_ave``,
``large_value``, and ``initial_vec_type``. The existing random generator
depends on MPI ownership and OpenMP threads, so a fixed seed alone does not
make samples identical across different process/thread counts.

Sector manifest
~~~~~~~~~~~~~~~

After constructing a nonempty symmetry basis and validating the sector options,
HPhi writes ``output/symmetry_sector.dat`` before starting the solver. When
``OutputDataHead=1``, the name is
``output/<CDataFileHead>_symmetry_sector.dat``. Ordinary runs without ``TransSym``
and definition-file generation with ``-sdry`` do not write this file.
An output error terminates the calculation on all MPI ranks.

The first line is ``format=HPhiSymmetrySector version=1``. Subsequent lines have
the form ``key=value`` and record the method, model, site count, fixed quantum
numbers, full canonical dimension (``full_dim``), sector dimension
(``sector_dim``), group order, optional ``momentum_index``, basis layout,
MPI ranks, OpenMP thread limit, term counts, and solver parameters.
The manifest describes the input sector; its presence does not certify solver
completion or convergence.

Three versioned fingerprints are included:

* ``group_digest``: ``hphi-group-fnv1a64-v1`` hashes the permutations and
  characters after sorting operations lexicographically by permutation.
  Renumbering operations does not change it. Character components are rounded
  to integer multiples of :math:`10^{-10}`; this quantization is not a general
  equivalence test for floating-point inputs near a rounding boundary.
* ``sector_digest``: ``hphi-sector-multiset-v1:count:xor:sum`` hashes each
  representative state, orbit size, and stabilizer size with FNV-1a 64 and
  combines the hashes without dependence on entry order or MPI ownership.
  Integers use fixed-width little-endian encoding. The sum is modulo
  :math:`2^{64}`. Norms and Hamiltonian diagonal values are excluded.
* ``hamiltonian_digest``: ``hphi-parsed-hamiltonian-fnv1a64-v2`` records the
  supported parsed Hamiltonian terms in their stored order, using the exact
  binary64 coefficient bits. Version 2 adds on-site potentials, pair hopping,
  and the split diagonal/off-diagonal InterAll arrays. Equivalent Hamiltonians expressed in different
  term orders or decompositions may have different fingerprints.

Sector identification uses the model, fixed quantum numbers, sector dimension,
``group_digest``, and ``sector_digest`` together. Conjugate representations may
share the same ``sector_digest``, so it must not be used alone. Hamiltonian
changes do not change the sector identity. These fingerprints are diagnostics,
not collision-free proofs of equivalence, and the order-independent sector
fingerprint alone does not validate vector-component ordering for checkpoint
input. Existing calculation outputs and ``CalcTimerRankStats.dat`` are unchanged.

.. raw:: latex

   \newpage
