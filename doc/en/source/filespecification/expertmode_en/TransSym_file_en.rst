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

*  The ``Lanczos``, ``CG``, ``TPQ`` (microcanonical TPQ), ``cTPQ``, ``FullDiag``, and ``TimeEvolution`` methods support ``Spin`` with
   :math:`S=1/2` and fixed ``2Sz``, ``SpinlessFermion`` with fixed ``Ncond``,
   and ``Hubbard`` / ``tJ`` with fixed ``Nup`` and ``Ndown``.
   Spin-one-half ``SpinGC`` without fixed Sz is also supported as described below.
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
   :math:`10^{-10}`. For these canonical models, ``PairLift``, ``NBodyInterAll``, and anomalous
   terms remain unsupported. SpinGC supports ``PairLift`` as described below. The raw spinless solver still rejects off-diagonal ``InterAll``;
   this extension applies to ``TransSym``. Standard-mode generation is unchanged.

   Correlation functions (``OneBodyG``, ``TwoBodyG``, ``ThreeBodyG``,
   ``FourBodyG``, ``SixBodyG``, and ``NBodyG``) are computed in the sector for
   ``Lanczos``, ``CG``, ``TPQ``, ``cTPQ``, and ``TimeEvolution``. The values are
   the exact expectation values in the sector; rows that are mapped onto each
   other by a symmetry operation have the same value, and an operator that
   leaves the sector gives zero. For canonical ``Spin``,
   ``ThreeBodyG``/``FourBodyG``/``SixBodyG`` are rejected at startup as in the
   raw basis; ``NBodyG`` expresses the same products. ``AnomalousG``, spectrum
   calculations, restart, and the input and output of Hamiltonians are not
   supported together with this file. Eigenvector I/O is available for CG and
   TimeEvolution through the sector checkpoint format below. Unsupported
   combinations terminate with an error.


SpinGC sectors (spin one-half, expert mode)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``CalcModel=4`` with ``TransSym`` selects a spatial-symmetry sector of the
complete spin-one-half space, of dimension :math:`2^{N_{\rm site}}` before
projection. It does **not** fix total :math:`S_z`. Omit ``2Sz``, ``Nup``,
``Ndown`` and ``Ncond`` entirely: even an explicit zero is rejected.
The site count must satisfy ``0 < Nsite < CHAR_BIT * sizeof(unsigned long)``.
Bit 0 denotes down, bit 1 up, and site 0 is the least significant bit.
Basis representatives are the smallest bit strings in their orbits; their
coefficients are positive real in the normalized projected states.
No full-space state vector is required by the sector solver.

Supported Hamiltonian families are on-site ``Trans`` (including transverse
and complex fields), ``Ising``, ``Exchange``, ``CoulombInter``, ``Hund``,
``PairLift``, and ``InterAll`` built from on-site spin matrix units, including
terms that change total :math:`S_z`. Hermiticity and group invariance are
still required. For :math:`E_i^{ab}=|a\rangle_i\langle b|`, one real
``PairLift`` row ``i j J`` means

.. math::

   J(E_i^{10}E_j^{10}+E_i^{01}E_j^{01}).

There is no additional factor of one-half. Reversed and duplicate rows add;
a same-site row vanishes. These rules do not relax the fixed-Sz restriction
of canonical ``Spin``.

.. list-table:: SpinGC sector method and feature support
   :header-rows: 1
   :widths: 17 24 22 19 18

   * - Method (CalcType)
     - Result
     - Correlations
     - Vector import/export
     - MPI layout
   * - Lanczos (0)
     - Low-energy states
     - All six formats
     - No / No
     - Distributed or replicated
   * - mTPQ (1)
     - Single-sector samples
     - All six formats
     - No / No
     - Distributed or replicated
   * - FullDiag (2)
     - All sector eigenvalues
     - No
     - No / No
     - Replicated metadata
   * - CG / LOBCG (3)
     - Low-energy states
     - All six formats
     - Yes / Yes
     - Distributed or replicated
   * - TimeEvolution (4)
     - Static or driven evolution
     - All six formats
     - Required / Optional
     - Distributed or replicated
   * - cTPQ (5)
     - Single-sector samples
     - All six formats
     - No / No
     - Distributed or replicated

The six correlation formats are ``OneBodyG``, ``TwoBodyG``, ``ThreeBodyG``,
``FourBodyG``, ``SixBodyG`` and ``NBodyG``. ``OneBodyG`` requires on-site
operators; off-site rows are rejected. Both aggregate and legacy output
formats are available. FullDiag supports LAPACK (Solver 0, one rank),
ScaLAPACK (1) and ELPA (3), subject to the restrictions below.

:math:`S_z=\sum_i S_i^z` and :math:`S_z^2` are evaluated with the actual
sector vector, not inferred from a fixed quantum number. Existing output
columns are unchanged: the CG energy file reports ``Sz``; TPQ/cTPQ/TE Flct
files report the first and second magnetization moments. CG does not add a
new ``Sz2`` column. mTPQ and cTPQ estimate a **single spatial-symmetry sector**,
not the full SpinGC thermal ensemble. Summing sectors requires their proper
statistical weights; a single sector sample is not such a sum.

CG vector import evaluates supplied states without restarting optimization.
TE import starts a new time grid; it can use a different Hamiltonian within
the same sector (quench). Same rank count, ownership, model, sector and phase
are required. Neither import is solver restart; ``ReStart`` remains rejected.
SpinGC checkpoints retain version 1 with ``model=4`` and zero ``nup/ndown/ne``
header fields. Their Hamiltonian fingerprint is
``hphi-parsed-hamiltonian-fnv1a64-v3``, which includes parsed PairLift rows.
Canonical models retain their v2 fingerprints and existing checkpoint format.
The manifest has ``fixed_quantities=none`` and ``full_dim=2^Nsite``.

``TEOneBody`` and ``TETwoBody`` may drive SpinGC using the right-endpoint
Taylor rule described below. Every used slice is checked before propagation.
``Laser`` is rejected; use invariant on-site ``TEOneBody`` or spin-product
``TETwoBody`` entries instead. General spin, Boost, Kondo, new Standard-mode
SpinGC momentum input, spin-axis rotations, global spin flip, antiunitary
operations, multidimensional irreducible representations, spectrum,
Hamiltonian I/O, solver restart and rank-changing checkpoint redistribution
are unsupported. ``CoulombIntra``, ``PairHop``, ``NBodyInterAll`` and
``AnomalousG`` are not supported in this SpinGC sector path. FullDiag also
rejects eigenvector/correlation output, distributed basis metadata, MAGMA,
and nonserial ``ExpecMode``.

Eight-site transverse-field CG example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In an empty directory, save the following standard-library Python code as
``make_input.py`` and run ``python3 make_input.py``. It writes portable expert
input files for :math:`H=-\sum_{i=0}^7 S_i^x`, translations with character 1
(momentum zero), and CG. The positive ``Trans`` coefficients implement the
minus sign in HPhi's transfer convention.

.. code-block:: python

   from pathlib import Path

   def definition(name, rows, count=None, keyword="NData"):
       rows = list(rows)
       header = "====\n{} {}\n====\n====\n====\n".format(
           keyword, len(rows) if count is None else count)
       Path(name).write_text(header + "".join(
           " ".join(map(str, row)) + "\n" for row in rows))

   Path("sym.def").write_text(
       "CalcMod calc.def\nModPara mod.def\nLocSpin loc.def\n"
       "TransSym group.def\nTrans trans.def\n")
   Path("calc.def").write_text(
       "CalcType 3\nCalcModel 4\nOutputMode 0\nOutputDataHead 1\n")
   Path("mod.def").write_text(
       "====\nModel_Parameters 0\n====\n====\n====\n"
       "CDataFileHead zvo\nCParaFileHead zqp\n====\n"
       "Nsite 8\nLanczos_max 400\ninitial_iv -1\nexct 1\n"
       "LanczosEps 18\nLanczosTarget 1\nLargeValue 100\nPreCG 0\n")
   definition("loc.def", ((i, 1) for i in range(8)))
   definition("trans.def", ((i, a, i, b, 0.5, 0)
              for i in range(8) for a, b in [(1, 0), (0, 1)]))
   definition("group.def", [(g, 1.0, 0.0) for g in range(8)] +
              [(g, i, (i+g) % 8, 1) for g in range(8) for i in range(8)],
              count=8, keyword="NQPTrans")

Run ``HPhi -e sym.def`` (or ``mpiexec -np 4 HPhi -e sym.def``) with
``OMP_NUM_THREADS=1``. The converged ground state has :math:`E=-4`,
:math:`\langle S_z\rangle=0` and :math:`\langle S_z^2\rangle=2`.
``output/zvo_energy.dat`` reports energy and Sz. To reconstruct Sz2 in this
CG example, request :math:`S_i^z S_j^z` via ``TwoBodyG`` and sum over all
:math:`i,j` (each :math:`S_i^z=(E_i^{11}-E_i^{00})/2`).
This example does not add or require a Standard-mode keyword.

Sector Lanczos basis layout
~~~~~~~~~~~~~~~~~~~~~~~~~~~

With ``TransSym``, sector Lanczos runs in the symmetry-reduced basis. The
default basis layout is ``distributed``: each MPI rank holds only its owned
rows and the ghost entries required for those rows. Matrix-vector products use
the existing distributed plan and halo exchange, while vector norms, inner
products, energy, and variance are combined with MPI reductions. When the
number of MPI ranks exceeds the sector dimension, some ranks can own zero rows;
these ranks still participate in the collective operations.

For developer comparison and diagnosis, set
``HPHI_SYMMETRY_BASIS_LAYOUT=replicated`` to select the rollback layout. This
setting is not recommended for normal operation. The layout choice does not
change output filenames or formats, or the sector manifest schema.

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
``OutputGreenFormat=1`` selects the usual aggregate files for SS/Norm/Flct and
for the correlation functions. Correlation functions are recomputed at the TPQ
evaluation steps selected by ``ExpecInterval``; requesting every site pair at
every step can be expensive for large sectors. Restart, vector I/O, and
spectrum remain unavailable with ``TransSym``.

These results estimate the trace within a **single symmetry sector**; they
are not a thermal average over the entire fixed-quantum-number space.
The manifest records ``ensemble=single_symmetry_sector``, ``num_ave``,
``large_value``, and ``initial_vec_type``. The existing random generator
depends on MPI ownership and OpenMP threads, so a fixed seed alone does not
make samples identical across different process/thread counts.

Sector cTPQ
~~~~~~~~~~~

``CalcType=5`` uses the same sector storage and outputs as mTPQ. It applies
:math:`\sum_{n=0}^{n_{\max}}(-\Delta\beta H_q/2)^n/n!` and normalizes after
each step. The Taylor truncation remains the user's convergence parameter.
Without ``InvTemp``, :math:`\Delta\beta=1/\mathrm{LargeValue}`;
``ExpandCoef`` defaults to 10 when omitted. An explicitly supplied order must
be a positive integer less than ``INT_MAX``. Zero or non-finite step norms
and matrix-vector failures terminate the calculation.

An ``InvTemp`` file uses rows ``beta nmax physcal eigen``. Beta must start
at zero and be finite and nondecreasing. Repeated beta values represent a zero step. ``nmax`` is a positive
integer less than ``INT_MAX``, and both flags are 0 or 1. The order on row
``i`` advances from beta ``i`` to beta ``i+1``; the final order is unused.
Sector cTPQ rejects any nonzero ``eigen`` flag, including at the initial or
final point. Vector I/O and restart remain unsupported.
For an explicit schedule, correlation functions are written at the initial
point and at rows whose ``physcal`` flag is 1. Without ``InvTemp``,
``ExpecInterval`` selects the correlation steps. ``OutputGreenFormat=1``
collects all samples and steps into the indexed aggregate correlation files.
The manifest records ``canonical_tpq_steps``, ``beta_schedule`` and either
the uniform step/order or all explicit beta/order rows. As for mTPQ, these
outputs describe a single sector, not the sum over all sectors.

Sector time evolution with a static Hamiltonian
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In expert mode, ``CalcType=4`` evolves a sector checkpoint with a static
Hamiltonian. Set ``InputEigenVec=1``, ``ReStart=0``, a positive integer
``ExpandCoef``, and ``SpectrumVec seed_eigenvec_0`` in the namelist. The
rank files ``output/seed_eigenvec_0_rank_<rank>.dat`` must use the sector
format below. A CG run with ``OutputEigenVec=1`` provides such files.
The basis, sector, MPI size and layout must match; the Hamiltonian may differ
for a quench. The default layout is distributed; replicated remains available
through ``HPHI_SYMMETRY_BASIS_LAYOUT=replicated``.

Before running TE, preserve all rank files of the CG seed under a separate
prefix, or change the TE ``CDataFileHead``. If the input and output prefixes
are identical, row-zero output ``<prefix>_eigenvec_0_rank_<rank>.dat``
overwrites the seed files. The separate sector final filename does not
prevent this row-zero overwrite.

Provide a ``TEOneBody`` time grid whose number of terms is zero at every row.
``Lanczos_max`` selects the number of rows. Times must be finite and
nondecreasing. The first row records the imported state without propagation;
subsequent rows apply the degree-``ExpandCoef`` Taylor polynomial of
:math:`\exp[-i H(t_j-t_{j-1})]`, followed by normalization. Repeated times are
allowed. The input checkpoint's solver step and time are recorded as
provenance, but do not resume its clock. For time-dependent interactions and Peierls driving, see the next section.

``SS``, ``Norm`` and ``Flct`` contain sector expectation values. ``Norm`` is
the norm before each step's normalization; Taylor truncation can make it
differ from one. Increase ``ExpandCoef`` or reduce the time spacing to check
convergence.
Correlation functions are evaluated after propagation at the steps selected
by ``ExpecInterval``. ``OutputGreenFormat=1`` collects them into files indexed
by the zero-based time-grid step. Requesting every site pair at every step can
be expensive for large sectors. ``ReStart`` remains unsupported.
The sector manifest records the actual time grid, Taylor order, and source
checkpoint's method, state, step, time and Hamiltonian digest.

With ``OutputEigenVec=1`` and positive ``OutputInterval``, periodic files are
``<prefix>_eigenvec_<step>_rank_<rank>.dat``; the final file is separately
named ``<prefix>_eigenvec_final_rank_<rank>.dat``. In these sector headers,
``step`` is the completed zero-based grid row and ``time`` is that row's
physical time. ``state_index`` retains the imported state's label.
The final file does not overwrite row zero. Any of these checkpoints can
seed a new same-sector run via ``SpectrumVec``; this is a new time grid,
not a restart. Raw TE binary files use their existing, different format.

Time-dependent sector Hamiltonians
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``TEOneBody`` or ``TETwoBody`` may drive the sector Hamiltonian.
``Laser`` is available only for the canonical models, not SpinGC.
Use one driving family per calculation. One-body and two-body entries are
added to the static Hamiltonian; Peierls driving changes the phases of its
parsed transfer coefficients. Diagonal and off-diagonal terms are supported
for the four canonical models above, including spinless fermions. The raw
spinless solver's diagonal-TE restriction is unchanged.

Before loading the initial vector or writing time-series data, HPhi checks
every used time slice for finite coefficients, conserved quantum numbers,
and invariance under the specified group. A violation at a later time stops
the run before any propagation. ``Laser`` uses the nine existing laser
parameters and times ``Tinit + step * TimeSlice``; finite nonnegative spacing
is required. The first row records the initial state at ``Tinit``. A
``TEOneBody``/``TETwoBody`` file supplies its own finite, nondecreasing grid.

For the interval ending at :math:`t_j`, sector TE freezes the Hamiltonian at
:math:`H(t_j)` and applies its Taylor polynomial to the previous state. All
powers, including the linear term, use this same matrix. This is a
right-endpoint piecewise-constant approximation to a continuously varying
Hamiltonian: increasing Taylor order alone does not remove the time-grid
error. Check both time spacing and Taylor order. This sector convention is
separate from the legacy raw dynamic-TE recurrence.

The current implementation retains the basis representatives, ordering,
phase and MPI ownership. At each time point it updates the diagonal values,
recreates the distributed lookup directory from the existing local basis,
and rebuilds the multiplication plan. This handles terms that appear or
disappear without enumerating the raw Hilbert space again. It does not yet
cache a common sparsity pattern or update only plan coefficients. The manifest records ``te_hamiltonian=time_dependent``,
``te_integrator=right_endpoint_taylor``, ``te_plan_update=rebuild_with_fixed_basis``, and
``te_hamiltonian_<step>`` for each effective parsed Hamiltonian. The same
digest is written into that row's checkpoint header.

Sector vector checkpoints
~~~~~~~~~~~~~~~~~~~~~~~~~

Expert-mode CG accepts ``OutputEigenVec=1`` and ``InputEigenVec=1``.
The files are ``output/<CDataFileHead>_eigenvec_<state>_rank_<rank>.dat``
with zero-based state and rank indices. All ranks write a file, including
empty owners. Reading requires the same model, fixed quantum numbers,
group/character, sector, MPI rank count, basis layout, local ownership,
global-index-to-representative mapping, and basis phase convention.
Legacy raw-basis vector files are rejected. ``InputEigenVec=2`` and
``ReStart`` remain unsupported. A CG input run evaluates the supplied states;
it does not solve for new eigenstates or resume iterations.

The source Hamiltonian digest is recorded separately and may differ from the
current Hamiltonian, allowing a same-sector quench. The load log reports this
change. Such a state import is distinct from a solver restart.
All metadata is checked collectively before loading vector data into the
active vector. Truncation, extra data, non-finite values, a squared global
norm differing from one by more than :math:`10^{-8}`, and checksum mismatches
are rejected. Files from different checkpoint sets cannot be mixed.

The version-1 binary format is portable across endianness. It consists of
28 little-endian unsigned 64-bit header words followed by ``local_dim``
pairs of IEEE binary64 real/imaginary values. The unused vector element zero
is not stored. Header words, in order, are:

.. code-block:: text

   magic version phase scalar model nsite nup ndown ne
   raw_dim sector_dim ranks rank layout offset local_dim group_digest
   sector_count sector_xor sector_sum order_digest hamiltonian_digest
   source_method state_index step time_bits payload_xor payload_sum

``magic`` is the eight bytes ``HPHISV1\n``; version and phase are 1;
scalar is 128; layout is 0 (replicated) or 1 (distributed). ``offset`` is
zero-based. ``source_method`` uses the CalcType number, and ``time_bits``
is an IEEE binary64 physical time (zero for CG). Phase 1 uses the normalized
:math:`\sum_g\overline{\chi(g)}T_g|r\rangle` with the smallest representative
and a positive real coefficient at that representative.

Group, sector and Hamiltonian fingerprints use the manifest algorithms.
The order digest is FNV-1a-64 over each owned entry's one-based global index,
representative integer, orbit size and stabilizer size, each encoded in eight
little-endian bytes. The payload digest starts with rank, offset and local
length (eight bytes each), followed by the serialized complex coefficients.
The header stores both XOR and unsigned-modulo-:math:`2^{64}` sum of those
per-rank hashes; these are integrity fingerprints, not cryptographic hashes.
Each rank writes a temporary ``.part`` file, and publishes it after all writes
have completed successfully. Publication errors fail the run; readers verify
that every rank belongs to the same complete checkpoint set.

Sector FullDiag
~~~~~~~~~~~~~~~

Expert-mode ``CalcType=2`` computes all eigenvalues within the chosen sector.
``Solver=0`` (LAPACK, one MPI rank), ``Solver=1`` (ScaLAPACK), and ``Solver=3``
(ELPA) are supported when compiled in. The basis metadata must be replicated;
``HPHI_SYMMETRY_BASIS_LAYOUT=distributed`` is rejected for FullDiag.
All ranks use the full sector dimension for solver work arrays. ELPA with
multiple ranks builds only its owned column panel. ScaLAPACK currently builds
a replicated Hamiltonian but does not allocate the unused replicated
eigenvector matrix. ELPA rejects a sector smaller than its process grid;
reduce the rank count in that case.

The output is ``output/<CDataFileHead>_energy_sector.dat`` with zero-based
eigenvalue indices and energies, accompanied by the sector manifest.
This filename always has the data-file prefix. Ordinary FullDiag still writes
``Eigenvalue.dat``. Sector FullDiag does not produce eigenstate observables,
correlation files, or eigenvectors. ``ExpecMode`` other than 0 and MAGMA are
unsupported. Hamiltonian/vector I/O and restart remain unavailable.
The manifest records the solver, matrix storage, and ``output_scope=eigenvalues``.

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
* ``hamiltonian_digest``: for canonical models, ``hphi-parsed-hamiltonian-fnv1a64-v2`` records the
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
