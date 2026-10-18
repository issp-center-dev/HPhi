.. _Subsec:restart_v3:

v3 solver restart files
=======================

CG, mTPQ, cTPQ and TimeEvolution use these files for ``ReStart`` continuation.
The saved and current MPI process counts may differ. HPhi reconstructs the
global vector from the saved contiguous intervals and reads their intersections
with the current ownership. A ``TransSym`` v3 continuation may also change
metadata layout 0/1 between ``replicated`` and ``distributed`` without changing
global vector order. Independently published TPQ sets may therefore have
different saved rank counts and layouts; all rank files within one set must
agree. Old sector v1/v2 state imports also permit rank-count and layout changes.

Names, units and manifest
-------------------------

Vectors are named
``<CDataFileHead>_restart_set<k>_comp<c>_gen<g>_rank_<r>.dat``. ``k`` is a CG
state, TPQ sample, or zero for TimeEvolution; component 0 is the state and CG
component 1 is its search direction. The generation is 0 or 1. The publication
unit is all states and both components for CG, one sample for TPQ/cTPQ, and set
0 for TimeEvolution.

``<CDataFileHead>_restart_manifest.dat`` is the sole pointer to published data::

   hphi_restart_manifest 1
   method <CG|TPQ|cTPQ|TimeEvolution>
   basis <raw|sector>
   set <k> ranks <P> gen <g> step <STEP> ncomp <1|2>

CG rows have one rank count, generation and step. TPQ rows are independent.
``ReStart=3`` never modifies these files. A new ``ReStart=1`` publication
replaces unrelated old rows. Writers create every ``.part`` file, rename the
complete unit, then atomically publish the manifest. Only after publication
are superseded files removed. The manifest therefore never names a partial
unit.

Binary header
-------------

Each vector contains 40 little-endian unsigned 64-bit words, followed by
``LOCAL_DIM`` pairs of IEEE binary64 real and imaginary values. Vector index
zero is omitted.

.. list-table:: v3 header words
   :header-rows: 1
   :widths: 7 22 71

   * - Index
     - Name
     - Meaning
   * - 0
     - MAGIC
     - ``0x0a31565349485048`` (``HPHISV1\n``)
   * - 1
     - VERSION
     - ``3``
   * - 2
     - PHASE
     - Basis phase convention; checked even for import
   * - 3
     - SCALAR
     - ``128`` for binary64 complex payloads
   * - 4
     - MODEL
     - Normalized input model identity
   * - 5
     - NSITE
     - Number of sites in the saved basis
   * - 6
     - NUP
     - Fixed up count, otherwise zero
   * - 7
     - NDOWN
     - Fixed down count, otherwise zero
   * - 8
     - NE
     - Fixed particle count, otherwise zero
   * - 9
     - FULL_DIM
     - Sector full dimension; zero for raw
   * - 10
     - DIM
     - Global vector dimension
   * - 11
     - RANKS
     - Saved MPI process count; equals the manifest row
   * - 12
     - RANK
     - Saved rank index
   * - 13
     - LAYOUT
     - 2 for raw; 0 replicated or 1 distributed for sector
   * - 14
     - OFFSET
     - Zero-based global start of this payload
   * - 15
     - LOCAL_DIM
     - Number of payload entries
   * - 16
     - GROUP
     - Sector group digest; zero for raw
   * - 17
     - SECTOR_COUNT
     - Partition-independent sector count
   * - 18
     - SECTOR_XOR
     - Partition-independent sector XOR digest
   * - 19
     - SECTOR_SUM
     - Partition-independent sector sum digest
   * - 20
     - ORDER
     - Sector global-order XOR digest
   * - 21
     - HAMILTONIAN
     - Informational Hamiltonian digest; zero/unknown for raw
   * - 22
     - METHOD
     - Source ``CalcType``
   * - 23
     - STATE
     - Set index; equals the file name and manifest
   * - 24
     - STEP
     - Next iteration/step, or terminal CG row as selected by FLAGS
   * - 25
     - TIME
     - TimeEvolution physical-time bits; zero otherwise
   * - 26
     - DATA_XOR
     - Partition-independent payload XOR digest
   * - 27
     - DATA_SUM
     - Partition-independent payload sum digest
   * - 28
     - LOCAL_SITE_MASK
     - Kondo localized-site mask; zero otherwise
   * - 29
     - FIXED_FLAGS
     - Raw fixed quantities or sector Kondo flags
   * - 30
     - BASIS_KIND
     - 0 raw, 1 sector
   * - 31
     - COMPONENT
     - 0 state, 1 CG direction
   * - 32
     - TOTAL2SZ
     - Fixed raw ``2Sz`` bit representation, otherwise zero
   * - 33
     - SITE_DIGEST
     - Raw site/local-spin identity digest
   * - 34
     - ORDER_SUM
     - Sector global-order sum digest
   * - 35
     - NSET
     - CG state count or source TPQ sample count; 1 for TE
   * - 36
     - SCALAR_A
     - CG eigenvalue or TPQ first norm
   * - 37
     - SCALAR_B
     - cTPQ accumulated inverse temperature; zero otherwise
   * - 38
     - SCALAR_C
     - cTPQ ``delta_tau``; zero otherwise
   * - 39
     - FLAGS
     - Bit 0 marks a terminal converged CG checkpoint; other bits are invalid

Comparison and digests
----------------------

Across ranks of one set/component, only ``RANK``, ``OFFSET`` and
``LOCAL_DIM`` differ. CG state and direction components may additionally
differ in ``COMPONENT`` and data digests. CG states may differ in ``STATE``,
data digests and ``SCALAR_A``; all share step, flags, set count, generation
and rank count. TPQ samples may have independent state, step, data/scalars
and manifest rows, but their basis identity and method agree.

The data digest hashes each global index and the exact real/imaginary bits,
then reduces both XOR and unsigned sum across ranks. Sector order digests hash
the global index, representative, orbit size and stabilizer size. The pair of
reductions detects corruption and ordering changes without depending on the
partition.

Rank redistribution
-------------------

v3 and sector v1/v2 files carry their saved rank count, offsets and local
dimensions. Keep every file from rank 0 through the saved count minus one,
including files with zero payload. HPhi validates the complete saved set and
then reads each overlap between a saved interval and the current interval.
``RestartSourceRanks`` is not required for these self-describing formats. v3
continuation ignores it and uses the count in each manifest row. A v3 or
sector v1/v2 state import treats a positive value as an optional equality
check against the saved count.

For sector v1/v2, HPhi first reproduces the old partition-dependent DATA
checksum from the saved files. Only after that check succeeds is the payload
redistributed. The partition-independent sector digest verifies the basis
set; the canonical representative-state sort fixes its global order. Source
and destination layouts may be either 0 or 1, but every source file in one
checkpoint must declare the same valid layout.

Legacy raw files have only a native 12-byte ``STEP``/local-dimension header.
They do not record the saved rank count or identity. A rank-changing legacy
read therefore requires ``RestartSourceRanks=<saved count>``. HPhi requires
all declared files, verifies their headers, lengths, common step and summed
dimension, and reports that only dimensions—not model or sector identity—were
checked.

The dedicated cross-rank regression coverage is:

.. list-table:: Rank-changing restart/import regression matrix
   :header-rows: 1
   :widths: 18 18 20 44

   * - Format
     - Purpose
     - Basis/method
     - Tested source to destination
   * - v3
     - continuation
     - raw CG
     - 4 to 1, 4 to 16, 16 to 4, and 1 to 4 ranks
   * - v3
     - continuation
     - raw mTPQ/cTPQ/TE
     - mTPQ 4 to 2, cTPQ 4 to 2, and TE 4 to 1 ranks
   * - v3
     - continuation
     - sector CG
     - 3 distributed to 6 replicated and 16 distributed ranks, including empty owners
   * - v3
     - state import
     - raw CG/TE
     - 4 to 1 rank
   * - sector v1/v2
     - state import
     - Spin/Kondo
     - Spin v1 4 to 1/2/16 and Kondo v2 4 to 1 rank; layout changes and checksum corruption are covered
   * - legacy raw
     - state import
     - CG
     - 4 to 1 and 4 to 16 ranks, including rejection of a missing empty source rank
   * - legacy raw
     - continuation
     - mTPQ
     - 4 to 2 ranks

There is not yet a dedicated rank-changing numerical regression for sector v3
mTPQ/cTPQ/TimeEvolution continuation or sector v3 state import, nor for legacy
raw CG/cTPQ/TimeEvolution continuation or TE state import. These paths use the
same interval reader, but the absence of a dedicated case should not be read as
additional numerical validation.

Recovery and legacy selection
-----------------------------

For every required set, HPhi first uses a manifest row. A missing or corrupt
referenced v3 file is fatal and never falls back to legacy data. Without a
row, complete legacy files are validated and used. ``ReStart=2`` starts a new
set only if neither source exists; ``ReStart=3`` reports an error. v1/v2
sector files are import-only. Legacy CG restarts additional iterations with
labels starting at one and do not recover the old CG history prefix. A legacy
rank change requires ``RestartSourceRanks`` as described above.

History recovery preserves rows before ``STEP``, removes later rows and an
incomplete EOF fragment, and rejects corrupt interior text. Terminal CG keeps
its already emitted terminal row only when it still satisfies the current
convergence test. Aggregate TPQ rows are appended in execution order; compare
multiple resumed sets after sorting by ``(set, step, observable)``.

Using a v3 vector as an imported state
--------------------------------------

Read the selected set and generation from the manifest. Copy every
``..._restart_set<k>_comp0_gen<g>_rank_<r>.dat`` file to the existing import
name, ``<SpectrumVec>_rank_<r>.dat`` for TimeEvolution or
``<CDataFileHead>_eigenvec_<state>_rank_<r>.dat`` for CG. Keep all source rank
files and set ``RestartSourceRanks`` to that count when specifying it. This
parameter is optional for v3, but must match the saved count when present.
The current rank count and, for a sector, layout may differ. Import starts a
new solver clock; it is not continuation.
