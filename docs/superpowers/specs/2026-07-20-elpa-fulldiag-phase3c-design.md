# ELPA FullDiag Phase 3c: CSR-based trace kernel for the energy family

Status: v3 — revised after design-review rounds 1 (10 must_fix) and 2
(5 must_fix)
Depends on: phase 3b (`2026-07-11-elpa-fulldiag-phase3-design.md` §3), the
phase-3b migration note, and the merged `feature/elpa-fulldiag` state as of
0f15aabf (develop merged, including the momentum/symmetry work).

## 1. Goal and scope

Remove the dominant remaining `ExpecMode 2` cost: the energy-family
evaluator (`expec_energy_flct()`), which still runs per state through the
Mode-1 fallback path and was measured to account for most of the Mode-2
observable wall time once the GF kernels landed (clavius N=4900: GF kernels
~1.1 s of a 26 s wall; kugui N=16384: `CalcPhys` 20.8 s, GF kernels ~2 s).

In scope:

* A per-rank sparse (CSR) representation of the full Hamiltonian H,
  collected once per run by re-running the audited makeHam element
  enumeration through a collector sink (two-pass count/fill; §3a).
* A streaming energy-family kernel: for each owned eigenstate x,
  y = H x (no MPI), `energy` = Re(x†y), `Phys.var` = y†y — **the existing
  field contract is preserved: `Phys.var` stores ⟨H²⟩, not the variance**
  (downstream consumers such as CalcByLanczos subtract the squared energy
  themselves; FullDiag's output path is unchanged). The fluctuation
  columns (doublon, num, Sz and their squares, per the model's existing
  column set) come from precomputed per-basis diagonal coefficient arrays.
* Plan integration as a new quantity slot `TRACE_Q_ENERGY` with an
  explicit build→finalize lifecycle (§3c) and the exclusive demotion
  taxonomy extended by two energy-specific reasons (memory gate;
  input-Hamiltonian source).
* Model coverage: **every model reachable through makeHam** (the
  Hubbard/tJ/Kondo group, HubbardGC, Spin, SpinGC — including general
  spin). This is broader than the GF kernels (no half-spin restriction and
  no per-model probe work), but it is NOT "every HPhi model": makeHam has
  no SpinlessFermion(GC) branch, so those models cannot reach FullDiag's
  matrix generation in the first place; the kernel claims exactly
  makeHam's coverage, nothing more.

Out of scope (unchanged Mode-1 behavior): S² (`expec_totalspin`), NBodyG,
AnomalousG, GF-kernel model coverage, tiled GF buffers, CalcTimer
integration (stdout timing lines only, as in 3b).

## 2. Inherited invariants

* **`Phys.var` semantics**: the field stores ⟨H²⟩ per state (see above);
  the kernel computes it genuinely as y†y — never substituted from the
  eigenvalue. Agreement with Modes 0/1 at the equivalence-test tolerance.
  No clamping: if downstream arithmetic produces a small negative
  variance from rounding, that behavior is identical to today's because
  the subtraction happens in unchanged downstream code.
* **Single writer**: when the energy kernel is active it produces ALL
  energy-family Phys fields for its states and `expec_energy_flct()` is
  skipped for them. The family is never split.
* **Single agreed plan**: all consumers observe one final plan (see the
  build→finalize lifecycle in §3c, which replaces 3b's single-call
  immutability with a two-phase construction that is immutable from
  finalize onward).
* **Guards**: the collector inherits the AddHamElem validity guards
  (`tmp_off > 0`, `dmv != 0.0`) unchanged — it hooks the same macro.
* **ExpecMode 0/1 byte-identity**: the all-fallback short-circuit is
  unchanged; no behavior change outside `ExpecMode 2`.

## 3. Design

### 3a. CSR collection (`src/expec_trace_ham.c`, new TU)

**Sink abstraction.** `hamstore.h` gains an explicit sink mode
`enum HamSink { HAM_SINK_DENSE_REPLICATED, HAM_SINK_DENSE_PANEL,
HAM_SINK_TRACE_COLLECT }` (a single global `iHamSinkMode`, superseding the
boolean `iHamPanelActive` checks inside the macro; the two dense modes map
onto today's branches so their behavior is untouched). makeHam is
refactored so that BOTH of the following are selected by the sink mode,
not by `iHamPanelActive` alone:

1. **Storage initialization**: the dense clearing loops (replicated
   `Ham[i][j]=0` / panel memset) are SKIPPED in collect mode — this is the
   round-1 must-fix: after ELPA freed `Ham_local`, re-running makeHam with
   the replicated flag would dereference the never-allocated `Ham`.
2. **Column range**: collect mode selects the FULL range (1..i_max) on
   every rank, independent of panel ownership.

makeHam remains non-reentrant (it already is, via its globals); collect
mode does not change that and is documented on the enum.

**Two-pass count/fill (replaces round-1's growable triplets).**

* **Pass 1 (count)**: run the enumeration with a counting sink — no
  triplet storage, just per-row counts in a `(N+1)`-sized array (this
  array is the future `rowptr`; see the lifetime graph). Its own
  allocation happens BEFORE the gate and its failure takes the same
  synchronized-demotion path as a gate failure (`local_ok = 0`).
* **Allocation lifetime graph and gate**: every size below is computed
  with `sizeof`-based checked arithmetic (3b `uintmax_t` guard style —
  overflow checked before every multiplication/conversion; no hard-coded
  byte constants):
  | array | size | lifetime |
  |---|---|---|
  | counts→rowptr | `(N+1)·sizeof(long int)` | pass 1 → kernel teardown (in-place prefix sum turns counts into rowptr) |
  | fill cursors | `(N+1)·sizeof(long int)` | pass 2 only (freed after merge) |
  | colidx | `nnz_raw·sizeof(long int)` | pass 2 → teardown |
  | val | `nnz_raw·sizeof(double complex)` | pass 2 → teardown |
  | y (kernel) | `N·sizeof(double complex)` | streaming → teardown |
  | diagonal coeff. arrays | `n_diag·N·sizeof(double)` | fill once → teardown |
  | sort workspace | `k_max·(sizeof(long int)+sizeof(double complex))`, `k_max` = max raw row length (known after pass 1) | pass 2 only (freed after compaction) |
  The gated peak is the SUM of all rows (they coexist during pass 2 /
  streaming; nothing reallocs, nothing shrinks — after merge
  `nnz_merged ≤ nnz_raw` simply leaves unused capacity). The gate
  compares this sum against the `HPHI_TRACE_BUF_MAX_MB` per-quantity cap
  (rank-0-parsed + Bcast, as today). The collector API takes the cap as
  an explicit parameter (testability; the env path supplies it in
  production).
* **Pass 2 (fill)**: entries go directly into their row segment via the
  cursors, then each row is sorted by column index with a STABLE in-row
  sort (hybrid: insertion sort below a small threshold, stable merge
  sort above it using the ONE reusable sort workspace from the lifetime
  graph — sized to the maximum raw row length, allocated with the other
  arrays under the gate, its failure taking the same synchronized
  demotion path — deterministic order at O(k log k) worst case, so a few
  operator-dense rows cannot cause quadratic blowup), and adjacent
  duplicates are summed in emission order — fully deterministic
  accumulation. Collection runs the enumeration loops SINGLE-THREADED
  (OpenMP disabled in collect mode; the diagonal's OpenMP loop writes
  exactly one entry per column j, routed in collect mode to an
  index-addressed per-j slot — no append, thread-safe even if its
  pragma stays). Row-count increments and the prefix-sum additions are
  overflow-checked just like the byte sizes.
* **Compaction (CSR validity after merge)**: merging shortens rows in
  place, leaving per-row tails of dead entries. A single forward
  compaction pass then moves each row's merged entries toward the start
  of `colidx`/`val` while rewriting `rowptr`; the ORIGINAL raw row
  boundaries remain available throughout because they live in the
  still-allocated cursor array (freed only after compaction). Final
  invariants, asserted in debug builds: `rowptr[0] == 0`, `rowptr`
  monotonically nondecreasing, `rowptr[N] == nnz_merged ≤ nnz_raw`.
  The result is a conventional compact CSR; the arrays keep their raw
  capacity (no realloc).
* **Canonicalization**: exact-zero entries (as collected, or zero after
  summation) are RETAINED, not pruned. Combined with the deterministic
  enumeration this makes nnz and the structure a deterministic,
  RANK-IDENTICAL function of the input (model, def files, coupling
  values). Note the claim is rank-determinism, not coupling-independence:
  the inherited `dmv != 0.0` guards mean coupling values do influence
  which entries are emitted — identically on every rank.
* **Sink-mode discipline**: the collector saves the current sink mode,
  sets `HAM_SINK_TRACE_COLLECT` around each enumeration pass, and
  RESTORES the saved mode on every exit path — success, gate/allocation
  demotion, and error — before any fallback evaluation can run.
* **Cost**: enumeration runs twice; measured MakeHam cost is ~1e-5 of the
  diagonalization, so this is negligible (re-verified by the timing line).

**Diagonal.** Audited: `makeHam.c` already routes `list_Diagonal[j]`
through `AddHamElem` (line ~116), so the collector receives the diagonal
through the same macro; the unit test pins the equality anyway.

**When / scratch.** Collection runs after diagonalization and eigenvector
redistribution, inside the ExpecMode-2 orchestration, before any state is
evaluated. Task 1 of the plan audits (against the phase-2 write inventory)
that no makeHam site reads v0/v1 as INPUT in collect mode; they are
treated as free scratch at that point. `list_Diagonal` remaining intact
from the generation phase is pinned as an INVARIANT (Task-1 audit
verifies nothing between generation and observables writes it; the unit
test pins the collected diagonal). If the audit ever finds a mutation,
the remedy is a side-effect-free recomputation helper — NOT a re-run of
`diagonalcalc()`, which writes check files and progress output and must
not be called from observable orchestration.

**Failure classes** (distinct, both synchronized in §3c):

* Gate exceeded / allocation failed → energy demotion (memory reason),
  synchronized via the §3c finalize reduction.
* makeHam returns nonzero in collect mode → this indicates an invalid
  re-enumeration, not a size problem; treated as a hard error, because
  the same enumeration succeeded during generation and a silent fallback
  would mask corruption. Termination is collective-safe: the failing
  rank prints the message and calls the project's global abort path
  (`exitMPI`, which reaches `MPI_Abort`), so surviving ranks cannot
  block in the finalize collective.

### 3b. Streaming kernel

Per owned eigenstate x (state-panel layout, no MPI):

1. **Frozen dispatch preconditions** (matching the ACTUAL current code:
   `expec_energy_flct()` consumes x from `v0`, computes the fluctuation
   fields from `v0`, and copies `v0 → v1` itself): the driver loads x
   into `v0` in BOTH branches — this keeps the ExpecMode-0/1 path
   byte-identical, since that path is simply the fallback branch. In the
   FALLBACK branch, `expec_energy_flct()` then runs unchanged (its own
   `v0 → v1` copy included). In the KERNEL branch, the driver performs
   the `v0 → v1` copy itself (reproducing the evaluator's postcondition
   `v1 == x` that all downstream evaluators — GF fallbacks, S², NBodyG,
   AnomalousG — rely on), then invokes the kernel.
2. The kernel reads x from `v1`, computes y = H x into its PRIVATE y
   buffer (part of the collector context, counted in the gate), and
   writes `energy = Re(x†y)`, `Phys.var = y†y`, and the fluctuation
   fields per the frozen table below. **Postcondition: the kernel writes
   neither v0 nor v1.** (The fallback leaves `v0` holding its mltply
   result; the Task-1 audit re-verifies that nothing downstream in the
   FullDiag per-state loop reads `v0` after the energy step, so the
   branches' `v0` difference is unobservable; the equivalence suites
   would catch a violation. Defined remedy if the audit ever finds a
   `v0` consumer: the kernel branch copies its y buffer into `v0` —
   y IS H x, exactly the fallback's `v0` postcondition — at the cost of
   one vector copy.)
3. SpMV is OpenMP-parallel over rows (deterministic per-row dot products;
   no cross-row accumulation).
4. **Frozen per-model fluctuation-field table** (audited against
   `expec_energy_flct.c` at 0f15aabf; the kernel must reproduce EXACTLY
   this, including the not-written cells, which today retain whatever
   value the field previously held):

   | Model case (dispatch) | doublon/doublon2 | num/num2 | Sz/Sz2 | num_up/num_down | n_diag arrays |
   |---|---|---|---|---|---|
   | Hubbard, tJ, tJGC, Kondo, KondoGC (`_Hubbard`) | basis-diag D(k), D²(k) | basis-diag N(k), N²(k) | basis-diag: 0.5·S(k), 0.25·S²(k) | 0.5·(num±Sz_raw) derived from the same sums | 3: D(k), N(k), S(k) (squares derived on the fly) |
   | HubbardGC (`_HubbardGC`) | same as above | same | same | same | 3 |
   | SpinGC, half (`_HalfSpinGC`) | constants 0, 0 | constants NsiteMPI, NsiteMPI² | basis-diag 0.5·S(k), 0.25·S²(k) | 0.5·(NsiteMPI±Sz_raw) | 1: S(k) |
   | SpinGC, general (`_GeneralSpinGC`) | constants 0, 0 | constants NsiteMPI, NsiteMPI² | basis-diag (SiteToBit) | 0.5·(NsiteMPI±Sz_raw) | 1: S(k) |
   | Spin, canonical (half AND general) | constants 0, 0 | constants NsiteMPI, NsiteMPI² | constants 0.5·Total2SzMPI, (0.5·Total2SzMPI)² | **NOT WRITTEN** (stale-preserving) | 0 |
   | SpinlessFermion, canonical | constants 0, 0 | constants NeMPI, NeMPI² | constants 0, 0 | **NOT WRITTEN** | 0 (but statically ineligible: no makeHam branch) |
   | SpinlessFermionGC (`_SpinlessFermionGC`) | constants 0, 0 | basis-diag N(k), N²(k) | constants 0, 0 | per current code | ineligible: no makeHam branch |

   Here S(k) is the raw 2·Sz bit sum and D(k)/N(k) the per-basis doublon
   and particle-number counts; each basis-diagonal array is filled once
   by extracting the per-k bit computations of `expec_energy_flct()`
   into shared per-k helpers (3b-style function extraction: the Mode-1
   evaluator and the filler call the SAME helper — no algebra is
   reimplemented). Constant fields are assigned once per state exactly
   as the current code does; NOT-WRITTEN cells stay not written.

Timing lines (rank 0, item-17 frozen format):
`ExpecMode 2 timing (rank 0): energy map=%.3fs stream=%.3fs output=%.3fs`
(map = both enumeration passes + merge; output = Phys-field stores).

### 3c. Plan lifecycle and integration

* `TRACE_Q_ENERGY` joins the quantity enum (`TRACE_Q_NQUANT` = 3).
* **Two-phase construction replacing single-call immutability**:
  1. `TraceBuildPlan()` (unchanged signature) marks the energy slot
     PROVISIONAL: statically eligible unless `iInputHam != 0` (see
     below); GF slots are final at this point exactly as in 3b.
  2. The orchestrator (`phys_distributed.c`, MPI-capable) runs the
     collector, then calls a new `TraceFinalizeEnergyPlan(plan, local_ok,
     nnz_raw)` which performs ONE `MPI_Allreduce(MPI_MIN)` over
     `long long buf[3]`. A SUCCESSFUL rank contributes
     `{1, (long long)nnz_raw, -(long long)nnz_raw}` (the conversion is
     checked: `nnz_raw` must be representable, else the rank treats
     itself as failed); a FAILED rank (gate exceeded, any allocation
     failure INCLUDING the pre-gate counts array — such a rank has no
     valid nnz) contributes `{0, LLONG_MAX, LLONG_MAX}`. Every rank
     applies the same verdict logic to the identical reduced values, in
     this order: (1) `reduced[0] == 0` → demote the energy slot on
     EVERY rank (memory reason) — this covers one-rank, several-rank,
     and all-rank failures, and the sentinel values keep failed ranks
     from polluting the consistency check; (2) otherwise (all ranks
     succeeded) `reduced[1] != -reduced[2]` → abort on ALL ranks with a
     clear message (non-deterministic enumeration is a correctness
     error, not a fallback case). No rank-0-only diagnosis: all ranks
     hold the same reduced buffer and take the same branch.
  3. From finalize onward the plan is immutable; `TraceReportPlan()`,
     the kernel, and `phys_stateparallel_local_loop()` all run AFTER
     finalize and observe the same final plan. (The local loop keeps its
     MPI-free property: finalize happens in the orchestrator TU.)
* **Input-Hamiltonian eligibility**: when the run's H comes from
  `InputHam` (`iInputHam != 0`, reachable e.g. on the replicated and
  ScaLAPACK paths), re-running makeHam would build a DIFFERENT matrix.
  The energy slot is statically demoted with its own reason/INFO line.
  (Conversion of an input matrix to CSR is a possible future extension —
  out of scope.)
* **Symmetry-basis runs**: if the Task-1 audit finds any FullDiag
  configuration (e.g. the new momentum/symmetry work merged from develop)
  in which `expec_energy_flct` evaluates in a basis other than the one
  makeHam enumerates, those configurations are statically demoted with
  the unsupported-configuration INFO line until separately validated.
  The audit outcome is recorded in the plan document.
* **CSR ownership**: the collector returns a context object
  (`TraceHamCsr`: rowptr/colidx/val/nnz + y buffer + diagonal arrays)
  owned by the orchestrator and passed BY POINTER to the local loop next
  to the plan; it is freed by the orchestrator after the loop. No
  process-global CSR state; the only global is the transient
  `iHamSinkMode` during the two enumeration passes.
* **INFO lines** (verbatim strings — the single source for code, tests,
  and documentation):
  | condition | line |
  |---|---|
  | kernel active | `INFO: ExpecMode 2: the energy/fluctuation family uses the trace kernel.` |
  | memory demotion | `INFO: ExpecMode 2: the energy/fluctuation family uses the ExpecMode-1 fallback (the Hamiltonian buffer would exceed HPHI_TRACE_BUF_MAX_MB).` |
  | InputHam | `INFO: ExpecMode 2: the energy/fluctuation family uses the ExpecMode-1 fallback (the Hamiltonian was read from InputHam).` |
  | unsupported configuration | `INFO: ExpecMode 2: the energy/fluctuation family uses the ExpecMode-1 fallback (unsupported configuration).` |

  The fixed line becomes `INFO: ExpecMode 2: S2, NBodyG, and AnomalousG
  always use the ExpecMode-1 path in this version.` As in 3b, the log
  output indents these lines; tests match the quoted substrings.

### 3d. Dispatch

`phys_distributed_local.c`: when `plan->kernel[TRACE_Q_ENERGY]==1`, the
per-state loop calls the energy kernel instead of `expec_energy_flct()`;
everything else keeps its 3b dispatch. The authoritative dispatch
sequence is §3b step 1: the driver loads x into `v0` in both branches;
the `v0 → v1` copy is performed by `expec_energy_flct()` itself in the
fallback branch (unchanged, preserving byte-identity) and by the DRIVER
only in the kernel branch.

Numerical accumulation contract: `x†y`, `y†y`, and the diagonal-field
sums use OpenMP `reduction` clauses with `double complex` /`double`
accumulators — the same reproducibility class as the Mode-1 evaluators
(which already use OpenMP reductions); run-to-run bitwise stability
across thread counts is NOT claimed by either path, and the equivalence
criterion remains the tolerance-based one.

## 4. Memory and performance model

* nnz counts merged structural entries (zeros retained). Hubbard chain
  L=10: ~21·N ≈ 1.3e6 → CSR ≈ 32 MB. Transverse-field SpinGC L=14
  (N=16384): ≈ (L+1)·N ≈ 2.5e5 entries → ≈ 6 MB. Both far under the
  1024 MiB default cap; the cap matters for operator-dense
  InterAll/NBodyInterAll inputs, where the collector demotes cleanly (the
  count pass sizes the problem before any big allocation). Note the CSR
  is replicated per rank — the gate is per-rank, consistent with 3b's
  per-rank budgets, and the INFO line makes the demotion visible.
* Per-state cost: one CSR SpMV (O(nnz) flat-array flops) replaces one
  `mltply()` traversal of equal flop order but much higher per-state
  overhead (bit manipulation, branches, function dispatch). Expected
  outcome: the energy-family share of Mode-2 `CalcPhys` shrinks toward
  the GF-kernel profile; measured gate in §5.

## 5. Testing

1. **Unit test `expec_trace_ham_check`** (noMPI + MPI, registered like
   `expec_trace_map_check`):
   * Dense-vs-CSR: for small instances of EVERY enabled family —
     Hubbard (canonical), HubbardGC, tJ, **tJGC**, Kondo, **KondoGC**,
     Spin(1/2), Spin (general, S=1), SpinGC(1/2), SpinGC (general) —
     build the CSR and compare against the replicated dense makeHam
     matrix with a scale-aware criterion (|a−b| ≤ atol + rtol·max|H|,
     atol 1e-13). The enabled-model claim and this tested set MUST
     match one-to-one.
   * Independent reference: analytically known matrices compared against
     hand-written dense matrices — 2-site Hubbard at half filling,
     2-site Heisenberg, and **one complex Hermitian case** (2-site
     Hubbard with a complex hopping phase) verifying Hx, x†Hx, and y†y —
     catches errors common to both enumeration paths, including
     conjugation/orientation/phase mistakes that real-valued references
     cannot see.
   * Kernel-level energy/⟨H²⟩ and field check: apply the kernel to
     NORMALIZED NON-eigenvector states and compare energy, `Phys.var`
     (⟨H²⟩), and every fluctuation field of §3b's frozen table against
     the reference computation — the direct validation of the variance
     path, which zvo_phys alone cannot provide (Phys.var is not among
     the gathered/printed columns). Includes a duplicate-cancellation
     case (large opposite-sign duplicate entries) exercising the stable
     merge. All Phys fields are initialized to distinctive sentinel
     values before each kernel call, and the NOT-WRITTEN cells of §3b's
     table are asserted to retain their sentinel bit patterns.
   * Memory-gate boundary: via the collector API's explicit byte-cap
     parameter (not the env), one case just above and one just below the
     exact predicted size — where the prediction is the full §3a lifetime
     sum INCLUDING the sort workspace; the env path is covered once (any
     value).
   * Synchronized-demotion injection (MPI variant): an injectable
     allocation-failure hook fails the collector on exactly one rank;
     assert the energy slot demotes on ALL ranks, the context is cleaned
     up, the sink mode is restored, and the Mode-1 fallback produces
     equivalent output.
2. **Equivalence script**: existing 5 cases assert the energy kernel INFO
   line and compare all zvo_phys columns across Modes 0/1/2 at 1e-8. New
   case: small tJ chain — energy kernel ACTIVE, both GF quantities print
   unsupported-model fallback (taxonomy orthogonality), zvo_phys at 1e-8.
   New negative case: `InputHam` run asserts the input-Hamiltonian INFO
   line and Mode-2 == Mode-0 outputs.
3. **Regression**: full noMPI suite + MPI np=2/np=3 rounds (ExpecMode 0/1
   byte-identity re-verified by existing suites).
4. **Benchmark gate** (item-17 methodology): the 3b Mode-2 case (N=16384
   SpinGC, 16 procs, kugui `CalcPhys` 20.8 s) re-measured; target Mode-2
   `CalcPhys` ≤ 0.5× the 3b value, energy map/stream breakdown recorded;
   if missed, record a break-even analysis (3b precedent).

## 6. Documentation

* CalcMod `ExpecMode` (ja/en): energy/fluctuation moves from the
  always-fallback list to the traced set with its reason list (memory,
  input Hamiltonian, unsupported configuration); INFO catalogue updated
  verbatim (§3c's table is the source). The `HPHI_TRACE_BUF_MAX_MB`
  entry is extended: default (1024 MiB), that it now also caps the
  energy family's Hamiltonian buffer, and the troubleshooting step
  (raise the cap, or accept the Mode-1 fallback the INFO line reports).
  A one-line terminology note clarifies that "trace kernel" refers to
  the precomputed-mapping streaming evaluation technique, not to the
  matrix trace Tr(·). Deliberately NOT added: a per-quantity bypass
  environment variable (e.g. HPHI_DISABLE_TRACE_ENERGY) — `ExpecMode 1`
  is the supported bypass for the whole trace layer, and per-quantity
  toggles would multiply the test matrix for no user need.
* Appendix `parallel_fulldiag_{ja,en}.rst`: ExpecMode-2 description and
  benchmark column refresh after implementation.
* PR migration note addendum.

## 7. Non-goals and follow-ups

S² tracing (3d candidate); GF model expansion; tiled GF buffers;
CalcTimer integration; InputHam→CSR conversion; GPU 2-stage (issue #286).

## 8. Review dispositions

**Round 1** (all ten must_fix addressed in the body): storage-init bypass
(§3a sink abstraction), OpenMP-safe collection (§3a single-thread +
indexed diagonal), v0/v1 contract (§3b), `Phys.var` = ⟨H²⟩ contract
(§1/§2), rank-synchronized finalize (§3c), explicit build→finalize
lifecycle (§3c), two-pass memory accounting (§3a), model-coverage claim
scoped to makeHam + independent reference tests (§1/§5), InputHam
eligibility (§3c), kernel-level var validation (§5). Open questions:
both distributed backends eligible (with the InputHam rule); CSR
ownership is a context object; buffer contract per §3b; symmetry-basis
reachability is a Task-1 audit with static demotion as the safe default;
rank nnz mismatch aborts (correctness), allocation failure demotes
synchronized (capacity); zero entries retained.

**Round 2** (all five must_fix addressed): implementable single
`MPI_Allreduce(MIN)` over `{ok, nnz, −nnz}` with all-rank identical
verdict logic (§3c); dispatch preconditions corrected to the actual
`v0`-consuming evaluator with the driver-side `v0→v1` copy in the kernel
branch (§3b); complete allocation lifetime graph with `sizeof`-based
sizes, counts→rowptr in-place reuse, explicit cursors, and pre-gate
count-array failure handling (§3a); the per-model fluctuation-field
table frozen IN THIS DESIGN including constant and NOT-WRITTEN cells
(§3b); tJGC/KondoGC added to the mandatory test matrix with the
claim-equals-tested-set rule (§5). Round-2 should_fix folded in: stable
in-row sort with a cancellation test (§3a/§5); nnz claim narrowed to
rank-determinism, acknowledging the coupling-dependent `dmv` guards
(§3a); sink-mode save/restore on every exit path (§3a);
`list_Diagonal` pinned as an invariant with a side-effect-free
recompute as the only remedy — never re-running `diagonalcalc()` (§3a);
complex Hermitian reference matrix (§5); single-rank failure-injection
test for synchronized demotion (§5). Remaining round-2 open questions
resolved by those sections: counts transform in place into rowptr with a
separate cursor array (§3a); the sink-mode refactor covers every
`AddHamElem`/`HAM_OWNED_COL` site (makeHam.c, nbody_interall.c,
anomalous_pair.c — the audited 33-site inventory from the blocker fix is
the checklist, Task-1); the collector's per-rank replication risk is
accepted and documented in §4 (per-rank gate + INFO visibility), with
column-major storage and dense-conversion alternatives recorded as
rejected (§3a's row-major choice keeps deterministic OpenMP SpMV).

**Round 3** (both must_fix addressed): failed-rank sentinel encoding
`{0, LLONG_MAX, LLONG_MAX}` with demote-before-consistency verdict order
and checked `long long` conversion (§3c); in-place compaction pass with
cursor-array raw boundaries and asserted CSR invariants (§3a).
Should_fix folded in: §3d now defers to §3b's single authoritative
dispatch sequence; overflow checks extended to count increments and
prefix sums (§3a); hybrid stable sort removes the quadratic worst case
(§3a); accumulation contract stated (OpenMP reductions, same
reproducibility class as Mode 1 — §3d); sentinel assertions for
NOT-WRITTEN fields (§5); collective-safe abort via `exitMPI` for
enumeration failure (§3a). Open questions: the symmetry/momentum
eligibility predicate is a Task-1 audit deliverable recorded in the plan
(static demotion is the default until validated); the `v0`-consumer
question has a defined remedy (copy y into `v0`, §3b) should the audit
find one.

**Round 4** (one must_fix addressed): the stable-merge sort workspace is
now an explicit lifetime-graph row (sized `k_max` from pass 1, gated,
synchronized-demotion on failure) and the boundary tests predict the
full sum including it (§3a/§5). Confirmed assumptions: the Task-1 audits
are mandatory prerequisites; every collector allocation happens before
`TraceFinalizeEnergyPlan`, so no post-finalize rank-local allocation can
diverge the plan.
