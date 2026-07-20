# ELPA FullDiag Phase 3c: CSR-based trace kernel for the energy family

Status: draft for design review (Phase 0)
Depends on: phase 3b (`2026-07-11-elpa-fulldiag-phase3-design.md` §3, trace
kernels for one-/two-body Green functions), the phase-3b migration note, and
the merged `feature/elpa-fulldiag` state as of 0f15aabf (develop merged).

## 1. Goal and scope

Remove the dominant remaining `ExpecMode 2` cost: the energy-family
evaluator (`expec_energy_flct()`), which still runs per state through the
Mode-1 fallback path and was measured to account for most of the Mode-2
observable wall time once the GF kernels landed (clavius N=4900: GF kernels
~1.1 s; kugui N=16384: `CalcPhys` 20.8 s, GF kernels ~2 s of it).

In scope:

* A per-rank sparse (CSR) representation of the full Hamiltonian H,
  collected once per run by re-running the audited makeHam element
  enumeration with a collecting sink.
* A streaming energy-family kernel: for each owned eigenstate x,
  y = H x (no MPI), energy = Re(x†y), ⟨H²⟩ = y†y, var = ⟨H²⟩ − ⟨H⟩²,
  plus the diagonal fluctuation quantities (doublon, num, Sz and their
  squares, per the model's existing column set) via precomputed per-basis
  coefficient arrays.
* Plan integration as a new quantity slot `TRACE_Q_ENERGY` with the
  existing exclusive demotion taxonomy (memory gate only).
* **Model coverage: every model FullDiag supports** — the collector reuses
  makeHam, so no per-model capability table applies to the energy family
  (unlike the GF kernels). This includes tJ/Kondo/general spin, which the
  GF kernels do not cover.

Out of scope (unchanged Mode-1 behavior): S² (`expec_totalspin`), NBodyG,
AnomalousG, the GF-kernel model coverage, tiled GF buffers, CalcTimer
integration of the trace timings (stdout timing lines only, as in 3b).

## 2. Inherited invariants (unchanged)

* **`var` semantics**: `var` is an eigenvector-quality check; it must be
  computed as a genuine ⟨H²⟩ − ⟨H⟩² per state (never substituted by the
  eigenvalue). The kernel satisfies this via y†y with y = H x. The Mode-0/1
  numerical agreement criterion for `var` stays the equivalence-test 1e-8.
* **Single writer**: when the energy kernel is active it produces ALL
  zvo_phys columns for its states; `phys_stateparallel_local_loop()` must
  then skip `expec_energy_flct()` for those states. No partial split of the
  energy family (it is one evaluator writing one row).
* **Single immutable plan**: `TraceBuildPlan()` remains the only
  kernel-vs-fallback decision point; all consumers (INFO reporting, kernel,
  fallback loop) read the same `TraceExecutionPlan` instance.
* **Cross-rank agreement**: the plan inputs for the energy slot must be
  rank-independent. The CSR size (nnz) is a deterministic function of the
  model/def files and is identical on every rank (every rank enumerates the
  full H), so the memory-gate outcome is identical without extra
  communication. The `HPHI_TRACE_BUF_MAX_MB` cap stays rank-0-parsed +
  Bcast, as in 3b.
* **Guards**: the collector inherits the AddHamElem validity guards from
  the blocker fix (`tmp_off > 0`, `dmv != 0.0`) unchanged, because it hooks
  the same macro the dense sinks use.

## 3. Design

### 3a. CSR collection (`src/expec_trace_ham.c`, new TU)

* **Sink hook**: `hamstore.h` gains a third sink. A new global flag
  `iHamTraceCollect` is checked FIRST inside `AddHamElem`; when set, the
  macro appends the triplet `(irow, jcol, val)` to the collector instead of
  writing any dense storage. The two dense branches are untouched. The
  collector pass runs with `iHamPanelActive == 0` so that every
  column-range site in makeHam (`hs_jb = 1, hs_je = i_max`) enumerates the
  FULL column range on every rank.
* **When**: after diagonalization and eigenvector redistribution, at plan
  build/execution time inside the ExpecMode-2 path — v0/v1 and the other
  generation scratch are reusable at that point. MakeHam's enumeration cost
  is negligible (measured ~1e-5 of the diagonalization time).
* **Diagonal**: the collected matrix must equal the dense H including the
  diagonal contribution (`list_Diagonal`). If the diagonal write site does
  not go through `AddHamElem`, the collector adds an explicit diagonal
  append pass; either way the unit test in §5 pins the equality.
* **Storage & growth**: triplets in growable arrays (doubling realloc);
  running total bytes checked against the energy slot's memory budget
  (§3c) during growth — on exceed, free everything and record the memory
  demotion (no partial state).
* **Merge**: sort triplets by (row, col), sum duplicates, build CSR
  (rowptr: (N+1)×8 B; colidx: nnz×8 B; val: nnz×16 B). Duplicate summation
  order may differ from dense accumulation order; agreement is therefore
  checked at 1e-13 (unit test) and 1e-8 (physics), not bitwise.
* **Failure handling**: any allocation failure inside the collector is
  treated exactly like the memory gate (demote to fallback, free, INFO
  line) — never a fatal error, because the Mode-1 path remains available.

### 3b. Streaming kernel

Per owned eigenstate x (state-panel layout, no MPI in the loop):

1. y = H x via CSR SpMV (OpenMP over rows, as the GF kernels do).
2. energy = Re(x†y); ⟨H²⟩ = y†y; var = ⟨H²⟩ − energy².
3. Diagonal quantities: for each model-defined column (doublon, num, Sz,
   and squares), precompute once per run a per-basis coefficient array
   d_q[k] (by extracting the same per-k bit computations
   `expec_energy_flct()` performs today into shared per-k helpers — a
   function-extraction refactor in the 3b style: the Mode-1 evaluator and
   the coefficient filler call the SAME extracted helper, no algebra is
   reimplemented), then ⟨D⟩ = Σ_k |x_k|² d_q[k] and ⟨D²⟩ = Σ_k |x_k|²
   d_q[k]² (diagonal operators square diagonally).
4. Write the same Phys fields the Mode-1 evaluator writes (identical
   column set per model, identical output formatting downstream).

Timing: two rank-0 stdout lines in the established format —
`ExpecMode 2 timing (rank 0): energy map=%.3fs stream=%.3fs output=%.3fs`
(map = collection+merge; output = the Phys-field write is folded into
stream if not separable; the line format freeze from item 17 applies).

### 3c. Plan integration

* `TRACE_Q_ENERGY` joins the quantity enum (`TRACE_Q_NQUANT` grows to 3).
* Static capability: always eligible (no model gate). `no_operators` and
  `demoted_shared_evaluator` are structurally 0 for the energy slot (H
  always exists; `expec_energy_flct` is not shared with any always-fallback
  quantity — S² lives in `expec_totalspin`, a separate evaluator).
* Memory gate: `gbuf_bytes[TRACE_Q_ENERGY]` = the verified CSR + diagonal
  coefficient allocation size: final CSR = 24·nnz + 8·(N+1); transient
  triplet peak during collection = 32·nnz (8 B row + 8 B col + 16 B val);
  the gate checks max(final, transient) + 8·N·(number of diagonal
  coefficient arrays).
  Checked against the same `HPHI_TRACE_BUF_MAX_MB` per-quantity cap.
  Because nnz is not known before enumeration, the gate is enforced DURING
  collection (§3a) rather than precomputed; the plan is still built first
  with `kernel[TRACE_Q_ENERGY]=1` and may be downgraded by the collector
  before any state is evaluated — this is the one deliberate relaxation of
  the "immutable plan" rule, bounded as follows: the downgrade happens at a
  single point (collector exit), before the INFO report and before any
  evaluation, on all ranks identically (nnz is rank-independent), so every
  consumer still observes one final, agreed plan. TraceReportPlan() prints
  AFTER the collector so the INFO lines reflect the final plan.
* INFO lines: the per-quantity line set gains `energy/fluctuation`
  variants: kernel line and memory-fallback line (no unsupported-model and
  no no-operators variant for this slot). The fixed line becomes
  `INFO: ExpecMode 2: S2, NBodyG, and AnomalousG always use the
  ExpecMode-1 path in this version.` (energy/fluctuation removed from it).

### 3d. Dispatch

* `phys_distributed.c` (orchestrator): collection runs once (all ranks)
  between plan build and the per-state loop; timing lines printed with the
  existing ones.
* `phys_distributed_local.c` (MPI-free): when
  `plan->kernel[TRACE_Q_ENERGY]==1`, the per-state loop calls the energy
  kernel's per-state evaluation instead of `expec_energy_flct()`; all other
  quantities keep their 3b dispatch. ExpecMode 0/1 behavior is
  byte-identical to before (all-fallback plan short-circuit, as in 3b).

## 4. Memory and performance model

* nnz(H) ≈ N · (off-diagonal terms per column + 1). Hubbard chain L=10:
  ~21·N ≈ 1.3e6 → CSR ≈ 32 MB. Transverse-field SpinGC L=14: ~(L+1)·N ≈
  2.5e5·... ≈ 6 MB. Well under the 1024 MiB default cap for every case in
  the benchmark suite; the cap exists for pathological operator-dense
  models (large InterAll sets).
* Per-state cost: 1 SpMV (O(nnz) flops, flat arrays) replaces one full
  `mltply()` traversal (same flop order but per-state bit-manipulation and
  branch overhead). Expected effect: the energy-family share of Mode-2
  `CalcPhys` drops toward the GF-kernel profile; benchmark gate in §5.

## 5. Testing

1. **Unit test `expec_trace_ham_check`** (noMPI + MPI variants, registered
   like `expec_trace_map_check`): for small instances of EVERY reachable
   model family — Hubbard, HubbardGC, Spin(1/2), SpinGC(1/2), general-spin
   Spin (S=1), tJ, Kondo, SpinlessFermion — build the CSR via the collector
   and compare against the dense replicated makeHam matrix elementwise
   (tol 1e-13). This directly certifies the model-complete claim, including
   models the GF kernels do not support. Also: a memory-gate boundary case
   (tiny `HPHI_TRACE_BUF_MAX_MB` forces the collector demotion path and the
   INFO line).
2. **Equivalence script** (`test/fulldiag_expecmode_equiv.sh`): existing 5
   cases now also assert the energy INFO line (kernel active) and continue
   to compare ALL zvo_phys columns (energy, var, doublon, num, Sz, S²…)
   across Modes 0/1/2 at 1e-8. One NEW case: a small tJ chain — asserts
   energy kernel ACTIVE while both GF quantities print the
   unsupported-model fallback line (taxonomy orthogonality), zvo_phys
   agreement at 1e-8.
3. **Regression**: full noMPI suite + MPI rounds np=2/np=3 (the ExpecMode
   0/1 byte-identity guarantee is re-verified by the existing suites).
4. **Benchmark gate** (kugui or clavius, the item-17 methodology): the 3b
   Mode-2 benchmark case (N=16384 SpinGC, 16 procs — kugui `CalcPhys`
   20.8 s) re-measured; target: Mode-2 `CalcPhys` ≤ 0.5× the 3b value with
   the energy map/stream breakdown recorded. If the target is missed,
   record a break-even analysis instead of dropping the item (the 3b
   precedent).

## 6. Documentation

* CalcMod `ExpecMode` entry (ja/en): move energy/fluctuation from the
  always-fallback list to the traced set (with its own fallback reason
  list: memory gate only, model-independent); update the INFO line
  catalogue verbatim.
* Appendix `parallel_fulldiag_{ja,en}.rst`: update the ExpecMode-2
  description (trace kernels now cover energy family for all models +
  one-/two-body GF for the four supported models) and re-run the
  ExpecMode benchmark column after the implementation lands.
* PR migration note addendum (Mode-2 semantics unchanged for outputs;
  var still genuine; new INFO lines).

## 7. Non-goals and follow-ups

* S² tracing (`expec_totalspin`) — 3d candidate.
* GF-kernel model expansion (tJ/Kondo `is_gc` grouping) — separate.
* Tiled GF buffers, CalcTimer integration — recorded follow-ups.
* GPU ELPA 2-stage memory profile — issue #286.

## 8. Risks / open questions (for the design review)

1. **Diagonal write site**: if `list_Diagonal` is applied outside
   `AddHamElem`, the collector needs the explicit diagonal pass (§3a);
   the unit test pins the total either way. Resolution planned at
   implementation time; no design impact.
2. **Plan downgrade point** (§3c): the bounded relaxation of plan
   immutability (single downgrade point, pre-report, rank-identical).
   Alternative rejected: a counting pre-pass to size nnz before the plan —
   doubles the enumeration for no observable benefit.
3. **Scratch reuse**: the collector runs post-diagonalization; v0/v1 are
   scratch there. If any model's makeHam path reads v0/v1 as INPUT (not
   scratch), the collector must snapshot them first — to be audited in
   Task 1 of the plan (the phase-2 write-transformation inventory lists
   every site).
4. **Float ordering**: duplicate-merge order differs from dense `+=`
   order; agreement is by tolerance (1e-13 matrix / 1e-8 physics), not
   bitwise — consistent with the 3b precedent for kernel-vs-evaluator
   comparisons.
