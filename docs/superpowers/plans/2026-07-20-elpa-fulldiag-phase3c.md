# Phase 3c: CSR-Based Energy-Family Trace Kernel — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the per-state `expec_energy_flct()` fallback in
`ExpecMode 2` with a CSR SpMV kernel whose matrix is collected once from the
audited makeHam enumeration, for every makeHam-reachable model.

**Architecture:** A sink-mode abstraction in `hamstore.h` lets the existing
makeHam enumeration run a third time (after diagonalization) into a
two-pass count/fill collector that builds a compact CSR per rank. A
streaming kernel computes energy = Re(x†y), Phys.var = y†y (⟨H²⟩ — the
existing field contract), and the fluctuation columns from precomputed
diagonal coefficient arrays. A new plan slot `TRACE_Q_ENERGY` follows a
build→finalize lifecycle with one MIN-Allreduce for rank agreement.

**Tech Stack:** C99, MPI, OpenMP, CMake/ctest; existing phase-3b trace
infrastructure (`src/expec_trace.c`, `src/phys_distributed*.c`).

**Spec:** `docs/superpowers/specs/2026-07-20-elpa-fulldiag-phase3c-design.md`
(v3, design-review converged). Where this plan and the spec disagree, the
spec governs; stop and escalate.

## Global Constraints

- `Phys.var` stores ⟨H²⟩ = y†y, NEVER the variance (spec §1/§2).
- ExpecMode 0/1 must remain byte-identical: the fallback branch runs
  `expec_energy_flct()` unchanged; the all-fallback short-circuit in
  `TraceBuildPlan()` is untouched (spec §2).
- Single writer: kernel active ⇒ `expec_energy_flct()` skipped for those
  states; the family is never split (spec §2).
- Canonical Spin and SpinlessFermion leave `Phys.num_up`/`Phys.num_down`
  NOT WRITTEN (stale-preserving). The kernel must reproduce this (spec §3b
  frozen table).
- INFO strings are the verbatim table in spec §3c — copy exactly; tests and
  docs match those substrings.
- All size arithmetic uses the 3b checked style: `uintmax_t` guards before
  every multiplication/addition/conversion; row-count increments and prefix
  sums included (spec §3a). No `realloc` anywhere in the collector.
- Collection runs single-threaded (OpenMP disabled in collect mode) except
  the index-addressed diagonal slot; the in-row sort is stable (insertion
  below threshold 16, stable merge above, using the one gated workspace).
- Statically ineligible for the energy kernel: `iInputHam != 0`;
  symmetry/momentum configurations per the Task-1 audit predicate;
  SpinlessFermion(GC) (no makeHam branch).
- Equivalence tolerances: dense-vs-CSR `|a−b| ≤ 1e-13 + 1e-13·max|H|`;
  physics 1e-8 (existing equiv-script tolerance).
- Commits follow the repository convention (imperative subject, body
  explaining why), each task's tests pass before its commit.

## File Structure

| File | Responsibility |
|---|---|
| `src/include/hamstore.h` (modify) | Sink-mode enum + collect branch in `AddHamElem` |
| `src/global.c`, `src/include/global.h` (modify) | `iHamSinkMode` global |
| `src/makeHam.c` (modify) | Skip dense clearing in collect mode; full column range in collect mode |
| `src/include/expec_trace_ham.h` (create) | `TraceHamCsr`, collector + kernel + finalize API |
| `src/expec_trace_ham.c` (create) | Two-pass collector, merge/compaction, gate, SpMV kernel, coefficient fillers |
| `src/expec_energy_flct.c` (modify) | Extract per-k bit computations into shared helpers (3b style) |
| `src/include/expec_energy_flct.h` (modify) | Helper declarations |
| `src/expec_trace.c`, `src/include/expec_trace.h` (modify) | `TRACE_Q_ENERGY` slot, provisional build, finalize, INFO lines |
| `src/phys_distributed.c` (modify) | Orchestration: collect → finalize → report; timing lines |
| `src/phys_distributed_local.c` (modify) | Kernel-vs-fallback dispatch for the energy step |
| `test/unit/expec_trace_ham_check.c` (create) | Dense-vs-CSR, references, kernel-level, gate, injection tests |
| `test/CMakeLists.txt` (modify) | Register the unit test (pattern: `expec_trace_map_check`) |
| `test/fulldiag_expecmode_equiv.sh` (modify) | Energy INFO asserts; tJ case; InputHam negative case |
| `doc/{ja,en}/.../CalcMod_file_*.rst`, appendix, migration note (modify) | Task 7 |

---

### Task 1: Audits and the sink-mode refactor (behavior-preserving)

**Files:**
- Modify: `src/include/hamstore.h`, `src/include/global.h` (~line 81),
  `src/global.c` (~line 171), `src/makeHam.c` (lines 96–112 init block;
  column-range sites at ~137, ~188, ~233 and every other
  `iHamPanelActive ? HamColBegin : 1` site)
- Create: `.superpowers/sdd/task-1-audit.md` (audit report, gitignored dir)

**Interfaces:**
- Consumes: current globals `iHamPanelActive`, `Ham`, `Ham_local`.
- Produces: `enum HamSinkMode { HAM_SINK_DENSE_REPLICATED = 0,
  HAM_SINK_DENSE_PANEL, HAM_SINK_TRACE_COLLECT };`
  `extern int iHamSinkMode;` (global.h) — dense modes are DERIVED from
  `iHamPanelActive` at makeHam entry so existing callers need no change;
  `HAM_SINK_TRACE_COLLECT` is set only by the Task-2 collector.
  Two function-pointer-free collect hooks (weak coupling to Task 2):
  `void (*HamCollectSinkFn)(long int irow, long int jcol,
  double complex val);` stored in a global
  `extern HamCollectSinkFn hamCollectSink;` called by the macro in collect
  mode. Produces the audit report consumed by Tasks 2–5.

- [ ] **Step 1: Write the audit report** `.superpowers/sdd/task-1-audit.md`
  answering, with file:line evidence, each item: (a) every
  `AddHamElem`/`HAM_OWNED_COL` site (start from the 33-site inventory in the
  2026-07-17 blocker-fix report: 28 makeHam.c + 4 nbody_interall.c + 1
  anomalous_pair.c) and whether nbody_interall.c/anomalous_pair.c sites are
  reachable during FullDiag Hamiltonian generation (they are GF/pair paths —
  document why they are or are not part of H); (b) does any makeHam site
  read `v0`/`v1` as INPUT (phase-2 write-transformation inventory says no —
  re-verify); (c) is `list_Diagonal` written anywhere between
  `diagonalcalc()` and the observable phase (grep all writers); (d) the
  symmetry/momentum eligibility predicate: which `X->Def` flags the
  develop-merged momentum work sets (inspect `src/mltplySpinSym.c`,
  `src/include/readdef.h` diff) and whether any FullDiag+ExpecMode-2 run
  can reach a non-makeHam basis — output: the exact predicate expression
  for Task 5, or "unreachable, no predicate needed"; (e) confirm the
  §3b frozen table against `expec_energy_flct.c` line by line; (f) confirm
  no downstream consumer reads `v0` after the energy step in
  `phys_distributed_local.c`'s loop (current order: energy → cisajs →
  cisajscktaltdc → nbodyg → anomalousg → totalspin, all called with `v1`).

- [ ] **Step 2: Write the failing compile check** — add to
  `src/include/hamstore.h` (after the includes):

```c
typedef void (*HamCollectSinkFn)(long int irow, long int jcol,
                                 double complex val);

enum HamSinkMode {
  HAM_SINK_DENSE_REPLICATED = 0,
  HAM_SINK_DENSE_PANEL,
  HAM_SINK_TRACE_COLLECT
};
```

and change the macro body to:

```c
#define AddHamElem(irow, jcol, val)                                    \
  do {                                                                 \
    long int hs_i_ = (long int)(irow);                                 \
    long int hs_j_ = (long int)(jcol);                                 \
    assert(hs_i_ >= 1);                                                \
    if (iHamSinkMode == HAM_SINK_TRACE_COLLECT) {                      \
      hamCollectSink(hs_i_, hs_j_, (val));                             \
    } else if (iHamPanelActive) {                                      \
      assert(hs_j_ >= HamColBegin && hs_j_ <= HamColEnd);              \
      assert(hs_i_ <= HamPanelLd);                                     \
      Ham_local[(hs_j_ - HamColBegin) * HamPanelLd + (hs_i_ - 1)]      \
        += (val);                                                      \
    } else {                                                           \
      Ham[hs_i_][hs_j_] += (val);                                      \
    }                                                                  \
  } while (0)
```

  Declarations in `global.h` next to `iHamPanelActive`:
  `extern int iHamSinkMode;` and `extern HamCollectSinkFn hamCollectSink;`
  (include `hamstore.h`'s typedef via a forward-compatible plain function
  pointer type `void (*)(long int, long int, double complex)` if include
  order forbids the typedef — keep global.h self-contained). Definitions in
  `global.c`: `int iHamSinkMode = 0; ... hamCollectSink = NULL;`

- [ ] **Step 3: makeHam collect-mode behavior** — in `src/makeHam.c`:
  (1) wrap the dense-initialization block (the `if (!iHamPanelActive) {...
  Ham clearing ...} else {... panel clearing ...}` at lines ~96–112) in
  `if (iHamSinkMode != HAM_SINK_TRACE_COLLECT) { ... }`;
  (2) at every column-range site replace
  `iHamPanelActive ? HamColBegin : 1` with
  `(iHamPanelActive && iHamSinkMode != HAM_SINK_TRACE_COLLECT) ? HamColBegin : 1`
  and the matching `HamColEnd`/`i_max` expression the same way (collect
  mode always enumerates the full range) — introduce local helper macros
  `HS_JB(X)`/`HS_JE(X)` at the top of makeHam.c to avoid 28 hand-edits
  drifting;
  (3) guard the diagonal OpenMP loop's body: in collect mode it must call
  `hamCollectSink(j, j, list_Diagonal[j])` (index-addressed, keeps its
  pragma — the Task-2 sink routes j==row==col diagonal calls to the
  preallocated per-j slot);
  (4) disable OpenMP for the off-diagonal enumeration in collect mode: the
  loops are not OpenMP today except the diagonal — VERIFY in the audit and
  note it; if any other pragma exists, gate it with
  `iHamSinkMode != HAM_SINK_TRACE_COLLECT`.

- [ ] **Step 4: Behavior-preservation regression** —

Run: `cmake --build build_noMPI -j 8 && ctest --test-dir build_noMPI`
Expected: 142/142 PASS (the refactor is a no-op for both dense modes:
`iHamSinkMode` is 0 and never set to COLLECT yet).

- [ ] **Step 5: Commit**

```bash
git add src/include/hamstore.h src/include/global.h src/global.c src/makeHam.c
git commit -m "Introduce the Hamiltonian sink-mode abstraction (phase 3c Task 1)"
```

---

### Task 2: The CSR collector

**Files:**
- Create: `src/include/expec_trace_ham.h`, `src/expec_trace_ham.c`
- Modify: `src/CMakeLists.txt` (add expec_trace_ham.c to the HPhi sources)
- Test: `test/unit/expec_trace_ham_check.c` (create; part 1),
  `test/CMakeLists.txt` (register — copy the `expec_trace_map_check` block
  at lines 442–474 verbatim, renaming target and source; keep the `rt` and
  MPI link conditions)

**Interfaces:**
- Consumes: Task 1's sink mode + `hamCollectSink` hook; `HPHI_TRACE_BUF_MAX_MB`
  cap value passed in (NOT read from env here).
- Produces (in `expec_trace_ham.h`):

```c
typedef struct {
  long int n;          /* matrix dimension (idim_max) */
  long int nnz;        /* merged entries; rowptr[n] == nnz */
  long int *rowptr;    /* size n+1 */
  long int *colidx;    /* size >= nnz (raw capacity retained) */
  double complex *val; /* size >= nnz */
  double complex *y;   /* kernel scratch, size n */
  int n_diag;          /* number of diagonal coefficient arrays (0..3) */
  double *diag[3];     /* D(k), N(k), S(k) as applicable, each size n */
} TraceHamCsr;

/* Returns 1 on success (csr populated), 0 on demotion (gate/allocation;
   csr fully freed, sink mode restored). cap_bytes: the per-quantity cap.
   fail_alloc_at: test hook, -1 in production; allocation #k fails when
   k == fail_alloc_at (0-based). */
int TraceHamCollect(struct BindStruct *X, size_t cap_bytes,
                    int fail_alloc_at, TraceHamCsr *csr);
void TraceHamFree(TraceHamCsr *csr);
/* Exact gated peak for given counts — exposed for the boundary test. */
size_t TraceHamGatedBytes(long int n, long int nnz_raw, long int k_max,
                          int n_diag);
```

- Collection algorithm (all in `expec_trace_ham.c`, static helpers):
  counting sink (increments per-row counts; overflow-checked), fill sink
  (writes into row segments via cursors), both installed into
  `hamCollectSink` around `makeHam(X)` calls with save/set/restore of
  `iHamSinkMode` on every exit path; in-place prefix sum counts→rowptr;
  stable hybrid in-row sort (insertion < 16, else stable merge with the one
  `k_max`-sized workspace); adjacent-duplicate sum; forward compaction with
  cursor-held raw boundaries; debug asserts `rowptr[0]==0`, monotonic,
  `rowptr[n]==nnz`. `n_diag`/`diag[]` filling is Task 3 (collector
  allocates per the model's n_diag from the frozen table; fills zeros
  until Task 3 lands — unit test part 1 checks matrix only).
- makeHam returning nonzero in collect mode → print
  `"ERROR: ExpecMode 2 Hamiltonian re-enumeration failed"` and call
  `exitMPI(-1)` (collective-safe abort; spec §3a).

- [ ] **Step 1: Write the failing unit test (part 1)** —
  `test/unit/expec_trace_ham_check.c` with a `main()` that (pattern:
  `expec_trace_map_check.c` — same stan-file-free direct `BindStruct`
  setup used there, or `-sdry`-generated defs under a scratch dir, matching
  whichever that test uses; read it first and follow it):
  - For each model fixture — Hubbard L=4 half-filled, HubbardGC L=4,
    tJ L=4 (2 holes), tJGC L=4, Kondo 2×2, KondoGC 2×2, Spin-1/2 L=6
    Sz=0, Spin S=1 L=4 Sz=0, SpinGC-1/2 L=6 Γ=0.5, SpinGC S=1 L=4 —
    (a) run legacy replicated makeHam into dense `Ham`; (b) run
    `TraceHamCollect` with `cap_bytes = SIZE_MAX/2`, `fail_alloc_at=-1`;
    (c) expand CSR to dense and assert
    `cabs(dense_csr[i][j] - Ham[i][j]) <= 1e-13 + 1e-13 * hmax` for all
    i,j, where `hmax` = max |Ham element|.
  - Gate boundary: compute `size_t want = TraceHamGatedBytes(...)` from
    the counting quantities of the Hubbard fixture (re-derive nnz_raw with
    an initial oversized-cap collect, `TraceHamFree`, then re-collect);
    assert collect SUCCEEDS with `cap_bytes = want` and FAILS (returns 0,
    csr zeroed) with `cap_bytes = want - 1`.
  - Injection: `fail_alloc_at = 0..5` each → returns 0, no leak (run under
    the MPI variant in Task 6's failure test; here assert clean return and
    that a subsequent normal collect succeeds — sink restored).

- [ ] **Step 2: Register and run to verify it fails**

Run: `cmake --build build_noMPI -j 8 --target expec_trace_ham_check 2>&1 | tail -5`
Expected: link failure (`TraceHamCollect` undefined).

- [ ] **Step 3: Implement `expec_trace_ham.c`** per the interface above.
  Checked-arithmetic helper (file-local, mirrors 3b's guards):

```c
static int checked_bytes(uintmax_t a, uintmax_t b, size_t *out) {
  if (a != 0 && b > UINTMAX_MAX / a) return 0;
  uintmax_t r = a * b;
  if (r > SIZE_MAX) return 0;
  *out = (size_t)r; return 1;
}
```

  `TraceHamGatedBytes` = sum of: `(n+1)*sizeof(long int)` (rowptr),
  `(n+1)*sizeof(long int)` (cursors), `nnz_raw*sizeof(long int)` (colidx),
  `nnz_raw*sizeof(double complex)` (val), `n*sizeof(double complex)` (y),
  `n_diag*n*sizeof(double)` (coefficients),
  `k_max*(sizeof(long int)+sizeof(double complex))` (sort workspace) —
  every term via `checked_bytes`, returning `SIZE_MAX` on overflow (which
  can never pass a real cap).

- [ ] **Step 4: Run the unit test**

Run: `ctest --test-dir build_noMPI -R expec_trace_ham_check --output-on-failure`
Expected: PASS. Then the full suite: `ctest --test-dir build_noMPI` → all PASS.

- [ ] **Step 5: Commit**

```bash
git add src/include/expec_trace_ham.h src/expec_trace_ham.c src/CMakeLists.txt test/unit/expec_trace_ham_check.c test/CMakeLists.txt
git commit -m "Add the two-pass CSR Hamiltonian collector (phase 3c Task 2)"
```

---

### Task 3: Diagonal coefficient extraction

**Files:**
- Modify: `src/expec_energy_flct.c`, `src/include/expec_energy_flct.h`,
  `src/expec_trace_ham.c` (fill `diag[]`)
- Test: extend `test/unit/expec_trace_ham_check.c` (part 2)

**Interfaces:**
- Produces per-k helpers extracted from the existing per-model loops
  (names frozen here; each returns the RAW quantity the loop computes
  today — the caller applies the same scaling the current code applies):

```c
/* expec_energy_flct.h */
void EnergyFlctCoeff_Hubbard(struct BindStruct *X, long int k,
                             double *D, double *N, double *S);   /* canonical list_1 basis */
void EnergyFlctCoeff_HubbardGC(struct BindStruct *X, long int k,
                               double *D, double *N, double *S); /* GC basis */
void EnergyFlctCoeff_HalfSpinGC(struct BindStruct *X, long int k, double *S);
void EnergyFlctCoeff_GeneralSpinGC(struct BindStruct *X, long int k, double *S);
```

- The four existing evaluator loops (`expec_energy_flct_Hubbard`,
  `_HubbardGC`, `_HalfSpinGC`, `_GeneralSpinGC`) are rewritten to call the
  SAME helpers inside their loops (3b-style extraction: move the loop-body
  bit arithmetic verbatim into the helper; the evaluator's accumulation and
  scaling lines stay in place). No numerical change.
- `expec_trace_ham.c`'s collector fills `csr->diag[]` per the frozen table:
  Hubbard family/HubbardGC → n_diag=3 (D,N,S); SpinGC (both) → n_diag=1
  (S); canonical Spin → n_diag=0.

- [ ] **Step 1: Write the failing test (part 2)** — for the Hubbard and
  SpinGC-1/2 fixtures: run the legacy `expec_energy_flct()` on a random
  normalized vector in `v0`, snapshot the 8 Phys fields; then compute the
  same fields from `csr->diag[]` arrays (⟨D⟩=Σ|x|²D(k) etc. with the frozen
  table's scalings: Sz=0.5·ΣS, Sz2=0.25·ΣS², num_up/down=0.5(num±ΣS));
  assert agreement ≤ 1e-12. For the canonical-Spin fixture: set
  `X->Phys.num_up = 4321.0` sentinel, run the constant-path computation,
  assert num_up is UNCHANGED (stale-preserving) and doublon==0,
  num==NsiteMPI, Sz==0.5·Total2SzMPI.

- [ ] **Step 2: Run to verify it fails** (helpers undefined → link error).

- [ ] **Step 3: Extract the helpers and fill `diag[]`.** Byte-identity
  check for the extraction itself:

Run: `ctest --test-dir build_noMPI` (full suite)
Expected: all PASS — the equivalence/regression tests pin the evaluator.

- [ ] **Step 4: Run the new test part** → PASS.

- [ ] **Step 5: Commit**

```bash
git add src/expec_energy_flct.c src/include/expec_energy_flct.h src/expec_trace_ham.c test/unit/expec_trace_ham_check.c
git commit -m "Extract per-basis fluctuation coefficients and fill the CSR context (phase 3c Task 3)"
```

---

### Task 4: The streaming kernel

**Files:**
- Modify: `src/expec_trace_ham.c`, `src/include/expec_trace_ham.h`
- Test: extend `test/unit/expec_trace_ham_check.c` (part 3)

**Interfaces:**
- Produces:

```c
/* Evaluate the energy family for the state in v1 (v1 == x, length n,
   1-based like the evaluators). Writes X->Phys.energy, X->Phys.var and
   the frozen-table fluctuation fields. NEVER writes v0 or v1. */
int TraceEnergyEvalState(struct BindStruct *X, const TraceHamCsr *csr);
```

- SpMV: `#pragma omp parallel for` over rows; per-row dot product into
  `csr->y`; then `energy`/`var` via `reduction(+:...)` loops (double /
  double complex accumulators — same reproducibility class as Mode 1).
- Fluctuation fields per the frozen table: basis-diagonal models use
  `diag[]` sums (squares from the same arrays on the fly); constant models
  assign the constants; canonical Spin/SpinlessFermion write NOTHING to
  num_up/num_down.

- [ ] **Step 1: Write the failing test (part 3)** —
  - Complex Hermitian reference: 2-site Hubbard with complex hopping
    `t = 0.3 + 0.4i` (hand-written 4×4 (per spin sector) dense reference
    in the test): assert CSR-vs-reference elementwise ≤ 1e-13; for three
    fixed non-eigenvector normalized states x: dense y=Hx vs kernel
    energy/var: `|energy - x†Hx| ≤ 1e-12`, `|var - |Hx|²| ≤ 1e-12`.
  - Cancellation: a fixture whose def file defines the same transfer twice
    with ±1e8 amplitudes (duplicates summing to 0) — assert the merged
    entry is exactly the sum the dense path produces (tolerance as in
    Task 2) and kernel energy matches dense.
  - Sentinels: set every Phys field to distinct sentinels before
    `TraceEnergyEvalState`; assert the NOT-WRITTEN cells retain them.

- [ ] **Step 2: Run to verify it fails.**
- [ ] **Step 3: Implement `TraceEnergyEvalState`.**
- [ ] **Step 4: Run unit test + full noMPI suite** → PASS.
- [ ] **Step 5: Commit**

```bash
git add src/expec_trace_ham.c src/include/expec_trace_ham.h test/unit/expec_trace_ham_check.c
git commit -m "Add the CSR streaming energy-family kernel (phase 3c Task 4)"
```

---

### Task 5: Plan slot, finalize, dispatch, INFO, timing

**Files:**
- Modify: `src/include/expec_trace.h` (enum + docs), `src/expec_trace.c`
  (build + report), `src/phys_distributed.c` (orchestration ~lines
  126–140), `src/phys_distributed_local.c` (dispatch ~line 89)

**Interfaces:**
- `TRACE_Q_ENERGY` added BEFORE `TRACE_Q_NQUANT`.
- `TraceBuildPlan()` energy slot: `kernel=1` provisionally unless
  `X->Def.iInputHam != 0` (new `demoted_input_ham[TRACE_Q_NQUANT]` field,
  energy-only) or the Task-1 symmetry predicate holds (new
  `demoted_unsupported_config[...]`, energy-only; omit the field entirely
  if the audit concluded "unreachable"). GF slots unchanged.
- New:

```c
/* expec_trace.h */
void TraceFinalizeEnergyPlan(TraceExecutionPlan *plan, int local_ok,
                             long int nnz_raw);
```

  implemented in `src/phys_distributed.c`-adjacent MPI code (NOT the
  MPI-free local TU): success rank contributes
  `{1, (long long)nnz_raw, -(long long)nnz_raw}` (checked conversion,
  else treat self as failed), failed rank `{0, LLONG_MAX, LLONG_MAX}`;
  one `MPI_Allreduce(MPI_IN_PLACE, buf, 3, MPI_LONG_LONG, MPI_MIN, ...)`;
  verdict order: `buf[0]==0` → demote everywhere (set
  `kernel[TRACE_Q_ENERGY]=0`, `demoted_memory[...]=1`); else
  `buf[1] != -buf[2]` → `fprintf` + `exitMPI(-1)`.
- Orchestration order in `phys_distributed.c`: build plan → if energy slot
  provisional: `TraceHamCollect` (cap from the SAME Bcast-ed `gbuf_max`) →
  `TraceFinalizeEnergyPlan` → `TraceReportPlan` (now prints final state) →
  loop. Timing: wrap collect in `clock_gettime` pair; print
  `ExpecMode 2 timing (rank 0): energy map=%.3fs stream=%.3fs output=%.3fs`
  after the loop (stream time accumulated inside the kernel via the
  existing TraceGetTimings pattern; output folded into stream is
  acceptable — keep the three-field format).
- `TraceReportPlan()`: use the verbatim INFO table from spec §3c; the
  fixed line drops "energy/fluctuation".
- `phys_distributed_local.c` dispatch (replacing the unconditional call):

```c
    if (!plan->kernel[TRACE_Q_ENERGY]) {
      if (expec_energy_flct(X) != 0) { rc = -1; break; }
    } else {
      for (j = 0; j < NN; j++) v1[j + 1] = v0[j + 1]; /* evaluator's v0->v1 postcondition */
      if (TraceEnergyEvalState(X, ham_csr) != 0) { rc = -1; break; }
    }
```

  (`ham_csr` threaded as a new `const TraceHamCsr *` parameter of
  `phys_stateparallel_local_loop()`, NULL when the slot is fallback.)

- [ ] **Step 1: Failing test** — extend the unit test (part 4, MPI
  variant): build a 2-rank plan where rank 1 collects with
  `fail_alloc_at=0`; assert after `TraceFinalizeEnergyPlan` BOTH ranks
  have `kernel[TRACE_Q_ENERGY]==0` and `demoted_memory==1`.
- [ ] **Step 2: Verify it fails** (function undefined).
- [ ] **Step 3: Implement** (enum, fields, finalize, report, orchestration,
  dispatch, timing).
- [ ] **Step 4: Run** noMPI suite (142+) AND MPI rounds:
  `MPIRUN="mpiexec -np 2" ctest --test-dir <mpi-build> -R "expec_trace|fulldiag_expecmode_equiv|elpa"` → PASS (equiv script still
  passes because Mode-2 energy output equals Mode-0 at 1e-8).
- [ ] **Step 5: Commit**

```bash
git add src/include/expec_trace.h src/expec_trace.c src/phys_distributed.c src/phys_distributed_local.c test/unit/expec_trace_ham_check.c
git commit -m "Wire the energy trace kernel into the ExpecMode-2 plan (phase 3c Task 5)"
```

---

### Task 6: Equivalence-script and MPI failure coverage

**Files:**
- Modify: `test/fulldiag_expecmode_equiv.sh`

**Interfaces:** Consumes the verbatim INFO strings (spec §3c) and the
existing case/assert helpers in the script (read lines 63–300 first).

- [ ] **Step 1: Extend the script (failing first):**
  - Every existing Mode-2 case additionally greps
    `"the energy/fluctuation family uses the trace kernel."`.
  - New case 6 (tJ chain, L=4, J=1, t=1, 2 holes, correlation): asserts the
    energy kernel line AND both GF `(unsupported model)` lines; zvo_phys
    Modes 0/1/2 agree at 1e-8.
  - New case 7 (InputHam negative): case-1 fixture with `OutputHam 1` run
    once (np=1) to produce the Ham file, then `InputHam` + `ExpecMode 2`:
    asserts `"(the Hamiltonian was read from InputHam)"` and Mode-2 ==
    Mode-0 zvo_phys at 1e-8.
- [ ] **Step 2: Run np=2 and np=3 rounds** — before Task 5's merge these
  asserts fail (that is the failing-first evidence if Task 6 is developed
  against Task 4's tree; when executed after Task 5, instead flip one
  expected string to a wrong value, observe the script fail, restore —
  record the non-vacuousness check in the report).
- [ ] **Step 3: Verify green**: `MPIRUN="mpiexec -np 2" ctest -R fulldiag_expecmode_equiv --test-dir <mpi-build>` and np=3 → PASS.
- [ ] **Step 4: Commit**

```bash
git add test/fulldiag_expecmode_equiv.sh
git commit -m "Cover the energy trace kernel in the ExpecMode equivalence script (phase 3c Task 6)"
```

---

### Task 7: Documentation

**Files:**
- Modify: `doc/ja/source/filespecification/expertmode_ja/CalcMod_file_ja.rst`,
  `doc/en/source/filespecification/expertmode_en/CalcMod_file_en.rst`
  (the `ExpecMode` entries), `doc/{ja,en}/source/technical/parallel_fulldiag_{ja,en}.rst`
  (ExpecMode-2 description), `docs/superpowers/specs/` migration-note
  addendum file `2026-07-20-phase3c-migration-note.md` (create)

- [ ] **Step 1: CalcMod (ja/en)**: move energy/fluctuation out of the
  always-fallback sentence; add its three fallback reasons and the four
  verbatim INFO lines (spec §3c table); extend the `HPHI_TRACE_BUF_MAX_MB`
  text (default 1024 MiB; now also caps the energy family's Hamiltonian
  buffer; troubleshooting: raise the cap or accept the reported fallback);
  add the one-line terminology note ("trace kernel" = precomputed-mapping
  streaming evaluation, not the matrix trace).
- [ ] **Step 2: Appendix (ja/en)**: update the ExpecMode-2 paragraph
  (energy family now traced for all makeHam models; S²/NBodyG/AnomalousG
  remain Mode-1) — do NOT touch the benchmark numbers yet (Task 8 refreshes
  them with measurements).
- [ ] **Step 3: Migration note**: Mode-2 outputs unchanged semantically;
  `var` still genuine ⟨H²⟩−⟨H⟩² downstream; new INFO lines listed; memory
  note (per-rank CSR ≈ 24·nnz bytes; cap shared with GF buffers).
- [ ] **Step 4: Build check**: both `sphinx-build -b html -q` runs, zero
  warnings. Commit:

```bash
git add doc docs/superpowers/specs/2026-07-20-phase3c-migration-note.md
git commit -m "Document the energy-family trace kernel (phase 3c Task 7)"
```

---

### Task 8: Hardware validation and benchmark gate (controller task — not a subagent)

- [ ] clavius: rsync + rebuild `build_elpa`; ELPA rounds np=2/np=3
  (`elpa|expec_trace|green_partial|fulldiag_expecmode_equiv|fulldiag_elpa`)
  → all PASS; Debug (assert) build one equiv round (CSR invariant asserts
  live).
- [ ] kugui: rebuild `build_kugui`; re-run the appendix ExpecMode benchmark
  (`test/manual/benchmarks/kugui_bench_cpu.pbs`, m0/m1/m2 rows only is
  fine): record `CalcPhys` and the energy map/stream breakdown. Gate:
  N=16384 Mode-2 `CalcPhys` ≤ 10.4 s (0.5× the 3b value 20.8 s); if
  missed, record the break-even analysis instead (spec §5.4).
- [ ] Update `test/manual/elpa_gpu_check.md` (phase-3c section) and the
  appendix benchmark table/figure (`make_bench_figs.py` rerun) with the
  new numbers; sphinx zero-warning check; commit.

---

## Self-Review

- Spec coverage: §3a → Tasks 1–2; §3b table/helpers → Tasks 3–4; §3c
  lifecycle/INFO → Task 5; §5.1 → Tasks 2–4 (parts 1–3) + Task 5 (part 4);
  §5.2 → Task 6; §5.3 → every task's regression step; §5.4 → Task 8; §6 →
  Task 7. InputHam/symmetry demotions → Task 5; audits → Task 1.
- Placeholders: none — every step names exact files, code, and commands;
  the two audit-dependent branches (symmetry predicate, v0 consumer) have
  their defined outcomes written in Tasks 1/5.
- Type consistency: `TraceHamCsr`, `TraceHamCollect`, `TraceHamGatedBytes`,
  `TraceEnergyEvalState`, `TraceFinalizeEnergyPlan`,
  `EnergyFlctCoeff_*` names are used identically across Tasks 2–6.
