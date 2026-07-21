# PR migration-note draft: phase 3c (energy-family trace kernel)

Status: draft, to be pasted into the PR description (the "Migration
notes" section) once phase 3c lands, appended after the phase 3b
migration-note bullets
(`docs/superpowers/specs/2026-07-12-phase3b-migration-note.md`). Written
to match that note's tone (keyword additions/deprecations, behavior
changes with attribution). Task 8's clavius benchmark gate is the
authority for any final numeric speedup claim in the PR text; this note
makes no speed promise beyond "expected to help."

---

## Migration notes (phase 3c)

- **`ExpecMode 2` now also traces the energy/fluctuation family**
  (previously: energy/fluctuation, including the `var` column,
  unconditionally used the `ExpecMode 1` path, same as `S2`, `NBodyG`,
  and `AnomalousG`). HPhi precomputes the Hamiltonian once per rank into
  a compact CSR (compressed sparse row) matrix, then streams every
  owned eigenstate through a CSR sparse matrix-vector product instead of
  a full `mltply`-style traversal, amortizing the per-state Hamiltonian
  traversal overhead the same way the phase 3b Green-function kernels
  amortize per-state operator dispatch. `S2`, `NBodyG`, and
  `AnomalousG` are **not** kernelized in this phase and remain
  unconditionally on the `ExpecMode 1` path.

- **Output semantics are unchanged.** `ExpecMode 2` still writes the
  same columns with the same values (within the existing floating-point
  tolerance) as `ExpecMode 0`/`1` — this phase only changes which code
  path computes the energy family, not what it computes or how it is
  reported. In particular `var` is still a genuine
  `⟨H²⟩ − ⟨H⟩²` downstream: `Phys.var` continues to store `⟨H²⟩` exactly
  as the `ExpecMode 1` evaluator (`expec_energy_flct`) always has, and
  downstream code subtracts `⟨H⟩²` to obtain the variance. This field
  contract does not change whether the trace kernel or the
  `ExpecMode 1` fallback produced the value.

- **Model coverage is much broader than the phase 3b Green-function
  kernels.** The phase 3b one-body/two-body trace kernels are validated
  and enabled for only four (model, spin-representation) rows
  (`Hubbard`, `HubbardGC`, half-integer `Spin`, half-integer `SpinGC`).
  The energy-family kernel, by contrast, is eligible for **every**
  model that `FullDiag` can run — every makeHam-reachable model,
  including `tJ`/`tJGC`, `Kondo`/`KondoGC`, and general-spin
  `Spin`/`SpinGC` (`S ≥ 1`) — because it re-collects the Hamiltonian
  generically rather than relying on a per-model capability table.
  Consequently the energy family has **no** "unsupported model"
  fallback case and **no** shared-evaluator fallback case (unlike the
  two-body Green function, it does not share its evaluator with any
  other always-fallback quantity).

- **Only two fallback reasons for the energy family**, checked in this
  order:
  1. **`InputHam`**: if this run's Hamiltonian was read from `InputHam`
     (`InputHam 1`), the trace kernel cannot rebuild a matching
     Hamiltonian by re-enumerating the model — doing so would re-derive
     the matrix from the model definition rather than reading back the
     one that was actually diagonalized — so the energy family
     unconditionally falls back to `ExpecMode 1`.
  2. **The `HPHI_TRACE_BUF_MAX_MB` memory gate**: this existing cap
     (introduced in phase 3b for the Green-function result buffers) now
     *also* bounds the energy family's own per-rank Hamiltonian buffer.
     If the projected CSR size would exceed the cap, the energy family
     falls back to `ExpecMode 1`.

  Each rank-synchronized outcome is reported by a dedicated rank-0
  `INFO` line at the start of the run, distinct from the existing
  one-body/two-body `INFO` lines:

  ```
  INFO: ExpecMode 2: the energy/fluctuation family uses the trace kernel.
  INFO: ExpecMode 2: the energy/fluctuation family uses the ExpecMode-1 fallback (the Hamiltonian buffer would exceed HPHI_TRACE_BUF_MAX_MB).
  INFO: ExpecMode 2: the energy/fluctuation family uses the ExpecMode-1 fallback (the Hamiltonian was read from InputHam).
  ```

  The fixed always-fallback line changes accordingly — energy/fluctuation
  is removed from it, since it is no longer unconditionally on the
  `ExpecMode 1` path:

  ```
  INFO: ExpecMode 2: S2, NBodyG, and AnomalousG always use the ExpecMode-1 path in this version.
  ```

  (Previously this line also listed "energy/fluctuation,"; see
  `src/expec_trace.c`'s `TraceReportPlan()` for the exact strings, also
  quoted verbatim in the `ExpecMode` CalcMod documentation.)

- **Memory note.** The energy family's per-rank CSR Hamiltonian buffer
  is approximately `24·nnz` bytes (`colidx`: 8 bytes/entry as `long
  int`, `val`: 16 bytes/entry as `double complex`; `nnz` here is the
  number of entries makeHam emits on this rank, which the buffer
  capacity and the memory gate are sized for — the merged count
  `csr->nnz = rowptr[N]` after summing duplicates can be smaller;
  row-pointer and diagonal-coefficient overhead is comparatively
  negligible). This
  buffer **shares the existing `HPHI_TRACE_BUF_MAX_MB` cap**
  (default 1024 MiB, integer range `[1, 1048576]`) with the phase 3b
  Green-function result buffers — it is not a separate, additional
  budget. As with the existing cap, it is parsed from the environment
  on rank 0 only (for `ExpecMode 2` runs) and broadcast to every rank,
  so it only needs to be set in the launch environment.

- **The CSR is replicated per rank, by design.** Each rank builds the CSR
  for the *full* column range, so the per-rank CSR does **not** shrink as
  the MPI-rank count grows. This is required, not incidental: the
  state-panel layout gives each rank *complete* eigenvectors to evaluate
  without per-state communication, and computing `y = H·x` for a
  locally-complete `x` needs every row of `H` locally; a distributed CSR
  would reintroduce the per-state collectives that `ExpecMode 1`/`2`
  exist to remove. The CSR is nonetheless compact: it is sparse, with
  `nnz ≈ T·N` emitted entries (`T` = Hamiltonian terms per column,
  roughly linear in `N` for a fixed model; `T` grows only slowly with
  size, e.g. `T ≈ 21` for the Hubbard chain at `L=10`). Its dominant
  `colidx`/`val` storage is
  `≈ 24·nnz` bytes (plus `O(N)` row-pointer/work/diagonal arrays), about
  32 MB at `N = 63504`. At moderate rank counts this is much smaller than
  the `O(N²/P)` eigenvector state panel, but because the CSR is
  `P`-independent while the panel shrinks with `P`, the CSR's fixed cost
  becomes relatively more significant as `P` grows and can dominate at
  large `P` or for operator-dense inputs; there the per-rank
  `HPHI_TRACE_BUF_MAX_MB` gate falls the family back to `ExpecMode 1` and
  reports it. (This corrects an earlier troubleshooting note that wrongly
  suggested adding MPI ranks shrinks the CSR.)

- **Terminology note.** "Trace kernel" (the name used throughout this
  feature, in `ExpecMode 2`'s description, and in the source comments)
  refers to the precomputed-mapping/precomputed-CSR streaming
  evaluation technique introduced in phases 3b/3c — it does not refer
  to the matrix trace `Tr(·)`.

- **Guarantee unchanged**: `ExpecMode` still changes only evaluation
  *speed*. `0`, `1`, and `2` give identical physics (all columns,
  including `var`) up to floating-point rounding (summation order
  differs between kernels, so results are not bit-identical).

- **Purely additive**: no existing `ExpecMode` value's behavior changes
  for the models/configurations that were already using `ExpecMode 1`
  or the phase 3b Green-function kernels; no keyword is renamed,
  deprecated, or removed; no output file format changes. The only
  behavior change for existing runs is that energy-family evaluation on
  `ExpecMode 2` may now go through the new trace kernel instead of the
  `ExpecMode 1` path whenever the model is not reading `InputHam` and
  the CSR fits under `HPHI_TRACE_BUF_MAX_MB` — output values are
  unaffected within tolerance.

## Reviewer pointers

- CalcMod documentation: `doc/en/source/filespecification/expertmode_en/CalcMod_file_en.rst`
  and `doc/ja/source/filespecification/expertmode_ja/CalcMod_file_ja.rst`
  (`ExpecMode` entry, item `2`, energy-family paragraphs immediately
  following the phase 3b Green-function paragraphs).
- Appendix: `doc/en/source/technical/parallel_fulldiag_en.rst` and
  `doc/ja/source/technical/parallel_fulldiag_ja.rst` ("Trace-kernel
  evaluation (ExpecMode 2)" section; benchmark tables are unchanged in
  this phase and will be refreshed by Task 8).
- Dispatch plan / `INFO` reporting for the energy slot:
  `src/expec_trace.c` (`TraceBuildPlan()`, `TraceFinalizeEnergyPlan()`,
  `TraceReportPlan()`), `src/include/expec_trace.h` (`TRACE_Q_ENERGY`,
  `demoted_input_ham[]`).
- CSR Hamiltonian collector: `src/expec_trace_ham.c`,
  `src/include/expec_trace_ham.h` (`TraceHamCsr`, `TraceHamCollect()`,
  `TraceHamGatedBytes()`).
- MPI orchestration (rank-synchronized finalize, CSR lifetime):
  `src/phys_distributed.c`.
- Unit tests: `test/unit/expec_trace_ham_check.c` (dense-vs-CSR
  equivalence across every enabled family, memory-gate boundary,
  synchronized-demotion injection).
- Equivalence tests (production path, Modes 0/1/2, including the
  energy-kernel INFO line, the tJ orthogonality case, and the
  `InputHam` negative case): `test/fulldiag_expecmode_equiv.sh`.
- Manual hardware-verification checklist (ELPA/GPU): `test/manual/elpa_gpu_check.md`.
