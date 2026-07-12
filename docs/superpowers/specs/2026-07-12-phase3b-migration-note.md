# PR migration-note draft: phase 3b (ExpecMode 2 trace kernels)

Status: draft, to be pasted into the PR description (the "Migration
notes" section) once phase 3b lands, appended after the phase 3a
migration-note bullets (`docs/superpowers/specs/2026-07-11-phase3a-migration-note.md`).
Written to match that note's tone (keyword additions/deprecations,
behavior changes with attribution). Usage guidance below is deliberately
measurement-agnostic; Task 8's clavius benchmark gate (Mode 0/1/2,
N=4900, np=4) is the authority for any final numeric speedup claim in
the PR text -- if that run shows Mode 2 slower than Mode 1 for some
case, update the guidance bullet accordingly before merging.

---

## Migration notes (phase 3b)

- **`ExpecMode 2` is now a real evaluation kernel** (previously: reserved,
  ran as `ExpecMode 1` with a downgrade `INFO` line). It applies **only**
  to the one-body (`expec_cisajs`-equivalent) and two-body
  (`expec_cisajscktaltdc`-equivalent) Green functions: for each operator,
  HPhi precomputes the basis-state mapping (destination state and
  amplitude) once, then streams every owned eigenstate through that
  mapping in a dense loop, amortizing the per-state operator-dispatch
  overhead that `ExpecMode 1` still pays for these two quantities.
  Energy/fluctuation (including the `var` column), `S2`, `NBodyG`, and
  `AnomalousG` are **not** kernelized in this phase -- they always run on
  the `ExpecMode 1` path, unconditionally. (Rationale, for reviewers:
  `var` is an eigenvector-quality check, so substituting the eigenvalue
  would violate the "`ExpecMode` changes only speed" guarantee; the
  energy family and diagonal quantities share one evaluator
  [`expec_energy_flct`], so partial kernelization would create a
  double-writer. Trace-kernelizing the energy family is left as a
  candidate for a future phase.)

- **Capability table (model support).** The trace kernel is validated
  and enabled for exactly four (model, spin-representation) rows:
  `Hubbard`, `HubbardGC`, `Spin` (half-integer only), `SpinGC`
  (half-integer only). Any other model (`tJ`/`tJGC`, `Kondo`/`KondoGC`,
  general-spin, Spinless) always falls back to `ExpecMode 1` for both
  Green-function quantities ("unsupported model").

- **Three further per-quantity runtime fallbacks** (checked in a fixed
  order — shared-evaluator, then no-operators, then the result-buffer
  cap — and mutually exclusive with each other and with the
  unsupported-model case) apply even on a supported model:
  - **No-operators rule**: a quantity with no operators of its kind
    defined in the input falls back trivially (INFO: "no operators of
    this kind are defined") — nothing is computed or written either way.
  - **Memory gate**: if a quantity's result buffer would exceed the new
    `HPHI_TRACE_BUF_MAX_MB` cap, that quantity falls back to
    `ExpecMode 1`.
  - **Shared-evaluator rule (two-body only)**: the two-body Green
    function shares its evaluator with the ThreeBodyG/FourBodyG/SixBodyG
    (N-body) Green functions. Whenever any N-body GF is requested, the
    two-body GF falls back together with it, so the N-body evaluator
    still runs (via the `ExpecMode 1` path) instead of being silently
    skipped, and the two-body output is never written twice. The
    one-body quantity is unaffected by this rule.

  Each quantity's outcome (trace kernel / one of the four fallback
  reasons) is reported by a rank-0 `INFO` line at the start of the run;
  see `src/expec_trace.c`'s `TraceReportPlan()` for the exact strings
  (also quoted verbatim in the `ExpecMode` CalcMod documentation).

- **New environment variable: `HPHI_TRACE_BUF_MAX_MB`** (integer MiB,
  range `[1, 1048576]`, default `1024`). Caps the trace kernel's
  per-rank, per-quantity result buffer only (not the state panel or any
  other allocation). Read from the environment on rank 0 only and
  broadcast, so it is sufficient to set it in the launch environment
  (e.g. the job script) rather than on every rank.

- **The `var` column is unchanged.** Because the energy/fluctuation
  family always runs on the `ExpecMode 1` path in this phase (see
  above), `var` is computed identically under `ExpecMode 0`, `1`, and
  `2` -- there is no new code path for it to diverge through.

- **Guarantee unchanged**: `ExpecMode` still changes only evaluation
  *speed*. `0`, `1`, and `2` give identical physics (all columns,
  including `var`) up to floating-point rounding (1e-8; summation order
  differs between kernels, so results are not bit-identical).

- **Usage guidance update**: `ExpecMode 2` is designed to further amortize
  per-state operator overhead for one-body/two-body Green-function-heavy
  FullDiag workloads on the supported models above, and is expected to be
  at least as fast as `ExpecMode 1` for such correlation-function-heavy
  workloads. See `test/manual/elpa_gpu_check.md` for measured results
  (the phase 3b benchmark gate compares Mode 0/1/2 wall time with a
  mapping-extraction/streaming/output breakdown). This phase does not
  promise a specific speedup number in the documentation; the benchmark
  results are the reference for anyone tuning `ExpecMode` choice.

- **Purely additive**: no existing `ExpecMode` value's behavior changes
  (`0` and `1` are untouched), no keyword is renamed, deprecated, or
  removed, and no output file format changes. The only new environment
  variable is `HPHI_TRACE_BUF_MAX_MB`, which defaults to a value (1024
  MiB) large enough that no existing FullDiag run should hit the gate
  unless it opts into a very large Green-function output.

## Reviewer pointers

- CalcMod documentation: `doc/en/source/filespecification/expertmode_en/CalcMod_file_en.rst`
  and `doc/ja/source/filespecification/expertmode_ja/CalcMod_file_ja.rst`
  (`ExpecMode` entry, item `2`, immediately after the phase 3a `1` bullet).
- Trace-kernel plan construction, capability table, and `INFO` reporting:
  `src/expec_trace.c` (`TraceBuildPlan()`, `kTraceCap`, `TraceReportPlan()`,
  `TraceGbufMaxBytesFromEnv()`).
- Trace-kernel internals (mapping-probe adapters, streaming, output):
  `src/expec_trace_internal.h`, `src/expec_trace.c`
  (`TraceStreamOneBody()`/`TraceStreamTwoBody()` and their output-phase
  counterparts).
- Capability-table evidence (why each of the four rows is `TRUE`, and
  which golden tests / unit tests verify each): the per-row comment block
  above `kTraceCap` in `src/expec_trace.c`, and
  `docs/superpowers/specs/2026-07-11-expec-call-inventory.md` §2c.
- `HPHI_TRACE_BUF_MAX_MB` parsing/broadcast: `src/expec_trace.c`
  (`TraceGbufMaxBytesFromEnv()`, rank-0-only `getenv`) and
  `src/phys_distributed.c` (the single MPI touch: broadcasting the parsed
  cap to every rank).
- Unit tests (mapping correctness, purity, one-body/two-body streaming,
  memory-gate boundaries): `test/unit/expec_trace_map_check.c`
  (`expec_trace_map_check` ctest).
- Equivalence tests (production path, mode 0 vs 1 vs 2, including the
  shared-evaluator and unsupported-model fallback assertions):
  `test/fulldiag_expecmode_equiv.sh` (`fulldiag_expecmode_equiv_np2`/
  `_np3` ctest cases; see `assert_kernel_plan()` and
  `assert_kernel_plan_shared_evaluator()`).
- Manual hardware-verification checklist (ELPA/GPU): `test/manual/elpa_gpu_check.md`,
  phase 3b items (map-check unit test, production-path equivalence,
  benchmark gate, optional GPU point).
