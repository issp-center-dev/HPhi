# PR migration-note draft: phase 3a (ExpecMode / state-parallel FullDiag observables)

Status: draft, to be pasted into PR #276's description (the "Migration
notes" section) once phase 3a lands. Written to match the tone of that
section's existing bullets (keyword additions/deprecations, behavior
changes with attribution).

---

## Migration notes (phase 3a)

- **New keyword: `ExpecMode` (CalcMod, default `0`).** Selects the
  evaluation kernel for full-diagonalization (`CalcType` = FullDiag)
  observables (energy, N, Sz, S2, doublon, Green functions):
  - `0`: conventional evaluation (existing behavior, unchanged).
  - `1`: state-task-parallel evaluation -- each MPI rank evaluates
    observables for its own contiguous block of eigenstates
    independently, with no communication during the per-state loop.
    Aggregate Green-function output is written as rank-local partial
    files and merged into the final aggregate files by rank 0 once every
    rank's manifest reports success.
  - `2`: reserved for a trace-kernel evaluation mode (planned for phase
    3b). Not yet implemented; currently runs as `ExpecMode 1` and prints
    an `INFO` line saying so.
  - Eligibility: nonzero `ExpecMode` requires `CalcType` = FullDiag with
    `Solver` 1 (ScaLAPACK) or 3 (ELPA); any other combination is rejected
    at startup. With exactly one MPI process, `ExpecMode` is
    automatically reverted to `0` (an `INFO` line is printed) since
    results are identical for a single process either way.
  - Guarantee: `ExpecMode` changes only evaluation *speed* -- `0`, `1`,
    and `2` give identical physics up to floating-point rounding
    (summation order differs between kernels, so results are not
    bit-identical).

- **Behavior fix (attributed to phase 3a): distributed FullDiag `S2`/`Sz`
  are no longer zero-filled.** Before this phase, FullDiag runs with
  `Solver` 1 (ScaLAPACK) or 3 (ELPA) and more than one MPI process
  (`ExpecMode 0`, the only mode that previously existed) skipped the S2/Sz
  calculation entirely -- `S2` and `Sz` were reported as `0`, and the
  stdout progress line used a shortened format that omitted the `S2`
  column. As of phase 3a, these distributed runs compute `S2`/`Sz` on rank
  0 (same as the serial path), and the stdout progress line always uses
  the single-process (serial) format, with the `S2` column present, for
  every `Solver`/`ExpecMode` combination. This is an intentional
  correctness fix, independent of `ExpecMode`'s value -- it also applies
  to the (unchanged, default) `ExpecMode 0`. Anyone diffing FullDiag
  stdout/`zvo_phys*` output against a pre-phase-3a distributed run should
  expect `S2`/`Sz` to change from `0` to their correct values, and the
  progress-line column layout to change from the shortened form to the
  serial form.

- No existing keyword is deprecated or removed in this phase. `ExpecMode`
  is purely additive and defaults to the pre-existing behavior (`0`)
  except for the `S2`/`Sz` fix above, which applies unconditionally to
  distributed FullDiag runs regardless of `ExpecMode`.

## Reviewer pointers

- CalcMod documentation: `doc/en/source/filespecification/expertmode_en/CalcMod_file_en.rst`
  and `doc/ja/source/filespecification/expertmode_ja/CalcMod_file_ja.rst`
  (`ExpecMode` entry, immediately after `NGPU`).
- Validation / demotion: `src/readdef.c` (`ReadcalcmodFile`, ExpecMode
  range + eligibility + nproc==1 demotion), `src/ErrorMessage.c`
  (`cErrExpecMode`).
- Mode dispatch / trace-kernel downgrade INFO: `src/phys.c` (top of
  `phys()`, `#ifdef _SCALAPACK` block).
- Mode 1 driver: `src/phys_distributed.c` (MPI orchestration) and
  `src/phys_distributed_local.c` (MPI-free per-rank state loop).
- Mode 0 S2/Sz unification: `src/phys.c` (the `use_scalapack` FullDiag
  branch, `ExpecLocalEnter()`/`expec_totalspin()`/`ExpecLocalLeave()`
  block replacing the old zero-fill).
- Equivalence tests: `test/fulldiag_expecmode_equiv.sh`
  (`fulldiag_expecmode_equiv_np2`/`_np3` ctest cases); merge semantics:
  `test/unit/green_partial_merge_check.c`.
- Manual hardware-verification checklist (ELPA/GPU): `test/manual/elpa_gpu_check.md`,
  items 10-14.
