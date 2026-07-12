# ExpecLocal call inventory (phase 3a, Task 1)

Status: frozen (Step 1-4 complete). This document is the interface contract
consumed by:
- Task 3 (ExpecLocal hooks + defensive guards): implements the local-mode
  behavior for every wrapperMPI/FileIO function listed below, inserts the
  `if (iExpecLocal) return ...;` defensive guards in `nbody_correlation.c` /
  `anomalous_pair.c`, wraps them in `EXPEC_LOCAL_GUARDED_BEGIN/END` markers,
  and empties `TEMP_UNGUARDED_FILES` in `test/check_expec_local_calls.sh`.
- Task 6 (Mode 1 driver split): adds `src/phys_distributed_local.c` to the
  `FILES` variable in `test/check_expec_local_calls.sh` once that file
  exists (NOT `src/phys_distributed.c`, which legitimately keeps raw MPI
  collectives and is permanently out of scope for this guard).

Regeneration: every call-site table below was produced by
`test/strip_c_comments.py <file>` (state-machine comment stripper, fails
loudly on empty output) piped through
`grep -noE "(SumMPI_[a-z]+|MaxMPI_[a-z]+|BcastMPI_[a-z]+|BarrierMPI|NormMPI_dc|VecProdMPI|MPI_[A-Za-z_]+|exitMPI|fopenMPI|childfopenMPI)\("`
against the current `feature/elpa-fulldiag` tree. This is the same
stripper/pattern used by `test/check_expec_local_calls.sh`, so the guard and
this document cannot silently drift apart.

## 1. Direct expec-reachable layer (guard-scanned files)

These six files are the "expec 到達層" -- the files FullDiag observable
evaluation calls directly, and the exact `FILES` set scanned by
`test/check_expec_local_calls.sh`.

| File | Calls found (count) | Classification | Notes |
|---|---|---|---|
| `src/expec_energy_flct.c` | `SumMPI_dc` x2, `SumMPI_d` x22 | **return-input** | All reductions; no raw MPI, no exitMPI, no file I/O. Local mode: every `SumMPI_*` becomes a no-op passthrough. |
| `src/expec_cisajs.c` | `SumMPI_dc` x7, `childfopenMPI` x1 | **return-input** / **local-open** | Reductions -> return-input. `childfopenMPI` (line 152, via `fopenMPI`) -> local-open (rank-0 gate removed, opens locally, NULL on failure). Also calls indirect-layer `child_*_MPI*` helpers (see §2) gated by `CheckPE()`. |
| `src/expec_cisajscktaltdc.c` | `SumMPI_dc` x16, `childfopenMPI` x4 | **return-input** / **local-open** | Same pattern as `expec_cisajs.c`. `childfopenMPI` at lines 198/204/212/220 (green-output partial channels, see design doc `GreenOutputSetPartialSuffix`). |
| `src/expec_totalspin.c` | `SumMPI_dc` x8 | **return-input** | Reductions only. Also calls indirect-layer `child_*_MPI*` helpers (see §2), gated by `isite > X->Def.Nsite`. |
| `src/nbody_correlation.c` | `SumMPI_dc` x24, `exitMPI` x18, `MPI_Sendrecv` x18, `childfopenMPI` x1 | **return-input** / **local-open** / **defensive-guard (required)** | `SumMPI_dc`/`childfopenMPI` as above. The 18 `MPI_Sendrecv` + 18 `exitMPI` (36 raw-MPI call sites) are the `*_partner_rank` site-swap paths audited in §3 below -- currently unreachable under replica FullDiag, but must get an explicit `if (iExpecLocal) return ...;` defensive guard (Task 3) since "should be unreachable" is not the same as "guaranteed never reachable by a future code path." Currently exempted from the guard script via `TEMP_UNGUARDED_FILES`. |
| `src/anomalous_pair.c` | `SumMPI_dc` x4, `exitMPI` x4, `MPI_Sendrecv` x4, `childfopenMPI` x1 | **return-input** / **local-open** / **defensive-guard (required)** | Same pattern as `nbody_correlation.c` (8 raw-MPI call sites: 4 `MPI_Sendrecv` + 4 `exitMPI`), same §3 audit, same `TEMP_UNGUARDED_FILES` exemption. |

**ExpecLocal allow-list (ties to `test/check_expec_local_calls.sh`
`ALLOW`)**: `SumMPI_dc SumMPI_d SumMPI_li SumMPI_i fopenMPI childfopenMPI
stdoutMPI`. `SumMPI_li` / `SumMPI_i` are not currently called from the six
files above but are included per the brief's frozen allow-list (harmless
headroom; the guard's detection regex only flags what actually appears).
`stdoutMPI` is used throughout as a bare `FILE*` macro argument to
`fprintf` (e.g. `expec_energy_flct.c:77`), never as `stdoutMPI(...)` --
it never matches the call-detection pattern and needs no local-mode
transformation beyond what the design doc §3 already specifies (rank-0-only
display, buffered elsewhere).

No `BcastMPI_*`, `MaxMPI_*`, `BarrierMPI`, `NormMPI_*`, or `VecProdMPI`
calls exist anywhere in the six direct files -- confirmed by the same
extraction. Per the design doc §3 matrix, since there is no current use,
Task 3 should implement these as "forbidden under `iExpecLocal` (debug
assert)" rather than inventing a no-communication semantics for an unused
path.

## 2. Indirect layer (child_*/GC_child_* definitions reached from §1)

Step 1 also requires extracting the same call inventory for every
`child_*`/`GC_child_*` function *definition* file reached from the six
direct files (not just the four `mltply*Core.c` files named in the task
brief's prose -- tracing the actual call graph surfaced two more files).
Method: every distinct `child_*`/`GC_child_*` identifier called from the
six files in §1 was located via
`grep -lE "^[A-Za-z_].*\b<name>\s*\(" src/*.c`, then each defining file was
run through the same stripper+grep extraction.

| Defining file | Calls found | Reachability from §1 | Classification |
|---|---|---|---|
| `src/mltplyHubbardCore.c` | (none: `SumMPI_*`/`MPI_*`/`exitMPI`/`fopenMPI` all absent) | N/A -- no MPI surface | **unreachable (no MPI to guard)** |
| `src/mltplySpinCore.c` | (none) | Defines `child_exchange_spin_element`, `child_Spin_CisAis`, `child_SpinGC_CisAis`, `child_GC_CisAit_spin_MPIdouble`, etc. called from `expec_cisajs.c`/`expec_cisajscktaltdc.c`/`expec_totalspin.c` | **unreachable (no MPI to guard)** |
| `src/mltplyMPIHubbardCore.c` | `MPI_Sendrecv` x18 (paired with `exitMPI` x18, i.e. 36 raw-MPI call sites), plus `CheckPE()` (site-boundary test, not itself an MPI call) | Defines `child_CisAisCjtAjt_Hubbard_MPI`, `child_CisAisCjtAku_Hubbard_MPI`, `child_CisAjtCkuAku_Hubbard_MPI`, `child_CisAjtCkuAlv_Hubbard_MPI` and their `GC_` variants, called from `expec_cisajscktaltdc.c` (e.g. line 860) | **unreachable under replica FullDiag** -- every call site is gated by `CheckPE(org_isiteN-1, X)==TRUE` for at least one site (`expec_cisajscktaltdc.c:853`); `CheckPE` (`mltplyMPIHubbardCore.c:38-46`) returns `TRUE` iff `org_isite+1 > X->Def.Nsite`, which is the identical `site >= Nsite` boundary audited in §3. Same conclusion, same "defensive guard is still required" caveat. |
| `src/mltplyMPISpinCore.c` | `MPI_Sendrecv` / `exitMPI` pairs throughout (spin analogue of the above) | Defines `child_GC_CisAis_spin_MPIdouble`, `child_GC_CisAisCjuAju_*_MPI{double,single}`, `child_GC_CisAitCiuAiv_spin_MPI*`, etc. called from `expec_cisajscktaltdc.c`/`expec_totalspin.c` | **unreachable under replica FullDiag** -- same `CheckPE`/site-boundary gating as `mltplyMPIHubbardCore.c`. |
| `src/mltplyMPIHubbard.c` | `MPI_Sendrecv` x10 / `exitMPI` x10 | Defines `child_general_hopp_MPIdouble`/`MPIsingle`, `child_GC_general_hopp_MPIdouble`/`MPIsingle`, called directly from `expec_cisajs.c` (lines 368-465) | **unreachable under replica FullDiag** -- call sites gated by `org_isite1 > X->Def.Nsite \|\| org_isite2 > X->Def.Nsite` (`expec_cisajs.c:443/460`); the callee itself (`mltplyMPIHubbard.c:104-136`, `child_GC_general_hopp_MPIdouble`) computes `origin = myrank ^ (mask1+mask2)` and always does `MPI_Sendrecv` with that `origin` once entered -- but under `Nsite==NsiteMPI` the caller-side gate means this function body is never entered at all. |
| `src/mltplyMPISpin.c` | `MPI_Sendrecv` x6 / `exitMPI` x6 | Defines `child_general_int_spin_MPIdouble`/`MPIsingle`, `child_general_int_spin_TotalS_MPIdouble`, called directly from `expec_totalspin.c` (lines 314/342/344) and `expec_cisajscktaltdc.c` (lines 1021/1041) | **unreachable under replica FullDiag** -- same `isite > X->Def.Nsite` caller-side gating (`expec_totalspin.c:315/341`). |

**Guard-scope decision for the indirect layer**: `test/check_expec_local_calls.sh`
`FILES` intentionally does **not** include the four `mltplyMPI*Core.c` /
`mltplyMPIHubbard.c` / `mltplyMPISpin.c` files, matching the design doc's
own framing of the check script's scope
(`docs/superpowers/specs/2026-07-11-elpa-fulldiag-phase3-design.md`, the
"将来の退行防止として" paragraph, which names only `expec_*.c` /
`nbody_correlation` / `anomalous_pair`). These indirect-layer files are
general-purpose MPI site-decomposition machinery shared with Lanczos/TPQ
(not FullDiag-specific), so guarding "any raw MPI call in this file"
would be nonsensical for their primary (non-FullDiag) callers. Their
FullDiag-reachability was instead fully audited above and rests on the
*same* `Nsite==NsiteMPI` invariant as §3 -- i.e. no separate proof burden.
If Task 3 later needs a defensive guard at this layer too (rather than
relying solely on the caller-side `CheckPE`/`isite>Nsite` gates plus the
`iExpecLocal` guards in `nbody_correlation.c`/`anomalous_pair.c`), that is
a Task 3 design decision informed by this table, not a gap this guard
script silently hides -- the reachability chain is fully written down here.

## 3. `partner_rank` audit (Step 2 -- proof obligation)

**Claim to verify**: in replicated FullDiag mode (`iFlgScaLAPACK=1`, no
MPI site-separation, i.e. `Nsite == NsiteMPI`), every `*_partner_rank`
function in `nbody_correlation.c` and `anomalous_pair.c` always returns
`partner_rank == myrank`, making the `partner != myrank` / `origin !=
myrank` raw-MPI branches (`MPI_Sendrecv` + `exitMPI`) unreachable.

**Structural precondition -- `Nsite == NsiteMPI` under `iFlgScaLAPACK=1`**:
- `src/readdef.c:310-320` (comment, contemporaneous with the flag's
  derivation): *"iFlgScaLAPACK now doubles as the internal
  'distributed-eigenvector FullDiag' flag ... check.c, and CheckMPI.c
  (iFlgScaLAPACK==1 -> replicated Hilbert-space treatment, NsiteMPI=Nsite,
  no site separation)."*
- `src/readdef.c:321`: `X->iFlgScaLAPACK = (X->iSolver == SOLVER_SCALAPACK
  || X->iSolver == SOLVER_ELPA) ? 1 : 0;`
- `src/check.c:99-109`: structural enforcement --
  ```c
  if(X->Def.iFlgScaLAPACK == 0) {
    if (CheckMPI(X) != TRUE) { return MPIFALSE; }   // may set Nsite < NsiteMPI
  }
  else{
    X->Def.NsiteMPI = X->Def.Nsite;                  // Nsite == NsiteMPI, no site split
    X->Def.Total2SzMPI = X->Def.Total2Sz;
  }
  ```
  When `iFlgScaLAPACK==1`, `CheckMPI()` (the function that assigns
  inter-process/local site counts and can make `Nsite < NsiteMPI`,
  `src/CheckMPI.c:69-105`) is **never called**; `NsiteMPI` is set equal to
  the already-local `Nsite` instead. So `Nsite == NsiteMPI` holds by
  construction for every replicated-FullDiag run, independent of `nproc`.
  Consequently every site index in the problem satisfies `site < Nsite`
  (no site can be `>= Nsite`, since `Nsite` equals the total site count).

**Per-function proof, `nbody_correlation.c`**:
- `apply_hubbardgc_rank_annihilate` / `apply_hubbardgc_rank_create`
  (defined at `nbody_correlation.c:569` / `583`; guard line at `:576` /
  `:590`, both `if (site < D->Nsite) return 1;`): under `Nsite==NsiteMPI`
  this guard is true for every site
  index that can occur, so `rank_state` is **never modified** by any call
  -> `apply_hubbardgc_rank_part` (line 619) returns `*partner_rank =
  current_rank` unchanged -> `nbodyg_hubbardgc_partner_rank` (line 653)
  always yields `partner_rank == current_rank == myrank`.
- `apply_spinless_rank_annihilate` / `apply_spinless_rank_create`
  (`nbody_correlation.c:721,734`, identical `if (site < D->Nsite) return
  1;` guard): same argument -> `nbodyg_spinless_partner_rank` (line 795)
  always yields `partner_rank == myrank`.
- `nbodyg_general_spin_partner_rank` (`nbody_correlation.c:840-899`): loop
  body guarded by `if (site < X->Def.Nsite) continue;` (line 862) for
  every factor -- under `Nsite==NsiteMPI` this is true for every factor,
  so the loop body that would modify `partner`/`side` never executes ->
  `*partner_rank = (int)partner` (line 898) returns the unmodified
  `current_rank == myrank`.
- Call sites (e.g. `nbody_correlation.c:1106-1135`): `int origin = myrank;
  ... if (origin == myrank) { <local path> } else { MPI_Sendrecv(...);
  ... exitMPI(...) on error }` -- since `origin` is proven `== myrank`
  above, the `else` (raw-MPI) branch is unreachable in replica FullDiag.

**Per-function proof, `anomalous_pair.c`**:
- `apply_anomalous_rank_annihilate` / `apply_anomalous_rank_create`
  (defined at `anomalous_pair.c:222` / `236`; guard line at `:229` /
  `:243`) use the same `if (site < D->Nsite) return 1;` guard, feeding
  `apply_anomalous_rank_term` (line 249) and then
  `anomalous_hubbardgc_partner_rank` (line 272) ->
  `*partner_rank = current_rank == myrank`.
- Call sites (`anomalous_pair.c:376-403`, `489-519`): `int partner =
  myrank; ... if (partner == myrank) { <local path> } else {
  MPI_Sendrecv(...); }` -- unreachable else-branch by the same argument.

**Conclusion**: the "cannot prove unreachable" fallback in the task brief
is **not needed** -- reachability is disproven with concrete, file/line
evidence above, for all three model families (Hubbard-type NBodyG,
spinless NBodyG, general-spin NBodyG) and both files. The invariant rests
on one structural fact (`check.c:99-109`, `Nsite==NsiteMPI` whenever
`iFlgScaLAPACK==1`) plus purely mechanical site-index guards in each
`*_rank_*` helper -- no floating-point or runtime-only reasoning involved.

**Even so, a defensive guard is mandatory (Task 3), not optional**: proof
of current unreachability is not a substitute for a runtime safety net,
because (a) the invariant depends on `check.c`'s branch structure staying
intact under future refactors, (b) `iExpecLocal` may in principle be
entered from a code path that does not go through `check.c`'s current
gating in the future, and (c) a single-rank `exitMPI`/`MPI_Sendrecv` call
that *does* fire under `iExpecLocal` would hang the other ranks with no
diagnostic. Task 3 must therefore add `if (iExpecLocal) return
<error-rc>;` immediately before each `partner != myrank` /
`origin != myrank` branch in both files (numbered call sites: 5 in
`nbody_correlation.c` -- `nbodyg_hubbardgc_partner_rank`/spinless/general-spin
x2 orientations each plus the shared helpers -- and 2 in
`anomalous_pair.c`), and wrap the guarded raw-MPI region with
`/* EXPEC_LOCAL_GUARDED_BEGIN */` / `/* EXPEC_LOCAL_GUARDED_END */` so
`test/check_expec_local_calls.sh` can drop `nbody_correlation.c` and
`anomalous_pair.c` from `TEMP_UNGUARDED_FILES`.

## 4. Current unguarded raw-MPI count (Step 4)

Running `test/check_expec_local_calls.sh` today (before Task 3) with an
**empty** `TEMP_UNGUARDED_FILES` reports the following disallowed
call sites (none are `EXPEC_LOCAL_GUARDED`-marked yet, since Task 3 has
not run):

| File | `MPI_Sendrecv(` | `exitMPI(` | Total raw-MPI sites |
|---|---:|---:|---:|
| `src/nbody_correlation.c` | 18 | 18 | 36 |
| `src/anomalous_pair.c` | 4 | 4 | 8 |
| **Total** | **22** | **22** | **44** |

Per the task brief, the exact count is informational (this table is the
source of truth, not an externally-fixed acceptance number). With the
current `TEMP_UNGUARDED_FILES="src/nbody_correlation.c
src/anomalous_pair.c"`, `check_expec_local_calls` skips these two files
entirely and **PASSES**; `ctest -R fulldiag` remains 17/17 (verified in
`build_noMPI`). Task 3 must empty `TEMP_UNGUARDED_FILES` once every one of
these 44 sites is behind an `if (iExpecLocal) return ...;` guard wrapped in
`EXPEC_LOCAL_GUARDED_BEGIN/END` markers (per §3's conclusion, all 44 are
provably unreachable in replica FullDiag mode today, but still require the
runtime safety net).

## 5. Guard mechanism summary

`test/check_expec_local_calls.sh`:
1. Strips comments from each file in `FILES` via `test/strip_c_comments.py`
   (fails loudly on a stripper error or empty output).
2. Pass 1 (on the **original**, pre-strip source): collects line-number
   ranges between `/* EXPEC_LOCAL_GUARDED_BEGIN */` and
   `/* EXPEC_LOCAL_GUARDED_END */` markers. The stripper preserves every
   newline (including inside `/* ... */` blocks), so stripped-file line
   numbers line up exactly with original-file line numbers -- no
   hardcoded line numbers anywhere in the script.
3. Pass 2 (on the **stripped** source): greps for the same token pattern
   used to build this document's tables. Any match whose line falls
   inside a guarded range, or whose captured name is on the `ALLOW` list,
   is not a violation. Anything else is.
4. Files listed in `TEMP_UNGUARDED_FILES` skip step 2/3 entirely (the
   stripper still runs, so a broken stripper is caught even for these
   files).
5. Exits 1 with a `VIOLATION: <file>:<line>: ...` message per finding if
   anything remains; otherwise prints `PASSED` and exits 0.

Verified by manual injection during Task 1 development (not committed):
inserting an unmarked `MPI_Barrier(MPI_COMM_WORLD)` call into
`src/expec_totalspin.c` made the guard fail with a `VIOLATION` pointing at
the injected line; wrapping the same call in
`EXPEC_LOCAL_GUARDED_BEGIN/END` markers made the guard pass again. The
injected code was reverted before committing (confirmed via `git diff`
showing no residual changes to `src/expec_totalspin.c`).

Registered as ctest `check_expec_local_calls` via
`add_hphi_test_with_srcdir(check_expec_local_calls)` in
`test/CMakeLists.txt`, next to the other validation-style
srcdir-parameterized tests (e.g. `green_output_format`).
