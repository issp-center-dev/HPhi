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

## 2b. `mltply()` subtree (addendum, final whole-branch review)

Step 1/2 above traced every `child_*`/`GC_child_*` identifier called
*directly* from the six §1 files. That tracing missed one call that goes
through a completely different entry point: `src/expec_energy_flct.c:195`
calls `mltply(X, v0, v1)` (`v0 += H*v1`, needed for the energy/fluctuation
observable), and `mltply()` is the same top-level Hamiltonian-multiply
dispatcher used by Lanczos/TPQ -- it fans out into a batched-MPI subtree
that §1/§2 never enumerated. This addendum closes that gap; it does not
change any guard behavior (`test/check_expec_local_calls.sh`'s `FILES` set
is unchanged -- `mltply.c`/`mltplyHubbard.c`/`mltplySpin.c`/
`mltplySpinless.c`/`mltplyMPIBatched.c` are general-purpose
site-decomposition machinery shared with Lanczos/TPQ, matching the same
"shared machinery, not `expec_*`-specific" reasoning §2 already applied to
`mltplyMPI*Core.c`), it only records the reachability argument for the
record.

| File | Role reached from `expec_energy_flct.c:195`'s `mltply()` call | MPI surface | Reachability under Mode 1 (replica FullDiag) |
|---|---|---|---|
| `src/mltply.c` | Top dispatcher: switches on `X->Def.iCalcModel` to `mltplyHubbardGC`/`mltplyHubbard`/`mltplySpin`/`mltplySpinGC`/`mltplySpinlessFermion`, then reduces `X->Large.prdct` | `SumMPI_dc` at `mltply.c:155` (`X->Large.prdct = SumMPI_dc(X->Large.prdct);`) | **reachable, but already covered**: `SumMPI_dc` is on the frozen ExpecLocal `ALLOW` list (§1) and is hooked by Task 3's ExpecLocal no-op-passthrough semantics like every other `SumMPI_dc` call site. |
| `src/mltplyHubbard.c` | `mltplyHubbardGC` (defined `:436`) calls the batched-InterAll init `InitializeMPIBatchedInterAll_HubbardGC` (`:567`) plus batched Transfer/DoubleTransfer inits (`:478`, `:495`); canonical `mltplyHubbard` (defined `:169`) calls the canonical-model batched inits (`:209`, `:226`). Both also keep a non-batched per-term `CheckPE`-gated fallback loop (e.g. canonical InterAll loop `:293-355`). | `MPI_Comm_size` (in the callees, no data exchange) plus the raw `MPI_Sendrecv` machinery reached transitively through `mltplyMPIBatched.c` (see below) | **unreachable under replica FullDiag** -- see reachability argument below. |
| `src/mltplySpin.c` | Canonical `Spin` Exchange path calls `InitializeMPIBatchedExchange_Spin` and its group loop (`:275`, `:283-286`); `SpinGC` Exchange path calls `InitializeMPIBatchedExchange_SpinGC` and its group loop (`:788`, similar shape). Both also keep a non-batched per-term fallback (`HPHI_MPI_NOBATCH` branch, e.g. `:314-330`). | same shape as `mltplyHubbard.c` | **unreachable under replica FullDiag** -- same argument. |
| `src/mltplySpinless.c` | `mltplySpinlessFermion` (defined `:115`) calls `InitializeMPIBatchedTransfers_SpinlessFermionGC` (`:162`) / `InitializeMPIBatchedTransfers_SpinlessFermion` (`:180`) depending on GC vs canonical model | same shape | **unreachable under replica FullDiag** -- same argument. |
| `src/mltplyMPIBatched.c` | Defines every `Initialize*`/`X_child_*_MPI*_batched` function named above. The real raw-MPI calls (`MPI_Sendrecv`, e.g. `:361-366`, `:495-505`, `:847-852`, and further pairs at the InterAll/Exchange variants) live *inside* the `X_child_*_batched` functions, which are only invoked from each caller's `for (g = 0; g < batched->num_groups; g++)` loop. | `MPI_Sendrecv` (paired sends/receives, batched-group communication) | **unreachable under replica FullDiag** -- every group loop above is driven by `batched->num_groups`, which the reachability argument below shows is always 0. |

**Reachability argument** (mirrors the `Nsite==NsiteMPI` invariant already
proved in §3, applied to this subtree instead of `nbody_correlation.c`/
`anomalous_pair.c`):

- Every `Initialize*` function in `mltplyMPIBatched.c` decides whether a
  term needs cross-rank batching using one of two equivalent site-locality
  tests: `CheckPE(site, X)` (InterAll-type terms -- e.g.
  `InitializeMPIBatchedInterAll_HubbardGC`'s `ComputeInterAllOrigin()`
  helper, `mltplyMPIBatched.c:1372-1489`, sets `any_interPE` only if
  `CheckPE()` is `TRUE` for at least one of the four sites) or an explicit
  `site + 1 > X->Def.Nsite` inequality (Transfer/Exchange-type terms -- e.g.
  `InitializeMPIBatchedTransfers_SpinlessFermionGC`'s `site1_local !=
  site2_local` check, `mltplyMPIBatched.c:90-99`, and
  `InitializeMPIBatchedExchange_Spin`'s `(site0+1 > Nsite) != (site1+1 >
  Nsite)` check, `mltplyMPIBatched.c:2483-2487`). `CheckPE()` itself
  (`mltplyMPIHubbardCore.c:38-46`, already audited in §2) returns `TRUE`
  iff `org_isite+1 > X->Def.Nsite` -- the identical condition.
- ExpecMode/Mode 1 is only reachable at all when Solver is 1 (ScaLAPACK) or
  3 (ELPA) (`src/readdef.c`'s `cErrExpecMode` gate), and both solvers set
  `X->Def.iFlgScaLAPACK = 1` (`readdef.c:321`), which forces
  `NsiteMPI = Nsite` with no site separation (`check.c:99-109`, the same
  structural fact §3 already establishes). So every site index satisfies
  `site < Nsite` whenever Mode 1 can run at all -- `CheckPE(site,X)` is
  `FALSE` and `site+1 > Nsite` is `FALSE` for every site, in every one of
  the functions above.
- Consequently: `any_interPE` never becomes `TRUE` (`ComputeInterAllOrigin`
  always takes the `return -1;` "all sites are local" path), and
  `site1_local != site2_local` / `(site0+1>Nsite) != (site1+1>Nsite)` are
  never true. Every `Initialize*` function's unique-origin counter
  (`num_unique` / `total_terms`) therefore stays `0`, so `batched->num_groups
  == 0` in every case, and the `for (g = 0; g < num_groups; g++)` loops that
  would call the `X_child_*_MPI*_batched` functions (where the actual
  `MPI_Sendrecv` calls live) always execute zero times. The non-batched
  per-term fallback branches (e.g. `mltplyHubbard.c`'s canonical InterAll
  loop) are gated by the same `CheckPE`/`site>Nsite` condition and are
  unreachable for the identical reason.
- This rests on the same structural invariant as §3 (`Nsite==NsiteMPI`
  whenever `iFlgScaLAPACK==1`), so it carries no separate proof burden --
  it is the same argument, applied to a different call path that Step 1/2
  did not originally walk.

**Guard-scope note**: no change to `test/check_expec_local_calls.sh`'s
`FILES` is implied or needed by this addendum -- these five files are
general-purpose Hamiltonian-multiply machinery shared with Lanczos/TPQ
(not `expec_*`-specific), exactly like the `mltplyMPI*Core.c` files §2
already excludes from the guard's scope for the same reason.

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

## 2c. Mapping-probe family audit (phase 3b, Task 2)

Per-family audit of every element-function branch the LOCAL (intra-process)
path of `expec_cisajs.c` / `expec_cisajscktaltdc.c` can reach for the four
trace-kernel candidate models (Hubbard, HubbardGC, Spin-half, SpinGC-half) x
two quantities (one-body / two-body). This is the precondition list for Task 5
capability TRUE-ing. Line numbers are as of this commit.

Common facts referenced below:

- **kprime / vector-index invariant.** Every extracted `*_map` core reports a
  0-based destination `kprime`; the pre-refactor original always read/wrote its
  result vector at `kprime+1`. Canonical GetOffComp yields a 1-based `off`, so
  `kprime = off-1`; the grand-canonical bare-bit `tmp_off` is 0-based, so
  `kprime = tmp_off`. Diagonal families use `kprime = j-1`. So `kprime+1` is the
  slot the Mode-1 code touched in every family.
- **M_CORR purity.** In M_CORR the `tmp_v0[...] += ...` write is gated by
  `X->Large.mode == M_MLTPLY || == M_CALCSPEC`, so it never fires; only the
  returned `dam_pr` (Hubbard core) / the `else` branch reading `tmp_v1`
  (SpinGC-half core) is live. The SpinGC-half element funcs additionally have an
  `H_CORR` branch reading `conj(tmp_v0[...])`; extraction uses M_CORR only, so
  that branch is dead too. The probes never pass a `tmp_v0` at all.
- **Reachable-helper write-set (column vi), used by the purity test's snapshot
  list.** The only mutable target any reachable helper writes is `X->Large`
  (via `general_hopp_GetInfo` / `general_int_GetInfo` and, for the
  Spin/SpinGC-half two-body, `Rearray_Interactions` which sets `X->Large.tmp_V`
  through GetInfo — but the two-body Spin/SpinGC path computes `isA_up/isB_up`
  directly from `X->Def.Tpow` and does NOT call `general_int_GetInfo`, so for
  those `X->Large` is touched only by the driver's own field writes
  `i_max/irght/ilft/ihfbit/mode`). `GetOffComp` reads `list_2_1`/`list_2_2`;
  `SgnBit` is pure; `child_*` read `list_1` / bare bits and `X->Large.is*`.
  **No reachable helper writes any global array** (`list_1`/`list_2_*` are
  read-only). The extraction driver additionally snapshots and restores the
  ENTIRE `X->Large` on entry/exit, so post-extraction `X` is byte-identical.
  Snapshot fields verified by the purity test: `X->Large.{mode, i_max, irght,
  ilft, ihfbit, is1_spin, is2_spin, is3_spin, is4_spin, isA_spin, isB_spin,
  A_spin, B_spin, is1_up, is1_down, is2_up, is2_down, tmp_V, tmp_J, isite1,
  isite2, isite3, isite4}`.

### 2c.1 One-body (cisajs)

| (i) function @ location | (ii) kprime source | (iii) amplitude | (iv) M_CORR purity evidence | (v) adapter | (vi) reachable-helper write-set |
|---|---|---|---|---|---|
| `GC_CisAis` @ mltplyHubbardCore.c (diagonal, HubbardGC) | diagonal `j-1` | occupation 0/1 (implicit coupling 1.0) | write gated `M_MLTPLY\|M_CALCSPEC` (in-func `if`) | `GC_CisAis_TraceProbe` | `general_hopp_GetInfo`->`X->Large.is*/A_spin`; probe pure |
| `GC_CisAjt` @ mltplyHubbardCore.c (off-diag, HubbardGC) | bare-bit `tmp_off` (`list_1_j^sum`) | Fermion sign `SgnBit`, coupling 1.0 | write gated `M_MLTPLY\|M_CALCSPEC` | `GC_CisAjt_TraceProbe` | `general_hopp_GetInfo`; `SgnBit` pure |
| `CisAjt` @ mltplyHubbardCore.c (off-diag, canonical Hubbard) | canonical `off-1` via `GetOffComp` | Fermion sign, coupling 1.0 | write gated `M_MLTPLY\|M_CALCSPEC` | `CisAjt_TraceProbe` | `general_hopp_GetInfo`; `GetOffComp` reads `list_2_*`; `SgnBit` pure |
| canonical diagonal (inline `list_1[j]&is`, expec_cisajs.c:482-487) | diagonal `j-1` | occupation 0/1 | pure read only (no element func) | driver-inline in `TraceMapExtractOneBody` | reads `list_1`,`X->Def.Tpow` |
| `child_Spin_CisAis` @ mltplySpinCore.c:210 (diagonal, Spin-half) | diagonal `j-1` | 0/1 spin match | pure function (no vector) | `child_Spin_CisAis_TraceProbe` | reads `list_1`; pure |
| `child_SpinGC_CisAis` @ mltplySpinCore.c:227 (diagonal, SpinGC-half) | diagonal `j-1` | 0/1 spin match | pure function | `child_SpinGC_CisAis_TraceProbe` | pure (bare bit) |
| `child_SpinGC_CisAit` @ mltplySpinCore.c:247 (transverse, SpinGC-half) | bare-bit `tmp_off` | flip sign +1 | pure function | `child_SpinGC_CisAit_TraceProbe` | pure (bare bit) |

Zero-result one-body (write an empty map -> GF 0, matching Mode 1): canonical
Hubbard cross-spin under `iFlgSzConserved` (expec_cisajs.c:427-433); Kondo
localized-vs-itinerant pair (:436-446); Spin/SpinGC off-diagonal
`org_isite1 != org_isite2` (expec_cisajs.c:567-569 / :721-724).

### 2c.2 Two-body (cisajscktaltdc)

Hubbard / HubbardGC do **not** call `Rearray_Interactions`; they call
`general_int_GetInfo` with a fixed `tmp_V = 1.0` (expec_cisajscktaltdc.c:857,916
canonical; :725,777 GC), so amplitude = `1.0 * tmp_sgn`. Branch selection is by
`isite1==isite2` / `isite3==isite4` (the same four-way `is*` test the Mode-1
dispatch uses).

| (i) function @ location | (ii) kprime source | (iii) amplitude | (iv) M_CORR purity | (v) adapter |
|---|---|---|---|---|
| `CisAisCisAis_element` (canonical, diag) | `j-1` | `tmp_V*tmp_sgn` | gated write | `CisAisCisAis_element_TraceProbe` |
| `CisAisCjtAku_element` (canonical) | `child_CisAjt` off `-1` | `tmp_V*tmp_sgn` | gated write | `CisAisCjtAku_element_TraceProbe` |
| `CisAjtCkuAku_element` (canonical) | `child_CisAjt` off `-1` | `tmp_V*tmp_sgn` | gated write | `CisAjtCkuAku_element_TraceProbe` |
| `CisAjtCkuAlv_element` (canonical) | `child_CisAjt` off `-1` (after `child_GC_CisAjt` intermediate) | `tmp_V*tmp_sgn` | gated write | `CisAjtCkuAlv_element_TraceProbe` |
| `GC_CisAisCisAis_element` (GC, diag) | `j-1` | `tmp_V*tmp_sgn` | gated write | `GC_CisAisCisAis_element_TraceProbe` |
| `GC_CisAisCjtAku_element` (GC) | bare-bit `tmp_off` | `tmp_V*tmp_sgn` | gated write | `GC_CisAisCjtAku_element_TraceProbe` |
| `GC_CisAjtCkuAku_element` (GC) | bare-bit `tmp_off` | `tmp_V*tmp_sgn` | gated write | `GC_CisAjtCkuAku_element_TraceProbe` |
| `GC_CisAjtCkuAlv_element` (GC) | bare-bit `tmp_off` | `tmp_V*tmp_sgn` | gated write | `GC_CisAjtCkuAlv_element_TraceProbe` |

All Hubbard/GC helper write-set (vi): `general_int_GetInfo` -> `X->Large.is*/
A_spin/B_spin/isA_spin/isB_spin/tmp_V/isite*`; `child_CisAjt` -> `GetOffComp`
reads `list_2_*`; `child_GC_CisAjt`/`child_CisAis`/`SgnBit` pure.

**SpinGC-half two-body.** Uses `Rearray_Interactions(...,2)` then, only if the
reordered pair is onsite in both factors (`org_isite1==org_isite2 &&
org_isite3==org_isite4`, expec_cisajscktaltdc.c:1995), computes
`isA_up=Tpow[isite2-1]`, `isB_up=Tpow[isite4-1]` and dispatches by
`(sigma1==sigma2?, sigma3==sigma4?)`. `tmp_V` from Rearray is folded into the
amplitude.

| (i) function @ location (:2002/:2008/:2014/:2020) | (ii) kprime | (iii) amp | (v) adapter |
|---|---|---|---|
| `GC_CisAisCisAis_spin_element` (diag) | `j-1` | `tmp_V*tmp_sgn` | `GC_CisAisCisAis_spin_element_TraceProbe` |
| `GC_CisAisCitAiu_spin_element` | bare-bit `tmp_off` | `tmp_V*tmp_sgn` | `GC_CisAisCitAiu_spin_element_TraceProbe` |
| `GC_CisAitCiuAiu_spin_element` | bare-bit `tmp_off` | `tmp_V*tmp_sgn` | `GC_CisAitCiuAiu_spin_element_TraceProbe` |
| `GC_CisAitCiuAiv_spin_element` | bare-bit `tmp_off` (after intermediate) | `tmp_V*tmp_sgn` | `GC_CisAitCiuAiv_spin_element_TraceProbe` |

Note: the SpinGC-half element funcs carry an `H_CORR` branch (reads
`conj(tmp_v0)`); extraction is M_CORR only so it is dead. The write-set (vi) for
this path is `X->Large.{i_max,irght,ilft,ihfbit,mode}` (driver-set) only — the
element funcs read `isA_up/isB_up` as explicit args, NOT from `X->Large`, and
`general_int_GetInfo` is NOT called here.

**Spin-half CANONICAL two-body (`expec_cisajscktalt_SpinHalf`,
expec_cisajscktaltdc.c:979-1101 — previously unmapped, enumerated here).** After
`Rearray_Interactions(...,2)`, the LOCAL branch
(`org_isite1<=Nsite && org_isite3<=Nsite`, :1064) has THREE reachable families:

| (i) branch @ location | (ii) kprime | (iii) amplitude | (v) adapter |
|---|---|---|---|
| density-density diagonal (`s1==s2 && s3==s4`, :1071) `CisAisCisAis_spin_element` | `j-1` | `tmp_V*tmp_sgn` | `CisAisCisAis_spin_element_TraceProbe` |
| same-index reduction (`o1==o3 && s1==s4 && s3==s2`, :1073-1079) inline `child_Spin_CisAis` | `j-1` | `tmp_V * occ` | driver reuses `child_Spin_CisAis_TraceProbe`, folds `tmp_V` |
| exchange (`s1==s4 && s2==s3`, :1081-1088) `child_exchange_spin_element` | canonical `tmp_off-1` | **`tmp_sgn` ONLY — NO `tmp_V`** | `child_exchange_spin_element_TraceProbe` |

**Amplitude anomaly (recorded per brief).** The exchange branch at
expec_cisajscktaltdc.c:1085-1087 computes `dmv = vec[j]*tmp_sgn;
dam_pr += conj(vec[tmp_off])*dmv` — it does **not** multiply by `tmp_V`, unlike
the diagonal/same-index branches. This is faithfully preserved: the exchange
adapter reports the bare sign and `TraceMapExtractTwoBody`'s Spin exchange branch
does not reintroduce `tmp_V`. (In practice this Sz-conserving spin-flip pair has
`tmp_V=+1` for the ordering that reaches this branch, but the code path omits it
regardless, and we match Mode 1 exactly.)

**Rearray-nonzero semantics (verified per brief).** When
`Rearray_Interactions` returns non-zero (an irregular pair that is neither
`i1==i2 & i3==i4` nor `i1==i4 & i3==i2`), the Mode-1 SpinHalf/SpinGCHalf paths
write a literal `0.0` correlation row and `continue`
(expec_cisajscktaltdc.c:1007-1012 SpinHalf, :1952-1957 SpinGCHalf). This is the
NORMAL "unconventional pair -> zero" behavior, **not** an error path (the
function still returns 0). The extraction driver mirrors it by returning the
`TraceMap.n == 0` sentinel; Task 4's kernel writes the 0.0 row for it.

### 2c.3 Capability conclusions (feeds Task 5)

- **HubbardGC / SpinGC-half (one-body + two-body):** all reachable families have
  clean `*_map` cores and adapters; locally unit-tested (GC bare-bit basis) in
  `test/unit/expec_trace_map_check.c`. Candidates for TRUE after Task 5 golden.
- **Hubbard / Spin-half canonical (one-body + two-body):** all reachable
  families have clean cores and adapters, but the canonical basis
  (`GetOffComp`/`list_1`) is not locally unit-testable; verification is deferred
  to the Task 3/4 clavius early checkpoints and the Task 5 canonical golden.
- **Exclusions:** general-spin (`iFlgGeneralSpin==1`) is out of scope — the
  driver returns `-1` for it. Kondo/tJ diagonal one-body num operators and the
  Kondo localized-site zero rows are handled as empty maps; their capability
  flip still awaits Task 5 coverage. Spinless/Kondo remain FALSE per the plan.
