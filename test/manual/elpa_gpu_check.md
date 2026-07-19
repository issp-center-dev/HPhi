# ELPA GPU manual verification protocol (no GPU CI)

Run before releases that touch the FullDiag/ELPA path.
Target: clavius (single node multi-GPU) and one multi-node GPU system.

## Build
    cmake -DUSE_ELPA=ON -DELPA_ROOT=<prefix> ..   # ELPA >= 2023.11.001 (CUDA build)
    make HPhi elpa_eigen_check

## Checklist
1. `ctest -R elpa_eigen_check` ... expect PASS (CPU path sanity)
2. `ctest -R fulldiag_elpa_hubbard_chain` ... expect PASS
3. GPU smoke (1 rank / 1 GPU):
   `Solver 3`, `NGPU 1`, 8-site Hubbard chain FullDiag (N=4900);
   compare zvo_phys energies against a `Solver 0` run (tol 1e-8);
   confirm "Using ELPA (GPU)" in stdout and GPU utilization in nvidia-smi.
4. GPU multi-rank (ranks = GPUs per node): same comparison, same 8-site case
   (N=4900; larger sizes are impractical for dense FullDiag).
5. Failure-path check: run with `NGPU 1` against a CPU-only ELPA build;
   expect a clear error mentioning NGPU 0 fallback instruction, no silent
   CPU execution.
6. `ctest -R elpa_redist_check` ... expect PASS (phase 2: panel
   redistribution correctness, any np >= 2).
7. `ctest -R fulldiag_spingc_gamma` ... expect PASS (phase 2: transverse-field
   regression, serial).
8. Startup rejection: `Solver 3` + `OutputHam 1` with nproc > 1 must abort
   at startup with the OutputHam/InputHam error (no run).
9. Memory spot check (phase 2): for a fixed N (e.g. 8-site Hubbard, N=4900),
   compare MaxRSS (`/usr/bin/time -v`) between np=1 and np=P; per-rank memory
   must drop roughly as 1/P.
10. `ctest -R fulldiag_expecmode_equiv` (np=2, np=3) ... expect PASS (phase 3a:
    `ExpecMode` 0 vs 1 vs 2 equivalence -- every `zvo_phys*`/Green aggregate
    output file must match within tolerance 1e-8). As of phase 3b (item 15
    below), the same script's case 1 (Hubbard) exercises the real
    `ExpecMode 2` trace kernel rather than a downgrade to Mode 1 -- see
    item 15 for the current, per-case kernel/fallback expectations.
11. `ctest -R green_partial_merge_check` ... expect PASS (phase 3a: partial
    (`.part<rank>`) file merge success/failure semantics, any np >= 2).
12. Mode 1 benchmark: same fixed-N FullDiag case as item 9, `Solver 3`,
    np >= 4; compare the observable-evaluation wall time (`CalcTimer.dat`)
    between `ExpecMode 0` and `ExpecMode 1`; expect `ExpecMode 1` faster
    (redundant per-rank re-evaluation of every eigenstate eliminated) with
    both giving matching physics (per item 10).
13. S2/Sz distributed-`ExpecMode 0` fix check: `Solver 1` or `Solver 3`,
    nproc > 1, `ExpecMode 0` (default); confirm stdout progress lines show
    a non-zero `S2=` column matching the `Solver 0` / serial reference
    (not the pre-phase-3a zero-filled `S2=0.000000`), and that the line
    format matches the single-process (serial) layout exactly (same
    fields, same order, `S2` column present).
14. `ExpecMode` eligibility / demotion checks: `ExpecMode 1` with `Solver 0`
    (or with a non-FullDiag `CalcType`) aborts at startup with the
    ExpecMode eligibility error (`cErrExpecMode`); `ExpecMode 1` at nproc=1
    prints the "reverts to 0" INFO and matches the `ExpecMode 0` reference.
    As of phase 3b, `ExpecMode 2` no longer unconditionally runs as
    `ExpecMode 1`; see item 15 for its current per-quantity kernel/fallback
    INFO lines.
15. `ctest -R expec_trace_map_check` ... expect PASS (phase 3b: ExpecMode-2
    trace-kernel unit tests -- mapping-probe correctness/purity, one-body
    and two-body streaming, and the `HPHI_TRACE_BUF_MAX_MB` memory-gate
    boundary, for the GC models only (HubbardGC/SpinGC-half); canonical
    Hubbard/Spin coverage comes from item 16's equivalence cases 3/5).
16. `ctest -R fulldiag_expecmode_equiv` (np=2, np=3), production path ...
    expect PASS (phase 3b: cases 1/2/3/5 -- Hubbard chain, SpinGC Gamma
    chain, canonical Spin chain, and the dedicated Hubbard one-body/
    two-body golden case, all of them supported-model rows -- must have
    `ExpecMode 2` print "ExpecMode 2: one-body Green functions use the
    trace kernel." and "ExpecMode 2: two-body Green functions use the
    trace kernel." (not a downgrade), with `zvo_phys*`/Green aggregate
    output, including the `var` column, matching `ExpecMode 0`/`1`
    within 1e-8. Case 4 (SpinGC honeycomb with ThreeBodyG/FourBodyG/
    SixBodyG defined) must instead print the two-body shared-evaluator
    fallback line ("... they share their evaluator with
    three-/four-/six-body Green functions.") while still selecting the
    one-body trace kernel. (HubbardGC and the unsupported-model fallback
    are not exercised by this script -- their coverage is
    `expec_trace_map_check`'s unit tests plus the per-row evidence
    comment above `kTraceCap` in `src/expec_trace.c`; item 15 above
    already covers HubbardGC.)
17. **Benchmark gate** (phase 3b, spec §6): fixed-N FullDiag case as item 9
    (8-site Hubbard chain, N=4900), `Solver 3`, np=4, one-body + two-body
    GF defined for all states; compare observable-evaluation wall time
    (`CalcTimer.dat`) across `ExpecMode 0`/`1`/`2`. Target: `ExpecMode 2`
    >= `ExpecMode 1`. Record the measured times **and** a breakdown of
    `ExpecMode 2`'s time into mapping-extraction vs streaming vs output
    phases (completion requires the breakdown, not just the totals); if
    the target is missed, record a break-even analysis instead of
    dropping the item, since the usage guidance in the `ExpecMode`
    documentation and the phase 3b migration note may need to be revised
    to match the measured result.
    Measurement method (final whole-branch review fix): the breakdown
    comes straight from the rank-0 stdout lines the `ExpecMode 2`
    orchestrator (`src/phys_distributed.c`, after `ExpecLocalLeave()`)
    prints for every quantity that ran as the trace kernel:
    `  ExpecMode 2 timing (rank 0): one-body map=%.3fs stream=%.3fs output=%.3fs`
    `  ExpecMode 2 timing (rank 0): two-body map=%.3fs stream=%.3fs output=%.3fs`
    (`map`/`stream`/`output` are exactly the mapping-extraction/streaming/
    output phases above); totals are wall-clock plus the `CalcTimer.dat`
    expec section, as in phase 3a.
18. Optional GPU data point (phase 3b, if clavius GPUs are free -- check
    `nvidia-smi` first and do not contend with other users' jobs): L=10
    Hubbard (N ~= 63504) FullDiag with `Solver 3`, `NGPU` > 0; record one
    wall-time comparison across `ExpecMode 0`/`1`/`2` for the GPU
    diagonalization + observable-evaluation path.
Record results (date, host, ELPA version, commit) at the bottom of this file.

---

## Results

### 2026-07-10 — clavius (s76), CPU phase — PASSED
- Commit: 8732e240 / ELPA 2025.06.001 (conda-forge, mpi_openmpi, non-threaded libelpa) / OpenMPI 5.0.8 / gcc 13.3 / conda env `hphi_elpa`
- Build: `cmake -DUSE_ELPA=ON -DELPA_INCLUDE_DIR=$CONDA_PREFIX/include/elpa -DELPA_LIBRARY=$CONDA_PREFIX/lib/libelpa.so -DSCALAPACK_LIBRARIES="-L$CONDA_PREFIX/lib -lscalapack"` — manual-override detection path exercised; `ELPA_HAVE_SETUP_GPU=1` detected (2025.06 API).
- Checklist 1 `elpa_eigen_check`: PASS at np=1,2,3,4,6,8,12,16 (residual, orthogonality, LAPACK agreement). ctest (exact:3) PASS.
- Checklist 2 `fulldiag_elpa_hubbard_chain` (np=2, multi-rank HPhi end-to-end incl. startup gates): PASS. `fulldiag_solver_keyword` on ELPA build (capability branches): PASS.
- Checklist 5 failure path: `Solver 3` + `NGPU 1` against CPU-only ELPA → clean abort (exit 14) with "the linked ELPA has no NVIDIA GPU support / Set NGPU 0" — no silent CPU fallback.
- Defects found on hardware and fixed in 8732e240: (a) nblk must be capped so every process row/col owns a block (ELPA_ERROR_SETUP otherwise); (b) ELPA2 2stage gives inaccurate eigenvectors with capped nblk (e.g. 24 on 4x2) → 1stage fallback when capped.
- Additional (Codex-recommended) checks, all PASS: N(=4) < process grid (5x5, np=25) aborts with the clear "reduce ranks" error on all ranks (no hang); uncapped nblk=64 multi-rank 2stage (L=6 chain, N=400, np=2) matches Solver 0 energies exactly (maxdiff 0.0 over 400 eigenvalues).
- Checklist 3/4: see GPU phase below.

### 2026-07-10 — clavius (s76), GPU phase — PASSED
- ELPA 2025.06.001 built from source with CUDA (`--enable-nvidia-gpu-kernels --with-NVIDIA-GPU-compute-capability=sm_89`, CUDA 12.9, prefix `~/opt/elpa-2025.06-cuda`, versioned include dir — exercised FindELPA's ELPA_ROOT glob path; `ELPA_HAVE_SETUP_GPU=1`).
- Checklist 3 (GPU smoke, 1 rank / 1 GPU): 8-site Hubbard chain FullDiag (N=4900), `Solver 3` + `NGPU 1`, `CUDA_VISIBLE_DEVICES=1`. "Using ELPA (GPU)" printed; HPhi observed as GPU compute app in nvidia-smi during the run; exit 0; all 4900 eigenvalues match the `Solver 0` LAPACK reference exactly (printed precision).
- Checklist 4 (multi-rank GPU): np=2, `NGPU 2`, both RTX 6000 Ada visible. Two distinct GPU UUIDs active during the run (ELPA round-robin = 1 rank per GPU as designed); exit 0; eigenvalues again match exactly.
- Checklist 5 (second variant): with the CPU conda libelpa shadowing the CUDA one via RPATH, the run aborted cleanly with "the linked ELPA has no NVIDIA GPU support / Set NGPU 0" — no silent CPU execution. **Troubleshooting note:** if a CPU-only libelpa with the same soname is on the loader path (e.g. conda env), it can shadow the CUDA build; ensure the CUDA ELPA lib dir wins (LD_LIBRARY_PATH/rpath) or remove the CPU copy.
- Protocol note: the checklist originally prescribed 12/14-site chains; those Hilbert dimensions are impractical for dense FullDiag, so the checklist now prescribes the validated 8-site case (N=4900).

### Benchmark gate (phase 1, N=4900, np=2, LapackDiag step from CalcTimer.dat)
| Backend | Diag time |
|---|---|
| Solver 1 ScaLAPACK pzheev (2 CPU ranks) | 149.8 s |
| Solver 3 ELPA CPU 2stage (2 CPU ranks) | 25.8 s (5.8x vs pzheev) |
| Solver 3 ELPA GPU 1stage (2 ranks x 2 RTX 6000 Ada) | 4.2 s (35x vs pzheev) |

### 2026-07-11 — clavius (s76), Phase 2 (distributed Hamiltonian generation) — PASSED
- Commit under test: a87c3e2f (phase-2 complete). Build: fresh `-DUSE_ELPA=ON -DELPA_ROOT=~/opt/elpa-2025.06-cuda` (CUDA ELPA 2025.06.001).
- `elpa_redist_check` (panel vs replicated fill, element-wise): PASS at np=1,2,3,4,6,8.
- `elpa_eigen_check`: PASS at np=1,2,3,4,8. ctest: `fulldiag_elpa_hubbard_chain` (incl. new OutputHam-rejection case), `elpa_redist_check`, `fulldiag_solver_keyword`, `fulldiag_spingc_gamma` all PASS.
- Distributed-generation physics (Solver 3 CPU, panel path) vs Solver 0, all eigenvalues at printed precision: Hubbard L=6 (N=400) np=4 maxdiff 0.0; Hubbard L=8 (N=4900) np=2 and np=4 maxdiff 0.0; SpinGC L=8 Gamma=0.5 (N=256, covers the reworked transverse-field hunk) np=2 maxdiff 0.0.
- Memory scaling (L=8, N=4900, MaxRSS of largest rank via /usr/bin/time -v): replicated np=1 = 1.93 GB -> distributed np=4 = 0.63 GB per rank (O(N^2/P) confirmed).
- GPU distributed generation: np=2 x 2 GPUs (`NGPU 2`, both RTX 6000 Ada free and used — 2 distinct GPU UUIDs observed), exit 0, all 4900 eigenvalues match Solver 0 exactly.

### 2026-07-12 — clavius (s76), Phase 3a (ExpecMode / state-parallel observables) — PASSED
- Commit under test: dd5abc84 + 920287a2 (test-def fix, see below). Build: `-DUSE_ELPA=ON -DELPA_ROOT=~/opt/elpa-2025.06-cuda` (Release) + a Debug (assert-enabled) build with explicit `-DSCALAPACK_LIBRARIES` for spot checks.
- **Environment workaround (record for future runs):** another user's `nvidia-cuda-mps-server` made `cuInit()` hang inside UCX's CUDA module constructor at process start (recvmsg on the MPS socket), so *any* HPhi run — CPU-only included — froze before `main`. Bypass: `export CUDA_MPS_PIPE_DIRECTORY=/tmp/nonexistent-mps-hphi` (cuInit then fails fast and UCX proceeds CPU-only). This also explains the initial `fulldiag_solver_keyword` 30-min ctest timeout, which passed (109 s) after the bypass.
- ctest (serial, ELPA build): `fulldiag_solver_keyword` PASS, `check_expec_local_calls` PASS.
- `elpa_statepanel_check` direct sweep: PASS at np=1,2,3,4,8. `elpa_statepanel_zero_owner` (exact:3, N=2): PASS. Other ELPA tests (`elpa_eigen_check`, `elpa_redist_check`, `fulldiag_elpa_hubbard_chain` incl. the strengthened 5-column reference check, `green_partial_merge_check`) PASS in the MPIRUN=np2/np3 ctest rounds.
- **Defect found by the equivalence test and fixed (920287a2):** case-4 `green6.def` referenced site 11 in the 8-site honeycomb cell. HPhi does not validate ThreeBody/FourBody/SixBodyG site indices against Nsite; an out-of-range site takes the inter-PE branch and crashes with an integer divide-by-zero (`Tpow[11]==0`) in `child_GC_CisAitCiuAiv_spin_MPIsingle` — in **both** ExpecMode 0 and 1 (pre-existing robustness gap, not a phase-3a regression; upstream validation fix filed as follow-up).
- `fulldiag_expecmode_equiv_np2` / `_np3`: PASS after the def fix (all 4 cases: Hubbard+NBodyG, SpinGC Gamma, canonical Spin, SpinGC honeycomb with 3/4/6-body GF; energies + all aggregate Green files equal within 1e-8, ExpecMode 2 downgrade INFO verified).
- **Benchmark gate** (L=8 Hubbard chain, N=4900, one-body+two-body GF for all states, Solver 3 CPU, np=4): `expec_energy_flct` timer section 3.96 s -> 0.98 s = **4.04x (~P=4)**; total wall 48.7 s -> 22.4 s. `zvo_phys` Mode 0 vs Mode 1: maxdiff **0.0** over 4900 states; per-eigen Green file sets identical (9809 files each).
- **Mode 0 S²/Sz unification** (behavior fix): distributed Mode 0 (np=4) vs serial Solver 0 — `<H>/<N>/<Sz>` maxdiff 0.0; `<S2>` maxdiff 0.0 over all 246 energy-isolated states (differences appear only inside degenerate subspaces, where the eigenbasis is not unique — expected). Header/stdout now identical to the serial format.
- **Solver 1 (ScaLAPACK) + ExpecMode 1** (L=6 chain, N=400, np=2): Mode 0 vs 1 maxdiff 0.0 (exercises the non-ELPA leg of the eligibility matrix).
- **Zero-owner ranks through the full driver** (SpinGC L=2, N=4, np=6 -> ranks 4-5 own zero states): Mode 0 vs 1 maxdiff 0.0, clean exit.
- **Assert-enabled Debug build** (ExpecLocal nesting/exit asserts, redistribution guards live): zero-owner Mode 1 run EXIT=0, maxdiff 0.0 vs Mode 0.

### 2026-07-12 — clavius (s76), Phase 3b (ExpecMode 2 trace kernels) — PASSED
- Commit under test: 93083ba5 (+ cae2b75f instrumentation). Build: `-DUSE_ELPA=ON -DELPA_ROOT=~/opt/elpa-2025.06-cuda` (Release). `CUDA_MPS_PIPE_DIRECTORY` bypass in effect (see Phase 3a note).
- Item 15 `expec_trace_map_check` (np=1): ALL PASS (mapping validity/purity, streaming, gate boundaries, no-operators, shared-evaluator plan cases).
- Item 16 production-path equivalence: `fulldiag_expecmode_equiv_np2`/`_np3` PASS (serial 2/2, np=2 suite 8/8, np=3 suite 8/8 incl. statepanel/zero-owner/merge/hubbard-chain). Cases 1/2/3/5 select both trace kernels; case 4 correctly demotes two-body via the shared-evaluator rule while keeping the one-body kernel. All output trees (incl. `var`) match within 1e-8.
  - Development-stage checkpoints (recorded for provenance): forced-kernel equiv np=2/3 after Task 4; the strengthened production run after Task 5 caught a REAL bug (mode 2 dropped `zvo_ThreeBody/FourBody/SixBody_eigen.dat` because `expec_cisajscktaltdc` shares its evaluator with the two-body GF) — fixed by the plan-level shared-evaluator demotion (dc6791ee).
- **Item 17 benchmark gate** (L=8 Hubbard chain, N=4900, one-body 16 + two-body 48 pairs for all states, Solver 3 CPU, np=4):
  | Mode | wall | expec_energy_flct (CalcTimer) |
  |---|---|---|
  | 0 | 64.2 s | 0.75 + 5.13 s |
  | 1 | 32.3 s | 0.18 + 1.27 s |
  | **2** | **26.0 s** | 0.12 + 0.86 s |
  Target `ExpecMode 2 >= ExpecMode 1` met (26.0 s vs 32.3 s; 2.5x vs Mode 0). Breakdown (rank-0 lines, per item 17's method): one-body map=0.000s stream=0.233s output=0.036s; two-body map=0.001s stream=0.781s output=0.052s — the kernels spend ~1.1 s where the Mode-1 evaluators spent ~7.4 s of the wall delta; the remaining Mode-2 time is the (deliberately un-kernelized) energy-family fallback (`mltply` per state). `zvo_phys` maxdiff 0.0 vs Mode 0 for both modes; per-eigen Green file sets identical (9809 files each).
- Fallback-reason INFO lines exercised END-TO-END on the production path: (a) shared-evaluator (case 4); (b) no-operators (TwoBodyG line removed from namelist → "no operators of this kind are defined", run exits 0); (c) memory gate (96 duplicated two-body pairs at `HPHI_TRACE_BUF_MAX_MB=1` → "result buffer would exceed HPHI_TRACE_BUF_MAX_MB", physics maxdiff 0.0 and the two-body per-eigen file BYTE-IDENTICAL to Mode 0 via the fallback). Note: the gate arithmetic means realistic small correlation defs (48 pairs × NC=1225 ≈ 0.9 MiB) fit under even the 1 MiB floor — demotion requires genuinely large nops×NC.
- Item 18 (optional L=10 GPU point): SKIPPED this round — N≈63504 dense diagonalization extrapolates to ~2.5 h occupying both GPUs (N³ scaling from 4.2 s at N=4900), and GPU 0 was running another user's MACE workload; per the GPU-sharing rules the run was not attempted. The item remains open for a quiet window.

### 2026-07-17 — clavius (s76), PR #276 blocker fix (invalid-row AddHamElem) — VERIFIED
- Maintainer report (tmisawa): ASan heap-buffer-overflow at makeHam.c:278 with `Solver 3` / nproc>1 — `AddHamElem(tmp_off, ...)` called with `tmp_off=0` on invalid/no-op transitions; the distributed panel indexes `(irow-1)` → underrun (the legacy replicated matrix silently absorbed row 0).
- Fix (commit 8729c326): validity guards at the 9 unsafe makeHam.c call sites (4 × `tmp_off > 0`, 5 × `dmv != 0.0` where the helper's out-param can be stale; skipping a `+= 0` is numerically exact) + `assert(hs_i_ >= 1)` (both branches) and `assert(hs_i_ <= HamPanelLd)` (panel branch) inside AddHamElem. 24 remaining sites audited as safe (diagonal, GC bare-bit-bounded, or already gated — 33 sites total: 28 in makeHam.c, 4 in nbody_interall.c, 1 in anomalous_pair.c; count corrected per the maintainer's recount, the original record said 34/25); inventory in the fix report.
- Follow-up hardening (same review, commit after 64d593ad): the diag failure path in `lapack_diag_elpa()` now frees `Z_vec` and NULLs it; the two remaining PR-introduced `descinit_` sites (the ELPA `descA`/`descZ_vec` pair in lapack_diag.c and the generation-panel descriptor in `RedistPanelToBlockCyclic`) got the same synchronized-verdict check as `RedistBlockCyclicToStatePanel`. (`diag_scalapack_cmp`'s descriptors predate this PR — legacy Solver-1 path, left as an upstream follow-up.)
- RED→GREEN on hardware (valgrind 3.27, conda env; ASan unusable under mpiexec here — shadow-range conflict, google/sanitizers#856): pre-fix binary → `Invalid read of size 16 ... 16 bytes before a block of size 10,368 alloc'd (setmem_large)` in makeHam on BOTH ranks; post-fix binary → **zero Invalid accesses** in the same window. (Both runs later stop inside ELPA's AVX-512 kernels with SIGILL — a valgrind instruction-decoder limitation, after generation completes; the pre/post comparison window is identical.)
- Debug (assert-enabled) build: full `fulldiag_elpa_hubbard_chain` run to completion, no assert trips. Release physics: chain + equiv np=2/3 + statepanel + merge 5/5+5/5 PASS.
- Remaining valgrind contexts are environment noise (OpenMPI/PMIx/libnvidia-ml) plus a PRE-EXISTING `GetFileName`/`ReadDefFileNInt` uninitialised-strlen in legacy def parsing (untouched by this branch; upstream follow-up candidate).

### 2026-07-18 — kugui (ISSP HPE Cray, AMD EPYC 7763 + Mellanox IB), Intel toolchain + MULTI-NODE — PASSED
- Commit under test: af858a0c. Toolchain: Intel oneAPI 2022.2.1 (classic icc/ifort) + Intel MPI 2021.7.1 + MKL 2022 ScaLAPACK; RHEL 8 / glibc 2.28; cmake 3.20.2; PBS. ELPA 2025.06.001 built from source (mpiicc/mpiifort, MKL, `--disable-openmp`).
- This closes the protocol's "one multi-node system" target (CPU): first true multi-node validation.
- **Two platform pitfalls found, diagnosed (gdb stacks + A/B tests), and worked around — now documented in the installation manual (ja/en)**:
  1. **ELPA AVX-512 kernel on a non-AVX-512 CPU**: ELPA compiles AVX-512 kernels regardless of build host and its runtime selection picked `single_hh_trafo_complex_AVX512_1hv_double` on EPYC → SIGILL in every uncapped-nblk 2-stage path (eigen np=1 N=97, equiv case4 N=256), while capped 1-stage paths worked — deceptive partial-pass pattern. Not a compiler-flag issue (persisted with plain `-O2`). Fix: configure ELPA with `--disable-avx512 --disable-avx512-kernels`.
  2. **Intel MPI gatherv deadlock (mlx/UCX provider)**: `ExpecMode 1/2`'s all_* `MPI_Gatherv` deadlocked (root in Waitall, sender in Ssend — gdb-confirmed) under the default tuned `linear_ssend` algorithm and `I_MPI_ADJUST_GATHERV=1/2`; completes with `I_MPI_ADJUST_GATHERV=3` or `I_MPI_FABRICS=shm` (same binary/call — judged an Intel MPI 2021.7 platform defect, not an HPhi protocol issue; OpenMPI on the workstation only ever exercised the shm path).
- Single node (with both workarounds): map_check, eigen/redist/statepanel np=1,2,4,8 (incl. all failure-injection phases), merge np=2, serial ctest, hubbard_chain, equiv np=2/3 (all 5 cases × 3 modes) — ALL PASS.
- **Multi-node (2 nodes × 4 ranks = 8, `I_MPI_ADJUST_GATHERV=3`)**: eigen/redist/statepanel/merge at 8 ranks over 2 nodes PASS (note: each unit run took ~15 min under the mlx provider — budget walltime accordingly). End-to-end L=8 Hubbard N=4900, Solver 3, all states one-+two-body GF:
  | Mode | wall (2 nodes / 8 ranks) |
  |---|---|
  | 0 | 11 m 52 s |
  | 1 | 52 s |
  | **2** | **49 s** |
  `zvo_phys` maxdiff 0.0 for Modes 1/2 vs Mode 0. The Mode-0 wall is dominated by per-state inter-node collectives (eigenvector gather + reductions × 4900 states) — exactly the cost the phase-3 state-task parallelism removes (**~14x faster**); Mode 2 ≥ Mode 1 holds here too (kernel breakdown: one-body map 0.000/stream 0.069/output 0.768 s; two-body 0.001/0.241/0.625 s).

### 2026-07-19 — kugui, manual-appendix benchmark sweep (CPU + A100 GPU) — COMPLETED
- Purpose: size-vs-time data for the new manual appendix (`doc/{ja,en}/source/technical/parallel_fulldiag_*.rst`). Binaries under test correspond EXACTLY to source commit dca3e750 (the working tree additionally carried documentation-only changes; `src/` and `test/` were byte-identical to dca3e750), rebuilt on kugui in `~/HPhi-elpa/build_kugui` (CPU) and `~/HPhi-elpa/build_kugui_gpu` (GPU).
- Reproduction: PBS scripts `~/kugui_bench_cpu.pbs` (sweep incl. run matrix and stan.in template), `~/kugui_build_gpu.pbs` (ELPA CUDA + HPhi build), `~/kugui_bench_gpu.pbs` (GPU sweep) on kugui. Env: `module load intel intel-mpi` (+ `cuda/12.4` for GPU), MKL 2022.2.1, `I_MPI_ADJUST_GATHERV=3`, CPU ELPA `~/opt/elpa-2025.06-cpu`, GPU ELPA `~/opt/elpa-2025.06-a100`. Timings = `LapackDiag` / `CalcPhys` sections of each run's `output/CalcTimer.dat`.
- System: SpinGC chain, `J=1`, `Gamma=0.5`, `L=8,10,12,14` (`N=2^L=256..16384`), standard mode `-sdry` + `calcmod.def` edits, `outputmode="none"` for the solver sweep / `"correlation"` + `OutputGreenFormat 1` for the ExpecMode sweep.
- CPU (PBS 858288, F1cpu, EPYC 7763, 16 cores per config — Solver 0: 1p x 16t; Solver 1/3: 16p x 1t; `I_MPI_ADJUST_GATHERV=3`): all 24 runs rc=0. `LapackDiag` at N=16384: S0 1929.6 s / S1 536.9 s / S3 172.6 s. `CalcPhys` (Solver 3, 16p) at N=16384: Mode 0 424.5 s / Mode 1 26.4 s / Mode 2 20.8 s (trace kernels confirmed active via INFO lines).
- GPU: ELPA 2025.06.001 rebuilt with CUDA 12.4 (`--enable-nvidia-gpu-kernels --enable-nvidia-sm80-gpu --with-NVIDIA-GPU-compute-capability=sm_80 --disable-avx512{,-kernels}`; needs `-lstdc++` appended to LIBS — first attempt failed linking `elpa2_print_kernels` with undefined `std::ios_base::Init`) as `~/opt/elpa-2025.06-a100`; HPhi rebuilt against it (`build_kugui_gpu`, PBS 858296). Sweep on F2acc (858297, EPYC 7763 + 2x A100-SXM4-40GB): all 8 runs rc=0, "Using ELPA (GPU)" confirmed. `LapackDiag` at N=16384: 1 GPU 48.2 s / 2 GPUs 35.5 s.
- Full tables/figures: manual appendix (both languages); raw dirs kugui `~/bench_appendix/{cpu,gpu}_L*_*`.
