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
    boundary, for Hubbard/HubbardGC/Spin-half/SpinGC-half).
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
