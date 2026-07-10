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
