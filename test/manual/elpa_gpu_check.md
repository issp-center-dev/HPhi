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
   `Solver 3`, `NGPU 1`, 12-site Hubbard chain FullDiag;
   compare zvo_phys energies against a `Solver 0` run (tol 1e-8);
   confirm "Using ELPA (GPU)" in stdout and GPU utilization in nvidia-smi.
4. GPU multi-rank (ranks = GPUs per node): same comparison at 14 sites.
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
- Checklist 3/4 (GPU smoke, CUDA ELPA >= 2023.11.001): PENDING — requires a CUDA build of ELPA (RTX 6000 Ada, CUDA 12.9 available).
