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
