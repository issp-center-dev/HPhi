# Manual-appendix benchmark scripts (kugui)

Reproduction scripts for the measurements in the manual appendix
"Parallel full diagonalization" (`doc/{ja,en}/source/technical/
parallel_fulldiag_{ja,en}.rst`) and the corresponding records in
`../elpa_gpu_check.md` (2026-07-19/20 entries). All PBS scripts target
the ISSP supercomputer kugui (PBS Pro; `module load intel intel-mpi`,
MKL 2022.x, plus `cuda/12.4` for GPU builds/runs) and assume:

- HPhi source tree at `~/HPhi-elpa`, CPU build in `build_kugui`,
  GPU build in `build_kugui_gpu`
- ELPA 2025.06.001 source tarball at `~/src/elpa-2025.06.001.tar.gz`
- CPU ELPA installed at `~/opt/elpa-2025.06-cpu`
  (see `../elpa_gpu_check.md` 2026-07-18 for its configure line;
  `--disable-avx512 --disable-avx512-kernels` is required on EPYC)
- results are written under `~/bench_appendix/`

| Script | What it measures | Appendix artifact |
|---|---|---|
| `kugui_bench_cpu.pbs` | SpinGC chain L=8..14 (N=2^L): Solver 0/1/3 diagonalization (16 cores per config) + ExpecMode 0/1/2 observables (Solver 3, 16 procs, aggregate GF) | CPU columns of both tables / both figures |
| `kugui_build_gpu.pbs` | Builds ELPA with CUDA sm_80 (`--enable-nvidia-gpu-kernels --enable-nvidia-sm80-gpu`; note the required `-lstdc++` in LIBS) and HPhi against it | prerequisite for all GPU runs |
| `kugui_bench_gpu.pbs` | Same sweep, ELPA GPU x1 (np=1) and x2 (np=2, one node) | GPU x1 / x2 columns |
| `kugui_bench_gpu4.pbs` | Same sweep, GPU x4 = 2 nodes x 2 A100 (np=4, `I_MPI_ADJUST_GATHERV=3`) | GPU x4 column |
| `kugui_gpu_probe.pbs` | Max-N bracketing on 4 GPUs: N=32,768 and N=48,620 | "Maximum feasible size" (GPU) |
| `kugui_l10_gpu.pbs` | L=10 Hubbard (N=63,504) on 4 GPUs — fails in ELPA device allocation (protocol item 18 closure) | "Maximum feasible size" (GPU boundary) |
| `kugui_l10_cpu.pbs` | Same N=63,504 on 4 CPU nodes x 32 ranks | "Maximum feasible size" (CPU) |
| `bench_times.csv` | Collected timers (`L,tag,LapackDiag,CalcPhys` from each run's `output/CalcTimer.dat`) | input to the figures |
| `make_bench_figs.py` | Renders `doc/figs/fulldiag_{solver,expecmode}_bench.png` and prints the RST table rows: `python3 make_bench_figs.py bench_times.csv doc/figs` | both figures / both tables |

Timing convention: diagonalization = `LapackDiag` section and
observables = `CalcPhys` section of `CalcTimer.dat`; each run directory
keeps `run.log`, `time.log` (`/usr/bin/time -v`), and `output/`.
