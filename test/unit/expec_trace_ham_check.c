/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 The University of Tokyo */

/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */

/*
 * Serial (no-MPI) unit test for the phase-3c Task-2 CSR Hamiltonian collector
 * (TraceHamCollect / TraceHamFree / TraceHamGatedBytes, src/expec_trace_ham.c).
 *
 * Part 1 (this task) checks the MATRIX only: rowptr / colidx / val. For each
 * model fixture it drives HPhi's real setup pipeline from a StdFace-generated
 * def set (StdFace_main -> ReadDef -> setmem -> sz -> diagonalcalc), builds the
 * legacy replicated dense Ham via makeHam(), then runs TraceHamCollect() (which
 * re-enumerates the SAME makeHam() in HAM_SINK_TRACE_COLLECT mode), expands the
 * CSR back to dense and asserts they agree to 1e-13 + 1e-13*hmax elementwise.
 *
 * The two energy-family models StdFace cannot emit directly (tJ / tJGC) are
 * produced by generating the corresponding Hubbard / HubbardGC def set and
 * patching CalcModel in calcmod.def -- the def-file structure is identical and
 * sz()/makeHam() build the constrained tJ basis from the same modpara; the
 * dense-vs-CSR identity is model-agnostic (both go through one makeHam()).
 *
 * Two further checks (also matrix-level) are on the Hubbard fixture:
 *   - Gate boundary: want = TraceHamGatedBytes(counting quantities). Collect
 *     SUCCEEDS at cap==want and DEMOTES (returns 0, csr zeroed) at cap==want-1.
 *   - Allocation-failure injection: fail_alloc_at = 0..5 each returns 0 with no
 *     leak and the sink restored, proven by a subsequent normal collect that
 *     succeeds.
 *
 * The coefficient arrays csr->diag[] (n_diag per the frozen model table) are
 * allocated and zeroed here; their filling and verification are Tasks 3-4.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/stat.h>
#ifdef MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#include "struct.h"
#include "global.h"
#include "DefCommon.h"
#include "StdFace_main.h"
#include "readdef.h"
#include "xsetmem.h"
#include "HPhiTrans.h"
#include "check.h"
#include "sz.h"
#include "diagonalcalc.h"
#include "makeHam.h"
#include "hamstore.h"
#include "wrapperMPI.h"
#include "expec_trace.h"
#include "expec_trace_ham.h"
#include "expec_energy_flct.h"
#include "CalcTime.h"

/* Globals provided by the linked src/global.c (Ham, list_*, v0/v1, sink hooks,
   myrank/nproc/stdoutMPI, iHamPanelActive/iHamSinkMode). */

static int g_failures = 0;
static char g_base[4096];

static void expect_true(const char *name, int cond) {
  if (!cond) { fprintf(stderr, "  FAIL %s\n", name); g_failures++; }
  else fprintf(stderr, "  ok   %s\n", name);
}

/* Overwrite calcmod.def's CalcModel line (used for tJ / tJGC). */
static void patch_calcmod(int model) {
  char lines[128][512];
  int nl = 0, i;
  FILE *f = fopen("calcmod.def", "r");
  if (f == NULL) return;
  while (nl < 128 && fgets(lines[nl], sizeof(lines[0]), f)) nl++;
  fclose(f);
  f = fopen("calcmod.def", "w");
  for (i = 0; i < nl; i++) {
    char key[64];
    if (sscanf(lines[i], "%63s", key) == 1 && strcmp(key, "CalcModel") == 0)
      fprintf(f, "CalcModel  %d\n", model);
    else
      fputs(lines[i], f);
  }
  fclose(f);
}

/* Overwrite trans.def with a caller-supplied body (used by part 3 to inject a
   complex hopping and a cancelling ±1e8 duplicate pair). `body` is written
   verbatim after StdFace generates the def set but before ReadDefFileIdxPara
   parses it, so makeHam builds H from exactly these transfer lines. */
static void write_trans(const char *body) {
  FILE *f = fopen("trans.def", "w");
  if (f == NULL) return;
  fputs(body, f);
  fclose(f);
}

/* Drive the full setup for one fixture (cwd already inside its scratch dir).
   If trans_override != NULL it replaces trans.def (complex/cancellation cases).
   Returns 0 on success with X populated + dense Ham built; -1 on setup error. */
static int setup_fixture_trans(struct BindStruct *X, const char *stan,
                               int patch_model, const char *trans_override) {
  FILE *f;
  memset(X, 0, sizeof(*X));
  iHamPanelActive = 0;
  iHamSinkMode = 0;
  hamCollectSink = NULL;

  f = fopen("stan.in", "w");
  if (f == NULL) return -1;
  fputs(stan, f);
  fclose(f);

  StdFace_main("stan.in");
  if (patch_model >= 0) patch_calcmod(patch_model);
  if (trans_override != NULL) write_trans(trans_override);

  mkdir("output", 0777); /* check() writes CHECK_*.dat here */

  setmem_HEAD(X);
  if (ReadDefFileNInt("namelist.def", &X->Def, &X->Boost) != 0) return -1;
  setmem_def(X, &X->Boost);
  if (ReadDefFileIdxPara(&X->Def, &X->Boost) != 0) return -1;
  SetConvergenceFactor(&X->Def);
  if (HPhiTrans(X) != 0) return -1;
  if (check(X) == MPIFALSE) return -1;
  if (setmem_large(X) != 0) return -1;
  if (sz(X, list_1, list_2_1, list_2_2) != 0) return -1;
  diagonalcalc(X);

  /* Legacy replicated dense pass. */
  iHamPanelActive = 0;
  iHamSinkMode = 0;
  if (makeHam(X) != 0) return -1;
  return 0;
}

/* Thin wrapper: the original (no trans override) signature parts 1-2 call. */
static int setup_fixture(struct BindStruct *X, const char *stan, int patch_model) {
  return setup_fixture_trans(X, stan, patch_model, NULL);
}

/* Expand a collected CSR to a fresh dense buffer (0-based, row-major n*n). */
static double complex *csr_to_dense(const TraceHamCsr *csr) {
  long int n = csr->n, i, k;
  double complex *dc = (double complex *)calloc((size_t)n * (size_t)n,
                                                sizeof(double complex));
  for (i = 0; i < n; i++) {
    for (k = csr->rowptr[i]; k < csr->rowptr[i + 1]; k++) {
      dc[i * n + csr->colidx[k]] += csr->val[k];
    }
  }
  return dc;
}

static void run_matrix_fixture(const char *label, const char *subdir,
                               const char *stan, int patch_model) {
  struct BindStruct X;
  TraceHamCsr csr;
  double complex *dc;
  double hmax = 0.0;
  long int n, i, j;
  int ok, mism = 0;

  fprintf(stderr, "[%s]\n", label);
  if (chdir(g_base) != 0) { expect_true(label, 0); return; }
  mkdir(subdir, 0777);
  if (chdir(subdir) != 0) { expect_true(label, 0); return; }

  if (setup_fixture(&X, stan, patch_model) != 0) {
    fprintf(stderr, "  FAIL %s: setup failed\n", label);
    g_failures++;
    return;
  }
  n = (long int)X.Check.idim_max;

  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++) {
      double a = cabs(Ham[i][j]);
      if (a > hmax) hmax = a;
    }

  ok = TraceHamCollect(&X, SIZE_MAX / 2, -1, &csr);
  {
    char nm[256];
    snprintf(nm, sizeof(nm), "%s: TraceHamCollect succeeds", label);
    expect_true(nm, ok == 1);
  }
  if (!ok) return;

  {
    char nm[256];
    snprintf(nm, sizeof(nm), "%s: csr->n == idim_max", label);
    expect_true(nm, csr.n == n);
    snprintf(nm, sizeof(nm), "%s: rowptr[0]==0 && rowptr[n]==nnz", label);
    expect_true(nm, csr.rowptr[0] == 0 && csr.rowptr[n] == csr.nnz);
  }

  dc = csr_to_dense(&csr);
  for (i = 0; i < n && mism == 0; i++) {
    for (j = 0; j < n; j++) {
      double d = cabs(dc[i * n + j] - Ham[i + 1][j + 1]);
      if (d > 1e-13 + 1e-13 * hmax) {
        fprintf(stderr, "  mismatch (%ld,%ld): csr=(% .3e,% .3e) Ham=(% .3e,% .3e) |d|=%.3e\n",
                i, j, creal(dc[i * n + j]), cimag(dc[i * n + j]),
                creal(Ham[i + 1][j + 1]), cimag(Ham[i + 1][j + 1]), d);
        mism = 1;
        break;
      }
    }
  }
  {
    char nm[256];
    snprintf(nm, sizeof(nm), "%s: dense CSR == replicated Ham (n=%ld, nnz=%ld)",
             label, n, csr.nnz);
    expect_true(nm, mism == 0);
  }
  free(dc);
  TraceHamFree(&csr);
}

/* Gate boundary + allocation-failure injection, on the Hubbard fixture. */
static void run_gate_and_injection(void) {
  struct BindStruct X;
  TraceHamCsr csr;
  const char *stan =
    "L = 4\n"
    "model = \"Hubbard\"\n"
    "method = \"FullDiag\"\n"
    "lattice = \"chain\"\n"
    "t = 1.0\n"
    "U = 4.0\n"
    "nelec = 4\n"
    "2Sz = 0\n";
  long int n, i, k_max = 0, nnz_raw;
  int n_diag, ok, fa;
  size_t want;

  fprintf(stderr, "[gate boundary + injection: Hubbard L=4]\n");
  if (chdir(g_base) != 0) { expect_true("gate: chdir base", 0); return; }
  mkdir("gate", 0777);
  if (chdir("gate") != 0) { expect_true("gate: chdir sub", 0); return; }

  if (setup_fixture(&X, stan, -1) != 0) {
    fprintf(stderr, "  FAIL gate: setup failed\n");
    g_failures++;
    return;
  }
  n = (long int)X.Check.idim_max;

  /* Oversized-cap collect to read back the counting quantities. Hubbard has no
     duplicate structural entries, so nnz_raw == csr.nnz and each row's raw
     count == its merged segment length. */
  ok = TraceHamCollect(&X, SIZE_MAX / 2, -1, &csr);
  expect_true("gate: initial oversized collect succeeds", ok == 1);
  if (!ok) return;
  nnz_raw = csr.nnz;
  n_diag = csr.n_diag;
  for (i = 0; i < n; i++) {
    long int len = csr.rowptr[i + 1] - csr.rowptr[i];
    if (len > k_max) k_max = len;
  }
  TraceHamFree(&csr);

  want = TraceHamGatedBytes(n, nnz_raw, k_max, n_diag);
  expect_true("gate: want is a real (non-overflow) size", want != SIZE_MAX);

  ok = TraceHamCollect(&X, want, -1, &csr);
  expect_true("gate: cap == want -> collect succeeds", ok == 1);
  if (ok) TraceHamFree(&csr);

  memset(&csr, 0x5A, sizeof(csr));
  ok = TraceHamCollect(&X, want - 1, -1, &csr);
  expect_true("gate: cap == want-1 -> collect demotes (returns 0)", ok == 0);
  expect_true("gate: cap == want-1 -> csr zeroed (rowptr NULL)", csr.rowptr == NULL);

  /* Allocation-failure injection: EVERY galloc TraceHamCollect performs must
     demote cleanly. TraceHamCollect (src/expec_trace_ham.c) does 6 fixed gallocs
     (0:rowptr 1:cursor 2:colidx 3:val 4:sort-ws 5:y) followed by n_diag
     coefficient-array gallocs (6..5+n_diag). The Hubbard fixture has n_diag==3,
     so the highest injectable index is 8 -- covering the diag[] allocations that
     the old 0..5 loop left untested. */
  {
    int max_fa = 5 + n_diag;
    expect_true("inject: fixture reaches the diag[] allocs (n_diag>0)", n_diag > 0);
    for (fa = 0; fa <= max_fa; fa++) {
      char nm[128];
      memset(&csr, 0x5A, sizeof(csr));
      ok = TraceHamCollect(&X, SIZE_MAX / 2, fa, &csr);
      snprintf(nm, sizeof(nm), "inject: fail_alloc_at=%d -> returns 0", fa);
      expect_true(nm, ok == 0);
      snprintf(nm, sizeof(nm), "inject: fail_alloc_at=%d -> csr zeroed", fa);
      expect_true(nm, csr.rowptr == NULL);
    }
  }

  /* Sink must be restored after every injected failure: a normal collect works. */
  ok = TraceHamCollect(&X, SIZE_MAX / 2, -1, &csr);
  expect_true("inject: subsequent normal collect succeeds (sink restored)", ok == 1);
  if (ok) TraceHamFree(&csr);
}

/* -----------------------------------------------------------------------
 * Part 2: energy-family diagonal coefficient arrays (csr->diag[]).
 *
 * For each eligible fixture we (a) load a fixed random normalized vector into
 * v0, (b) run the legacy expec_energy_flct() and snapshot its eight Phys
 * fluctuation fields, (c) run TraceHamCollect() and rebuild the SAME eight
 * fields from csr->diag[] applying the frozen-table scalings, and assert they
 * agree to 1e-12. The canonical-Spin fixture instead pins the stale-preserving
 * constant path (num_up/num_down untouched).
 * --------------------------------------------------------------------------*/

static double frand_pm1(void) { return 2.0 * ((double)rand() / (double)RAND_MAX) - 1.0; }

static void expect_close(const char *label, const char *field,
                         double a, double b, double tol) {
  double d = fabs(a - b);
  char nm[320];
  snprintf(nm, sizeof(nm), "%s: %s (csr=%.15g legacy=%.15g |d|=%.3e)",
           label, field, a, b, d);
  expect_true(nm, d <= tol);
}

static void run_coeff_fixture(const char *label, const char *subdir,
                              const char *stan, int patch_model) {
  struct BindStruct X;
  TraceHamCsr csr;
  long int n, i;
  double complex *xs;
  double nrm;
  int ok;
  double e_doublon, e_doublon2, e_num, e_num2, e_Sz, e_Sz2, e_num_up, e_num_down;
  double sumD = 0, sumD2 = 0, sumN = 0, sumN2 = 0, sumS = 0, sumS2 = 0;
  double c_doublon, c_doublon2, c_num, c_num2, c_Sz, c_Sz2, c_num_up, c_num_down;

  fprintf(stderr, "[coeff %s]\n", label);
  if (chdir(g_base) != 0) { expect_true(label, 0); return; }
  mkdir(subdir, 0777);
  if (chdir(subdir) != 0) { expect_true(label, 0); return; }

  if (setup_fixture(&X, stan, patch_model) != 0) {
    fprintf(stderr, "  FAIL %s: setup failed\n", label);
    g_failures++;
    return;
  }
  n = (long int)X.Check.idim_max;

  /* Fixed random normalized state into v0; keep an independent snapshot xs. */
  srand(20260721u);
  xs = (double complex *)malloc((size_t)(n + 1) * sizeof(double complex));
  nrm = 0.0;
  for (i = 1; i <= n; i++) {
    v0[i] = frand_pm1() + frand_pm1() * I;
    nrm += creal(conj(v0[i]) * v0[i]);
  }
  nrm = sqrt(nrm);
  for (i = 1; i <= n; i++) { v0[i] /= nrm; xs[i] = v0[i]; }

  /* Legacy evaluator (destroys v0); snapshot its eight fluctuation fields. */
  expec_energy_flct(&X);
  e_doublon  = X.Phys.doublon;
  e_doublon2 = X.Phys.doublon2;
  e_num      = X.Phys.num;
  e_num2     = X.Phys.num2;
  e_Sz       = X.Phys.Sz;
  e_Sz2      = X.Phys.Sz2;
  e_num_up   = X.Phys.num_up;
  e_num_down = X.Phys.num_down;

  ok = TraceHamCollect(&X, SIZE_MAX / 2, -1, &csr);
  {
    char nm[256];
    snprintf(nm, sizeof(nm), "%s: TraceHamCollect succeeds", label);
    expect_true(nm, ok == 1);
  }
  if (!ok) { free(xs); return; }

  for (i = 0; i < n; i++) {
    double w = creal(conj(xs[i + 1]) * xs[i + 1]);
    if (csr.n_diag == 3) {
      double D = csr.diag[0][i], N = csr.diag[1][i], S = csr.diag[2][i];
      sumD += w * D; sumD2 += w * D * D;
      sumN += w * N; sumN2 += w * N * N;
      sumS += w * S; sumS2 += w * S * S;
    } else if (csr.n_diag == 1) {
      double S = csr.diag[0][i];
      sumS += w * S; sumS2 += w * S * S;
    }
  }

  if (csr.n_diag == 3) {           /* Hubbard family / HubbardGC */
    c_doublon  = sumD;
    c_doublon2 = sumD2;
    c_num      = sumN;
    c_num2     = sumN2;
    c_Sz       = 0.5 * sumS;
    c_Sz2      = 0.25 * sumS2;
    c_num_up   = 0.5 * (sumN + sumS);
    c_num_down = 0.5 * (sumN - sumS);
  } else {                         /* SpinGC (n_diag == 1) */
    double Ns = (double)X.Def.NsiteMPI;
    c_doublon  = 0.0;
    c_doublon2 = 0.0;
    c_num      = Ns;
    c_num2     = Ns * Ns;
    c_Sz       = 0.5 * sumS;
    c_Sz2      = 0.25 * sumS2;
    c_num_up   = 0.5 * (Ns + sumS);
    c_num_down = 0.5 * (Ns - sumS);
  }

  expect_close(label, "doublon",  c_doublon,  e_doublon,  1e-12);
  expect_close(label, "doublon2", c_doublon2, e_doublon2, 1e-12);
  expect_close(label, "num",      c_num,      e_num,      1e-12);
  expect_close(label, "num2",     c_num2,     e_num2,     1e-12);
  expect_close(label, "Sz",       c_Sz,       e_Sz,       1e-12);
  expect_close(label, "Sz2",      c_Sz2,      e_Sz2,      1e-12);
  expect_close(label, "num_up",   c_num_up,   e_num_up,   1e-12);
  expect_close(label, "num_down", c_num_down, e_num_down, 1e-12);

  TraceHamFree(&csr);
  free(xs);
}

/* Canonical Spin: constant path must leave num_up/num_down UNWRITTEN. */
static void run_spin_sentinel(void) {
  struct BindStruct X;
  const char *stan =
    "L = 6\nmodel = \"Spin\"\nmethod = \"FullDiag\"\nlattice = \"chain\"\n"
    "J = 1.0\n2Sz = 0\n";
  long int n, i;

  fprintf(stderr, "[coeff Spin sentinel L=6]\n");
  if (chdir(g_base) != 0) { expect_true("spin-sentinel: chdir base", 0); return; }
  mkdir("spin_sentinel", 0777);
  if (chdir("spin_sentinel") != 0) { expect_true("spin-sentinel: chdir sub", 0); return; }

  if (setup_fixture(&X, stan, -1) != 0) {
    fprintf(stderr, "  FAIL spin-sentinel: setup failed\n");
    g_failures++;
    return;
  }
  n = (long int)X.Check.idim_max;

  srand(777u);
  for (i = 1; i <= n; i++) v0[i] = frand_pm1() + frand_pm1() * I;

  X.Phys.num_up   = 4321.0;   /* sentinels: the Spin path must not touch these */
  X.Phys.num_down = 8765.0;

  expec_energy_flct(&X);

  expect_true("spin-sentinel: num_up stale-preserved (==4321)",   X.Phys.num_up == 4321.0);
  expect_true("spin-sentinel: num_down stale-preserved (==8765)", X.Phys.num_down == 8765.0);
  expect_true("spin-sentinel: doublon == 0",  X.Phys.doublon == 0.0);
  expect_true("spin-sentinel: doublon2 == 0", X.Phys.doublon2 == 0.0);
  expect_true("spin-sentinel: num == NsiteMPI", X.Phys.num == (double)X.Def.NsiteMPI);
  expect_close("spin-sentinel", "Sz", X.Phys.Sz, 0.5 * (double)X.Def.Total2SzMPI, 1e-12);
}

/* -----------------------------------------------------------------------
 * Part 3: the streaming energy-family kernel TraceEnergyEvalState().
 *
 * The kernel READS the state from the global v1 and writes only Phys, so every
 * check here loads v1 (never v0) and reads back X.Phys. The independent
 * reference for the complex-Hermitian matrix is the legacy replicated dense
 * Ham[][] that setup_fixture built through the ordinary makeHam() path -- a
 * source of truth wholly separate from the kernel's CSR SpMV. Computing
 * x^dagger H x and |H x|^2 straight from Ham[][] (with explicit conj) is what
 * independently pins the kernel's complex conj/creal handling; transcribing
 * HPhi's internal transfer sign/index convention into a literal 4x4 would add a
 * second, fragile source of truth, so Ham[][] is used as the reference matrix.
 * --------------------------------------------------------------------------*/

/* Dense y = H x, energy = Re(x^dagger y), var = |H x|^2, straight from the
   1-based replicated Ham[][] -- fully independent of the CSR/kernel path. */
static void dense_energy_var(long int n, const double complex *x,
                             double *energy, double *var) {
  long int i, j;
  double complex e = 0.0;
  double v = 0.0;
  for (i = 1; i <= n; i++) {
    double complex yi = 0.0;
    for (j = 1; j <= n; j++) yi += Ham[i][j] * x[j];
    e += conj(x[i]) * yi;
    v += creal(conj(yi) * yi);
  }
  *energy = creal(e);
  *var = v;
}

/* Elementwise CSR-vs-Ham check (independent reference == makeHam dense). */
static int csr_matches_ham(const TraceHamCsr *csr, double *hmax_out) {
  long int n = csr->n, i, j;
  double complex *dc = csr_to_dense(csr);
  double hmax = 0.0;
  int mism = 0;
  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++) {
      double a = cabs(Ham[i][j]);
      if (a > hmax) hmax = a;
    }
  for (i = 0; i < n && mism == 0; i++)
    for (j = 0; j < n; j++)
      if (cabs(dc[i * n + j] - Ham[i + 1][j + 1]) > 1e-13 + 1e-13 * hmax) {
        mism = 1;
        break;
      }
  free(dc);
  if (hmax_out) *hmax_out = hmax;
  return mism == 0;
}

/* Complex-Hermitian reference: 2-site Hubbard, hopping t = 0.3 + 0.4i.
   Transfer term stored as value * c^dagger_{i,si} c_{j,sj}; each direction is
   listed with its Hermitian conjugate so makeHam builds a complex-Hermitian H
   (up spin = sigma 0, down = sigma 1). */
static void run_complex_energy(void) {
  struct BindStruct X;
  TraceHamCsr csr;
  long int n, i, s;
  int ok;
  const char *stan =
    "L = 2\nmodel = \"Hubbard\"\nmethod = \"FullDiag\"\nlattice = \"chain\"\n"
    "t = 1.0\nU = 4.0\nnelec = 2\n2Sz = 0\n";
  const char *trans =
    "======================== \n"
    "NTransfer       4  \n"
    "======================== \n"
    "========i_j_s_tijs====== \n"
    "======================== \n"
    "    1     0     0     0     0.300000000000000     0.400000000000000\n"
    "    0     0     1     0     0.300000000000000    -0.400000000000000\n"
    "    1     1     0     1     0.300000000000000     0.400000000000000\n"
    "    0     1     1     1     0.300000000000000    -0.400000000000000\n";

  fprintf(stderr, "[kernel complex-Hermitian energy: Hubbard L=2 t=0.3+0.4i]\n");
  if (chdir(g_base) != 0) { expect_true("kernel-cplx: chdir base", 0); return; }
  mkdir("kernel_cplx", 0777);
  if (chdir("kernel_cplx") != 0) { expect_true("kernel-cplx: chdir sub", 0); return; }

  if (setup_fixture_trans(&X, stan, -1, trans) != 0) {
    fprintf(stderr, "  FAIL kernel-cplx: setup failed\n"); g_failures++; return;
  }
  n = (long int)X.Check.idim_max;

  ok = TraceHamCollect(&X, SIZE_MAX / 2, -1, &csr);
  expect_true("kernel-cplx: TraceHamCollect succeeds", ok == 1);
  if (!ok) return;

  expect_true("kernel-cplx: CSR == dense Ham (<=1e-13, complex)",
              csr_matches_ham(&csr, NULL));
  /* the matrix must actually carry an imaginary part (guards a silent
     real-only regression that would make the conj test vacuous). */
  {
    double imax = 0.0;
    for (i = 1; i <= n; i++)
      for (s = 1; s <= n; s++)
        if (fabs(cimag(Ham[i][s])) > imax) imax = fabs(cimag(Ham[i][s]));
    expect_true("kernel-cplx: Ham has a non-zero imaginary part", imax > 0.1);
  }

  /* Analytic spot-check of makeHam's complex-element generation (spec 5.1):
     the CSR==Ham identity above shares its source (makeHam) with the reference,
     so a conjugation/phase bug common to BOTH would be invisible there. Here we
     pin the elements against the ANALYTIC value: every non-zero off-diagonal of
     this 2-site fixture is a single hop carrying the bare transfer t=0.3+0.4i up
     to a fermion sign, so it must equal +/-(0.3 +/- 0.4i) -- i.e. |Re|==0.3 and
     |Im|==0.4 exactly. Both imaginary signs must appear (the Hermitian
     conjugate partners Ham[i][j]=conj(Ham[j][i])), which catches a dropped or
     mis-signed imaginary part that the makeHam-vs-CSR check cannot see. */
  {
    int n_off = 0, seen_pos_im = 0, seen_neg_im = 0, bad_val = 0;
    for (i = 1; i <= n; i++)
      for (s = 1; s <= n; s++) {
        double re, im;
        if (i == s) continue;
        if (cabs(Ham[i][s]) < 1e-12) continue;
        re = creal(Ham[i][s]);
        im = cimag(Ham[i][s]);
        n_off++;
        if (fabs(fabs(re) - 0.3) > 1e-12 || fabs(fabs(im) - 0.4) > 1e-12) bad_val = 1;
        if (im > 0.0) seen_pos_im = 1;
        if (im < 0.0) seen_neg_im = 1;
      }
    expect_true("kernel-cplx: every off-diagonal == +/-(0.3+/-0.4i) analytically",
                n_off > 0 && bad_val == 0);
    expect_true("kernel-cplx: both conjugate imag signs present (+0.4 and -0.4)",
                seen_pos_im && seen_neg_im);
  }

  /* Three fixed non-eigenvector normalized states loaded into v1. */
  for (s = 0; s < 3; s++) {
    double complex xs[64];
    double nrm = 0.0, e_dense, var_dense;
    char nm[128];
    for (i = 1; i <= n; i++) {
      double re = cos(0.7 * (double)(i + 3 * s)) + 0.3 * (double)s;
      double im = sin(1.1 * (double)(i + 2 * s)) - 0.2 * (double)s;
      xs[i] = re + im * I;
      nrm += creal(conj(xs[i]) * xs[i]);
    }
    nrm = sqrt(nrm);
    for (i = 1; i <= n; i++) { xs[i] /= nrm; v1[i] = xs[i]; v0[i] = -777.0; }

    ok = (TraceEnergyEvalState(&X, &csr) == 0);
    snprintf(nm, sizeof(nm), "kernel-cplx[state %ld]: EvalState returns 0", s);
    expect_true(nm, ok);

    dense_energy_var(n, xs, &e_dense, &var_dense);
    snprintf(nm, sizeof(nm), "kernel-cplx[state %ld]: energy == x_dag H x", s);
    expect_close("kernel-cplx", nm, X.Phys.energy, e_dense, 1e-12);
    snprintf(nm, sizeof(nm), "kernel-cplx[state %ld]: var == |H x|^2", s);
    expect_close("kernel-cplx", nm, X.Phys.var, var_dense, 1e-12);

    /* v0 must be untouched by the kernel. */
    {
      int v0_ok = 1;
      for (i = 1; i <= n; i++) if (v0[i] != -777.0) v0_ok = 0;
      snprintf(nm, sizeof(nm), "kernel-cplx[state %ld]: v0 not written", s);
      expect_true(nm, v0_ok);
    }
    /* v1 must still hold the input state. */
    {
      int v1_ok = 1;
      for (i = 1; i <= n; i++) if (v1[i] != xs[i]) v1_ok = 0;
      snprintf(nm, sizeof(nm), "kernel-cplx[state %ld]: v1 preserved", s);
      expect_true(nm, v1_ok);
    }
  }
  TraceHamFree(&csr);
}

/* Cancellation: a duplicate +/-1e8 transfer pair (summing to 0) plus a normal
   hopping. The merged CSR entry must equal the dense-path sum exactly, i.e. the
   1e8's cancel in BOTH paths -- proven by the elementwise CSR==Ham identity and
   an O(1) post-cancellation |Ham|. */
static void run_cancellation(void) {
  struct BindStruct X;
  TraceHamCsr csr;
  long int n, i;
  int ok;
  double hmax = 0.0, e_dense, var_dense;
  const char *stan =
    "L = 2\nmodel = \"Hubbard\"\nmethod = \"FullDiag\"\nlattice = \"chain\"\n"
    "t = 1.0\nU = 4.0\nnelec = 2\n2Sz = 0\n";
  /* up hopping 0<->1 given as +1e8 and -1e8 (cancel to 0); down hopping = 1. */
  const char *trans =
    "======================== \n"
    "NTransfer       6  \n"
    "======================== \n"
    "========i_j_s_tijs====== \n"
    "======================== \n"
    "    1     0     0     0    100000000.000000000     0.000000000\n"
    "    0     0     1     0    100000000.000000000     0.000000000\n"
    "    1     0     0     0   -100000000.000000000     0.000000000\n"
    "    0     0     1     0   -100000000.000000000     0.000000000\n"
    "    1     1     0     1     1.000000000000000     0.000000000\n"
    "    0     1     1     1     1.000000000000000     0.000000000\n";

  fprintf(stderr, "[kernel cancellation: duplicate +/-1e8 transfer]\n");
  if (chdir(g_base) != 0) { expect_true("kernel-cancel: chdir base", 0); return; }
  mkdir("kernel_cancel", 0777);
  if (chdir("kernel_cancel") != 0) { expect_true("kernel-cancel: chdir sub", 0); return; }

  if (setup_fixture_trans(&X, stan, -1, trans) != 0) {
    fprintf(stderr, "  FAIL kernel-cancel: setup failed\n"); g_failures++; return;
  }
  n = (long int)X.Check.idim_max;

  ok = TraceHamCollect(&X, SIZE_MAX / 2, -1, &csr);
  expect_true("kernel-cancel: TraceHamCollect succeeds", ok == 1);
  if (!ok) return;

  /* Merged CSR entry == dense sum (both cancel the 1e8 pair to 0). */
  expect_true("kernel-cancel: CSR == dense Ham (merged == dense sum)",
              csr_matches_ham(&csr, &hmax));
  /* The 1e8 magnitudes really cancelled: nothing O(1e8) survives. */
  expect_true("kernel-cancel: post-cancellation |Ham| is O(1) (<100)", hmax < 100.0);

  {
    double complex xs[64];
    double nrm = 0.0;
    for (i = 1; i <= n; i++) {
      xs[i] = (0.5 + 0.1 * (double)i) + (0.2 - 0.05 * (double)i) * I;
      nrm += creal(conj(xs[i]) * xs[i]);
    }
    nrm = sqrt(nrm);
    for (i = 1; i <= n; i++) { xs[i] /= nrm; v1[i] = xs[i]; }
    ok = (TraceEnergyEvalState(&X, &csr) == 0);
    expect_true("kernel-cancel: EvalState returns 0", ok);
    dense_energy_var(n, xs, &e_dense, &var_dense);
    expect_close("kernel-cancel", "energy == dense", X.Phys.energy, e_dense, 1e-12);
    expect_close("kernel-cancel", "var == dense", X.Phys.var, var_dense, 1e-12);
  }
  TraceHamFree(&csr);
}

/* Sentinels: canonical Spin's constant path must leave num_up/num_down (the
   frozen NOT-WRITTEN cells) untouched, while energy/var and the constant
   fluctuation fields are overwritten. */
static void run_energy_sentinel(void) {
  struct BindStruct X;
  TraceHamCsr csr;
  long int n, i;
  int ok;
  double e_dense, var_dense;
  const char *stan =
    "L = 6\nmodel = \"Spin\"\nmethod = \"FullDiag\"\nlattice = \"chain\"\n"
    "J = 1.0\n2Sz = 0\n";

  fprintf(stderr, "[kernel sentinels: canonical Spin L=6]\n");
  if (chdir(g_base) != 0) { expect_true("kernel-sentinel: chdir base", 0); return; }
  mkdir("kernel_sentinel", 0777);
  if (chdir("kernel_sentinel") != 0) { expect_true("kernel-sentinel: chdir sub", 0); return; }

  if (setup_fixture(&X, stan, -1) != 0) {
    fprintf(stderr, "  FAIL kernel-sentinel: setup failed\n"); g_failures++; return;
  }
  n = (long int)X.Check.idim_max;

  ok = TraceHamCollect(&X, SIZE_MAX / 2, -1, &csr);
  expect_true("kernel-sentinel: TraceHamCollect succeeds", ok == 1);
  if (!ok) return;
  expect_true("kernel-sentinel: canonical Spin has n_diag == 0", csr.n_diag == 0);

  /* Fixed normalized real-ish state into v1. */
  {
    double nrm = 0.0;
    srand(31337u);
    for (i = 1; i <= n; i++) v1[i] = frand_pm1() + frand_pm1() * I;
    for (i = 1; i <= n; i++) nrm += creal(conj(v1[i]) * v1[i]);
    nrm = sqrt(nrm);
    for (i = 1; i <= n; i++) v1[i] /= nrm;
  }

  /* Distinct sentinels in every Phys field the kernel might touch. */
  X.Phys.energy   = 1111.0;
  X.Phys.var      = 2222.0;
  X.Phys.doublon  = 3333.0;
  X.Phys.doublon2 = 4444.0;
  X.Phys.num      = 5555.0;
  X.Phys.num2     = 6666.0;
  X.Phys.Sz       = 7777.0;
  X.Phys.Sz2      = 8888.0;
  X.Phys.num_up   = 4321.0;
  X.Phys.num_down = 8765.0;

  ok = (TraceEnergyEvalState(&X, &csr) == 0);
  expect_true("kernel-sentinel: EvalState returns 0", ok);

  /* Frozen NOT-WRITTEN cells retain their sentinels. */
  expect_true("kernel-sentinel: num_up UNWRITTEN (==4321)",   X.Phys.num_up == 4321.0);
  expect_true("kernel-sentinel: num_down UNWRITTEN (==8765)", X.Phys.num_down == 8765.0);
  /* Written constant cells took their frozen values (not the sentinels). */
  expect_true("kernel-sentinel: doublon written 0",  X.Phys.doublon == 0.0);
  expect_true("kernel-sentinel: doublon2 written 0", X.Phys.doublon2 == 0.0);
  expect_true("kernel-sentinel: num == NsiteMPI", X.Phys.num == (double)X.Def.NsiteMPI);
  expect_true("kernel-sentinel: num2 == NsiteMPI^2",
              X.Phys.num2 == (double)X.Def.NsiteMPI * (double)X.Def.NsiteMPI);
  expect_close("kernel-sentinel", "Sz == 0.5*Total2SzMPI",
               X.Phys.Sz, 0.5 * (double)X.Def.Total2SzMPI, 1e-12);
  expect_true("kernel-sentinel: energy overwritten (!= sentinel)", X.Phys.energy != 1111.0);
  expect_true("kernel-sentinel: var overwritten (!= sentinel)",    X.Phys.var != 2222.0);
  /* var == <H^2> checked against the independent dense reference over the SAME
     state (v1), not merely >= 0 (which is tautological for a sum of squares). */
  dense_energy_var(n, v1, &e_dense, &var_dense);
  expect_close("kernel-sentinel", "energy == dense <H>", X.Phys.energy, e_dense, 1e-12);
  expect_close("kernel-sentinel", "var == dense <H^2>",  X.Phys.var,    var_dense, 1e-12);

  TraceHamFree(&csr);
}

/* Kernel field-write scalings vs the production evaluator (n_diag==3).
 *
 * Part 2 checked the diag[] contents; this checks the KERNEL's own eight
 * fluctuation-field writes -- the scaling logic duplicated from
 * expec_energy_flct() into TraceEnergyEvalState (0.25 vs 0.5, the +/- in
 * num_up/num_down). Run the real expec_energy_flct() on a random normalized
 * state in v0 and snapshot its eight fields (plus energy/var); run the kernel on
 * the SAME state in v1; assert they agree. A swapped factor or transposed sign
 * would surface here directly, ahead of Task 6's end-to-end 1e-8 script. */
static void run_kernel_flct_equiv(void) {
  struct BindStruct X;
  TraceHamCsr csr;
  long int n, i;
  double complex *xs;
  double nrm;
  int ok;
  double e_doublon, e_doublon2, e_num, e_num2, e_Sz, e_Sz2, e_num_up, e_num_down;
  double e_energy, e_var;
  /* HubbardGC (Sz NOT conserved) so the basis carries a spread of S(k): a
     canonical 2Sz=0 sector would give sumS==0, making Sz2's 0.25 factor and the
     +/- in num_up/num_down vacuous (num_up==num_down==0.5*num). Here sumS!=0, so
     a swapped 0.25<->0.5 or a transposed num_up/down sign is caught. */
  const char *stan =
    "L = 2\nmodel = \"HubbardGC\"\nmethod = \"FullDiag\"\nlattice = \"chain\"\n"
    "t = 1.0\nU = 4.0\n";

  fprintf(stderr, "[kernel flct-equiv: HubbardGC L=2 vs expec_energy_flct]\n");
  if (chdir(g_base) != 0) { expect_true("kernel-flct: chdir base", 0); return; }
  mkdir("kernel_flct", 0777);
  if (chdir("kernel_flct") != 0) { expect_true("kernel-flct: chdir sub", 0); return; }

  if (setup_fixture(&X, stan, -1) != 0) {
    fprintf(stderr, "  FAIL kernel-flct: setup failed\n"); g_failures++; return;
  }
  n = (long int)X.Check.idim_max;

  /* Fixed random normalized state into v0; snapshot xs (v0 is destroyed by
     expec_energy_flct, which moves v0->v1 before computing energy). */
  srand(24680u);
  xs = (double complex *)malloc((size_t)(n + 1) * sizeof(double complex));
  nrm = 0.0;
  for (i = 1; i <= n; i++) {
    v0[i] = frand_pm1() + frand_pm1() * I;
    nrm += creal(conj(v0[i]) * v0[i]);
  }
  nrm = sqrt(nrm);
  for (i = 1; i <= n; i++) { v0[i] /= nrm; xs[i] = v0[i]; }

  expec_energy_flct(&X);
  e_doublon  = X.Phys.doublon;   e_doublon2 = X.Phys.doublon2;
  e_num      = X.Phys.num;       e_num2     = X.Phys.num2;
  e_Sz       = X.Phys.Sz;        e_Sz2      = X.Phys.Sz2;
  e_num_up   = X.Phys.num_up;    e_num_down = X.Phys.num_down;
  e_energy   = X.Phys.energy;    e_var      = X.Phys.var;

  ok = TraceHamCollect(&X, SIZE_MAX / 2, -1, &csr);
  expect_true("kernel-flct: TraceHamCollect succeeds", ok == 1);
  if (!ok) { free(xs); return; }
  expect_true("kernel-flct: HubbardGC has n_diag == 3", csr.n_diag == 3);
  /* Non-vacuous guard: the +/- sign path is only exercised when sumS!=0, i.e.
     num_up != num_down (their difference IS sumS). */
  expect_true("kernel-flct: fixture exercises Sz (num_up != num_down)",
              fabs(e_num_up - e_num_down) > 1e-6);
  expect_true("kernel-flct: fixture exercises Sz2 (Sz2 > 0)", e_Sz2 > 1e-9);

  for (i = 1; i <= n; i++) v1[i] = xs[i];
  ok = (TraceEnergyEvalState(&X, &csr) == 0);
  expect_true("kernel-flct: EvalState returns 0", ok);

  /* The kernel's own eight field writes must match the legacy evaluator. */
  expect_close("kernel-flct", "doublon",  X.Phys.doublon,  e_doublon,  1e-12);
  expect_close("kernel-flct", "doublon2", X.Phys.doublon2, e_doublon2, 1e-12);
  expect_close("kernel-flct", "num",      X.Phys.num,      e_num,      1e-12);
  expect_close("kernel-flct", "num2",     X.Phys.num2,     e_num2,     1e-12);
  expect_close("kernel-flct", "Sz",       X.Phys.Sz,       e_Sz,       1e-12);
  expect_close("kernel-flct", "Sz2",      X.Phys.Sz2,      e_Sz2,      1e-12);
  expect_close("kernel-flct", "num_up",   X.Phys.num_up,   e_num_up,   1e-12);
  expect_close("kernel-flct", "num_down", X.Phys.num_down, e_num_down, 1e-12);
  /* Bonus: the kernel's SpMV energy/var also match the mltply evaluator. */
  expect_close("kernel-flct", "energy",   X.Phys.energy,   e_energy,   1e-9);
  expect_close("kernel-flct", "var",      X.Phys.var,      e_var,      1e-9);

  TraceHamFree(&csr);
  free(xs);
}

/* -----------------------------------------------------------------------
 * Part 4: TraceFinalizeEnergyPlan() -- the build->finalize energy-plan
 * lifecycle (phase 3c Task 5). The function lives in the MPI-capable
 * orchestration TU src/phys_distributed.c (linked into this test); it performs
 * ONE MPI_Allreduce(MIN) over {ok, nnz, -nnz} and applies the same verdict on
 * every rank.
 *
 * nproc == 1 (serial / noMPI build, or an MPI singleton): the single-rank
 * buffer trivially passes, so we pin BOTH pass-through outcomes -- a success
 * keeps the energy kernel, a failure demotes it (memory reason).
 *
 * nproc >= 2 (launched under mpiexec): rank 1 injects an allocation failure in
 * its collect (fail_alloc_at=0), every other rank collects normally; after
 * TraceFinalizeEnergyPlan() EVERY rank must end demoted
 * (kernel[TRACE_Q_ENERGY]==0, demoted_memory==1). The collect runs with a
 * temporarily-replicated basis (nproc faked to 1 across setup+collect only) so
 * each rank builds an identical, non-distributed fixture and setup's
 * nproc-guarded MPI calls stay no-ops; the true nproc is restored before the
 * collective finalize, whose MPI_Allreduce spans the real ranks via
 * MPI_COMM_WORLD regardless of the global. The nnz values are irrelevant here
 * (any failed rank makes reduced[0]==0, short-circuiting the nnz-agreement
 * check before it is reached). */
static void run_finalize_check(void) {
  fprintf(stderr, "[finalize: energy build->finalize lifecycle]\n");

  if (nproc == 1) {
    TraceExecutionPlan plan;

    memset(&plan, 0, sizeof(plan));
    plan.kernel[TRACE_Q_ENERGY] = 1;
    TraceFinalizeEnergyPlan(&plan, 1, 12345L);
    expect_true("finalize np=1: local_ok=1 keeps the energy kernel",
                plan.kernel[TRACE_Q_ENERGY] == 1);
    expect_true("finalize np=1: local_ok=1 leaves demoted_memory clear",
                plan.demoted_memory[TRACE_Q_ENERGY] == 0);

    memset(&plan, 0, sizeof(plan));
    plan.kernel[TRACE_Q_ENERGY] = 1;
    TraceFinalizeEnergyPlan(&plan, 0, 0L);
    expect_true("finalize np=1: local_ok=0 demotes the energy kernel",
                plan.kernel[TRACE_Q_ENERGY] == 0);
    expect_true("finalize np=1: local_ok=0 sets demoted_memory",
                plan.demoted_memory[TRACE_Q_ENERGY] == 1);
    return;
  }

  {
    struct BindStruct X;
    TraceHamCsr csr;
    TraceExecutionPlan plan;
    char sub[64];
    int local_ok, saved_nproc, saved_myrank;
    long int nnz_raw;
    const char *stan =
      "L = 4\nmodel = \"Hubbard\"\nmethod = \"FullDiag\"\nlattice = \"chain\"\n"
      "t = 1.0\nU = 4.0\nnelec = 4\n2Sz = 0\n";

    snprintf(sub, sizeof(sub), "finalize_rank%d", myrank);
    if (chdir(g_base) != 0) { expect_true("finalize: chdir base", 0); return; }
    mkdir(sub, 0777);
    if (chdir(sub) != 0) { expect_true("finalize: chdir sub", 0); return; }

    memset(&csr, 0, sizeof(csr));

    /* Build+collect with a replicated basis (fake nproc=1); restore before the
       collective finalize. */
    saved_nproc = nproc;
    saved_myrank = myrank;
    nproc = 1;
    myrank = 0;
    if (setup_fixture(&X, stan, -1) != 0) {
      fprintf(stderr, "  FAIL finalize: setup failed (rank %d)\n", saved_myrank);
      g_failures++;
      nproc = saved_nproc;
      myrank = saved_myrank;
      return;
    }
    if (saved_myrank == 1) {
      local_ok = TraceHamCollect(&X, SIZE_MAX / 2, 0, &csr); /* forced failure */
      nnz_raw = 0;
    } else {
      local_ok = TraceHamCollect(&X, SIZE_MAX / 2, -1, &csr);
      nnz_raw = local_ok ? csr.nnz : 0;
    }
    nproc = saved_nproc;
    myrank = saved_myrank;

    memset(&plan, 0, sizeof(plan));
    plan.kernel[TRACE_Q_ENERGY] = 1; /* provisional, as TraceBuildPlan() leaves it */
    TraceFinalizeEnergyPlan(&plan, local_ok, nnz_raw);

    expect_true("finalize np>=2: rank1-fail demotes energy on THIS rank",
                plan.kernel[TRACE_Q_ENERGY] == 0);
    expect_true("finalize np>=2: rank1-fail sets demoted_memory on THIS rank",
                plan.demoted_memory[TRACE_Q_ENERGY] == 1);
    if (local_ok) TraceHamFree(&csr);
  }
}

/* -----------------------------------------------------------------------
 * Part 5: TraceModelEnergySupported() -- the positive supported-model
 * predicate that gates the energy trace kernel (src/expec_trace_ham.c).
 *
 * Pure function of the CalcModel enum: no MPI, no makeHam, no fixture. The
 * energy kernel implements the per-basis fluctuation semantics of exactly the
 * Hubbard/Kondo/tJ family (canonical + GC + N-conserved), SpinGC, and
 * canonical Spin; SpinlessFermion(GC) and any unknown model must be REJECTED
 * (their num/Sz semantics differ and are not implemented in the kernel, so a
 * true here would silently activate a wrong-value kernel -- the blocker this
 * predicate fixes). Runs in every build (MPI and noMPI).
 * --------------------------------------------------------------------------*/
static void run_energy_supported_predicate(void) {
  int i;
  int supported[] = { Hubbard, HubbardNConserved, Kondo, KondoNConserved,
                      tJ, tJNConserved, HubbardGC, KondoGC, tJGC, SpinGC, Spin };
  int unsupported[] = { SpinlessFermion, SpinlessFermionGC, 999 };

  fprintf(stderr, "[energy-supported predicate: TraceModelEnergySupported]\n");
  for (i = 0; i < (int)(sizeof(supported) / sizeof(supported[0])); i++) {
    char nm[96];
    snprintf(nm, sizeof(nm), "energy-supported: model %d supported (==1)", supported[i]);
    expect_true(nm, TraceModelEnergySupported(supported[i]) == 1);
  }
  for (i = 0; i < (int)(sizeof(unsupported) / sizeof(unsupported[0])); i++) {
    char nm[96];
    snprintf(nm, sizeof(nm), "energy-supported: model %d NOT supported (==0)", unsupported[i]);
    expect_true(nm, TraceModelEnergySupported(unsupported[i]) == 0);
  }
}

/* -----------------------------------------------------------------------
 * Part 6: end-to-end plan build + INFO reporting for the energy slot
 * (TraceBuildPlan / TraceReportPlan, src/expec_trace.c). Pure: no MPI, no
 * makeHam, no fixture -- a minimal BindStruct with only the Def fields
 * TraceBuildPlan reads. This is the local non-vacuousness proof for the
 * blocker-1 INFO wording: a SUPPORTED model keeps the energy kernel and prints
 * "the energy/fluctuation family uses the trace kernel.", while
 * SpinlessFermion is DEMOTED (demoted_unsupported_model set, kernel cleared)
 * and prints "... uses the ExpecMode-1 fallback (unsupported model).". The two
 * INFO strings are byte-exact from the shell equivalence script's greps.
 * --------------------------------------------------------------------------*/
static void build_plan_for_model(int model, TraceExecutionPlan *plan) {
  struct BindStruct X;
  memset(&X, 0, sizeof(X));
  X.Def.iExpecMode = EXPECMODE_TRACE;
  X.Def.iCalcModel = model;
  X.Def.iFlgGeneralSpin = FALSE;
  X.Def.iInputHam = 0;
  X.Def.NCisAjt = 0;
  X.Def.NCisAjtCkuAlvDC = 0;
  TraceBuildPlan(&X, 1, ((size_t)1) << 30, plan);
}

static int report_contains(const TraceExecutionPlan *plan, const char *needle) {
  char buf[4096];
  size_t nread;
  FILE *fp = tmpfile();
  int found;
  if (fp == NULL) return 0;
  TraceReportPlan(plan, fp);
  rewind(fp);
  nread = fread(buf, 1, sizeof(buf) - 1, fp);
  buf[nread] = '\0';
  fclose(fp);
  found = (strstr(buf, needle) != NULL);
  return found;
}

static void run_energy_plan_report(void) {
  TraceExecutionPlan plan;
  /* Match the FULL INFO line (the "INFO: ExpecMode 2: " prefix through the
     terminal period), not a loose substring, so a reworded prefix/suffix cannot
     slip past report_contains()'s strstr. Byte-exact from TraceReportPlan(). */
  const char *kernel_line =
    "INFO: ExpecMode 2: the energy/fluctuation family uses the trace kernel.";
  const char *unsupported_line =
    "INFO: ExpecMode 2: the energy/fluctuation family uses the ExpecMode-1 fallback (unsupported model).";

  fprintf(stderr, "[energy plan+report: TraceBuildPlan/TraceReportPlan]\n");

  /* Supported model (Hubbard): energy kernel stays active. */
  build_plan_for_model(Hubbard, &plan);
  expect_true("plan Hubbard: energy kernel active",
              plan.kernel[TRACE_Q_ENERGY] == 1);
  expect_true("plan Hubbard: not demoted-unsupported",
              plan.demoted_unsupported_model[TRACE_Q_ENERGY] == 0);
  expect_true("report Hubbard: prints the energy trace-kernel line",
              report_contains(&plan, kernel_line));

  /* Unsupported model (canonical SpinlessFermion): energy kernel demoted. */
  build_plan_for_model(SpinlessFermion, &plan);
  expect_true("plan SpinlessFermion: energy kernel demoted (kernel==0)",
              plan.kernel[TRACE_Q_ENERGY] == 0);
  expect_true("plan SpinlessFermion: demoted_unsupported_model set",
              plan.demoted_unsupported_model[TRACE_Q_ENERGY] == 1);
  expect_true("plan SpinlessFermion: NOT demoted_input_ham",
              plan.demoted_input_ham[TRACE_Q_ENERGY] == 0);
  expect_true("report SpinlessFermion: prints the unsupported-model fallback line",
              report_contains(&plan, unsupported_line));
  /* Non-vacuous: the demoted plan must NOT print the trace-kernel line. */
  expect_true("report SpinlessFermion: does NOT print the trace-kernel line",
              !report_contains(&plan, kernel_line));

  /* SpinlessFermionGC: same demotion. */
  build_plan_for_model(SpinlessFermionGC, &plan);
  expect_true("plan SpinlessFermionGC: energy kernel demoted (kernel==0)",
              plan.kernel[TRACE_Q_ENERGY] == 0);
  expect_true("plan SpinlessFermionGC: demoted_unsupported_model set",
              plan.demoted_unsupported_model[TRACE_Q_ENERGY] == 1);
}

int main(int argc, char **argv) {
  stdoutMPI = stdout;
  myrank = 0;
  nproc = 1;
#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  /* Allocate the timer arrays that the full binary sets up in InitTimer()
     (called from HPhiMain before any StartTimer). In an MPI build StartTimer()/
     StopTimer() dereference the global Timer[]/TimerStart[] pointers (NULL until
     InitTimer runs); the parts 1-3 fixtures reach them via
     expec_energy_flct(). InitTimer's body is a no-op in the noMPI build. */
  InitTimer();
#else
  (void)argc; (void)argv;
#endif

#ifdef _OPENMP
  /* The fixtures are tiny (<= 4 sites); with the default thread count on a
     large node (e.g. 64), sz()'s OpenMP team barriers livelock/stall when the
     directly-executed singleton binary shares the machine with the BLAS
     thread pool (observed on a 64-core host: main thread stuck spinning in
     gomp_team_barrier_wait_end inside calculate_jb_* with num_threads=64,
     alongside 128 idle OpenBLAS workers). One thread is ample here and keeps
     every OpenMP code path exercised (serially). */
  omp_set_num_threads(1);
#endif

  if (getcwd(g_base, sizeof(g_base)) == NULL) {
    fprintf(stderr, "getcwd failed\n");
    return 1;
  }
  strncat(g_base, "/expec_trace_ham_scratch", sizeof(g_base) - strlen(g_base) - 1);
  mkdir(g_base, 0777);

  fprintf(stderr, "== expec_trace_ham_check == (rank %d / %d)\n", myrank, nproc);

  /* Parts 1-3 are the single-process coverage: single-process
     StdFace/ReadDef/makeHam fixtures plus the streaming-kernel checks. They run
     whenever nproc == 1 -- in the noMPI build, and also in an MPI build invoked
     as a singleton (ctest runs the binary directly; OpenMPI/MPICH give it a
     one-rank MPI_COMM_WORLD). The one thing the singleton needs that the full
     HPhi binary does via InitTimer() -- allocating Timer[]/TimerStart[] before
     the StartTimer() calls reached through expec_energy_flct() -- is done in
     main() above. When nproc >= 2 (the mpiexec launcher) these single-process
     fixtures are skipped and only part 4 (the rank-synchronized finalize) plus
     the pure parts 5-6 run. */
  if (nproc == 1) {
  fprintf(stderr, "== Parts 1-3 (single-process fixtures) running (nproc==1) ==\n");

  run_matrix_fixture("Hubbard L=4 half-filled", "hubbard",
    "L = 4\nmodel = \"Hubbard\"\nmethod = \"FullDiag\"\nlattice = \"chain\"\n"
    "t = 1.0\nU = 4.0\nnelec = 4\n2Sz = 0\n", -1);

  run_matrix_fixture("HubbardGC L=4", "hubbardgc",
    "L = 4\nmodel = \"HubbardGC\"\nmethod = \"FullDiag\"\nlattice = \"chain\"\n"
    "t = 1.0\nU = 4.0\n", -1);

  run_matrix_fixture("tJ L=4 (2 holes)", "tj",
    "L = 4\nmodel = \"Hubbard\"\nmethod = \"FullDiag\"\nlattice = \"chain\"\n"
    "t = 1.0\nU = 4.0\nnelec = 2\n2Sz = 0\n", tJ);

  run_matrix_fixture("tJGC L=4", "tjgc",
    "L = 4\nmodel = \"HubbardGC\"\nmethod = \"FullDiag\"\nlattice = \"chain\"\n"
    "t = 1.0\nU = 4.0\n", tJGC);

  run_matrix_fixture("Kondo chain L=2", "kondo",
    "L = 2\nmodel = \"Kondo\"\nmethod = \"FullDiag\"\nlattice = \"chain\"\n"
    "t = 1.0\nJ = 4.0\nnelec = 2\n2Sz = 0\n", -1);

  run_matrix_fixture("KondoGC chain L=2", "kondogc",
    "L = 2\nmodel = \"KondoGC\"\nmethod = \"FullDiag\"\nlattice = \"chain\"\n"
    "t = 1.0\nJ = 4.0\n", -1);

  run_matrix_fixture("Spin-1/2 L=6 Sz=0", "spin",
    "L = 6\nmodel = \"Spin\"\nmethod = \"FullDiag\"\nlattice = \"chain\"\n"
    "J = 1.0\n2Sz = 0\n", -1);

  run_matrix_fixture("Spin S=1 L=4 Sz=0", "genspin",
    "L = 4\nmodel = \"Spin\"\nmethod = \"FullDiag\"\nlattice = \"chain\"\n"
    "J = 1.0\n2S = 2\n2Sz = 0\n", -1);

  run_matrix_fixture("SpinGC-1/2 L=6 Gamma=0.5", "spingc",
    "L = 6\nmodel = \"SpinGC\"\nmethod = \"FullDiag\"\nlattice = \"chain\"\n"
    "J = 1.0\nGamma = 0.5\n", -1);

  run_matrix_fixture("SpinGC S=1 L=4", "genspingc",
    "L = 4\nmodel = \"SpinGC\"\nmethod = \"FullDiag\"\nlattice = \"chain\"\n"
    "J = 1.0\n2S = 2\n", -1);

  run_gate_and_injection();

  /* ---- Part 2: diagonal coefficient arrays (csr->diag[]). ---- */
  run_coeff_fixture("Hubbard L=4 half-filled", "coeff_hubbard",
    "L = 4\nmodel = \"Hubbard\"\nmethod = \"FullDiag\"\nlattice = \"chain\"\n"
    "t = 1.0\nU = 4.0\nnelec = 4\n2Sz = 0\n", -1);

  run_coeff_fixture("HubbardGC L=4", "coeff_hubbardgc",
    "L = 4\nmodel = \"HubbardGC\"\nmethod = \"FullDiag\"\nlattice = \"chain\"\n"
    "t = 1.0\nU = 4.0\n", -1);

  run_coeff_fixture("tJGC L=4", "coeff_tjgc",
    "L = 4\nmodel = \"HubbardGC\"\nmethod = \"FullDiag\"\nlattice = \"chain\"\n"
    "t = 1.0\nU = 4.0\n", tJGC);

  run_coeff_fixture("KondoGC chain L=2", "coeff_kondogc",
    "L = 2\nmodel = \"KondoGC\"\nmethod = \"FullDiag\"\nlattice = \"chain\"\n"
    "t = 1.0\nJ = 4.0\n", -1);

  run_coeff_fixture("SpinGC-1/2 L=6 Gamma=0.5", "coeff_spingc",
    "L = 6\nmodel = \"SpinGC\"\nmethod = \"FullDiag\"\nlattice = \"chain\"\n"
    "J = 1.0\nGamma = 0.5\n", -1);

  run_coeff_fixture("SpinGC S=1 L=4", "coeff_genspingc",
    "L = 4\nmodel = \"SpinGC\"\nmethod = \"FullDiag\"\nlattice = \"chain\"\n"
    "J = 1.0\n2S = 2\n", -1);

  run_spin_sentinel();

  /* ---- Part 3: streaming energy-family kernel (TraceEnergyEvalState). ---- */
  run_complex_energy();
  run_cancellation();
  run_energy_sentinel();
  run_kernel_flct_equiv();

  } else {
    fprintf(stderr, "== Parts 1-3 SKIPPED (nproc==%d > 1; single-process only) ==\n",
            nproc);
  } /* if (nproc == 1) : parts 1-3 are single-process coverage */

  /* ---- Part 4: TraceFinalizeEnergyPlan (single-rank at nproc==1, the
     rank-synchronized 2-rank demotion under mpiexec). Collective at nproc>=2,
     so EVERY rank must reach it. ---- */
  run_finalize_check();

  /* ---- Part 5: pure supported-model predicate (both MPI and noMPI). ---- */
  run_energy_supported_predicate();

  /* ---- Part 6: plan build + INFO reporting for the energy slot. ---- */
  run_energy_plan_report();

  if (chdir(g_base) == 0) chdir("..");

  {
    int failed = (g_failures != 0);
#ifdef MPI
    /* A failure on ANY rank fails the whole test (rank 0 returns the verdict). */
    int any_failed = failed;
    MPI_Allreduce(MPI_IN_PLACE, &any_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    failed = any_failed;
    MPI_Finalize();
    if (myrank != 0) return failed ? 1 : 0;
#endif
    if (!failed) {
      fprintf(stderr, "ALL PASS\n");
      return 0;
    }
    fprintf(stderr, "FAILURES: %d\n", g_failures);
    return 1;
  }
}
