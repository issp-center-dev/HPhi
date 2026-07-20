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
#include "expec_trace_ham.h"

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

/* Drive the full setup for one fixture (cwd already inside its scratch dir).
   Returns 0 on success with X populated + dense Ham built; -1 on setup error. */
static int setup_fixture(struct BindStruct *X, const char *stan, int patch_model) {
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

  /* Allocation-failure injection: each of the first six allocations. */
  for (fa = 0; fa <= 5; fa++) {
    char nm[128];
    memset(&csr, 0x5A, sizeof(csr));
    ok = TraceHamCollect(&X, SIZE_MAX / 2, fa, &csr);
    snprintf(nm, sizeof(nm), "inject: fail_alloc_at=%d -> returns 0", fa);
    expect_true(nm, ok == 0);
    snprintf(nm, sizeof(nm), "inject: fail_alloc_at=%d -> csr zeroed", fa);
    expect_true(nm, csr.rowptr == NULL);
  }

  /* Sink must be restored after every injected failure: a normal collect works. */
  ok = TraceHamCollect(&X, SIZE_MAX / 2, -1, &csr);
  expect_true("inject: subsequent normal collect succeeds (sink restored)", ok == 1);
  if (ok) TraceHamFree(&csr);
}

int main(void) {
  stdoutMPI = stdout;
  myrank = 0;
  nproc = 1;

  if (getcwd(g_base, sizeof(g_base)) == NULL) {
    fprintf(stderr, "getcwd failed\n");
    return 1;
  }
  strncat(g_base, "/expec_trace_ham_scratch", sizeof(g_base) - strlen(g_base) - 1);
  mkdir(g_base, 0777);

  fprintf(stderr, "== expec_trace_ham_check ==\n");

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

  if (chdir(g_base) == 0) chdir("..");

  if (g_failures == 0) {
    fprintf(stderr, "ALL PASS\n");
    return 0;
  }
  fprintf(stderr, "FAILURES: %d\n", g_failures);
  return 1;
}
