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
 * Serial (no-MPI) unit test for the phase-3b Task-2 mapping-probe extraction
 * drivers (TraceMapExtractOneBody/TwoBody, src/expec_trace.c) and the
 * co-located `*_TraceProbe` adapters (src/mltplyHubbardCore.c /
 * src/mltplySpinCore.c). GC models only (HubbardGC L=4 -> n=256, SpinGC-half
 * L=4 -> n=16), so the grand-canonical bare-bit basis (index == bit pattern)
 * needs no list_1/list_2 and the whole test runs in one process.
 *
 * The three checks required by the plan (docs/.../2026-07-12-...-phase3b.md
 * Task 2 Step 3):
 *   (1) mapping validity  -- for a fixed random complex vector z, the streamed
 *       sum  sum_{k:kprime>=0} conj(z[kprime]) * amp * z[k]  equals the direct
 *       sum of the ORIGINAL element functions' dam_pr (M_CORR, vec=z) to 1e-13.
 *       One-body: diagonal, off-diagonal, and a zero-result operator.
 *       Two-body: all four element-family branches + same-index + zero-result.
 *   (2) purity -- double extraction is bit-identical (kprime and amp), and the
 *       snapshot fields the drivers may touch (all of X->Large; the extraction
 *       write-set per inventory 2c) are byte-restored. GC paths read no global
 *       array (list_1/list_2 unused), so X->Large is the whole write-set.
 *   (3) boundary -- a zero operator streams to 0, and the Rearray-irregular
 *       SpinGC pair returns the n==0 sentinel (Mode-1's "write a 0.0 row").
 *
 * Reference model: the "direct" value calls the UNMODIFIED original element
 * functions the Mode-1 path calls, so the test pins the driver's dispatch,
 * amplitude assembly (tmp_V folding, sign), and streaming index/conjugation
 * against the ground truth that the 18/18 regression already validates.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "struct.h"
#include "DefCommon.h"
#include "mltplyCommon.h"
#include "mltplyHubbardCore.h"
#include "mltplySpinCore.h"
#include "bitcalc.h"
#include "global.h"
#include "expec_trace.h"
#include "expec_trace_internal.h"
#include "rearray_interactions.h"
#include <stdlib.h> /* setenv/unsetenv (POSIX; matches this test's other libc use) */

/* `myrank` and stdoutMPI are provided by the linked src/global.c;
   Rearray_Interactions is the REAL definition from the linked
   src/rearray_interactions.c (same TU the production driver uses). */

/* --------------------------------------------------------------------------
 * Test harness
 * -------------------------------------------------------------------------- */
static int g_failures = 0;
#define TOL 1e-13

static void expect_close(const char *name, double complex a, double complex b) {
  double d = cabs(a - b);
  if (d > TOL) {
    fprintf(stderr, "  FAIL %-40s direct=(% .12e,% .12e) streamed=(% .12e,% .12e) |diff|=%.3e\n",
            name, creal(a), cimag(a), creal(b), cimag(b), d);
    g_failures++;
  } else {
    fprintf(stderr, "  ok   %-40s |diff|=%.3e\n", name, d);
  }
}
static void expect_true(const char *name, int cond) {
  if (!cond) { fprintf(stderr, "  FAIL %s\n", name); g_failures++; }
  else fprintf(stderr, "  ok   %s\n", name);
}

/* Deterministic pseudo-random complex vector, 1-based [1..n], vec[0]=0. */
static void fill_random(double complex *vec, long int n, unsigned int seed) {
  long int j;
  unsigned long s = seed;
  vec[0] = 0.0;
  for (j = 1; j <= n; j++) {
    double re, im;
    s = s * 6364136223846793005UL + 1442695040888963407UL;
    re = ((double)((s >> 11) & 0xFFFFF) / (double)0x100000) - 0.5;
    s = s * 6364136223846793005UL + 1442695040888963407UL;
    im = ((double)((s >> 11) & 0xFFFFF) / (double)0x100000) - 0.5;
    vec[j] = re + im * I;
  }
}

/* Streamed evaluation: sum over source k0 (0-based) with kprime>=0 of
   conj(z[kprime]) * amp * z[k0], z 0-based == vec[k0+1]. */
static double complex streamed_gf(const TraceMap *map, const double complex *vec) {
  long int k0;
  double complex acc = 0.0;
  for (k0 = 0; k0 < map->n; k0++) {
    if (map->kprime[k0] >= 0)
      acc += conj(vec[map->kprime[k0] + 1]) * map->amp[k0] * vec[k0 + 1];
  }
  return acc;
}

/* ---- direct references (unmodified original element functions) ---- */

static void set_large_common(struct BindStruct *X, long int n) {
  long unsigned int irght, ilft, ihfbit;
  GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit);
  X->Large.i_max = n;
  X->Large.irght = irght;
  X->Large.ilft = ilft;
  X->Large.ihfbit = ihfbit;
  X->Large.mode = M_CORR;
}

static double complex direct_onebody(struct BindStruct *X, int ip,
                                     double complex *vec, long int n) {
  long unsigned int o1 = X->Def.CisAjt[ip][0] + 1, s1 = X->Def.CisAjt[ip][1];
  long unsigned int o2 = X->Def.CisAjt[ip][2] + 1, s2 = X->Def.CisAjt[ip][3];
  long int j;
  double complex dam = 0.0;
  set_large_common(X, n);
  if (X->Def.iCalcModel == HubbardGC) {
    /* inline of GC_general_hopp (avoids linking the MPI-heavy mltplyHubbard.c);
       trans == 1.0, so the extra "* 1.0" matches GC_general_hopp verbatim. */
    long unsigned int is1, is2, Asum, Adiff, off = 0;
    general_hopp_GetInfo(X, o1, o2, s1, s2);
    is1 = X->Large.is1_spin; is2 = X->Large.is2_spin;
    Asum = X->Large.isA_spin; Adiff = X->Large.A_spin;
    if (is1 == is2)
      for (j = 1; j <= n; j++) dam += GC_CisAis(j, vec, vec, X, is1, 1.0) * 1.0;
    else
      for (j = 1; j <= n; j++) dam += GC_CisAjt(j, vec, vec, X, is1, is2, Asum, Adiff, 1.0, &off) * 1.0;
  } else { /* SpinGC half */
    if (o1 == o2) {
      long unsigned int is1 = X->Def.Tpow[o1 - 1];
      if (s1 == s2) {
        for (j = 1; j <= n; j++)
          dam += child_SpinGC_CisAis(j, X, is1, s1) * conj(vec[j]) * vec[j];
      } else {
        long unsigned int off = 0;
        for (j = 1; j <= n; j++) {
          int sgn = child_SpinGC_CisAit(j, X, is1, s2, &off);
          if (sgn != 0) dam += sgn * conj(vec[off + 1]) * vec[j];
        }
      }
    }
  }
  return dam;
}

static double complex direct_twobody(struct BindStruct *X, int ip,
                                     double complex *vec, long int n) {
  long int j;
  double complex dam = 0.0;
  set_large_common(X, n);
  if (X->Def.iCalcModel == HubbardGC) {
    long unsigned int o1 = X->Def.CisAjtCkuAlvDC[ip][0] + 1, s1 = X->Def.CisAjtCkuAlvDC[ip][1];
    long unsigned int o2 = X->Def.CisAjtCkuAlvDC[ip][2] + 1, s2 = X->Def.CisAjtCkuAlvDC[ip][3];
    long unsigned int o3 = X->Def.CisAjtCkuAlvDC[ip][4] + 1, s3 = X->Def.CisAjtCkuAlvDC[ip][5];
    long unsigned int o4 = X->Def.CisAjtCkuAlvDC[ip][6] + 1, s4 = X->Def.CisAjtCkuAlvDC[ip][7];
    long unsigned int off = 0, is1, is2, is3, is4, Asum, Adiff, Bsum, Bdiff;
    general_int_GetInfo(ip, X, o1, o2, o3, o4, s1, s2, s3, s4, 1.0);
    is1 = X->Large.is1_spin; is2 = X->Large.is2_spin; Asum = X->Large.isA_spin; Adiff = X->Large.A_spin;
    is3 = X->Large.is3_spin; is4 = X->Large.is4_spin; Bsum = X->Large.isB_spin; Bdiff = X->Large.B_spin;
    if (is1 == is2 && is3 == is4)
      for (j = 1; j <= n; j++) dam += GC_CisAisCisAis_element(j, is1, is3, 1.0, vec, vec, X, &off);
    else if (is1 == is2 && is3 != is4)
      for (j = 1; j <= n; j++) dam += GC_CisAisCjtAku_element(j, is1, is3, is4, Bsum, Bdiff, 1.0, vec, vec, X, &off);
    else if (is1 != is2 && is3 == is4)
      for (j = 1; j <= n; j++) dam += GC_CisAjtCkuAku_element(j, is1, is2, is3, Asum, Adiff, 1.0, vec, vec, X, &off);
    else
      for (j = 1; j <= n; j++) dam += GC_CisAjtCkuAlv_element(j, is1, is2, is3, is4, Asum, Adiff, Bsum, Bdiff, 1.0, vec, vec, X, &off);
  } else { /* SpinGC half */
    long unsigned int o1, o2, o3, o4, s1, s2, s3, s4;
    double complex tmp_V;
    long unsigned int off = 0;
    if (Rearray_Interactions(ip, &o1, &o2, &o3, &o4, &s1, &s2, &s3, &s4, &tmp_V, X, 2) != 0)
      return 0.0;
    if (o1 == o2 && o3 == o4) {
      long unsigned int isA = X->Def.Tpow[o2 - 1], isB = X->Def.Tpow[o4 - 1];
      if (s1 == s2 && s3 == s4)
        for (j = 1; j <= n; j++) dam += GC_CisAisCisAis_spin_element(j, isA, isB, s2, s4, tmp_V, vec, vec, X);
      else if (s1 == s2 && s3 != s4)
        for (j = 1; j <= n; j++) dam += GC_CisAisCitAiu_spin_element(j, s2, s4, isA, isB, tmp_V, vec, vec, X, &off);
      else if (s1 != s2 && s3 == s4)
        for (j = 1; j <= n; j++) dam += GC_CisAitCiuAiu_spin_element(j, s2, s4, isA, isB, tmp_V, vec, vec, X, &off);
      else
        for (j = 1; j <= n; j++) dam += GC_CisAitCiuAiv_spin_element(j, s2, s4, isA, isB, tmp_V, vec, vec, X, &off);
    }
  }
  return dam;
}

/* ---- BindStruct minimal setup ---- */

static int **alloc_ops(int nrow, int ncol) {
  int **p = (int **)malloc(sizeof(int *) * (size_t)nrow);
  int r;
  for (r = 0; r < nrow; r++) p[r] = (int *)calloc((size_t)ncol, sizeof(int));
  return p;
}
static void set_ob(int **CisAjt, int i, int i1, int s1, int i2, int s2) {
  CisAjt[i][0] = i1; CisAjt[i][1] = s1; CisAjt[i][2] = i2; CisAjt[i][3] = s2;
}
static void set_tb(int **T, int i, int i1, int s1, int i2, int s2, int i3, int s3, int i4, int s4) {
  T[i][0] = i1; T[i][1] = s1; T[i][2] = i2; T[i][3] = s2;
  T[i][4] = i3; T[i][5] = s3; T[i][6] = i4; T[i][7] = s4;
}

static void init_bind(struct BindStruct *X, int model, int nsite, long int n) {
  int nb = 2 * nsite; /* enough Tpow entries for Hubbard (2 per site) */
  long unsigned int *Tpow;
  int m;
  memset(X, 0, sizeof(*X));
  X->Def.Nsite = nsite;
  X->Def.iCalcModel = model;
  X->Def.iFlgGeneralSpin = 0;
  X->Def.iFlgSzConserved = 0;
  X->Check.idim_max = (unsigned long)n;
  Tpow = (long unsigned int *)malloc(sizeof(long unsigned int) * (size_t)(nb + 1));
  for (m = 0; m <= nb; m++) Tpow[m] = 1UL << m;
  X->Def.Tpow = Tpow;
}

/* Fields the extraction drivers may touch: the full X->Large write-set
   (inventory 2c column vi). Compared field-by-field (memcmp is padding-unsafe). */
static int large_equal(const struct LargeList *a, const struct LargeList *b) {
  return a->i_max == b->i_max && a->mode == b->mode &&
         a->irght == b->irght && a->ilft == b->ilft && a->ihfbit == b->ihfbit &&
         a->is1_spin == b->is1_spin && a->is2_spin == b->is2_spin &&
         a->is3_spin == b->is3_spin && a->is4_spin == b->is4_spin &&
         a->isA_spin == b->isA_spin && a->isB_spin == b->isB_spin &&
         a->A_spin == b->A_spin && a->B_spin == b->B_spin &&
         a->is1_up == b->is1_up && a->is1_down == b->is1_down &&
         a->is2_up == b->is2_up && a->is2_down == b->is2_down &&
         a->tmp_V == b->tmp_V && a->tmp_J == b->tmp_J &&
         a->isite1 == b->isite1 && a->isite2 == b->isite2 &&
         a->isite3 == b->isite3 && a->isite4 == b->isite4;
}

static void run_onebody_case(struct BindStruct *X, int ip, double complex *vec,
                             long int n, const char *name) {
  TraceMap map;
  double complex direct, streamed;
  int rc = TraceMapExtractOneBody(X, ip, &map);
  expect_true(name, rc == 0);
  if (rc != 0) return;
  streamed = streamed_gf(&map, vec);
  direct = direct_onebody(X, ip, vec, n);
  expect_close(name, direct, streamed);
  TraceMapFree(&map);
}
static void run_twobody_case(struct BindStruct *X, int ip, double complex *vec,
                             long int n, const char *name) {
  TraceMap map;
  double complex direct, streamed;
  int rc = TraceMapExtractTwoBody(X, ip, &map);
  expect_true(name, rc == 0);
  if (rc != 0) return;
  streamed = streamed_gf(&map, vec);
  direct = direct_twobody(X, ip, vec, n);
  expect_close(name, direct, streamed);
  TraceMapFree(&map);
}

/* ---- purity: double extraction bit-identical + X->Large restored ---- */
static void run_purity(struct BindStruct *X, int twobody, int ip, long int n,
                       const char *name) {
  TraceMap m1, m2;
  struct LargeList before;
  long int k;
  int ok = 1;
  /* poison X->Large with a recognizable pattern to prove full restore */
  memset(&X->Large, 0x5A, sizeof(X->Large));
  before = X->Large;
  if (twobody) { TraceMapExtractTwoBody(X, ip, &m1); TraceMapExtractTwoBody(X, ip, &m2); }
  else         { TraceMapExtractOneBody(X, ip, &m1); TraceMapExtractOneBody(X, ip, &m2); }
  ok = ok && large_equal(&before, &X->Large);
  ok = ok && (m1.n == m2.n);
  if (m1.n == m2.n && m1.n > 0) {
    for (k = 0; k < m1.n; k++) {
      if (m1.kprime[k] != m2.kprime[k]) ok = 0;
      if (m1.amp[k] != m2.amp[k]) ok = 0;
    }
  }
  expect_true(name, ok);
  TraceMapFree(&m1);
  TraceMapFree(&m2);
}

/* ---- Task 3 Step 2: TraceStreamOneBody() gbuf contents vs direct
 * expec_cisajs_HubbardGC / expec_cisajs_SpinGCHalf execution (via the same
 * direct_onebody() ground-truth helper the Step-3 mapping-validity checks
 * above already use), for a random 3-state panel. 1e-13 tolerance. ---- */
static void run_stream_onebody_case(struct BindStruct *X, long int n,
                                    unsigned int seed_base, const char *label) {
  const long int jb = 1, je = 3, ncols = 3;
  long int nops = (long int)X->Def.NCisAjt;
  double complex *panel = (double complex *)malloc(sizeof(double complex) * (size_t)(ncols * n));
  double complex *cols[3];
  double complex *gbuf;
  long int c, p, k;
  int rc;
  char nm[192];

  for (c = 0; c < ncols; c++) {
    cols[c] = (double complex *)malloc(sizeof(double complex) * (size_t)(n + 1));
    fill_random(cols[c], n, seed_base + (unsigned int)c * 7919u);
    /* panel[(col)*NN + k] = cols[c][k+1]: 0-based slice of the 1-based
       reference vector, exactly the layout TraceStreamOneBody()'s doc
       comment specifies (== phys_distributed_local.c:70-81's convention). */
    for (k = 0; k < n; k++) panel[c * n + k] = cols[c][k + 1];
  }

  gbuf = (double complex *)malloc(sizeof(double complex) * (size_t)(nops * ncols));
  rc = TraceStreamOneBody(X, panel, jb, je, n, ncols, gbuf);
  snprintf(nm, sizeof(nm), "%s TraceStreamOneBody rc==0", label);
  expect_true(nm, rc == 0);

  if (rc == 0) {
    for (p = 0; p < nops; p++) {
      for (c = 0; c < ncols; c++) {
        double complex direct = direct_onebody(X, (int)p, cols[c], n);
        snprintf(nm, sizeof(nm), "%s gbuf[pair=%ld,state=%ld] vs direct", label, p, c);
        expect_close(nm, direct, gbuf[p * ncols + c]);
      }
    }
  }

  free(gbuf);
  for (c = 0; c < ncols; c++) free(cols[c]);
  free(panel);
}

/* ---- Task 4 Step 2: TraceStreamTwoBody() gbuf contents vs direct
 * expec_cisajscktalt_HubbardGC / expec_cisajscktalt_SpinGCHalf execution (via
 * the same direct_twobody() ground-truth helper the two-body mapping-
 * validity checks above already use), for a random 3-state panel. 1e-13
 * tolerance. Mirrors run_stream_onebody_case() exactly, swapping
 * NCisAjt/CisAjt/TraceStreamOneBody/direct_onebody for their two-body
 * counterparts. ---- */
static void run_stream_twobody_case(struct BindStruct *X, long int n,
                                    unsigned int seed_base, const char *label) {
  const long int jb = 1, je = 3, ncols = 3;
  long int nops = (long int)X->Def.NCisAjtCkuAlvDC;
  double complex *panel = (double complex *)malloc(sizeof(double complex) * (size_t)(ncols * n));
  double complex *cols[3];
  double complex *gbuf;
  long int c, p, k;
  int rc;
  char nm[192];

  for (c = 0; c < ncols; c++) {
    cols[c] = (double complex *)malloc(sizeof(double complex) * (size_t)(n + 1));
    fill_random(cols[c], n, seed_base + (unsigned int)c * 7919u);
    for (k = 0; k < n; k++) panel[c * n + k] = cols[c][k + 1];
  }

  gbuf = (double complex *)malloc(sizeof(double complex) * (size_t)(nops * ncols));
  rc = TraceStreamTwoBody(X, panel, jb, je, n, ncols, gbuf);
  snprintf(nm, sizeof(nm), "%s TraceStreamTwoBody rc==0", label);
  expect_true(nm, rc == 0);

  if (rc == 0) {
    for (p = 0; p < nops; p++) {
      for (c = 0; c < ncols; c++) {
        double complex direct = direct_twobody(X, (int)p, cols[c], n);
        snprintf(nm, sizeof(nm), "%s gbuf[pair=%ld,state=%ld] vs direct", label, p, c);
        expect_close(nm, direct, gbuf[p * ncols + c]);
      }
    }
  }

  free(gbuf);
  for (c = 0; c < ncols; c++) free(cols[c]);
  free(panel);
}

/* ---- Task 3 Step 2: memory-gate boundary. With the capability table
 * forced on via the HPHI_TRACE_FORCE dev hook (kTraceCap itself stays all
 * FALSE until plan Task 5; this hook is the documented way to exercise the
 * gate without it -- see TraceParseForceEnv()'s doc comment in
 * expec_trace.c), TraceBuildPlan() is called directly (not through
 * TraceGbufMaxBytesFromEnv()/the MPI orchestrator -- the plan doc says the
 * env var is read by phys_distributed.c, not needed here) with a byte cap
 * exactly at, and one byte under, nops*nc_uniform*sizeof(double complex).
 * At the cap: fits -> kernel[ONEBODY]==1, demoted_memory[ONEBODY]==0.
 * One byte short: demoted -> kernel[ONEBODY]==0, demoted_memory[ONEBODY]==1,
 * gbuf_bytes[ONEBODY]==0. ---- */
static void test_memory_gate_boundary(void) {
  struct BindStruct X;
  TraceExecutionPlan plan;
  long int nc_uniform = 4;
  size_t exact_bytes;
  int **ob = alloc_ops(1, 4);

  fprintf(stderr, "[memory-gate boundary]\n");
  memset(&X, 0, sizeof(X));
  X.Def.iCalcModel = HubbardGC;
  X.Def.iFlgGeneralSpin = 0;
  X.Def.iExpecMode = EXPECMODE_TRACE;
  X.Def.CisAjt = ob;
  X.Def.NCisAjt = 5;           /* nops */
  X.Def.NCisAjtCkuAlvDC = 0;   /* keep TWOBODY out of this boundary check */

  setenv("HPHI_TRACE_FORCE", "onebody", 1);

  exact_bytes = (size_t)X.Def.NCisAjt * (size_t)nc_uniform * sizeof(double complex);

  TraceBuildPlan(&X, nc_uniform, exact_bytes, &plan);
  expect_true("gate boundary: cap==exact -> kernel[ONEBODY]==1",
             plan.kernel[TRACE_Q_ONEBODY] == 1);
  expect_true("gate boundary: cap==exact -> demoted_memory[ONEBODY]==0",
             plan.demoted_memory[TRACE_Q_ONEBODY] == 0);
  expect_true("gate boundary: cap==exact -> gbuf_bytes[ONEBODY]==nops*nc_uniform*16",
             plan.gbuf_bytes[TRACE_Q_ONEBODY] == exact_bytes);

  TraceBuildPlan(&X, nc_uniform, exact_bytes - 1, &plan);
  expect_true("gate boundary: cap==exact-1 -> kernel[ONEBODY]==0",
             plan.kernel[TRACE_Q_ONEBODY] == 0);
  expect_true("gate boundary: cap==exact-1 -> demoted_memory[ONEBODY]==1",
             plan.demoted_memory[TRACE_Q_ONEBODY] == 1);
  expect_true("gate boundary: cap==exact-1 -> gbuf_bytes[ONEBODY]==0",
             plan.gbuf_bytes[TRACE_Q_ONEBODY] == 0);

  unsetenv("HPHI_TRACE_FORCE");
  free(ob[0]);
  free(ob);
}

/* ---- Task 4 Step 2: TWOBODY memory-gate boundary. Mirrors
 * test_memory_gate_boundary() exactly, but forces "twobody" and sizes the
 * cap off X.Def.NCisAjtCkuAlvDC (NCisAjt stays 0 so ONEBODY never enters the
 * plan, keeping this check isolated to TRACE_Q_TWOBODY -- the plan's per-
 * quantity gate is independent by construction: one quantity's cap can
 * demote without affecting the other, see TraceBuildPlan()'s per-q loop). ---- */
static void test_memory_gate_boundary_twobody(void) {
  struct BindStruct X;
  TraceExecutionPlan plan;
  long int nc_uniform = 4;
  size_t exact_bytes;
  int **tb = alloc_ops(1, 8);

  fprintf(stderr, "[memory-gate boundary: TWOBODY]\n");
  memset(&X, 0, sizeof(X));
  X.Def.iCalcModel = HubbardGC;
  X.Def.iFlgGeneralSpin = 0;
  X.Def.iExpecMode = EXPECMODE_TRACE;
  X.Def.NCisAjt = 0;              /* keep ONEBODY out of this boundary check */
  X.Def.CisAjtCkuAlvDC = tb;
  X.Def.NCisAjtCkuAlvDC = 7;       /* nops */

  setenv("HPHI_TRACE_FORCE", "twobody", 1);

  exact_bytes = (size_t)X.Def.NCisAjtCkuAlvDC * (size_t)nc_uniform * sizeof(double complex);

  TraceBuildPlan(&X, nc_uniform, exact_bytes, &plan);
  expect_true("gate boundary: cap==exact -> kernel[TWOBODY]==1",
             plan.kernel[TRACE_Q_TWOBODY] == 1);
  expect_true("gate boundary: cap==exact -> demoted_memory[TWOBODY]==0",
             plan.demoted_memory[TRACE_Q_TWOBODY] == 0);
  expect_true("gate boundary: cap==exact -> gbuf_bytes[TWOBODY]==nops*nc_uniform*16",
             plan.gbuf_bytes[TRACE_Q_TWOBODY] == exact_bytes);

  TraceBuildPlan(&X, nc_uniform, exact_bytes - 1, &plan);
  expect_true("gate boundary: cap==exact-1 -> kernel[TWOBODY]==0",
             plan.kernel[TRACE_Q_TWOBODY] == 0);
  expect_true("gate boundary: cap==exact-1 -> demoted_memory[TWOBODY]==1",
             plan.demoted_memory[TRACE_Q_TWOBODY] == 1);
  expect_true("gate boundary: cap==exact-1 -> gbuf_bytes[TWOBODY]==0",
             plan.gbuf_bytes[TRACE_Q_TWOBODY] == 0);

  unsetenv("HPHI_TRACE_FORCE");
  free(tb[0]);
  free(tb);
}

static void test_hubbardgc(void) {
  struct BindStruct X;
  long int n = 256;
  double complex *vec = (double complex *)malloc(sizeof(double complex) * (size_t)(n + 1));
  int **ob = alloc_ops(3, 4);
  int **tb = alloc_ops(6, 8);
  fprintf(stderr, "[HubbardGC L=4, n=256]\n");
  init_bind(&X, HubbardGC, 4, n);
  fill_random(vec, n, 0x1234u);

  /* one-body */
  X.Def.CisAjt = ob; X.Def.NCisAjt = 3;
  set_ob(ob, 0, 0, 0, 0, 0);   /* diagonal n_{0up} */
  set_ob(ob, 1, 0, 0, 1, 0);   /* off-diagonal c^+_{0up} c_{1up} */
  set_ob(ob, 2, 0, 0, 0, 1);   /* cross-spin same-site c^+_{0up} c_{0down} */
  run_onebody_case(&X, 0, vec, n, "HubbardGC 1B diagonal");
  run_onebody_case(&X, 1, vec, n, "HubbardGC 1B off-diagonal");
  run_onebody_case(&X, 2, vec, n, "HubbardGC 1B cross-spin");
  run_purity(&X, 0, 1, n, "HubbardGC 1B purity");
  run_stream_onebody_case(&X, n, 0x1B57u, "HubbardGC 1B stream");

  /* two-body: four element-family branches + same-index */
  X.Def.CisAjtCkuAlvDC = tb; X.Def.NCisAjtCkuAlvDC = 6;
  set_tb(tb, 0, 0,0, 0,0, 1,0, 1,0);   /* CisAisCisAis: n_{0up} n_{1up} */
  set_tb(tb, 1, 0,0, 0,0, 1,0, 2,0);   /* CisAisCjtAku */
  set_tb(tb, 2, 0,0, 1,0, 2,0, 2,0);   /* CisAjtCkuAku */
  set_tb(tb, 3, 0,0, 1,0, 2,0, 3,0);   /* CisAjtCkuAlv */
  set_tb(tb, 4, 0,0, 0,0, 0,0, 0,0);   /* same-index (all identical) */
  set_tb(tb, 5, 1,0, 0,0, 2,0, 2,0);   /* CisAjtCkuAku, empty create-at-occupied path variety */
  run_twobody_case(&X, 0, vec, n, "HubbardGC 2B CisAisCisAis");
  run_twobody_case(&X, 1, vec, n, "HubbardGC 2B CisAisCjtAku");
  run_twobody_case(&X, 2, vec, n, "HubbardGC 2B CisAjtCkuAku");
  run_twobody_case(&X, 3, vec, n, "HubbardGC 2B CisAjtCkuAlv");
  run_twobody_case(&X, 4, vec, n, "HubbardGC 2B same-index");
  run_twobody_case(&X, 5, vec, n, "HubbardGC 2B off-diag zero-check");
  run_purity(&X, 1, 3, n, "HubbardGC 2B purity");
  run_stream_twobody_case(&X, n, 0x2C68u, "HubbardGC 2B stream");

  free(vec);
}

static void test_spingchalf(void) {
  struct BindStruct X;
  long int n = 16;
  double complex *vec = (double complex *)malloc(sizeof(double complex) * (size_t)(n + 1));
  int **ob = alloc_ops(3, 4);
  int **tb = alloc_ops(6, 8);
  TraceMap map;
  fprintf(stderr, "[SpinGC-half L=4, n=16]\n");
  init_bind(&X, SpinGC, 4, n);
  fill_random(vec, n, 0x9abcu);

  /* one-body */
  X.Def.CisAjt = ob; X.Def.NCisAjt = 3;
  set_ob(ob, 0, 0, 0, 0, 0);   /* diagonal Sz-like (child_SpinGC_CisAis) */
  set_ob(ob, 1, 0, 0, 0, 1);   /* transverse S^+/S^- (child_SpinGC_CisAit) */
  set_ob(ob, 2, 0, 0, 1, 0);   /* different-site -> zero-result (empty map) */
  run_onebody_case(&X, 0, vec, n, "SpinGC 1B diagonal");
  run_onebody_case(&X, 1, vec, n, "SpinGC 1B transverse");
  run_onebody_case(&X, 2, vec, n, "SpinGC 1B zero-result");
  run_purity(&X, 0, 1, n, "SpinGC 1B purity");
  run_stream_onebody_case(&X, n, 0x59A6u, "SpinGC 1B stream");

  /* two-body: four spin-family branches */
  X.Def.CisAjtCkuAlvDC = tb; X.Def.NCisAjtCkuAlvDC = 6;
  set_tb(tb, 0, 0,0, 0,0, 1,0, 1,0);   /* SzSz diagonal (CisAisCisAis_spin) */
  set_tb(tb, 1, 0,0, 0,0, 1,0, 1,1);   /* CisAisCitAiu */
  set_tb(tb, 2, 0,0, 0,1, 1,0, 1,0);   /* CisAitCiuAiu */
  set_tb(tb, 3, 0,0, 0,1, 1,0, 1,1);   /* CisAitCiuAiv (S^+ S^-) */
  set_tb(tb, 4, 0,0, 0,0, 0,0, 0,0);   /* same-index (onsite, all site 0) */
  set_tb(tb, 5, 0,0, 1,0, 2,0, 3,0);   /* 4 distinct sites -> Rearray-irregular */
  run_twobody_case(&X, 0, vec, n, "SpinGC 2B CisAisCisAis_spin");
  run_twobody_case(&X, 1, vec, n, "SpinGC 2B CisAisCitAiu");
  run_twobody_case(&X, 2, vec, n, "SpinGC 2B CisAitCiuAiu");
  run_twobody_case(&X, 3, vec, n, "SpinGC 2B CisAitCiuAiv");
  run_twobody_case(&X, 4, vec, n, "SpinGC 2B same-index");
  run_purity(&X, 1, 0, n, "SpinGC 2B purity");

  /* boundary: Rearray-irregular pair -> n==0 sentinel; streamed==0 */
  {
    int rc = TraceMapExtractTwoBody(&X, 5, &map);
    expect_true("SpinGC 2B Rearray-irregular rc==0", rc == 0);
    expect_true("SpinGC 2B Rearray-irregular n==0 sentinel", map.n == 0);
    expect_close("SpinGC 2B irregular streamed==0", 0.0, streamed_gf(&map, vec));
    TraceMapFree(&map);
  }

  /* stream test: covers all 6 pairs, including the Rearray-irregular one
     (index 5) above -- TraceStreamTwoBody()'s inner k-loop over map.n==0
     naturally contributes 0.0, matching direct_twobody()'s explicit
     "Rearray fails -> return 0.0" branch. */
  run_stream_twobody_case(&X, n, 0x5DE1u, "SpinGC 2B stream");

  free(vec);
}

int main(void) {
  fprintf(stderr, "== expec_trace_map_check ==\n");
  test_hubbardgc();
  test_spingchalf();
  test_memory_gate_boundary();
  test_memory_gate_boundary_twobody();
  if (g_failures == 0) {
    fprintf(stderr, "ALL PASS\n");
    return 0;
  }
  fprintf(stderr, "FAILURES: %d\n", g_failures);
  return 1;
}
