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
#include "expec_trace_internal.h"

/* (`myrank` and stdoutMPI are provided by the linked src/global.c.) */

/* --------------------------------------------------------------------------
 * Link-local verbatim copy of Rearray_Interactions (src/expec_cisajscktaltdc.c).
 * The production two-body driver calls the REAL one; linking that whole
 * translation unit would drag in the entire MPI two-body element library. This
 * copy is byte-for-byte the type==2 (two-body) reordering logic and MUST be
 * kept in sync with the source; the Task-5 clavius golden test exercises the
 * real function end-to-end.
 * -------------------------------------------------------------------------- */
int Rearray_Interactions(
    int i,
    long unsigned int *org_isite1, long unsigned int *org_isite2,
    long unsigned int *org_isite3, long unsigned int *org_isite4,
    long unsigned int *org_sigma1, long unsigned int *org_sigma2,
    long unsigned int *org_sigma3, long unsigned int *org_sigma4,
    double complex *tmp_V, struct BindStruct *X, int type) {
  long unsigned int t1, t2, t3, t4, u1, u2, u3, u4;
  (void)type; /* the test only uses type==2 (CisAjtCkuAlvDC) */
  t1 = (long unsigned int)X->Def.CisAjtCkuAlvDC[i][0] + 1;
  u1 = (long unsigned int)X->Def.CisAjtCkuAlvDC[i][1];
  t2 = (long unsigned int)X->Def.CisAjtCkuAlvDC[i][2] + 1;
  u2 = (long unsigned int)X->Def.CisAjtCkuAlvDC[i][3];
  t3 = (long unsigned int)X->Def.CisAjtCkuAlvDC[i][4] + 1;
  u3 = (long unsigned int)X->Def.CisAjtCkuAlvDC[i][5];
  t4 = (long unsigned int)X->Def.CisAjtCkuAlvDC[i][6] + 1;
  u4 = (long unsigned int)X->Def.CisAjtCkuAlvDC[i][7];
  if (t1 == t2 && t3 == t4) {
    if (t1 > t3) {
      *org_isite1 = t3; *org_sigma1 = u3; *org_isite2 = t4; *org_sigma2 = u4;
      *org_isite3 = t1; *org_sigma3 = u1; *org_isite4 = t2; *org_sigma4 = u2;
    } else {
      *org_isite1 = t1; *org_sigma1 = u1; *org_isite2 = t2; *org_sigma2 = u2;
      *org_isite3 = t3; *org_sigma3 = u3; *org_isite4 = t4; *org_sigma4 = u4;
    }
    *tmp_V = 1.0;
  } else if (t1 == t4 && t3 == t2) {
    if (t1 > t3) {
      *org_isite1 = t3; *org_sigma1 = u3; *org_isite2 = t2; *org_sigma2 = u2;
      *org_isite3 = t1; *org_sigma3 = u1; *org_isite4 = t4; *org_sigma4 = u4;
    } else {
      *org_isite1 = t1; *org_sigma1 = u1; *org_isite2 = t4; *org_sigma2 = u4;
      *org_isite3 = t3; *org_sigma3 = u3; *org_isite4 = t2; *org_sigma4 = u2;
    }
    *tmp_V = -1.0;
  } else {
    return -1;
  }
  return 0;
}

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

  free(vec);
}

int main(void) {
  fprintf(stderr, "== expec_trace_map_check ==\n");
  test_hubbardgc();
  test_spingchalf();
  if (g_failures == 0) {
    fprintf(stderr, "ALL PASS\n");
    return 0;
  }
  fprintf(stderr, "FAILURES: %d\n", g_failures);
  return 1;
}
