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
/**
 * @file expec_trace.c
 *
 * @brief ExpecMode 2 (trace-kernel) plan construction, INFO reporting, and
 * kernel dispatch skeleton.
 *
 * This translation unit is MPI-free by construction (no mpi.h, no raw MPI_*
 * call, no exitMPI; the only permitted wrapperMPI helpers are
 * SumMPI_dc/d/li/i, fopenMPI, childfopenMPI, stdoutMPI) and compiles in
 * EVERY build, including ENABLE_MPI=OFF -- it is added to
 * test/check_expec_local_calls.sh's FILES in the same commit that creates
 * it. The single MPI touch this feature needs (broadcasting the
 * HPHI_TRACE_BUF_MAX_MB cap so every rank's plan agrees) lives in the MPI
 * orchestration layer, src/phys_distributed.c, which calls
 * TraceGbufMaxBytesFromEnv() on rank 0 only and passes the Bcast result in
 * here as a plain size_t.
 *
 * Task 1 (this commit) ships the dispatch skeleton: kTraceCap below is all
 * FALSE, so TraceBuildPlan() always returns an all-fallback plan for every
 * (model, quantity) combination and expec_trace_owned_states() is a no-op.
 * ExpecMode 2 is therefore observably identical to ExpecMode 1 except for
 * the new TraceReportPlan() INFO lines. Later 3b tasks add the mapping-probe
 * adapters and the real one-body/two-body kernels; only Task 5's golden
 * cross-checks are permitted to flip a kTraceCap row to TRUE.
 */
#include "expec_trace.h"
#include "expec_trace_internal.h"
#include "DefCommon.h"
#include "mltplyCommon.h"
#include "mltplyHubbardCore.h"
#include "mltplySpinCore.h"
#include "bitcalc.h"
#include "global.h"
#include "rearray_interactions.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

/* Rearray_Interactions() (src/rearray_interactions.c): the two-body extraction
   driver reuses the SAME +-1 tmp_V-folding reordering pass the Mode-1
   Spin/SpinGC-half two-body path uses. It performs no MPI. */

/* Static capability table: (calc_model, flg_general_spin) -> per-quantity
   trace-kernel readiness. Rows are flipped to 1 ONLY by plan Task 5, after
   the golden cross-checks for EVERY reachable operator-family branch of that
   (model, quantity) pass. See the branch-coverage table in
   docs/superpowers/specs/2026-07-11-expec-call-inventory.md §2c. */
typedef struct { int calc_model; int flg_general_spin; int q[TRACE_Q_NQUANT]; } TraceCap;
static const TraceCap kTraceCap[] = {
  { Hubbard,   0, {0, 0} },
  { HubbardGC, 0, {0, 0} },
  { Spin,      0, {0, 0} },   /* half only; iFlgGeneralSpin==1 is deliberately
                                 not matched by this row */
  { SpinGC,    0, {0, 0} },
};
#define N_TRACE_CAP ((int)(sizeof(kTraceCap) / sizeof(kTraceCap[0])))

/**
 * @brief Checked result-buffer size calculator, shared by TraceBuildPlan()'s
 * memory-gate decision and the Task 3/4 kernel malloc()s -- one function, so
 * the "does it fit" question can never be answered two different ways.
 *
 * @return 0 if the buffer does not fit / is not representable; otherwise the
 * exact byte count nops*nc_uniform*sizeof(double complex).
 */
static size_t TraceGbufBytes(long int nops, long int nc_uniform, size_t max_bytes) {
  size_t sn, sc;
  if (nops <= 0 || nc_uniform <= 0) return 0;
  /* Rule out truncation on platforms where long int is wider than size_t,
     via uintmax_t (do not assume a 64-bit long/size_t). */
  if ((uintmax_t)nops > (uintmax_t)SIZE_MAX ||
      (uintmax_t)nc_uniform > (uintmax_t)SIZE_MAX) return 0;
  sn = (size_t)nops; sc = (size_t)nc_uniform;
  if (sc > SIZE_MAX / sizeof(double complex)) return 0;      /* ncols*16 overflow */
  if (sn > SIZE_MAX / (sc * sizeof(double complex))) return 0; /* nops*(ncols*16) overflow */
  if (sn * sc * sizeof(double complex) > max_bytes) return 0; /* cap exceeded */
  return sn * sc * sizeof(double complex);
}

size_t TraceGbufMaxBytesFromEnv(void) {
  const size_t default_bytes = ((size_t)1024) << 20; /* 1024 MiB */
  const char *e = getenv("HPHI_TRACE_BUF_MAX_MB");
  long int val;
  char *endptr;

  if (e == NULL || e[0] == '\0') return default_bytes;

  errno = 0;
  val = strtol(e, &endptr, 10);
  if (errno != 0 || endptr == e || *endptr != '\0' || val < 1 || val > 1048576) {
    fprintf(stderr,
            "  Warning: HPHI_TRACE_BUF_MAX_MB='%s' is not an integer in "
            "[1,1048576]; using the default 1024 MiB.\n", e);
    return default_bytes;
  }
  /* MiB->byte conversion, checked (value_mb > SIZE_MAX>>20 falls back). */
  if ((uintmax_t)val > (uintmax_t)(SIZE_MAX >> 20)) {
    fprintf(stderr,
            "  Warning: HPHI_TRACE_BUF_MAX_MB=%ld overflows size_t after "
            "MiB->byte conversion; using the default 1024 MiB.\n", val);
    return default_bytes;
  }
  return ((size_t)val) << 20;
}

/**
 * @brief development-only hook: plan Task 5 REMOVES this
 *
 * Parses HPHI_TRACE_FORCE as a comma-separated list of quantity names
 * (onebody, twobody). Only the quantities named are forced to
 * kernel[q]=1-eligible regardless of kTraceCap (the memory gate in
 * TraceBuildPlan() still applies on top); every quantity NOT named is never
 * forced, so an unimplemented kernel can never be silently selected for a
 * quantity the caller didn't explicitly ask for. Unrecognized tokens print
 * one stderr warning each and are otherwise ignored.
 */
static void TraceParseForceEnv(int force_q[TRACE_Q_NQUANT]) {
  const char *e;
  char buf[256];
  char *tok;
  size_t len;

  force_q[TRACE_Q_ONEBODY] = 0;
  force_q[TRACE_Q_TWOBODY] = 0;

  e = getenv("HPHI_TRACE_FORCE");
  if (e == NULL || e[0] == '\0') return;

  len = strlen(e);
  if (len >= sizeof(buf)) len = sizeof(buf) - 1;
  memcpy(buf, e, len);
  buf[len] = '\0';

  tok = strtok(buf, ",");
  while (tok != NULL) {
    if (strcmp(tok, "onebody") == 0) {
      force_q[TRACE_Q_ONEBODY] = 1;
    } else if (strcmp(tok, "twobody") == 0) {
      force_q[TRACE_Q_TWOBODY] = 1;
    } else {
      fprintf(stderr,
              "  Warning: HPHI_TRACE_FORCE token '%s' is not a recognized "
              "quantity name (onebody, twobody); ignored.\n", tok);
    }
    tok = strtok(NULL, ",");
  }
}

void TraceBuildPlan(const struct BindStruct *X, long int nc_uniform,
                    size_t gbuf_max_bytes, TraceExecutionPlan *plan) {
  int q, i;
  int cap_q[TRACE_Q_NQUANT] = {0, 0};
  int force_q[TRACE_Q_NQUANT];

  memset(plan, 0, sizeof(*plan));
  plan->nc_uniform = nc_uniform;

  if (X == NULL || X->Def.iExpecMode != EXPECMODE_TRACE) {
    /* ExpecMode!=2: all-fallback plan (every field already zeroed above). */
    return;
  }

  for (i = 0; i < N_TRACE_CAP; i++) {
    if (kTraceCap[i].calc_model == X->Def.iCalcModel &&
        kTraceCap[i].flg_general_spin == X->Def.iFlgGeneralSpin) {
      cap_q[TRACE_Q_ONEBODY] = kTraceCap[i].q[TRACE_Q_ONEBODY];
      cap_q[TRACE_Q_TWOBODY] = kTraceCap[i].q[TRACE_Q_TWOBODY];
      break;
    }
  }

  TraceParseForceEnv(force_q); /* development-only hook: plan Task 5 REMOVES this */

  for (q = 0; q < TRACE_Q_NQUANT; q++) {
    long int nops;

    if (!(cap_q[q] || force_q[q])) {
      /* unsupported model (and not dev-forced): stays kernel=0,
         demoted_memory=0 -- TraceReportPlan() reports this as
         "unsupported model". */
      continue;
    }

    nops = (q == TRACE_Q_ONEBODY) ? (long int)X->Def.NCisAjt
                                   : (long int)X->Def.NCisAjtCkuAlvDC;

    plan->gbuf_bytes[q] = TraceGbufBytes(nops, nc_uniform, gbuf_max_bytes);
    if (plan->gbuf_bytes[q] == 0) {
      plan->demoted_memory[q] = 1;
      plan->kernel[q] = 0;
    } else {
      plan->kernel[q] = 1;
    }
  }
}

void TraceReportPlan(const TraceExecutionPlan *plan, FILE *fp) {
  static const char *kQuantityName[TRACE_Q_NQUANT] = { "one-body", "two-body" };
  int q;

  for (q = 0; q < TRACE_Q_NQUANT; q++) {
    if (plan->kernel[q]) {
      fprintf(fp, "  INFO: ExpecMode 2: %s Green functions use the trace kernel.\n",
              kQuantityName[q]);
    } else if (plan->demoted_memory[q]) {
      fprintf(fp,
              "  INFO: ExpecMode 2: %s Green functions use the ExpecMode-1 "
              "fallback (result buffer would exceed HPHI_TRACE_BUF_MAX_MB).\n",
              kQuantityName[q]);
    } else {
      fprintf(fp,
              "  INFO: ExpecMode 2: %s Green functions use the ExpecMode-1 "
              "fallback (unsupported model).\n",
              kQuantityName[q]);
    }
  }
  fprintf(fp,
          "  INFO: ExpecMode 2: energy/fluctuation, S2, NBodyG, and "
          "AnomalousG always use the ExpecMode-1 path in this version.\n");
}

int expec_trace_owned_states(struct BindStruct *X, const TraceExecutionPlan *plan,
                             const double complex *panel,
                             long int jb, long int je, long int NN) {
  /* Task 1: kTraceCap ships all FALSE, so in production plan->kernel[q] is
     0 for every quantity and this is a pure no-op -- every quantity is left
     for phys_stateparallel_local_loop()'s ExpecMode-1 fallback to handle, as
     it always has. (The HPHI_TRACE_FORCE dev hook can flip kernel[q] on
     before a real kernel exists, for early Task 3/4 checkpoint testing
     only; that is an explicit opt-in outside the production dispatch path,
     see TraceParseForceEnv().) */
  (void)X;
  (void)plan;
  (void)panel;
  (void)jb;
  (void)je;
  (void)NN;
  return 0;
}

/* ================================================================== */
/* Basis-mapping extraction (phase 3b Task 2).                        */
/*                                                                    */
/* These drivers replicate the LOCAL (intra-process) dispatch of      */
/* expec_cisajs.c / expec_cisajscktaltdc.c with X->Large.mode=M_CORR  */
/* and, per operator pair, stream the co-located `*_TraceProbe`        */
/* adapters over k=1..idim_max to fill a TraceMap. In the replicated  */
/* FullDiag layout the trace kernel runs in (iFlgScaLAPACK, no MPI    */
/* site decomposition) every operator site is intra-process, so the   */
/* inter-PE branches of the Mode-1 dispatch are unreachable and are   */
/* intentionally not mirrored here.                                    */
/*                                                                    */
/* Purity: X->Large is snapshot on entry and restored on exit, so the */
/* GetInfo scratch writes leave X byte-identical (see the purity      */
/* unit test). The only writes any reachable helper makes are to      */
/* X->Large; list_1/list_2_* and the Def/Tpow tables are read-only.   */
/* ================================================================== */

/** @brief Allocate/zero TraceMap buffers of length n (kprime=-1, amp=0). */
static int trace_map_alloc(TraceMap *map, long int n) {
  long int k;
  map->n = n;
  map->is_diagonal = 0;
  map->kprime = (long int *)malloc(sizeof(long int) * (size_t)n);
  map->amp = (double complex *)malloc(sizeof(double complex) * (size_t)n);
  if (map->kprime == NULL || map->amp == NULL) {
    TraceMapFree(map);
    return -1;
  }
  for (k = 0; k < n; k++) {
    map->kprime[k] = -1;
    map->amp[k] = 0.0;
  }
  return 0;
}

void TraceMapFree(TraceMap *map) {
  if (map == NULL) return;
  free(map->kprime);
  free(map->amp);
  map->kprime = NULL;
  map->amp = NULL;
  map->n = 0;
  map->is_diagonal = 0;
}

/** @brief Post-extraction range assertion: every kprime in [-1, n-1]. */
static void trace_map_assert_range(const TraceMap *map) {
  long int k;
  if (map->n <= 0) return;
  for (k = 0; k < map->n; k++) {
    assert(map->kprime[k] >= -1 && map->kprime[k] < map->n);
  }
}

int TraceMapExtractOneBody(struct BindStruct *X, int ipair, TraceMap *map) {
  long int n = (long int)X->Check.idim_max;
  struct LargeList saved = X->Large; /* purity snapshot */
  long unsigned int irght, ilft, ihfbit;
  long unsigned int org_isite1, org_isite2, org_sigma1, org_sigma2;
  long int k;
  int rc = 0;

  memset(map, 0, sizeof(*map));
  if (trace_map_alloc(map, n) != 0) return -1;

  if (GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit) != 0) {
    TraceMapFree(map);
    return -1;
  }
  X->Large.i_max = n;
  X->Large.irght = irght;
  X->Large.ilft = ilft;
  X->Large.ihfbit = ihfbit;
  X->Large.mode = M_CORR;

  org_isite1 = (long unsigned int)X->Def.CisAjt[ipair][0] + 1;
  org_sigma1 = (long unsigned int)X->Def.CisAjt[ipair][1];
  org_isite2 = (long unsigned int)X->Def.CisAjt[ipair][2] + 1;
  org_sigma2 = (long unsigned int)X->Def.CisAjt[ipair][3];

  switch (X->Def.iCalcModel) {
  case HubbardGC: {
    long unsigned int isite1, isite2, Asum, Adiff;
    general_hopp_GetInfo(X, org_isite1, org_isite2, org_sigma1, org_sigma2);
    isite1 = X->Large.is1_spin;
    isite2 = X->Large.is2_spin;
    Asum = X->Large.isA_spin;
    Adiff = X->Large.A_spin;
    if (isite1 == isite2) {
      map->is_diagonal = 1;
      for (k = 1; k <= n; k++)
        GC_CisAis_TraceProbe(k, X, isite1, &map->kprime[k - 1], &map->amp[k - 1]);
    } else {
      for (k = 1; k <= n; k++)
        GC_CisAjt_TraceProbe(k, X, isite1, isite2, Asum, Adiff,
                             &map->kprime[k - 1], &map->amp[k - 1]);
    }
    break;
  }
  case Hubbard:
  case tJ:
  case tJGC:
  case Kondo:
  case KondoGC: {
    /* Sz-conserved cross-spin one-body, and Kondo localized-vs-itinerant
       pairs, yield a 0.0 GF: leave the (already all -1) empty map. */
    if (X->Def.iFlgSzConserved == TRUE && org_sigma1 != org_sigma2) break;
    if (X->Def.iCalcModel == Kondo || X->Def.iCalcModel == KondoGC) {
      if ((X->Def.LocSpn[org_isite1 - 1] == 1 && X->Def.LocSpn[org_isite2 - 1] == 0) ||
          (X->Def.LocSpn[org_isite1 - 1] == 0 && X->Def.LocSpn[org_isite2 - 1] == 1))
        break;
    }
    general_hopp_GetInfo(X, org_isite1, org_isite2, org_sigma1, org_sigma2);
    if (org_isite1 == org_isite2 && org_sigma1 == org_sigma2) {
      long unsigned int is = X->Def.Tpow[2 * org_isite1 - 2 + org_sigma1];
      map->is_diagonal = 1;
      for (k = 1; k <= n; k++) {
        map->kprime[k - 1] = k - 1;
        map->amp[k - 1] = (double complex)((list_1[k] & is) / is);
      }
    } else {
      long unsigned int isite1 = X->Large.is1_spin, isite2 = X->Large.is2_spin;
      long unsigned int Asum = X->Large.isA_spin, Adiff = X->Large.A_spin;
      for (k = 1; k <= n; k++)
        CisAjt_TraceProbe(k, X, isite1, isite2, Asum, Adiff,
                          &map->kprime[k - 1], &map->amp[k - 1]);
    }
    break;
  }
  case Spin: {
    if (X->Def.iFlgGeneralSpin != FALSE) { rc = -1; break; } /* general spin: unsupported */
    if (org_sigma1 == org_sigma2 && org_isite1 == org_isite2) {
      long unsigned int isite1 = X->Def.Tpow[org_isite1 - 1];
      map->is_diagonal = 1;
      for (k = 1; k <= n; k++)
        child_Spin_CisAis_TraceProbe(k, X, isite1, org_sigma1,
                                     &map->kprime[k - 1], &map->amp[k - 1]);
    } /* else off-diagonal spin hopping -> empty map (GF = 0) */
    break;
  }
  case SpinGC: {
    if (X->Def.iFlgGeneralSpin != FALSE) { rc = -1; break; }
    if (org_isite1 == org_isite2) {
      long unsigned int isite1 = X->Def.Tpow[org_isite1 - 1];
      if (org_sigma1 == org_sigma2) {
        map->is_diagonal = 1;
        for (k = 1; k <= n; k++)
          child_SpinGC_CisAis_TraceProbe(k, X, isite1, org_sigma1,
                                         &map->kprime[k - 1], &map->amp[k - 1]);
      } else {
        for (k = 1; k <= n; k++)
          child_SpinGC_CisAit_TraceProbe(k, X, isite1, org_sigma2,
                                         &map->kprime[k - 1], &map->amp[k - 1]);
      }
    } /* else empty map (GF = 0) */
    break;
  }
  default:
    rc = -1;
    break;
  }

  if (rc == 0) trace_map_assert_range(map);
  X->Large = saved; /* purity restore */
  if (rc != 0) TraceMapFree(map);
  return rc;
}

int TraceMapExtractTwoBody(struct BindStruct *X, int ipair, TraceMap *map) {
  long int n = (long int)X->Check.idim_max;
  struct LargeList saved = X->Large; /* purity snapshot */
  long unsigned int irght, ilft, ihfbit;
  long int k;
  int rc = 0;
  double complex tmp_V = 1.0;
  long unsigned int oi1, oi2, oi3, oi4, os1, os2, os3, os4;

  memset(map, 0, sizeof(*map));

  if (GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit) != 0)
    return -1;
  X->Large.i_max = n;
  X->Large.irght = irght;
  X->Large.ilft = ilft;
  X->Large.ihfbit = ihfbit;
  X->Large.mode = M_CORR;

  switch (X->Def.iCalcModel) {
  case HubbardGC:
  case Hubbard:
  case tJ:
  case tJGC:
  case Kondo:
  case KondoGC: {
    long unsigned int isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff;
    int is_gc = (X->Def.iCalcModel == HubbardGC || X->Def.iCalcModel == tJGC ||
                 X->Def.iCalcModel == KondoGC);
    oi1 = (long unsigned int)X->Def.CisAjtCkuAlvDC[ipair][0] + 1;
    os1 = (long unsigned int)X->Def.CisAjtCkuAlvDC[ipair][1];
    oi2 = (long unsigned int)X->Def.CisAjtCkuAlvDC[ipair][2] + 1;
    os2 = (long unsigned int)X->Def.CisAjtCkuAlvDC[ipair][3];
    oi3 = (long unsigned int)X->Def.CisAjtCkuAlvDC[ipair][4] + 1;
    os3 = (long unsigned int)X->Def.CisAjtCkuAlvDC[ipair][5];
    oi4 = (long unsigned int)X->Def.CisAjtCkuAlvDC[ipair][6] + 1;
    os4 = (long unsigned int)X->Def.CisAjtCkuAlvDC[ipair][7];
    /* Canonical Hubbard applies the Sz-conserved 0.0-row shortcut. */
    if (!is_gc && X->Def.iFlgSzConserved == TRUE && (os1 + os3 != os2 + os4)) {
      map->n = 0; /* irregular/forbidden -> Mode-1 writes a 0.0 row */
      break;
    }
    if (trace_map_alloc(map, n) != 0) { rc = -1; break; }
    tmp_V = 1.0;
    general_int_GetInfo(ipair, X, oi1, oi2, oi3, oi4, os1, os2, os3, os4, tmp_V);
    isite1 = X->Large.is1_spin;
    isite2 = X->Large.is2_spin;
    Asum = X->Large.isA_spin;
    Adiff = X->Large.A_spin;
    isite3 = X->Large.is3_spin;
    isite4 = X->Large.is4_spin;
    Bsum = X->Large.isB_spin;
    Bdiff = X->Large.B_spin;
    if (isite1 == isite2 && isite3 == isite4) {
      map->is_diagonal = 1;
      for (k = 1; k <= n; k++)
        if (is_gc)
          GC_CisAisCisAis_element_TraceProbe(k, isite1, isite3, tmp_V, X, &map->kprime[k - 1], &map->amp[k - 1]);
        else
          CisAisCisAis_element_TraceProbe(k, isite1, isite3, tmp_V, X, &map->kprime[k - 1], &map->amp[k - 1]);
    } else if (isite1 == isite2 && isite3 != isite4) {
      for (k = 1; k <= n; k++)
        if (is_gc)
          GC_CisAisCjtAku_element_TraceProbe(k, isite1, isite3, isite4, Bsum, Bdiff, tmp_V, X, &map->kprime[k - 1], &map->amp[k - 1]);
        else
          CisAisCjtAku_element_TraceProbe(k, isite1, isite3, isite4, Bsum, Bdiff, tmp_V, X, &map->kprime[k - 1], &map->amp[k - 1]);
    } else if (isite1 != isite2 && isite3 == isite4) {
      for (k = 1; k <= n; k++)
        if (is_gc)
          GC_CisAjtCkuAku_element_TraceProbe(k, isite1, isite2, isite3, Asum, Adiff, tmp_V, X, &map->kprime[k - 1], &map->amp[k - 1]);
        else
          CisAjtCkuAku_element_TraceProbe(k, isite1, isite2, isite3, Asum, Adiff, tmp_V, X, &map->kprime[k - 1], &map->amp[k - 1]);
    } else {
      for (k = 1; k <= n; k++)
        if (is_gc)
          GC_CisAjtCkuAlv_element_TraceProbe(k, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff, tmp_V, X, &map->kprime[k - 1], &map->amp[k - 1]);
        else
          CisAjtCkuAlv_element_TraceProbe(k, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff, tmp_V, X, &map->kprime[k - 1], &map->amp[k - 1]);
    }
    break;
  }
  case Spin: {
    /* Canonical Spin-half two-body: Rearray folds the +-1 into tmp_V. */
    if (X->Def.iFlgGeneralSpin != FALSE) { rc = -1; break; }
    if (Rearray_Interactions(ipair, &oi1, &oi2, &oi3, &oi4, &os1, &os2, &os3, &os4, &tmp_V, X, 2) != 0) {
      map->n = 0; /* irregular pair -> Mode-1 writes a 0.0 row */
      break;
    }
    if (trace_map_alloc(map, n) != 0) { rc = -1; break; }
    {
      long unsigned int isA_up = X->Def.Tpow[oi1 - 1];
      long unsigned int isB_up = X->Def.Tpow[oi3 - 1];
      if (os1 == os2 && os3 == os4) { /* density-density diagonal */
        map->is_diagonal = 1;
        for (k = 1; k <= n; k++)
          CisAisCisAis_spin_element_TraceProbe(k, isA_up, isB_up, os2, os4, tmp_V, X, &map->kprime[k - 1], &map->amp[k - 1]);
      } else if (oi1 == oi3 && os1 == os4 && os3 == os2) { /* same-index reduction */
        map->is_diagonal = 1;
        for (k = 1; k <= n; k++) {
          long int kp;
          double complex a;
          child_Spin_CisAis_TraceProbe(k, X, isA_up, os1, &kp, &a);
          map->kprime[k - 1] = kp;
          map->amp[k - 1] = tmp_V * a;
        }
      } else if (os1 == os4 && os2 == os3) { /* exchange: amp = sign, NO tmp_V */
        for (k = 1; k <= n; k++)
          child_exchange_spin_element_TraceProbe(k, X, isA_up, isB_up, os2, os4, &map->kprime[k - 1], &map->amp[k - 1]);
      } /* else empty map (GF = 0) */
    }
    break;
  }
  case SpinGC: {
    if (X->Def.iFlgGeneralSpin != FALSE) { rc = -1; break; }
    if (Rearray_Interactions(ipair, &oi1, &oi2, &oi3, &oi4, &os1, &os2, &os3, &os4, &tmp_V, X, 2) != 0) {
      map->n = 0;
      break;
    }
    if (trace_map_alloc(map, n) != 0) { rc = -1; break; }
    if (oi1 == oi2 && oi3 == oi4) {
      long unsigned int isA_up = X->Def.Tpow[oi2 - 1];
      long unsigned int isB_up = X->Def.Tpow[oi4 - 1];
      if (os1 == os2 && os3 == os4) {
        map->is_diagonal = 1;
        for (k = 1; k <= n; k++)
          GC_CisAisCisAis_spin_element_TraceProbe(k, isA_up, isB_up, os2, os4, tmp_V, X, &map->kprime[k - 1], &map->amp[k - 1]);
      } else if (os1 == os2 && os3 != os4) {
        for (k = 1; k <= n; k++)
          GC_CisAisCitAiu_spin_element_TraceProbe(k, os2, os4, isA_up, isB_up, tmp_V, X, &map->kprime[k - 1], &map->amp[k - 1]);
      } else if (os1 != os2 && os3 == os4) {
        for (k = 1; k <= n; k++)
          GC_CisAitCiuAiu_spin_element_TraceProbe(k, os2, os4, isA_up, isB_up, tmp_V, X, &map->kprime[k - 1], &map->amp[k - 1]);
      } else {
        for (k = 1; k <= n; k++)
          GC_CisAitCiuAiv_spin_element_TraceProbe(k, os2, os4, isA_up, isB_up, tmp_V, X, &map->kprime[k - 1], &map->amp[k - 1]);
      }
    } /* else (non-onsite pairing) -> empty map (GF = 0), matching Mode 1 */
    break;
  }
  default:
    rc = -1;
    break;
  }

  if (rc == 0) trace_map_assert_range(map);
  X->Large = saved; /* purity restore */
  if (rc != 0) TraceMapFree(map);
  return rc;
}
