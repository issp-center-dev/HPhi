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
#include "DefCommon.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

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
