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
 * Task 1 shipped the dispatch skeleton with kTraceCap all FALSE and a
 * development-only HPHI_TRACE_FORCE hook so Tasks 3/4's clavius early
 * checkpoints could exercise the ONEBODY/TWOBODY kernels ahead of the
 * capability table being trustworthy. Task 5's golden cross-checks (unit
 * tests in test/unit/expec_trace_map_check.c plus the equiv np=2/3
 * checkpoints in test/fulldiag_expecmode_equiv.sh) have now passed for
 * Hubbard, HubbardGC, Spin (half), and SpinGC (half) -- see kTraceCap's
 * per-row evidence comments below -- so those four rows are TRUE and the
 * HPHI_TRACE_FORCE hook has been REMOVED (it was documented from its
 * introduction as development-only, not shipped): the capability table is
 * now the only enablement path for ExpecMode 2.
 *
 * Task 3 adds the ONEBODY kernel body: TraceStreamOneBody() (declared in
 * expec_trace_internal.h, so it's directly unit-testable) streams every
 * X->Def.CisAjt pair, operator-outer / state-inner, into a gbuf sized
 * EXACTLY plan->gbuf_bytes[TRACE_Q_ONEBODY] (never recomputed here -- see
 * TraceBuildPlan()'s doc comment), and expec_trace_owned_states() then
 * writes gbuf out state-major in a separate output phase, reusing
 * GreenOutputKindUsesAggregate/OpenAggregate/CloseAggregate/WriteIndexPrefix
 * and the GREEN_ONEBODY_ROW_FORMAT row format exactly as
 * src/expec_cisajs.c's FullDiag path does, so the two paths can never
 * silently drift into byte-different output for the same values.
 *
 * Task 4 adds the TWOBODY kernel body: TraceStreamTwoBody() and
 * expec_trace_twobody_output() mirror the ONEBODY pair exactly (same
 * operator-outer/state-inner streaming formula, same separate-buffer/
 * separate-output-phase structure), sized by plan->gbuf_bytes[TRACE_Q_TWOBODY]
 * ONLY -- ONEBODY and TWOBODY are demoted to the ExpecMode-1 fallback
 * independently, so one quantity's memory-gate outcome never affects the
 * other's. The two-body row format is NOT a single macro like ONEBODY's:
 * src/expec_cisajscktaltdc.c's PRE-EXISTING call sites already split into
 * GREEN_TWOBODY_ROW_FORMAT and GREEN_TWOBODY_ROW_FORMAT_SP (see
 * green_row_format.h's doc comment for the exact split), so
 * expec_trace_twobody_output() selects between them per pair using only
 * that pair's raw operator indices and X->Def.iFlgSzConserved -- see its own
 * doc comment for the exact selection rule.
 *
 * For both quantities, the streaming phase completes (or fails) entirely
 * before the output phase opens its first file, so a mapping/allocation
 * failure never produces partial output for that quantity. Cross-quantity
 * atomicity is NOT guaranteed: if ONEBODY's output phase succeeds but
 * TWOBODY fails (or vice versa were the dispatch order reversed), the
 * succeeding quantity's part file(s) remain on disk, but the collective
 * rc=-1 that failure produces means GreenOutputMergePartials() never
 * publishes ANY aggregate this run (the manifest is all-or-nothing) -- the
 * same recovery model (rerun) Mode 1's mid-loop failures already have.
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
#include "green_output.h"
#include "green_row_format.h"
#include "FileIO.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

/* green_output.h / FileIO.h (phase 3b Task 3): the ONEBODY output phase
   below reuses GreenOutputKindUsesAggregate/OpenAggregate/CloseAggregate/
   WriteIndexPrefix and childfopenMPI exactly as src/expec_cisajs.c does --
   these are all on the frozen ExpecLocal wrapperMPI allow-list
   (childfopenMPI ultimately calls fopenMPI; see the file header above) or
   are themselves MPI-free (green_output.c's manifest bookkeeping performs no
   MPI_* call outside GreenOutputMergePartials(), which this TU never
   calls). */

/* Rearray_Interactions() (src/rearray_interactions.c): the two-body extraction
   driver reuses the SAME +-1 tmp_V-folding reordering pass the Mode-1
   Spin/SpinGC-half two-body path uses. It performs no MPI. */

/* Static capability table: (calc_model, flg_general_spin) -> per-quantity
   trace-kernel readiness. Rows are flipped to 1 ONLY by plan Task 5, after
   the golden cross-checks for EVERY reachable operator-family branch of that
   (model, quantity) pass. See the branch-coverage table in
   docs/superpowers/specs/2026-07-11-expec-call-inventory.md §2c.

   Phase 3b Task 5 evidence (why each row below is TRUE, machine-readable so
   a future dispatch-branch change can be checked against what was actually
   verified):

   - Hubbard   {1,1}: §2c.1 rows "CisAjt @ mltplyHubbardCore.c" + "canonical
     diagonal (inline ..., expec_cisajs.c:482-487)" (one-body); §2c.2 "Hubbard
     / HubbardGC do not call Rearray_Interactions" table, all 4 branches
     (CisAisCisAis_element / CisAisCjtAku_element / CisAjtCkuAku_element /
     CisAjtCkuAlv_element) + the Sz-violating 0.0-shortcut row (§2c.2
     "Rearray-nonzero semantics" paragraph, mirrored for canonical Hubbard by
     the `iFlgSzConserved` shortcut at expec_cisajscktaltdc.c:728-734) (two-
     body). Canonical basis (GetOffComp/list_1) is not locally unit-testable
     (§2c.3), so verification is the clavius forced-kernel checkpoint
     np=2/3 2026-07-12 (Task 3 onebody-only, Task 4 onebody+twobody) PLUS
     the Task 5 golden case fulldiag_expecmode_equiv.sh's case1_hubbard_nbodyg
     (ExpecMode-2 sanity sub-case) and case5_hubbard_onebody_twobody_golden
     (dedicated branch-coverage golden: up/down hops, forward/reversed pairs,
     the boundary site pair, diagonal density both spins, and all 4 two-body
     element families + the Sz-violating zero-result row, all sites < Nsite).
   - HubbardGC {1,1}: §2c.1 rows "GC_CisAis"/"GC_CisAjt" (one-body); §2c.2 GC
     branches of the same 4-way table (two-body). GC bare-bit basis is
     locally unit-testable: test/unit/expec_trace_map_check.c's
     test_hubbardgc() -- "HubbardGC 1B diagonal", "HubbardGC 1B off-diagonal",
     "HubbardGC 1B cross-spin", "HubbardGC 1B purity", "HubbardGC 1B stream",
     "HubbardGC 2B CisAisCisAis", "HubbardGC 2B CisAisCjtAku",
     "HubbardGC 2B CisAjtCkuAku", "HubbardGC 2B CisAjtCkuAlv",
     "HubbardGC 2B same-index", "HubbardGC 2B off-diag zero-check",
     "HubbardGC 2B purity", "HubbardGC 2B stream", all PASS; plus
     fulldiag_expecmode_equiv.sh's case1_hubbard_nbodyg ExpecMode-2 sub-case
     (NBodyG exercises the always-fallback path in the same run) and the
     clavius forced-kernel checkpoint np=2/3 2026-07-12.
   - Spin      {1,1}: half only (iFlgGeneralSpin==1 is deliberately not
     matched by this row -- general spin stays out of scope, §2c.3
     "Exclusions"). §2c.1 row "child_Spin_CisAis" (one-body); §2c.2
     "Spin-half CANONICAL two-body" table, all 3 reachable families
     (density-density diagonal / same-index reduction / exchange) plus the
     Rearray-irregular 0.0 row (§2c.2 "Rearray-nonzero semantics"). Canonical
     basis is not locally unit-testable (§2c.3), so verification is the
     clavius forced-kernel checkpoint np=2/3 2026-07-12 (case3_spin_chain,
     canonical Spin L=8) plus fulldiag_expecmode_equiv.sh's
     case3_spin_chain ExpecMode-2 sub-case (Task 5).
   - SpinGC    {1,1}: half only, same general-spin exclusion as Spin. §2c.1
     rows "child_SpinGC_CisAis"/"child_SpinGC_CisAit" (one-body); §2c.2
     "SpinGC-half two-body" table, all 4 families (GC_CisAisCisAis_spin /
     GC_CisAisCitAiu / GC_CisAitCiuAiu / GC_CisAitCiuAiv). GC bare-bit basis
     is locally unit-testable: test/unit/expec_trace_map_check.c's
     test_spingchalf() -- "SpinGC 1B diagonal", "SpinGC 1B transverse",
     "SpinGC 1B zero-result", "SpinGC 1B purity", "SpinGC 1B stream",
     "SpinGC 2B CisAisCisAis_spin", "SpinGC 2B CisAisCitAiu",
     "SpinGC 2B CisAitCiuAiu", "SpinGC 2B CisAitCiuAiv",
     "SpinGC 2B same-index", "SpinGC 2B purity",
     "SpinGC 2B Rearray-irregular rc==0/n==0 sentinel/irregular streamed==0",
     "SpinGC 2B stream", all PASS; plus fulldiag_expecmode_equiv.sh's
     case2_spingc_gamma, case4_spingc_honeycomb_manybody (ThreeBodyG/
     FourBodyG/SixBodyG stay on the always-fallback path in the same run)
     ExpecMode-2 sub-cases, and the clavius forced-kernel checkpoint
     np=2/3 2026-07-12.

   NOT flipped (stay at the {0,0} fallback default, so no row is listed
   below): tJ/tJGC/Kondo/KondoGC -- known is_gc grouping mismatch in
   TraceMapExtractTwoBody/expec_trace_twobody_output's `is_gc` test (tJGC/
   KondoGC are grouped with HubbardGC's is_gc dispatch, which the §2c audit
   never verified for the tJ/Kondo-specific localized/itinerant exclusions);
   general spin (iFlgGeneralSpin==1) and Spinless remain out of scope per the
   plan's Global Constraints. */
typedef struct { int calc_model; int flg_general_spin; int q[TRACE_Q_NQUANT]; } TraceCap;
static const TraceCap kTraceCap[] = {
  { Hubbard,   0, {1, 1} },
  { HubbardGC, 0, {1, 1} },
  { Spin,      0, {1, 1} },   /* half only; iFlgGeneralSpin==1 is deliberately
                                 not matched by this row */
  { SpinGC,    0, {1, 1} },
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

void TraceBuildPlan(const struct BindStruct *X, long int nc_uniform,
                    size_t gbuf_max_bytes, TraceExecutionPlan *plan) {
  int q, i;
  int cap_q[TRACE_Q_NQUANT] = {0, 0};

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

  for (q = 0; q < TRACE_Q_NQUANT; q++) {
    long int nops;

    if (!cap_q[q]) {
      /* unsupported model: stays kernel=0, demoted_memory=0 --
         TraceReportPlan() reports this as "unsupported model". */
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

/**
 * @brief Task 3 output phase: write gbuf's ONEBODY results state-major,
 * exactly mirroring src/expec_cisajs.c's FullDiag branching and row format.
 *
 * For each owned state n (jb..je), sets X->Phys.eigen_num = n-1 (same
 * per-state convention the ExpecMode-1 fallback loop uses immediately after
 * this -- see phys_stateparallel_local_loop(), which re-sets it every
 * iteration regardless, so there is no cross-talk), opens exactly the file
 * expec_cisajs() would open for that state (GreenOutputOpenAggregate() in
 * aggregate/partial mode, or a direct childfopenMPI() of
 * cFileName1BGreen_FullDiag otherwise), writes one GREEN_ONEBODY_ROW_FORMAT
 * row per pair (in X->Def.CisAjt order, i.e. the same order gbuf was filled
 * in), then closes it. The open failure is checked and propagated (matching
 * expec_cisajs()); the aggregate close's return value is deliberately NOT
 * checked here either, matching expec_cisajs() verbatim -- a close failure
 * still surfaces because GreenOutputCloseAggregate() sets the session's
 * sticky closed_ok=0 for this kind, which GreenOutputMergePartials() (called
 * by the orchestrator after this whole ExpecLocal session ends) checks
 * before publishing anything.
 *
 * @return 0 on success, -1 on the first open failure (Mode-1 parity: no
 * retry, no fallback -- the sticky manifest/collective rc=-1 path takes over
 * from here, same as a Mode-1 write failure would).
 */
static int expec_trace_onebody_output(struct BindStruct *X, long int jb, long int je,
                                      long int ncols, const double complex *gbuf) {
  long int nops = (long int)X->Def.NCisAjt;
  long int n, p;

  for (n = jb; n <= je; n++) {
    FILE *fp = NULL;
    char sdt[D_FileNameMax];

    X->Phys.eigen_num = (int)(n - 1); /* 0-based, same convention as Mode 1 */

    if (GreenOutputKindUsesAggregate(X, GreenOutputOneBody)) {
      if (GreenOutputOpenAggregate(X, GreenOutputOneBody, &fp) != 0) return -1;
    } else {
      sprintf(sdt, cFileName1BGreen_FullDiag, X->Def.CDataFileHead, X->Phys.eigen_num);
      if (childfopenMPI(sdt, "w", &fp) != 0) return -1;
    }

    for (p = 0; p < nops; p++) {
      long unsigned int i1 = (long unsigned int)X->Def.CisAjt[p][0];
      long unsigned int s1 = (long unsigned int)X->Def.CisAjt[p][1];
      long unsigned int i2 = (long unsigned int)X->Def.CisAjt[p][2];
      long unsigned int s2 = (long unsigned int)X->Def.CisAjt[p][3];
      double complex val = gbuf[p * ncols + (n - jb)];

      GreenOutputWriteIndexPrefix(fp, X);
      fprintf(fp, GREEN_ONEBODY_ROW_FORMAT, i1, s1, i2, s2, creal(val), cimag(val));
    }

    if (GreenOutputKindUsesAggregate(X, GreenOutputOneBody)) {
      GreenOutputCloseAggregate(GreenOutputOneBody, fp); /* return value not
          checked -- Mode-1 parity, see the function doc comment above */
    } else {
      fclose(fp);
    }
  }
  return 0;
}

/**
 * @brief Task 4 output phase: write gbuf's TWOBODY results state-major,
 * exactly mirroring src/expec_cisajscktaltdc.c's FullDiag branching and row
 * formats for the models this quantity's capability table (kTraceCap) can
 * select: HubbardGC, and the canonical Hubbard-family group (Hubbard/tJ/
 * tJGC/Kondo/KondoGC, all of which Mode 1 routes to
 * expec_cisajscktalt_Hubbard()), and Spin/SpinGC (half only). Only Hubbard/
 * HubbardGC/Spin/SpinGC are TRUE as of Task 5 (tJ/tJGC/Kondo/KondoGC stay on
 * the fallback path -- see kTraceCap's "NOT flipped" note), so this format-
 * selection logic is written for the full Hubbard family even though only
 * a subset of it is currently reachable via the trace kernel.
 *
 * Row-format selection replicates Mode 1's per-model dispatch WITHOUT
 * needing any state carried over from the streaming phase, because it only
 * depends on the pair's raw operator indices (X->Def.CisAjtCkuAlvDC[p]) and
 * X->Def.iFlgSzConserved -- both already available here, same as the
 * one-body output phase reads X->Def.CisAjt[p] directly:
 *   - HubbardGC/tJGC/KondoGC (is_gc): GREEN_TWOBODY_ROW_FORMAT, always
 *     (Mode 1's GC path never writes the Sz-conserved 0.0-shortcut row).
 *   - Hubbard/tJ/Kondo (canonical, not is_gc): GREEN_TWOBODY_ROW_FORMAT_SP
 *     for a pair with iFlgSzConserved==TRUE && sigma1+sigma3!=sigma2+sigma4
 *     (Mode-1's Sz-conserved-violation 0.0 shortcut row --
 *     TraceMapExtractTwoBody()'s map.n==0 sentinel for exactly this same
 *     condition, see its doc comment), else GREEN_TWOBODY_ROW_FORMAT (the
 *     normally-computed row).
 *   - Spin/SpinGC (half): GREEN_TWOBODY_ROW_FORMAT_SP always -- Mode 1's
 *     expec_cisajscktalt_SpinHalf/SpinGCHalf use this format for BOTH their
 *     Rearray-irregular 0.0 row and their normally-computed row.
 *
 * @return 0 on success, -1 on the first open failure (Mode-1 parity, same
 * as expec_trace_onebody_output() above).
 */
static int expec_trace_twobody_output(struct BindStruct *X, long int jb, long int je,
                                      long int ncols, const double complex *gbuf) {
  long int nops = (long int)X->Def.NCisAjtCkuAlvDC;
  long int n, p;
  int model = X->Def.iCalcModel;
  int is_gc = (model == HubbardGC || model == tJGC || model == KondoGC);
  int is_hubbard_family = (model == Hubbard || model == HubbardGC || model == tJ ||
                           model == tJGC || model == Kondo || model == KondoGC);

  for (n = jb; n <= je; n++) {
    FILE *fp = NULL;
    char sdt[D_FileNameMax];

    X->Phys.eigen_num = (int)(n - 1); /* 0-based, same convention as Mode 1 */

    if (GreenOutputKindUsesAggregate(X, GreenOutputTwoBody)) {
      if (GreenOutputOpenAggregate(X, GreenOutputTwoBody, &fp) != 0) return -1;
    } else {
      sprintf(sdt, cFileName2BGreen_FullDiag, X->Def.CDataFileHead, X->Phys.eigen_num);
      if (childfopenMPI(sdt, "w", &fp) != 0) return -1;
    }

    for (p = 0; p < nops; p++) {
      long unsigned int i1 = (long unsigned int)X->Def.CisAjtCkuAlvDC[p][0];
      long unsigned int s1 = (long unsigned int)X->Def.CisAjtCkuAlvDC[p][1];
      long unsigned int i2 = (long unsigned int)X->Def.CisAjtCkuAlvDC[p][2];
      long unsigned int s2 = (long unsigned int)X->Def.CisAjtCkuAlvDC[p][3];
      long unsigned int i3 = (long unsigned int)X->Def.CisAjtCkuAlvDC[p][4];
      long unsigned int s3 = (long unsigned int)X->Def.CisAjtCkuAlvDC[p][5];
      long unsigned int i4 = (long unsigned int)X->Def.CisAjtCkuAlvDC[p][6];
      long unsigned int s4 = (long unsigned int)X->Def.CisAjtCkuAlvDC[p][7];
      double complex val = gbuf[p * ncols + (n - jb)];
      const char *fmt;

      if (is_hubbard_family) {
        fmt = (!is_gc && X->Def.iFlgSzConserved == TRUE && (s1 + s3 != s2 + s4))
                  ? GREEN_TWOBODY_ROW_FORMAT_SP : GREEN_TWOBODY_ROW_FORMAT;
      } else {
        fmt = GREEN_TWOBODY_ROW_FORMAT_SP;
      }

      GreenOutputWriteIndexPrefix(fp, X);
      fprintf(fp, fmt, i1, s1, i2, s2, i3, s3, i4, s4, creal(val), cimag(val));
    }

    if (GreenOutputKindUsesAggregate(X, GreenOutputTwoBody)) {
      GreenOutputCloseAggregate(GreenOutputTwoBody, fp); /* return value not
          checked -- Mode-1 parity, see expec_trace_onebody_output()'s doc
          comment above */
    } else {
      fclose(fp);
    }
  }
  return 0;
}

int expec_trace_owned_states(struct BindStruct *X, const TraceExecutionPlan *plan,
                             const double complex *panel,
                             long int jb, long int je, long int NN) {
  long int ncols = (je >= jb) ? (je - jb + 1) : 0;
  double complex *gbuf;
  int rc;

  /* Zero-owner rank: nothing to stream, nothing to write (Mode-1's loop
     over an empty [jb,je] range is likewise a no-op) -- return immediately,
     before even looking at plan->kernel[], so a zero-owner rank's manifest
     attempted count stays 0 rather than recording an empty attempt. */
  if (ncols <= 0) return 0;

  if (plan->kernel[TRACE_Q_ONEBODY]) {
    /* plan->gbuf_bytes[TRACE_Q_ONEBODY] is the ONLY size this malloc may
       use (see TraceExecutionPlan's doc comment in expec_trace.h) -- no
       re-reading HPHI_TRACE_BUF_MAX_MB, no re-deriving nops*ncols here.
       kernel[q]==1 implies TraceBuildPlan() already verified this is
       nonzero, but check anyway rather than trust that invariant blindly. */
    if (plan->gbuf_bytes[TRACE_Q_ONEBODY] == 0) return -1;

    gbuf = (double complex *)malloc(plan->gbuf_bytes[TRACE_Q_ONEBODY]);
    if (gbuf == NULL) return -1; /* before any write: no partial output */

    if (TraceStreamOneBody(X, panel, jb, je, NN, ncols, gbuf) != 0) {
      free(gbuf); /* mapping extraction failed -- nothing was written yet */
      return -1;
    }

    rc = expec_trace_onebody_output(X, jb, je, ncols, gbuf);
    free(gbuf);
    if (rc != 0) return rc;
  }

  if (plan->kernel[TRACE_Q_TWOBODY]) {
    /* plan->gbuf_bytes[TRACE_Q_TWOBODY] is the ONLY size this malloc may
       use, exactly mirroring the TRACE_Q_ONEBODY block above -- see
       TraceExecutionPlan's doc comment in expec_trace.h. This is a SEPARATE
       buffer and a SEPARATE gate decision from ONEBODY's: one quantity can
       be demoted to the ExpecMode-1 fallback while the other still uses the
       trace kernel this run (TraceBuildPlan() decides each plan->kernel[q]
       independently). */
    if (plan->gbuf_bytes[TRACE_Q_TWOBODY] == 0) return -1;

    gbuf = (double complex *)malloc(plan->gbuf_bytes[TRACE_Q_TWOBODY]);
    if (gbuf == NULL) return -1; /* before any write: no partial output */

    if (TraceStreamTwoBody(X, panel, jb, je, NN, ncols, gbuf) != 0) {
      free(gbuf); /* mapping extraction failed -- nothing was written yet */
      return -1;
    }

    rc = expec_trace_twobody_output(X, jb, je, ncols, gbuf);
    free(gbuf);
    if (rc != 0) return rc;
  }

  return 0;
}

int TraceStreamOneBody(struct BindStruct *X, const double complex *panel,
                       long int jb, long int je, long int NN,
                       long int ncols, double complex *gbuf) {
  long int nops = (long int)X->Def.NCisAjt;
  long int p, n, k;

  /* The panel's stride NN must equal the Hilbert-space dimension every
     TraceMap is built over -- see src/phys_distributed_local.c:70-81's
     `v0[j + 1] = panel[(n - jb) * NN + j]` (NN doubles as both the
     FullDiag state count and, here, the per-state vector length; the two
     coincide because neig == idim_max in the replicated FullDiag driver). */
  assert((long int)X->Check.idim_max == NN);

  for (p = 0; p < nops; p++) {
    TraceMap map;

    if (TraceMapExtractOneBody(X, (int)p, &map) != 0) return -1;

    for (n = jb; n <= je; n++) {
      const double complex *z = panel + (n - jb) * NN;
      double complex acc = 0.0;

      for (k = 0; k < map.n; k++) {
        if (map.kprime[k] >= 0)
          acc += conj(z[map.kprime[k]]) * map.amp[k] * z[k];
      }
      gbuf[p * ncols + (n - jb)] = acc;
    }

    TraceMapFree(&map); /* only one TraceMap alive at a time (spec Sec.3.1) */
  }
  return 0;
}

int TraceStreamTwoBody(struct BindStruct *X, const double complex *panel,
                       long int jb, long int je, long int NN,
                       long int ncols, double complex *gbuf) {
  long int nops = (long int)X->Def.NCisAjtCkuAlvDC;
  long int p, n, k;

  /* Same panel-stride invariant as TraceStreamOneBody() -- see that
     function's doc comment. */
  assert((long int)X->Check.idim_max == NN);

  for (p = 0; p < nops; p++) {
    TraceMap map;

    if (TraceMapExtractTwoBody(X, (int)p, &map) != 0) return -1;

    for (n = jb; n <= je; n++) {
      const double complex *z = panel + (n - jb) * NN;
      double complex acc = 0.0;

      for (k = 0; k < map.n; k++) {
        if (map.kprime[k] >= 0)
          acc += conj(z[map.kprime[k]]) * map.amp[k] * z[k];
      }
      gbuf[p * ncols + (n - jb)] = acc;
    }

    TraceMapFree(&map); /* only one TraceMap alive at a time (spec Sec.3.1) */
  }
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
