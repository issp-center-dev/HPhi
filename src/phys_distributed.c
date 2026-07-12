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
 * @file phys_distributed.c
 *
 * @brief ExpecMode 1 (state-task-parallel) FullDiag observables driver:
 * MPI orchestration layer.
 *
 * This translation unit owns every genuine MPI collective of the
 * state-parallel path: the state-panel redistribution, the single
 * MIN-Allreduce failure rendezvous, the all_* Gatherv to rank 0, and the
 * collective GreenOutputMergePartials(). The MPI-free per-rank observable
 * loop lives in src/phys_distributed_local.c (which the ExpecLocal guard
 * scans); this file is permanently NOT scanned. The whole TU is compiled
 * only under _SCALAPACK (an empty TU otherwise, mirroring matrixscalapack.c),
 * because it references the distributed eigenvector globals Z_vec/descZ_vec
 * and the ScaLAPACK redistribution primitive.
 *
 * Phase 3b Task 1 also makes this the ExpecMode 2 (trace-kernel) dispatch
 * point: it builds the TraceExecutionPlan (identically on every rank -- see
 * src/expec_trace.c), reports the plan on rank 0, then runs
 * expec_trace_owned_states() before the ExpecMode-1 fallback loop so the
 * two never race on the same quantity. The HPHI_TRACE_BUF_MAX_MB cap is
 * parsed (rank 0 only) and Bcast ONLY when iExpecMode==EXPECMODE_TRACE, so
 * ExpecMode 0/1 runs never touch that environment variable (final
 * whole-branch review fix); TraceBuildPlan() never reads gbuf_max_bytes on
 * the ExpecMode!=TRACE path anyway (it returns the all-fallback plan before
 * that argument is used), so passing 0 there is safe.
 *
 * As of Task 5, the capability table (src/expec_trace.c's kTraceCap) gates
 * Hubbard, HubbardGC, half-integer Spin, and half-integer SpinGC -- every
 * other model, plus any of those four models' quantities that a runtime
 * check (shared-evaluator, no-operators, or the memory gate; see
 * TraceBuildPlan()) demotes, still falls back to the ExpecMode-1 path this
 * same call dispatches to right afterward. ExpecMode 0/1 always take the
 * all-fallback plan regardless (TraceBuildPlan() returns it unconditionally
 * whenever X->Def.iExpecMode != EXPECMODE_TRACE).
 *
 * Task 8 adds one more rank-0-only, post-Leave step: printing the trace
 * kernel's per-quantity phase-timing breakdown (map-extraction/streaming/
 * output) via TraceGetTimings() -- see that call site below for the exact
 * line format and src/include/expec_trace_internal.h for the accessor's
 * contract.
 */
#include "phys_distributed.h"
#ifdef _SCALAPACK
#include "matrixscalapack.h"   /* mpi.h, global.h, Z_vec/descZ_vec, use_scalapack,
                                  RedistBlockCyclicToStatePanel */
#include "green_output.h"
#include "wrapperMPI.h"
#include "DefCommon.h"
#include "expec_trace.h"
#include "expec_trace_internal.h"   /* TraceGetTimings() only -- see the file
                                       header above and that function's doc
                                       comment for why the orchestrator, not
                                       a kernel-internal caller, uses it */
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <stdint.h>

int phys_stateparallel(struct BindStruct *X, unsigned long int neig) {
  long int NN = (long int)neig;
  int P = nproc;
  long int NC = (NN + P - 1) / P;
  long int jb = (long int)myrank * NC + 1;
  long int je = ((long int)myrank + 1) * NC;
  long int ncols;
  int rc_local = 0;
  double complex *panel = NULL;
  int ok, gok, g;
  long int r;

  if (je > NN) je = NN;
  ncols = (je >= jb) ? (je - jb + 1) : 0;

  /* --- allocate the state panel; sync the malloc verdict with a MIN-Allreduce
     BEFORE any collective touches it (same pattern as lapack_diag_elpa). --- */
  panel = malloc((((NN * ncols) > 0) ? (size_t)(NN * ncols) : 1)
                 * sizeof(double complex));
  ok = (panel != NULL) ? 0 : -1;
  if (ok != 0)
    fprintf(stdout, "  Error: malloc failed for state panel (rank %d).\n", myrank);
  MPI_Allreduce(&ok, &gok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (gok != 0) {
    free(panel);
    return -1;
  }

  /* Redistribute the 2D block-cyclic eigenvectors into this rank's state
     column panel, then release Z_vec (peak memory drops back to O(N^2/P)). */
  if (RedistBlockCyclicToStatePanel(NN, Z_vec, descZ_vec, jb, ncols, NN, panel)
      != 0) {
    free(panel);
    return -1;
  }
  free(Z_vec);
  Z_vec = NULL;

  /* --- ExpecMode 2 plan: build once, identically on every rank (the
     capability table gates Hubbard/HubbardGC/half-Spin/half-SpinGC as of
     Task 5; every other model, and any runtime-demoted quantity of those
     four, gets the all-fallback plan -- see src/expec_trace.c). ---- */
  {
    TraceExecutionPlan tplan;
    /* overflow-free ceiling division (NN>=0, nproc>0 are the caller's
       precondition) -- distinct from NC above, which the plan API commits
       to calling nc_uniform to make explicit it is NOT a rank-local ncols. */
    long int nc_uniform = NN / (long int)nproc + ((NN % (long int)nproc) != 0);
    uint64_t gbuf_max = 0;
    /* Final whole-branch review fix: only ExpecMode 2 runs parse (and can
       warn about) HPHI_TRACE_BUF_MAX_MB. ExpecMode 0/1 leave gbuf_max at 0
       and skip the Bcast too -- TraceBuildPlan() never reads this argument
       on the ExpecMode!=TRACE path (see its own early return), so 0 is
       never divided by or otherwise used there. */
    if (X->Def.iExpecMode == EXPECMODE_TRACE) {
      if (myrank == 0) gbuf_max = (uint64_t)TraceGbufMaxBytesFromEnv(); /* rank 0 only */
      MPI_Bcast(&gbuf_max, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    }
    /* checked narrowing in case size_t is narrower than 64 bits here */
    TraceBuildPlan(X, nc_uniform,
                   (gbuf_max > (uint64_t)SIZE_MAX) ? (size_t)SIZE_MAX : (size_t)gbuf_max,
                   &tplan); /* ExpecMode!=2 => all-fallback plan */
    if (X->Def.iExpecMode == EXPECMODE_TRACE && myrank == 0)
      TraceReportPlan(&tplan, stdoutMPI);

    /* --- MPI-free per-rank observable session (partial green_output
       session). Single-exit invariant: once ExpecLocalEnter() has run,
       GreenOutputClearPartialSuffix() and ExpecLocalLeave() below MUST run
       exactly once before this block ends, regardless of whether the
       kernel dispatch or the fallback loop reports failure -- do NOT add an
       early return between Enter() and Leave(). --- */
    ExpecLocalEnter();                 /* also clears the sticky ExpecLocal error */
    GreenOutputSetPartialSuffix(myrank);
    rc_local = expec_trace_owned_states(X, &tplan, panel, jb, je, NN);
    if (rc_local == 0)
      rc_local = phys_stateparallel_local_loop(X, panel, jb, je, NN, &tplan);
    GreenOutputClearPartialSuffix();
    ExpecLocalLeave();

    /* Task 8 benchmark-breakdown line (final whole-branch review fix):
       rank 0 only, after Leave (so ExpecLocal's guard window is already
       closed -- this print is not itself scanned/restricted the way the
       per-rank observable session above is), and only when this rank's
       local trace-kernel dispatch actually completed
       (rc_local==0 -- a failed call may have left the accumulators
       mid-phase). One line per quantity this run selected the trace kernel
       for; a quantity that fell back is not printed (nothing was timed for
       it -- see TraceGetTimings()'s doc comment). */
    if (X->Def.iExpecMode == EXPECMODE_TRACE && myrank == 0 && rc_local == 0) {
      double timings[TRACE_Q_NQUANT][3];
      TraceGetTimings(timings);
      if (tplan.kernel[TRACE_Q_ONEBODY])
        fprintf(stdoutMPI,
                "  ExpecMode 2 timing (rank 0): one-body map=%.3fs stream=%.3fs output=%.3fs\n",
                timings[TRACE_Q_ONEBODY][0], timings[TRACE_Q_ONEBODY][1],
                timings[TRACE_Q_ONEBODY][2]);
      if (tplan.kernel[TRACE_Q_TWOBODY])
        fprintf(stdoutMPI,
                "  ExpecMode 2 timing (rank 0): two-body map=%.3fs stream=%.3fs output=%.3fs\n",
                timings[TRACE_Q_TWOBODY][0], timings[TRACE_Q_TWOBODY][1],
                timings[TRACE_Q_TWOBODY][2]);
    }
  }
  free(panel);

  /* --- single rendezvous: share the failure verdict across all ranks. --- */
  MPI_Allreduce(&rc_local, &g, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (g != 0)
    return -1;

  /* --- Gatherv the all_* summary arrays to rank 0. Each rank owns the
     block-contiguous segment [jb-1, je-1]; recvcounts/displs are derived from
     the deterministic block ownership. All all_* arrays are double, so the
     same recvcounts/displs are reused for every array. --- */
  {
    int *recvcounts = NULL, *displs = NULL;
    double *tmp = NULL;
    int sendcount;

    /* long->int width check: the largest int-valued quantity fed to
       MPI_Gatherv is a displ/count bounded by NN. */
    ok = 0;
    if (NN > (long int)INT_MAX) {
      ok = -1;
      if (myrank == 0)
        fprintf(stdout,
                "  Error: N=%ld exceeds INT_MAX; cannot Gatherv all_* arrays.\n",
                NN);
    }
    /* rank-0 receive buffers; sync malloc success before the collective. */
    if (ok == 0 && myrank == 0) {
      recvcounts = malloc((size_t)P * sizeof(int));
      displs = malloc((size_t)P * sizeof(int));
      tmp = malloc((size_t)NN * sizeof(double));
      if (recvcounts == NULL || displs == NULL || tmp == NULL) {
        ok = -1;
        fprintf(stdout,
                "  Error: malloc failed for all_* Gatherv buffers (rank 0).\n");
      }
    }
    MPI_Allreduce(&ok, &gok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (gok != 0) {
      free(recvcounts);
      free(displs);
      free(tmp);
      return -1;
    }

    if (myrank == 0) {
      for (r = 0; r < P; r++) {
        long int jb_r = r * NC + 1;
        long int je_r = (r + 1) * NC;
        long int nc_r;
        if (je_r > NN) je_r = NN;
        nc_r = (je_r >= jb_r) ? (je_r - jb_r + 1) : 0;
        recvcounts[r] = (int)nc_r;
        displs[r] = (int)(r * NC);  /* 0-based start; unused when nc_r==0 */
      }
    }
    sendcount = (int)ncols;

    /* rank 0 receives into a separate temp buffer, then copies into
       X->Phys.all_* (avoids MPI_IN_PLACE aliasing with its own segment).
       Senders pass &all_xxx[jb-1] with sendcount=ncols (0 for zero-owner). */
#define PHYS_GATHER_ONE(arr)                                                   \
    MPI_Gatherv((ncols > 0) ? &X->Phys.arr[jb - 1] : X->Phys.arr, sendcount,   \
                MPI_DOUBLE, tmp, recvcounts, displs, MPI_DOUBLE, 0,            \
                MPI_COMM_WORLD);                                               \
    if (myrank == 0) {                                                         \
      long int k;                                                             \
      for (k = 0; k < NN; k++) X->Phys.arr[k] = tmp[k];                        \
    }
    PHYS_GATHER_ONE(all_energy)
    PHYS_GATHER_ONE(all_doublon)
    PHYS_GATHER_ONE(all_sz)
    PHYS_GATHER_ONE(all_s2)
    PHYS_GATHER_ONE(all_num_up)
    PHYS_GATHER_ONE(all_num_down)
#undef PHYS_GATHER_ONE

    free(recvcounts);
    free(displs);
    free(tmp);
  }

  /* --- rank 0 re-renders the progress lines in state order, in the serial
     format (S2 column included). Mode 1 is only ever entered for FullDiag. --- */
  if (myrank == 0 && X->Def.iCalcType == FullDiag) {
    long int i;
    double tmp_N;
    for (i = 0; i < NN; i++) {
      if (X->Def.iCalcModel == Spin || X->Def.iCalcModel == SpinGC) {
        tmp_N = X->Def.NsiteMPI;
      } else {
        tmp_N = X->Phys.all_num_up[i] + X->Phys.all_num_down[i];
      }
      fprintf(stdoutMPI,
              "i=%5ld Energy=%10lf N=%10lf Sz=%10lf S2=%10lf Doublon=%10lf \n",
              i, X->Phys.all_energy[i], tmp_N, X->Phys.all_sz[i],
              X->Phys.all_s2[i], X->Phys.all_doublon[i]);
    }
  }

  /* --- publish the aggregate green_output files: rank 0 concatenates the
     rank-ordered (= state-ordered) partials. Collective; identical rc on all
     ranks; publishes/deletes only on global success. --- */
  if (GreenOutputMergePartials(X) != 0)
    return -1;

  return 0;
}

#endif /* _SCALAPACK */
