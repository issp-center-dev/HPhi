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
 * @file expec_trace.h
 *
 * @brief ExpecMode 2 (trace-kernel) dispatch plan: the orchestration-boundary
 * API between src/phys_distributed.c (MPI orchestration) and the trace
 * kernels.
 *
 * This header is FROZEN as of phase 3b Task 1: later 3b tasks (mapping-probe
 * adapters, the one-body/two-body kernels, the golden capability-table
 * flip) do not add to it. Kernel-internal types/functions live in a separate
 * header, src/expec_trace_internal.h (Task 2) -- that split marks the
 * orchestration boundary explicitly, even though HPhi never installs either
 * header externally.
 *
 * Task 1 ships the dispatch skeleton only: the capability table
 * (src/expec_trace.c's kTraceCap) is all FALSE, so TraceBuildPlan() always
 * returns an all-fallback plan and expec_trace_owned_states() is a no-op.
 * ExpecMode 2 therefore behaves exactly like ExpecMode 1 until later tasks
 * flip capability rows to TRUE (only after the golden cross-checks in Task
 * 5 pass for every reachable operator-family branch).
 */
#pragma once
#include <stddef.h>
#include <stdio.h>
#include <complex.h>
#include "struct.h"

/** @brief Trace-kernel-eligible quantities (ExpecMode 2). */
typedef enum {
  TRACE_Q_ONEBODY = 0,   /* expec_cisajs equivalent */
  TRACE_Q_TWOBODY,       /* expec_cisajscktaltdc equivalent */
  TRACE_Q_NQUANT
} TraceQuantity;

/**
 * @brief The single source of truth for which quantities the ExpecMode 2
 * trace kernel handles for this run, and what it needed to decide that.
 *
 * Built once per run by TraceBuildPlan() and never mutated afterwards. The
 * INFO reporting (TraceReportPlan()), the kernel dispatcher
 * (expec_trace_owned_states()), and the ExpecMode-1 fallback loop
 * (phys_stateparallel_local_loop()) all consume this same instance -- no
 * other code path is permitted to make an independent kernel-vs-fallback
 * decision for a quantity.
 */
typedef struct {
  /** kernel[q]==1: the trace kernel is responsible for quantity q this run.
      kernel[q]==0: quantity q falls back to the ExpecMode-1 path. Immutable
      after TraceBuildPlan() returns. */
  int kernel[TRACE_Q_NQUANT];
  /** demoted_memory[q]==1: quantity q is statically capability-table
      eligible (or HPHI_TRACE_FORCE-forced) but was demoted to fallback by
      the runtime memory gate (TraceGbufBytes() returned 0). */
  int demoted_memory[TRACE_Q_NQUANT];
  /** When kernel[q]==1, the verified allocation size (bytes) for that
      quantity's result buffer, i.e. TraceGbufBytes()'s return value. Task
      3/4's malloc() must use ONLY this value -- re-reading the environment
      variable or recomputing the size expression anywhere else is
      structurally excluded by this field being the sole source. 0 when
      kernel[q]==0. */
  size_t gbuf_bytes[TRACE_Q_NQUANT];
  /** The uniform block width NC=ceil(neig/nproc) used to build this plan.
      This is NOT any particular rank's local owned-column count (ncols);
      zero-ownership is decided by the caller via je<jb, never via this
      field, so that the plan is bit-identical across every rank. */
  long int nc_uniform;
} TraceExecutionPlan;

/**
 * @brief Build the ExpecMode 2 execution plan: static capability table AND
 * runtime memory gate (checked size computation), evaluated once.
 *
 * If X->Def.iExpecMode != EXPECMODE_TRACE, returns the all-fallback plan
 * (every kernel[q]==0) unconditionally -- this is how ExpecMode 0/1 share
 * this same dispatch code path with byte-identical behavior to before phase
 * 3b.
 *
 * @param[in] X calculation parameters (model, general-spin flag, GF
 *               operator counts, ExpecMode)
 * @param[in] nc_uniform the uniform block width NC=ceil(neig/nproc) --
 *               identical on every rank (the caller must NOT pass a
 *               rank-local ncols; doing so would break the plan's
 *               cross-rank agreement)
 * @param[in] gbuf_max_bytes the per-quantity, per-rank result-buffer byte
 *               cap already resolved and MPI_Bcast-ed by the caller (see
 *               TraceGbufMaxBytesFromEnv())
 * @param[out] plan the plan to populate
 */
void TraceBuildPlan(const struct BindStruct *X, long int nc_uniform,
                    size_t gbuf_max_bytes, TraceExecutionPlan *plan);

/**
 * @brief Parse HPHI_TRACE_BUF_MAX_MB (getenv only; no MPI). Call on rank 0
 * only -- the caller (phys_distributed.c) MPI_Bcasts the result so every
 * rank agrees even if the environment differs node-to-node.
 *
 * @return the byte cap: HPHI_TRACE_BUF_MAX_MB (1..1048576, MiB) if valid,
 * else the default of 1024 MiB. Invalid values fall back to the default and
 * print one stderr warning.
 */
size_t TraceGbufMaxBytesFromEnv(void);

/**
 * @brief rank-0-only: print the plan, one INFO line per quantity (with the
 * demotion reason when applicable) plus one fixed line noting that the
 * energy/fluctuation, S2, NBodyG, and AnomalousG families always use the
 * ExpecMode-1 path in this version.
 */
void TraceReportPlan(const TraceExecutionPlan *plan, FILE *fp);

/**
 * @brief Evaluate the quantities this plan marks kernel[q]==1 for every
 * state this rank owns, and write them to the GF output.
 *
 * Completeness is guaranteed per-quantity, decided before any output is
 * written for that quantity: if mapping extraction or buffer allocation
 * fails for a kernel-owned quantity, that quantity is not written at all
 * (rc=-1, no partial output for it).
 *
 * @return 0 on success, -1 on a local (rank-local) failure.
 */
int expec_trace_owned_states(struct BindStruct *X, const TraceExecutionPlan *plan,
                             const double complex *panel,
                             long int jb, long int je, long int NN);
