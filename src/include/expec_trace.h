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
 * Task 1 shipped the dispatch skeleton with the capability table
 * (src/expec_trace.c's kTraceCap) all FALSE, so TraceBuildPlan() returned an
 * all-fallback plan and expec_trace_owned_states() was a no-op. Task 5's
 * golden cross-checks (unit tests + equiv np=2/3 checkpoints) have since
 * passed for every reachable operator-family branch of Hubbard, HubbardGC,
 * Spin (half), and SpinGC (half), so those four rows are now TRUE for the
 * one-body and two-body quantities; tJ/tJGC/Kondo/KondoGC, general spin, and
 * Spinless remain on the ExpecMode-1 fallback path (see kTraceCap's
 * per-row comments for exactly which evidence backs each TRUE row).
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
  TRACE_Q_ENERGY,        /* expec_energy_flct equivalent (phase 3c) -- the CSR
                            energy/fluctuation family kernel. Unlike the two GF
                            slots (whose eligibility comes from the per-model
                            capability table kTraceCap), the energy slot is
                            eligible for EVERY model that reaches makeHam (see
                            TraceBuildPlan()), and its final kernel[]/demotion
                            state is decided by the two-phase build->finalize
                            lifecycle (TraceBuildPlan() then
                            TraceFinalizeEnergyPlan()), NOT by TraceBuildPlan()
                            alone -- the runtime CSR collector may demote it
                            rank-synchronously. */
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
      eligible but was demoted to fallback by the runtime memory gate
      (TraceGbufBytes() returned 0). */
  int demoted_memory[TRACE_Q_NQUANT];
  /** demoted_shared_evaluator[q]==1: quantity q is capability-table eligible
      (and would have passed the memory gate) but was demoted to fallback
      because its Mode-1 evaluator is SHARED with always-fallback quantities
      that must still run: expec_cisajscktaltdc() evaluates the two-body GF
      AND the ThreeBody/FourBody/SixBody GFs in one pass (it runs when ANY
      of NCisAjtCkuAlvDC/NTBody/NFBody/NSBody > 0 -- see
      expec_cisajscktaltdc.c:115), so when any multibody GF is defined the
      two-body GF must fall back together with them -- skipping the whole
      evaluator would silently drop the multibody outputs, and calling it
      anyway would double-write the two-body files (two writers, forbidden).
      Only TRACE_Q_TWOBODY can currently be demoted this way (expec_cisajs()
      handles one-body only -- it contains no NTBody/NFBody/NSBody
      references). NOTE: this field was added AFTER the Task-1 header freeze
      as a design-gap fix uncovered by the Task-5 golden tests (equiv case4:
      mode2 lost zvo_ThreeBody/FourBody/SixBody_eigen.dat) -- the freeze
      exception is deliberate and this comment is its record. */
  int demoted_shared_evaluator[TRACE_Q_NQUANT];
  /** no_operators[q]==1: quantity q is statically capability-table eligible
      (and is not demoted by the shared-evaluator rule) but this run defines
      ZERO operators of that kind (X->Def.NCisAjt==0 for TRACE_Q_ONEBODY,
      X->Def.NCisAjtCkuAlvDC==0 for TRACE_Q_TWOBODY), so there is nothing for
      the trace kernel to stream. Checked BEFORE the memory gate -- an empty
      operator table always makes TraceGbufBytes() return 0 too (nops<=0),
      so without this field the memory gate would set demoted_memory[q]=1
      and TraceReportPlan() would wrongly claim the result buffer exceeded
      HPHI_TRACE_BUF_MAX_MB for a quantity that has no result buffer to
      exceed anything with. Mutually exclusive with demoted_memory[q] and
      demoted_shared_evaluator[q] by construction (each quantity's outcome
      is decided by exactly one of these `if`/`else if` branches in
      TraceBuildPlan()'s per-q loop). NOTE: like demoted_shared_evaluator
      above, this field was added AFTER the Task-1 header freeze (final
      whole-branch review of phase 3b) -- the freeze exception is
      deliberate and this comment is its record. */
  int no_operators[TRACE_Q_NQUANT];
  /** demoted_input_ham[q]==1: quantity q was statically demoted to the
      ExpecMode-1 fallback because this run's Hamiltonian was read from
      InputHam (X->Def.iInputHam != 0). Re-running makeHam in the CSR collector
      would build a DIFFERENT matrix than the one that was diagonalized, so the
      energy family cannot use the trace kernel. Energy-only in practice
      (TRACE_Q_ENERGY): the GF slots never set this. Set by TraceBuildPlan()
      (static, pre-collect) and, being a build-time verdict, is exclusive with
      the finalize-time demoted_memory[TRACE_Q_ENERGY] -- TraceBuildPlan()
      leaves kernel[TRACE_Q_ENERGY]=0 when this is set, so
      TraceFinalizeEnergyPlan()'s collect/reduce is skipped entirely for an
      InputHam run. (Phase 3c Task 5; the Task-1 audit established that a
      FullDiag symmetry basis is unreachable, so NO unsupported-configuration
      predicate/field exists -- InputHam is the only static energy demotion.) */
  int demoted_input_ham[TRACE_Q_NQUANT];
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
 * @brief Phase 3c Task 5: finalize the PROVISIONAL energy slot with a single
 * rank-synchronized reduction, so every rank ends with the identical energy
 * verdict (build->finalize lifecycle, design spec 3c).
 *
 * After TraceBuildPlan() marks kernel[TRACE_Q_ENERGY]=1 provisionally, the
 * orchestrator runs the CSR collector (TraceHamCollect) on this rank, then
 * calls this to reconcile every rank's local outcome:
 *
 *   - @p local_ok is TraceHamCollect()'s return (1 success / 0 demotion or any
 *     allocation failure INCLUDING the pre-gate counts array).
 *   - @p nnz_raw is this rank's collected raw nnz (ham_csr.nnz) on success;
 *     ignored when local_ok==0 (pass 0).
 *
 * Protocol: a SUCCESS rank contributes {1, nnz, -nnz} (the long int -> long
 * long conversion is checked; if nnz_raw is negative or not representable this
 * rank treats itself as FAILED); a FAILED rank contributes
 * {0, LLONG_MAX, LLONG_MAX}. One MPI_Allreduce(MPI_MIN) over long long buf[3]
 * (a no-op that trivially passes at nproc==1 / non-MPI builds). Every rank then
 * applies the SAME verdict to the identical reduced buffer, IN THIS ORDER:
 *   (1) reduced[0]==0  -> at least one rank failed: demote energy on EVERY
 *       rank (kernel[TRACE_Q_ENERGY]=0, demoted_memory[TRACE_Q_ENERGY]=1);
 *   (2) else reduced[1] != -reduced[2] -> all ranks succeeded but their nnz
 *       disagree: the enumeration is nondeterministic across ranks, a
 *       CORRECTNESS error, not a fallback case -> fprintf(stderr,...) + a
 *       collective-safe exitMPI(-1) on ALL ranks.
 *   (3) else: energy stays kernel[TRACE_Q_ENERGY]=1 (final).
 *
 * From this call onward the plan is immutable: TraceReportPlan(), the kernel
 * dispatch, and phys_stateparallel_local_loop() all run afterward and observe
 * the same final plan. Skip this call entirely when TraceBuildPlan() already
 * left kernel[TRACE_Q_ENERGY]=0 (InputHam static demotion, or ExpecMode!=2).
 *
 * Implemented in the MPI-capable orchestration TU (src/phys_distributed.c),
 * NOT the MPI-free local loop.
 */
void TraceFinalizeEnergyPlan(TraceExecutionPlan *plan, int local_ok,
                             long int nnz_raw);

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
 * demotion reason when applicable) plus one fixed line noting that the S2,
 * NBodyG, and AnomalousG families always use the ExpecMode-1 path in this
 * version. As of phase 3c the energy/fluctuation family (TRACE_Q_ENERGY) is
 * itself a per-quantity INFO line -- trace kernel active, or the ExpecMode-1
 * fallback with its reason (HPHI_TRACE_BUF_MAX_MB memory gate, or InputHam) --
 * and is therefore DROPPED from the fixed always-fallback line. Call AFTER
 * TraceFinalizeEnergyPlan() so the energy line reflects the final verdict.
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
