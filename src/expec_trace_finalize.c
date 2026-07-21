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
 * @file expec_trace_finalize.c
 *
 * @brief Phase 3c Task 5: TraceFinalizeEnergyPlan() -- the ONE rank-
 * synchronized MPI step of the ExpecMode 2 energy build->finalize lifecycle,
 * carved into its own MPI-capable translation unit.
 *
 * Why a separate TU (not expec_trace.c, not phys_distributed.c):
 *
 *  - It CANNOT live in src/expec_trace.c: that TU is MPI-FREE by construction
 *    (see its file header) and is in test/check_expec_local_calls.sh's FILES
 *    list, so a raw MPI_Allreduce()/exitMPI() there would break both the
 *    documented invariant and the CI guard.
 *
 *  - It should NOT force the serial unit test to link src/phys_distributed.c:
 *    under -D_SCALAPACK that TU also compiles the phys_stateparallel()
 *    orchestrator, which references RedistBlockCyclicToStatePanel()
 *    (matrixscalapack.c) and the phys_distributed_local.c symbols -- none of
 *    which test/unit/expec_trace_ham_check.c links, so it fails to link on any
 *    _SCALAPACK/ELPA host (the defect this file fixes).
 *
 * So the finalize -- the only plan step that must issue a genuine MPI
 * collective -- gets its own tiny TU with NO ScaLAPACK dependency. It is linked
 * into the HPhi build (src/CMakeLists.txt) AND directly into the unit test
 * (test/CMakeLists.txt), which drives its single-rank pass-through serially and
 * the two-rank demotion under mpiexec. phys_distributed.c keeps calling it via
 * the declaration in expec_trace.h. This TU is NOT scanned by the ExpecLocal
 * guard (it is genuine MPI orchestration code, like phys_distributed.c).
 *
 * The full protocol contract is documented on the declaration in
 * src/include/expec_trace.h.
 */
#include "expec_trace.h"
#include "wrapperMPI.h"   /* exitMPI (collective-safe abort) */
#ifdef MPI
#include <mpi.h>
#endif
#include <limits.h>
#include <stdio.h>

void TraceFinalizeEnergyPlan(TraceExecutionPlan *plan, int local_ok,
                             long int nnz_merged) {
  long long buf[3];

  if (local_ok && nnz_merged >= 0 && (long long)nnz_merged != LLONG_MAX) {
    /* success rank with a representable merged nnz distinct from the failure
       sentinel (long int always fits long long; nnz_merged>=0 makes the
       negation safe). */
    buf[0] = 1;
    buf[1] = (long long)nnz_merged;
    buf[2] = -(long long)nnz_merged;
  } else {
    /* failed rank (gate exceeded, any allocation failure), or an nnz that
       cannot be encoded distinctly from the sentinel -> treat self as failed. */
    buf[0] = 0;
    buf[1] = LLONG_MAX;
    buf[2] = LLONG_MAX;
  }

#ifdef MPI
  MPI_Allreduce(MPI_IN_PLACE, buf, 3, MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
#endif

  /* Verdict order (identical reduced buffer on every rank -> identical branch),
     design spec 3c: demote-before-consistency so failed-rank sentinels never
     pollute the nnz agreement check. */
  if (buf[0] == 0) {
    plan->kernel[TRACE_Q_ENERGY] = 0;
    plan->demoted_memory[TRACE_Q_ENERGY] = 1;
  } else if (buf[1] != -buf[2]) {
    fprintf(stderr,
            "  Error: ExpecMode 2 energy trace kernel: the collected "
            "Hamiltonian nnz differs across MPI ranks (min nnz=%lld, "
            "max nnz=%lld) -- a nondeterministic makeHam enumeration, which is "
            "a correctness error, not a fallback case. Aborting.\n",
            buf[1], -buf[2]);
    exitMPI(-1);
  }
  /* else: every rank succeeded with an agreeing nnz -> energy stays
     kernel[TRACE_Q_ENERGY]=1 (final). */
}
