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
 * @file phys_distributed_local.c
 *
 * @brief ExpecMode 1 (state-task-parallel) MPI-free local observable loop.
 *
 * This translation unit is the "local loop" layer of the state-parallel
 * FullDiag observables driver. It must contain ZERO raw MPI_* calls, no
 * exitMPI, and no wrapperMPI helper outside the frozen ExpecLocal allow-list
 * (SumMPI_dc/d/li/i, fopenMPI, childfopenMPI, stdoutMPI) -- it does not even
 * include mpi.h. test/check_expec_local_calls.sh scans this file (it is in
 * that guard's FILES list) so the invariant is CI-enforced. All genuine MPI
 * collectives (redistribution, the MIN-Allreduce rendezvous, the all_*
 * Gatherv, and GreenOutputMergePartials) live in the companion orchestration
 * TU src/phys_distributed.c, which is permanently NOT scanned.
 *
 * Because it references only build-independent globals (v0/v1, X->Phys.*) and
 * always-available APIs, this file compiles fully in every build, including
 * the non-MPI / non-_SCALAPACK default build; it is simply never called there
 * (its only caller, phys_stateparallel(), is compiled out without _SCALAPACK).
 */
#include "phys_distributed.h"
#include "expec_energy_flct.h"
#include "expec_totalspin.h"
#include "expec_cisajs.h"
#include "expec_cisajscktaltdc.h"
#include "nbody_correlation.h"
#include "anomalous_pair.h"
#include "green_output.h"
#include "wrapperMPI.h"
#include "DefCommon.h"
#include "global.h"

/**
 * @brief MPI-free per-rank observable state loop. See phys_distributed.h.
 *
 * The expec_* call sequence, signatures, and the all_* assignment group are
 * copied verbatim from the serial path in src/phys.c (the source of truth):
 * the eigenvector is loaded into v0 (1-based) exactly as the serial loop
 * does, then expec_energy_flct() moves it into v1 (and puts H*v1 into v0), so
 * the subsequent evaluators -- including expec_totalspin() -- read the
 * eigenvector through their v1 argument. X->Phys.eigen_num stays 0-based.
 */
int phys_stateparallel_local_loop(struct BindStruct *X,
                                  double complex *panel,
                                  long int jb, long int je, long int NN) {
  int rc = 0;
  long int n, j;

  ExpecLocalEnter();                 /* also clears the sticky ExpecLocal error */
  GreenOutputSetPartialSuffix(myrank);

  for (n = jb; n <= je && rc == 0; n++) {
    X->Phys.eigen_num = (int)(n - 1); /* 0-based, as in the serial phys.c */

    for (j = 0; j < NN; j++) {
      v0[j + 1] = panel[(n - jb) * NN + j];
    }

    if (expec_energy_flct(X) != 0) { rc = -1; break; }
    if (expec_cisajs(X, v1) != 0) { rc = -1; break; }
    if (expec_cisajscktaltdc(X, v1) != 0) { rc = -1; break; }
    if (expec_nbodyg(X, v1) != 0) { rc = -1; break; }
    if (expec_anomalousg(X, v1) != 0) { rc = -1; break; }
    if (X->Def.iCalcType == FullDiag) {
      if (expec_totalspin(X, v1) != 0) { rc = -1; break; }
    }
    if (ExpecLocalError()) { rc = -1; break; }

    /* Record this rank's owned segment (index n-1), same assignment group as
       the tail of the serial phys.c loop. */
    X->Phys.all_energy[n - 1] = X->Phys.energy;
    X->Phys.all_doublon[n - 1] = X->Phys.doublon;
    X->Phys.all_sz[n - 1] = X->Phys.Sz;
    X->Phys.all_s2[n - 1] = X->Phys.s2;
    X->Phys.all_num_up[n - 1] = X->Phys.num_up;
    X->Phys.all_num_down[n - 1] = X->Phys.num_down;
  }

  GreenOutputClearPartialSuffix();   /* reached on both success and break paths */
  ExpecLocalLeave();
  return rc;
}
