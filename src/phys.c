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
 * @file phys.c
 *
 * @brief Calculate physical quantities for full diagonalization method
 *
 * After full diagonalization obtains all eigenvalues and eigenvectors,
 * this module computes physical observables for each eigenstate.
 *
 * Computed quantities (for each eigenstate i):
 * - Energy: E_i (from diagonalization)
 * - Number of particles: \f$\langle n\rangle\f$ (for Hubbard models)
 * - Total spin: \f$\langle S^2\rangle\f$, \f$\langle S_z\rangle\f$
 * - One-body Green's functions: \f$\langle c^\dagger_i c_j\rangle\f$
 * - Two-body Green's functions: \f$\langle c^\dagger_i c_j c^\dagger_k c_l\rangle\f$
 *
 * ScaLAPACK support:
 * - When _SCALAPACK is defined, eigenvectors are distributed across ranks
 * - GetEigenVector() gathers eigenvector to rank 0 for observable calculation
 *
 * Output:
 * - Results written to output files (one per eigenvalue)
 * - Energy spectrum to zvo_energy.dat
 *
 * @version 0.1, 0.2
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
#include "phys.h"
#include "expec_energy_flct.h"
#include "expec_totalspin.h"
#include "expec_cisajs.h"
#include "expec_cisajscktaltdc.h"
#include "nbody_correlation.h"
#include "anomalous_pair.h"
#include "green_output.h"
#include "wrapperMPI.h"
#include "FileIO.h"
#include "DefCommon.h"
#include "phys_distributed.h"
#include <assert.h>
#include <limits.h>
#ifdef _SCALAPACK
#include "matrixscalapack.h"
#endif

/**
 * @brief Write one global FullDiag eigenvector using the established binary
 * eigenvector-file layout.
 *
 * Iterative solvers write one local-basis slice per MPI rank. FullDiag uses
 * MPI to distribute dense-matrix columns or eigenstates instead, while each
 * eigenvector itself still spans the complete Hilbert space. Consequently a
 * FullDiag state has one file, with the decomposition suffix fixed to rank 0.
 */
int FullDiagOutputEigenvector(const struct BindStruct *X,
                              unsigned long int state,
                              const double complex *eigenvector) {
  char sdt[D_FileNameMax];
  FILE *fp = NULL;
  int iteration = 0;
  unsigned long int dimension;
  double complex unused = 0.0;
  int nfn;
  int failed = 0;

  if (X == NULL || eigenvector == NULL || state > (unsigned long int)INT_MAX) {
    fprintf(stderr, "Error: invalid FullDiag eigenvector output request.\n");
    return -1;
  }
  dimension = X->Check.idim_max;
  nfn = snprintf(sdt, sizeof(sdt), cFileNameOutputEigen,
                 X->Def.CDataFileHead, (int)state, 0);
  if (nfn < 0 || nfn >= (int)sizeof(sdt)) {
    fprintf(stderr,
            "Error: FullDiag eigenvector file name is too long (state=%lu).\n",
            state);
    return -1;
  }
  if (childfopenALL(sdt, "wb", &fp) != 0) return -1;

  if (fwrite(&iteration, sizeof(iteration), 1, fp) != 1 ||
      fwrite(&dimension, sizeof(dimension), 1, fp) != 1 ||
      fwrite(&unused, sizeof(unused), 1, fp) != 1 ||
      fwrite(eigenvector, sizeof(*eigenvector), dimension, fp) != dimension) {
    fprintf(stderr,
            "Error: failed to write FullDiag eigenvector file (%s).\n", sdt);
    failed = 1;
  }
  if (fclose(fp) != 0) {
    fprintf(stderr,
            "Error: failed to close FullDiag eigenvector file (%s).\n", sdt);
    failed = 1;
  }
  return failed ? -1 : 0;
}

/**
 * @brief Compute physical quantities for all eigenstates from full diagonalization
 *
 * Iterates over neig eigenstates and computes observables for each.
 * For ScaLAPACK, eigenvectors are gathered to rank 0 before computation.
 *
 * For each eigenstate:
 * 1. Copy eigenvector to v0 (working vector)
 * 2. Call expec_energy_flct() for energy (already known, but variance check)
 * 3. Call expec_totalspin() for \f$\langle S^2\rangle\f$, \f$\langle S_z\rangle\f$
 * 4. Call expec_cisajs() for one-body Green's functions
 * 5. Call expec_cisajscktaltdc() for two-body Green's functions
 *
 * @param X BindStruct with model and calculation parameters [inout]
 * @param neig Number of eigenvalues/eigenvectors to process [in]
 *
 * @version 0.2 Added general spin output
 * @version 0.1
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void phys(struct BindStruct *X, //!<[inout]
          unsigned long int neig //!<[in]
) {
  long unsigned int i, j, i_max;
  double tmp_N;
  int eigen_output_failed = 0;
  i_max = X->Check.idim_max;
#ifdef _SCALAPACK
  double complex *vec_tmp;
  int rank;
  if (use_scalapack && X->Def.iExpecMode != EXPECMODE_SERIAL) {
    int phys_rc = phys_stateparallel(X, neig);
    FreeDistributedEigenvectors(&Z_vec, descZ_vec, &use_scalapack);
    if (phys_rc != 0) exitMPI(-1);
    assert(!ExpecLocalActive());
    return;
  }
  if(use_scalapack){
  /* S2/Sz are now computed on rank 0 in the state loop below (ExpecLocal-
     wrapped), so the former "total spin is not calculated" notice is gone. */
  vec_tmp = malloc(i_max*sizeof(double complex));
  }
#endif
  if (GreenOutputInitializeAggregateFiles(X) != 0) {
    fprintf(stderr, "Error: initialize aggregate Green output files.\n");
    exitMPI(-1);
  }
  for (i = 0; i < neig; i++) {
#ifdef _SCALAPACK
    for (j = 0; j < i_max; j++) {
      v0[j + 1] = 0.0;
    }
    if(use_scalapack){
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#ifdef _ELPA
      if (X->Def.iSolver == SOLVER_ELPA) {
        GetEigenVectorBlock(i, i_max, Z_vec, descZ_vec, vec_tmp);
      } else {
        GetEigenVector(i, i_max, Z_vec, descZ_vec, vec_tmp);
      }
#else
      GetEigenVector(i, i_max, Z_vec, descZ_vec, vec_tmp);
#endif
      if(rank == 0) {
        for (j = 0; j < i_max; j++) {
          v0[j + 1] = vec_tmp[j];
        }
      }
      else{
        for (j = 0; j < i_max; j++) {
          v0[j + 1] = 0.0;
        }
      }
    } else {
        if (X->Def.iCalcType == FullDiag) {
            if (myrank == 0) {
                for (j = 0; j < i_max; j++) {
                    v0[j + 1] = L_vec[i][j];
                }
            }
        } else {
            for (j = 0; j < i_max; j++) {
                v0[j + 1] = L_vec[i][j];
            }

        }
    }
#else
    for (j = 0; j < i_max; j++) {
      v0[j + 1] = L_vec[i][j];
    }    
#endif

    if (X->Def.iCalcType == FullDiag &&
        X->Def.iOutputEigenVec == TRUE && myrank == 0 &&
        !eigen_output_failed) {
#ifdef _SCALAPACK
      const double complex *eigenvector =
          use_scalapack ? vec_tmp : L_vec[i];
#else
      const double complex *eigenvector = L_vec[i];
#endif
      if (FullDiagOutputEigenvector(X, i, eigenvector) != 0)
        eigen_output_failed = 1;
    }

    X->Phys.eigen_num = i;
    if (expec_energy_flct(X) != 0) {
      fprintf(stderr, "Error: calc expec_energy.\n");
      exitMPI(-1);
    }
    if (expec_cisajs(X, v1) != 0) {
      fprintf(stderr, "Error: calc OneBodyG.\n");
      exitMPI(-1);
    }
    if (expec_cisajscktaltdc(X, v1) != 0) {
      fprintf(stderr, "Error: calc TwoBodyG.\n");
      exitMPI(-1);
    }
    if (expec_nbodyg(X, v1) != 0) {
      fprintf(stderr, "Error: calc NBodyG.\n");
      exitMPI(-1);
    }
    if (expec_anomalousg(X, v1) != 0) {
      fprintf(stderr, "Error: calc AnomalousG.\n");
      exitMPI(-1);
    }
    
#ifdef _SCALAPACK
    if(use_scalapack){
      if (X->Def.iCalcType == FullDiag) {
        /* Mode 0 (distributed) S2/Sz unification: rank 0 holds the gathered
           eigenvector in v1 (moved there from v0 by expec_energy_flct), so
           compute totalspin locally under an ExpecLocal wrapper (reductions
           become no-communication pass-throughs). Other ranks hold a zeroed
           vector; their values are not used (display/all_* take rank 0's). */
        ExpecLocalEnter();
        if (myrank == 0) {
          if (expec_totalspin(X, v1) != 0) { ExpecLocalLeave(); exitMPI(-1); }
        } else {
          X->Phys.s2 = 0.0;
          X->Phys.Sz = 0.0;
        }
        ExpecLocalLeave();
      }
    }else{
      if (X->Def.iCalcType == FullDiag) {
        if (expec_totalspin(X, v1) != 0) {
          fprintf(stderr, "Error: calc TotalSpin.\n");
          exitMPI(-1);
        }
      }
    } 
#else
    if (X->Def.iCalcType == FullDiag) {
      if (expec_totalspin(X, v1) != 0) {
        fprintf(stderr, "Error: calc TotalSpin.\n");
        exitMPI(-1);
      }
    }
#endif
    
    if (X->Def.iCalcModel == Spin || X->Def.iCalcModel == SpinGC) {
      tmp_N = X->Def.NsiteMPI;
    } else {
      tmp_N = X->Phys.num_up + X->Phys.num_down;
    }

    if (X->Def.iCalcType == FullDiag){
      /* Unified serial format (S2 column) for both distributed and
         non-distributed paths: distributed Mode 0 now computes S2/Sz on
         rank 0 above rather than zero-filling. */
      fprintf(stdoutMPI, "i=%5ld Energy=%10lf N=%10lf Sz=%10lf S2=%10lf Doublon=%10lf \n", i, X->Phys.energy, tmp_N,
              X->Phys.Sz, X->Phys.s2, X->Phys.doublon);
    }
    else if (X->Def.iCalcType == CG)
      fprintf(stdoutMPI, "i=%5ld Energy=%10lf N=%10lf Sz=%10lf Doublon=%10lf \n", i, X->Phys.energy, tmp_N,
              X->Phys.Sz, X->Phys.doublon);
    X->Phys.all_energy[i] = X->Phys.energy;
    X->Phys.all_doublon[i] = X->Phys.doublon;
    X->Phys.all_sz[i] = X->Phys.Sz;
    X->Phys.all_s2[i] = X->Phys.s2;
    X->Phys.all_num_up[i] = X->Phys.num_up;
    X->Phys.all_num_down[i] = X->Phys.num_down;
  }
#ifdef MPI
  if (X->Def.iCalcType == FullDiag && X->Def.iOutputEigenVec == TRUE) {
    int any_eigen_output_failed = 0;
    MPI_Allreduce(&eigen_output_failed, &any_eigen_output_failed, 1, MPI_INT,
                  MPI_MAX, MPI_COMM_WORLD);
    eigen_output_failed = any_eigen_output_failed;
  }
#endif
#ifdef _SCALAPACK
  if(use_scalapack) {
    free(vec_tmp);
    FreeDistributedEigenvectors(&Z_vec, descZ_vec, &use_scalapack);
  }
#endif
  if (eigen_output_failed) exitMPI(-1);
  assert(!ExpecLocalActive());
}
