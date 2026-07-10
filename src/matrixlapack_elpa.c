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

#ifdef _ELPA
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <mpi.h>
#include <elpa/elpa.h>
#include "matrixlapack_elpa.h"

/* GPU enable option name; single point of change for future amd-gpu etc.
   (design doc section 3, GPU policy) */
static const char *ELPA_GPU_OPTION = "nvidia-gpu";

/* Share a local error across all ranks so every rank takes the same
   branch before each collective ELPA call (design doc section 4). */
static int SyncError(int ierr) {
  int gerr = 0;
  MPI_Allreduce(&ierr, &gerr, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  return gerr;
}

/* elpa_set for int values with mandatory error check. */
static int SetElpaInt(elpa_t handle, const char *name, int value) {
  int error = ELPA_OK;
  elpa_set(handle, name, value, &error);
  if (error != ELPA_OK) {
    fprintf(stdout, "  Error: elpa_set(\"%s\", %d) failed: %s\n",
            name, value, elpa_strerr(error));
    return -1;
  }
  return 0;
}

/**
 * @brief Diagonalize a 2D block-cyclic distributed Hermitian matrix with
 * ELPA (design doc section 3). Eigenvalues are returned on all ranks in w;
 * eigenvectors stay distributed in Z_distr. A_distr is destroyed.
 * @return 0 on success, -1 on failure (same value on all ranks).
 */
int diag_elpa_cmp(int xNsize, double complex *A_distr,
                  double complex *Z_distr, double *w,
                  int local_nrows, int local_ncols,
                  int myrow, int mycol, int nblk, int ngpu) {
  elpa_t handle = NULL;
  int error = ELPA_OK;
  int ierr = 0;

  if (elpa_init(20211125) != ELPA_OK) {
    fprintf(stdout, "  Error: the linked ELPA is older than API 20211125 (2021.11).\n");
    ierr = -1;
  }
  if (SyncError(ierr) != 0) return -1;

  handle = elpa_allocate(&error);
  if (error != ELPA_OK) {
    fprintf(stdout, "  Error: elpa_allocate failed: %s\n", elpa_strerr(error));
    ierr = -1;
  }

  /* Mandatory parameters: set BEFORE elpa_setup (ELPA manual sec. 2). */
  if (ierr == 0) {
    if (SetElpaInt(handle, "na", xNsize) != 0 ||
        SetElpaInt(handle, "nev", xNsize) != 0 ||
        SetElpaInt(handle, "local_nrows", local_nrows) != 0 ||
        SetElpaInt(handle, "local_ncols", local_ncols) != 0 ||
        SetElpaInt(handle, "nblk", nblk) != 0 ||
        SetElpaInt(handle, "mpi_comm_parent",
                   (int)MPI_Comm_c2f(MPI_COMM_WORLD)) != 0 ||
        SetElpaInt(handle, "process_row", myrow) != 0 ||
        SetElpaInt(handle, "process_col", mycol) != 0) {
      ierr = -1;
    }
  }
  if (SyncError(ierr) != 0) goto cleanup_fail;

  error = elpa_setup(handle);
  if (error != ELPA_OK) {
    fprintf(stdout, "  Error: elpa_setup failed: %s\n", elpa_strerr(error));
    ierr = -1;
  }
  if (SyncError(ierr) != 0) goto cleanup_fail;

  /* Tunable runtime options: set AFTER elpa_setup (ELPA manual sec. 2).
     GPU: 1stage is usually faster on GPU; 2stage on CPU. */
  if (ngpu >= 1) {
    if (SetElpaInt(handle, "solver", ELPA_SOLVER_1STAGE) != 0) ierr = -1;
    if (ierr == 0 && SetElpaInt(handle, ELPA_GPU_OPTION, 1) != 0) {
      fprintf(stdout,
              "  Error: the linked ELPA has no NVIDIA GPU support.\n"
              "         Set \"NGPU 0\" in calcmod.def for CPU execution.\n");
      ierr = -1;
    }
  } else {
    /* 2stage is preferred on CPU, but ELPA's 2stage solver returned
       deterministically inaccurate eigenvectors (residual ~1e-6) for
       capped block sizes (observed with nblk=24 on a 4x2 grid, ELPA
       2025.06.001). Capping only happens for matrices small relative to
       the process grid, where 1stage is safe and fast enough, so fall
       back to 1stage whenever the block size was capped below ELPA_NBLK. */
    int solver = (nblk < ELPA_NBLK) ? ELPA_SOLVER_1STAGE : ELPA_SOLVER_2STAGE;
    if (SetElpaInt(handle, "solver", solver) != 0) ierr = -1;
  }
  if (SyncError(ierr) != 0) goto cleanup_fail;

#ifdef _ELPA_GPU
  if (ngpu >= 1) {
    error = elpa_setup_gpu(handle);
    if (error != ELPA_OK) {
      fprintf(stdout,
              "  Error: elpa_setup_gpu failed: %s\n"
              "         Check the GPU environment, or set \"NGPU 0\" for CPU execution.\n",
              elpa_strerr(error));
      ierr = -1;
    }
    if (SyncError(ierr) != 0) goto cleanup_fail;
  }
#else
  if (ngpu >= 1) {
    /* readdef.c rejects this combination at startup; defense in depth. */
    fprintf(stdout, "  Error: this HPhi build has no ELPA GPU API (_ELPA_GPU).\n");
    ierr = -1;
    if (SyncError(ierr) != 0) goto cleanup_fail;
  }
#endif

  /* Type-generic macro resolves to the double-complex solver. */
  elpa_eigenvectors(handle, A_distr, w, Z_distr, &error);
  if (error != ELPA_OK) {
    fprintf(stdout, "  Error: elpa_eigenvectors failed: %s\n", elpa_strerr(error));
    ierr = -1;
  }
  if (SyncError(ierr) != 0) goto cleanup_fail;

  elpa_deallocate(handle, &error);
  elpa_uninit(&error);
  return 0;

cleanup_fail:
  if (handle != NULL) {
    elpa_deallocate(handle, &error);
  }
  elpa_uninit(&error);
  return -1;
}
#endif /* _ELPA */
