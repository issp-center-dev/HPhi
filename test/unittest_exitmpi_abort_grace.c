/*
HPhi  -  Quantum Lattice Model Simulator
Copyright (C) 2015 The University of Tokyo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
 * @file unittest_exitmpi_abort_grace.c
 *
 * @brief Regression driver for exitMPI(): a diagnostic that rank 0 prints
 *        right before every rank aborts must reach the output even when
 *        rank 0 is the last rank to call exitMPI().
 *
 * Every rank "detects" the same fatal input error. Rank 0 is delayed by
 * HPHI_TEST_RANK0_DELAY_MS milliseconds (default 1000) before it prints the
 * marker line to stdoutMPI and calls exitMPI(-1); the other ranks call
 * exitMPI(-1) immediately. Without a grace period in exitMPI(), the first
 * MPI_Abort() makes the launcher kill rank 0 while it is still sleeping, so
 * the marker never appears in the output. This is the mechanism behind the
 * intermittent loss of input-error messages observed in MPI test runs.
 *
 * Driven by test/exitmpi_abort_grace_mpi.sh with at least two ranks.
 * With a single rank the program exits with 77 (skip).
 */
#define _POSIX_C_SOURCE 200809L

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "global.h"
#include "wrapperMPI.h"

static void sleep_ms(long ms)
{
  struct timespec req, rem;
  if (ms <= 0) return;
  req.tv_sec = ms / 1000;
  req.tv_nsec = (ms % 1000) * 1000000L;
  while (nanosleep(&req, &rem) != 0 && errno == EINTR) req = rem;
}

int main(int argc, char *argv[])
{
  const char *env = getenv("HPHI_TEST_RANK0_DELAY_MS");
  long delay_ms = (env != NULL) ? atol(env) : 1000;

  InitializeMPI(argc, argv);
  if (nproc < 2) {
    fprintf(stdoutMPI,
            "Skipping: the exitMPI abort grace test needs at least two MPI ranks.\n");
    FinalizeMPI();
    return 77;
  }

  /* Emulate rank 0 lagging behind the other ranks when all of them
     detect the same input error at (almost) the same point. */
  if (myrank == 0) sleep_ms(delay_ms);
  fprintf(stdoutMPI, "Error: exitMPI abort grace marker (rank 0 diagnostic).\n");
  exitMPI(-1);
  return 0; /* not reached */
}
