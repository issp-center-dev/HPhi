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
 * @file wrapperMPI.c
 *
 * @brief MPI wrapper functions for portable parallel code
 *
 * Provides abstraction layer for MPI operations, allowing the same code
 * to compile with or without MPI support. When MPI is disabled, these
 * functions provide serial equivalents.
 *
 * Key global variables set by InitializeMPI():
 * - myrank: MPI process ID (0 for serial)
 * - nproc: Number of MPI processes (1 for serial)
 * - nthreads: Number of OpenMP threads
 * - stdoutMPI: Stdout for rank 0, /dev/null for others (avoids duplicate output)
 *
 * Communication patterns:
 * - MPI_Sendrecv: Primary pattern for inter-process data exchange
 * - Broadcasts: Parameters read by rank 0, broadcast to all
 * - Reductions: Sum/Max across processes for global quantities
 *
 * File I/O wrappers (fopenMPI, fgetsMPI):
 * - Ensure only rank 0 reads/writes parameter files
 * - Broadcast file contents to other ranks
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
#ifdef MPI
#include <mpi.h>
#endif
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wrapperMPI.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <math.h>
#include <complex.h>
#include "splash.h"
#include "global.h"

/**
@brief ExpecLocal mode flag (phase 3a). ON/OFF only via ExpecLocalEnter()/
ExpecLocalLeave(); every other function in this file must use the
ExpecLocal*() accessors below, never this variable directly, so the
invariant (no nesting, no stray writes) stays enforceable in one place.
*/
static int iExpecLocal = 0;
/**
@brief Sticky ExpecLocal error flag (phase 3a). Set by ExpecLocalSetError()
when a defensive guard fires during the local per-rank evaluation loop;
read back by the Mode 1 driver (ExpecLocalError()) to fail that state
without ever calling exitMPI() from inside the loop. Cleared on
ExpecLocalEnter().
*/
static int iExpecLocalError = 0;

/**
@brief Enter ExpecLocal mode. Nesting is a programming error and aborts
(assert). Also clears the sticky error flag for the new local-evaluation
session.
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void ExpecLocalEnter(void) {
  assert(!iExpecLocal);
  iExpecLocal = 1;
  iExpecLocalError = 0;
}/*void ExpecLocalEnter*/
/**
@brief Leave ExpecLocal mode. Asserts the mode was actually active (leaving
without entering is a programming error).
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void ExpecLocalLeave(void) {
  assert(iExpecLocal);
  iExpecLocal = 0;
}/*void ExpecLocalLeave*/
/**
@brief Read whether ExpecLocal mode is currently active.
@return Non-zero iff ExpecLocal mode is active.
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int ExpecLocalActive(void) {
  return iExpecLocal;
}/*int ExpecLocalActive*/
/**
@brief Record a deferred ExpecLocal error (e.g. a defensive raw-MPI guard
fired). Always defined regardless of ExpecLocalActive(), so guard call
sites never need to special-case it.
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
void ExpecLocalSetError(void) {
  iExpecLocalError = 1;
}/*void ExpecLocalSetError*/
/**
@brief Read the sticky ExpecLocal error flag accumulated since the last
ExpecLocalEnter(). Always defined regardless of ExpecLocalActive().
@return Non-zero iff ExpecLocalSetError() was called since ExpecLocalEnter().
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int ExpecLocalError(void) {
  return iExpecLocalError;
}/*int ExpecLocalError*/

/**
@brief MPI initialization wrapper
Process ID (::myrank), Number of processes (::nproc),
Number of threads (::nthreads), and pointer to the standard output
(::stdoutMPI) are specified here.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void InitializeMPI(int argc, char *argv[]){
  int ierr;

#ifdef MPI
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if(ierr != 0) exitMPI(ierr);
  /* Parse HPHI_MPI_NOBATCH on rank 0 and broadcast so every rank agrees
     (a per-rank divergence would deadlock the collective MPI_Sendrecv pattern). */
  if (myrank == 0) {
    const char *e = getenv("HPHI_MPI_NOBATCH");
    iFlgMPIBatch = (e != NULL && atoi(e) != 0) ? 0 : 1;
  }
  MPI_Bcast(&iFlgMPIBatch, 1, MPI_INT, 0, MPI_COMM_WORLD);
#else
  nproc = 1;
  myrank = 0;
  {
    const char *e = getenv("HPHI_MPI_NOBATCH");
    iFlgMPIBatch = (e != NULL && atoi(e) != 0) ? 0 : 1;
  }
#endif
  if (myrank == 0) stdoutMPI = stdout;
  else stdoutMPI = fopen("/dev/null", "w");
  splash();

#pragma omp parallel default(none) shared(nthreads)
#pragma omp master
#ifdef _OPENMP
  nthreads = omp_get_num_threads();
#else
  nthreads=1;
#endif
  fprintf(stdoutMPI, "\n\n#####  Parallelization Info.  #####\n\n");
  fprintf(stdoutMPI, "  OpenMP threads : %d\n", nthreads);
  fprintf(stdoutMPI, "  MPI PEs : %d \n\n", nproc);
  fprintf(stdoutMPI, "  MPI batching : %s\n\n", iFlgMPIBatch ? "ON" : "OFF");
}/*void InitializeMPI(int argc, char *argv[])*/
/**
@brief MPI Finitialization wrapper
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void FinalizeMPI(){
  int ierr;
#ifdef MPI
  ierr = MPI_Finalize();
  if (ierr != 0) fprintf(stderr, "\n  MPI_Finalize() = %d\n\n", ierr);
#endif
  if (myrank != 0) fclose(stdoutMPI);
}
/**
@brief MPI Abortation wrapper
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void exitMPI(
  int errorcode//!<[in] Error-code to be returned as that of this program
)
{
  int ierr;
  fflush(stdout);
#ifdef MPI
  fprintf(stdout,"\n\n #######  [HPhi] You DO NOT have to WORRY about the following MPI-ERROR MESSAGE.  #######\n\n");
  ierr = MPI_Abort(MPI_COMM_WORLD, errorcode);
  ierr = MPI_Finalize();
  if (ierr != 0) fprintf(stderr, "\n  MPI_Finalize() = %d\n\n", ierr);
#endif
  exit(errorcode);
}/*void exitMPI*/
/**
@brief MPI file I/O (open) wrapper.
Only the root node (::myrank = 0) should be open/read/write (small) parameter files.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
FILE* fopenMPI(
  const char* FileName,//!<[in] Input/output file
  const char* mode//!<[in] "w", "r", etc.
){
  FILE* fp;

  if (iExpecLocal) {
    /* ExpecLocal mode: every rank evaluates its own eigenstates
       independently, so every rank opens the file itself at the calling
       rank (no rank-0 gate, no /dev/null fallback for non-zero ranks).
       Callers that need distinct per-rank paths (e.g. green_output partial
       channels) are responsible for making FileName rank-unique; fopenMPI
       itself just stops assuming rank 0 is the only writer/reader. */
    fp = fopen(FileName, mode);
    return fp;
  }

  if (myrank == 0) fp = fopen(FileName, mode);
  else fp = fopen("/dev/null", "w");

  return fp;
}/*FILE* fopenMPI*/
/**
@brief MPI file I/O (get a line, fgets) wrapper.
Only the root node (::myrank = 0) reads and broadcast string.
@return The same as that of fgets
@author Mitsuaki Kawamura (The University of Tokyo)
*/
char* fgetsMPI(
  char* InputString,//!<[out] read line.
  int maxcount,//!<[in] Length of string
  FILE* fp//!<[in] file pointer
){
  int inull;
  char *ctmp;

  /* Not reachable from the FullDiag expec_* evaluation layer (see
     docs/superpowers/specs/2026-07-11-expec-call-inventory.md §1); a debug
     assert catches an accidental future call during ExpecLocal mode
     instead of inventing untested no-communication semantics for it (same
     style as the other assert-forbidden wrappers below, e.g. BarrierMPI). */
  assert(!iExpecLocal);
  ctmp = InputString;
  inull = 0;
  if (myrank == 0) {
    ctmp = fgets(InputString, maxcount, fp);
    if (ctmp == NULL){
      inull = 1;
    }
    
    while(*InputString == '\n' || strncmp(InputString, "#", 1)==0){
      ctmp = fgets(InputString, maxcount, fp);
      if (ctmp == NULL){
        inull=1;
        break;
      }
    }
  }
#ifdef MPI
  MPI_Bcast(InputString, maxcount, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&inull, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  if (myrank != 0 && inull == 1) {
    ctmp = NULL;
  }

  return ctmp;
}/*char* fgetsMPI*/
/**
@brief MPI barrier wrapper.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void BarrierMPI(){
  /* Not reachable from the FullDiag expec_* evaluation layer (see
     docs/superpowers/specs/2026-07-11-expec-call-inventory.md §1); a debug
     assert catches an accidental future call during ExpecLocal mode
     instead of inventing untested no-communication semantics for it. */
  assert(!iExpecLocal);
#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}/*void BarrierMPI()*/
/**
@brief MPI wrapper function to obtain maximum unsigned
long integer across processes.
@return Maximum value across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
unsigned long int MaxMPI_li(
  unsigned long int idim//!<[in] Value to be maximized
){
  /* Not reachable from the FullDiag expec_* evaluation layer (see
     docs/superpowers/specs/2026-07-11-expec-call-inventory.md §1); a debug
     assert catches an accidental future call during ExpecLocal mode
     instead of inventing untested no-communication semantics for it. */
  assert(!iExpecLocal);
#ifdef MPI
  int ierr;
  ierr = MPI_Allreduce(MPI_IN_PLACE, &idim, 1,
    MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  if(ierr != 0) exitMPI(-1);
#endif
  return(idim);
}/*unsigned long int MaxMPI_li*/
/**
@brief MPI wrapper function to obtain maximum Double
across processes.
@return Maximum value across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
double MaxMPI_d(
  double dvalue//!<[in] Value to be maximized
){
  /* Not reachable from the FullDiag expec_* evaluation layer (see
     docs/superpowers/specs/2026-07-11-expec-call-inventory.md §1); a debug
     assert catches an accidental future call during ExpecLocal mode
     instead of inventing untested no-communication semantics for it. */
  assert(!iExpecLocal);
#ifdef MPI
  int ierr;
  ierr = MPI_Allreduce(MPI_IN_PLACE, &dvalue, 1,
    MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if(ierr != 0) exitMPI(-1);
#endif
  return(dvalue);
}/*double MaxMPI_d*/
/**
@brief MPI wrapper function to obtain sum of Double
complex across processes.
@return Sumed value across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
double complex SumMPI_dc(
  double complex norm//!<[in] Value to be summed
){
  /* ExpecLocal mode: each rank evaluates its own eigenstate independently,
     so this reduction must not touch other ranks -- return the local
     value unchanged (per docs/superpowers/specs/
     2026-07-11-expec-call-inventory.md §1, the expec_* layer's dominant
     call). */
  if (iExpecLocal) return norm;
#ifdef MPI
  int ierr;
  ierr = MPI_Allreduce(MPI_IN_PLACE, &norm, 1,
    MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
  if(ierr != 0) exitMPI(-1);
#endif
  return(norm);
}/*double complex SumMPI_dc*/
/**
@brief MPI wrapper function to obtain sum of Double
across processes.
@return Sumed value across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
double SumMPI_d(
  double norm//!<[in] Value to be summed
){
  /* ExpecLocal mode: no-communication pass-through (see SumMPI_dc above). */
  if (iExpecLocal) return norm;
#ifdef MPI
  int ierr;
  ierr = MPI_Allreduce(MPI_IN_PLACE, &norm, 1,
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
  if(ierr != 0) exitMPI(-1);
#endif
  return(norm);
}/*double SumMPI_d*/
/**
@brief MPI wrapper function to obtain sum of unsigned
long integer across processes.
@return Sumed value across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
unsigned long int SumMPI_li(
  unsigned long int idim//!<[in] Value to be summed
){
  /* ExpecLocal mode: no-communication pass-through (see SumMPI_dc above).
     Not currently called from the expec_* layer, but kept on the frozen
     ExpecLocal allow-list as harmless headroom. */
  if (iExpecLocal) return idim;
#ifdef MPI
  int ierr;
  ierr = MPI_Allreduce(MPI_IN_PLACE, &idim, 1,
    MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  if(ierr != 0) exitMPI(-1);
#endif
  return(idim);
}/*unsigned long int SumMPI_li*/
/**
@brief MPI wrapper function to obtain sum of
integer across processes.
@return Sumed value across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
int SumMPI_i(
  int idim//!<[in] Value to be summed
) {
  /* ExpecLocal mode: no-communication pass-through (see SumMPI_dc above).
     Not currently called from the expec_* layer, but kept on the frozen
     ExpecLocal allow-list as harmless headroom. */
  if (iExpecLocal) return idim;
#ifdef MPI
  int ierr;
  ierr = MPI_Allreduce(MPI_IN_PLACE, &idim, 1,
                       MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(ierr != 0) exitMPI(-1);
#endif
  return(idim);
}/*int SumMPI_i*/
/**
@brief MPI wrapper function to broadcast an integer across processes.
@return Broadcasted value across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
int BcastMPI_i(
  int root,//!<[in] The source process of the broadcast
  int idim//!<[in] Value to be broadcasted
) {
  int idim0;
  idim0 = idim;
#ifdef MPI
  if (MPI_Bcast(&idim0, 1, MPI_INT, root, MPI_COMM_WORLD) != MPI_SUCCESS) {
    exitMPI(-1);
  }
#endif
  return(idim0);
}/*int BcastMPI_i*/
/**
@brief MPI wrapper function to broadcast unsigned long
integer across processes.
@return Broadcasted value across processes.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
unsigned long int BcastMPI_li(
  int root,//!<[in] The source process of the broadcast
  unsigned long int idim//!<[in] Value to be broadcasted
) {
  unsigned long int idim0;
  /* Not reachable from the FullDiag expec_* evaluation layer (see
     docs/superpowers/specs/2026-07-11-expec-call-inventory.md §1); a debug
     assert catches an accidental future call during ExpecLocal mode
     instead of inventing untested no-communication semantics for it. */
  assert(!iExpecLocal);
  idim0 = idim;
#ifdef MPI
    MPI_Bcast(&idim0, 1, MPI_UNSIGNED_LONG, root, MPI_COMM_WORLD);
#endif
  return(idim0);
}/*unsigned long int BcastMPI_li*/
/**
@brief Compute norm of process-distributed vector
@f$|{\bf v}_1|^2@f$
@return Norm @f$|{\bf v}_1|^2@f$
*/
double NormMPI_dc(
  unsigned long int idim,//!<[in] Local dimension of vector
  double complex *_v1//!<[in] [idim] vector to be producted
){
  double complex cdnorm=0;
  double dnorm =0;
  unsigned long int i;
  /* Not reachable from the FullDiag expec_* evaluation layer (see
     docs/superpowers/specs/2026-07-11-expec-call-inventory.md §1); a debug
     assert catches an accidental future call during ExpecLocal mode
     instead of inventing untested no-communication semantics for it.
     (The inner SumMPI_dc call below is already a no-op under
     ExpecLocal, so if this guard is ever relaxed, the reduction itself
     needs no further change.) */
  assert(!iExpecLocal);
  //DEBUG
#pragma omp parallel for default(none) private(i) firstprivate(myrank) shared(_v1, idim) reduction(+: cdnorm)
  for(i=1;i<=idim;i++){
    cdnorm += conj(_v1[i])*_v1[i];
  }
#ifdef MPI
  cdnorm = SumMPI_dc(cdnorm);
#endif
  dnorm=creal(cdnorm);
  dnorm=sqrt(dnorm);

  return dnorm;
}/*double NormMPI_dc*/
/**
@brief Compute conjugate scaler product of process-distributed vector
@f${\bf v}_1^* \cdot {\bf v}_2@f$
@return Conjugate scaler product @f${\bf v}_1^* \cdot {\bf v}_2@f$
*/
double complex VecProdMPI(
  long unsigned int ndim,//!<[in] Local dimension of vector
  double complex *v1,//!<[in] [ndim] vector to be producted
  double complex *v2//!<[in] [ndim] vector to be producted
){
  long unsigned int idim;
  double complex prod;

  /* Not reachable from the FullDiag expec_* evaluation layer (see
     docs/superpowers/specs/2026-07-11-expec-call-inventory.md §1); a debug
     assert catches an accidental future call during ExpecLocal mode
     instead of inventing untested no-communication semantics for it.
     (The inner SumMPI_dc call below is already a no-op under
     ExpecLocal, so if this guard is ever relaxed, the reduction itself
     needs no further change.) */
  assert(!iExpecLocal);
  prod = 0.0;
#pragma omp parallel for default(none) shared(v1,v2,ndim) private(idim) reduction(+: prod)
  for (idim = 1; idim <= ndim; idim++) prod += conj(v1[idim]) * v2[idim];
  prod = SumMPI_dc(prod);

  return(prod);
}/*double complex VecProdMPI*/
