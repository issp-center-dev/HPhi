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

#ifndef HPHI_EXPEC_TRACE_HAM_H
#define HPHI_EXPEC_TRACE_HAM_H

#include <complex.h>
#include <stddef.h>
#include "struct.h"

/**
 * @file expec_trace_ham.h
 *
 * Phase-3c Task 2: the per-rank CSR Hamiltonian collector for the ExpecMode 2
 * energy-family trace kernel.
 *
 * TraceHamCollect() re-runs makeHam()'s full enumeration in the
 * HAM_SINK_TRACE_COLLECT sink mode (Task 1) into a two-pass count/fill
 * collector, producing a compact, merged, per-row-sorted CSR for this rank.
 * An exact pre-allocation memory gate (TraceHamGatedBytes) demotes the whole
 * kernel back to the Mode-1 fallback when the projected peak exceeds the cap.
 *
 * Row/column indices are 0-based in the CSR: matrix element (i, j) is stored
 * for row i (== makeHam's 1-based irow minus 1) at column j (== jcol minus 1).
 * The matrix diagonal H(j,j) is a normal CSR entry; csr->diag[] holds the
 * SEPARATE energy-family diagonal coefficient arrays (D/N/S), which Task 3
 * fills -- Task 2 allocates n_diag of them (frozen per-model table) and zeroes
 * them.
 */
typedef struct {
  long int n;          /* matrix dimension (idim_max) */
  long int nnz;        /* merged entries; rowptr[n] == nnz */
  long int *rowptr;    /* size n+1 */
  long int *colidx;    /* size >= nnz (raw capacity retained) */
  double complex *val; /* size >= nnz */
  double complex *y;   /* kernel scratch, size n */
  int n_diag;          /* number of diagonal coefficient arrays (0..3) */
  double *diag[3];     /* D(k), N(k), S(k) as applicable, each size n */
} TraceHamCsr;

/**
 * Collect the Hamiltonian of X into a compact CSR (this rank's column range,
 * or the full matrix at nproc==1).
 *
 * Returns 1 on success (csr populated), 0 on demotion (gate exceeded or an
 * allocation failed; csr fully freed and zeroed, sink mode restored).
 *
 * cap_bytes:     the per-quantity byte cap the projected peak must fit under.
 * fail_alloc_at: test hook; -1 in production. Allocation number k (0-based, in
 *                the fixed order rowptr, cursors, colidx, val, sort-workspace,
 *                y, then the n_diag coefficient arrays) is forced to fail when
 *                k == fail_alloc_at, exercising the demotion/cleanup path.
 */
int TraceHamCollect(struct BindStruct *X, size_t cap_bytes,
                    int fail_alloc_at, TraceHamCsr *csr);

/** Free every buffer csr owns and zero the struct. Safe on a zeroed csr. */
void TraceHamFree(TraceHamCsr *csr);

/**
 * Exact projected peak (bytes) for the given counting quantities -- the value
 * the memory gate compares against cap_bytes, exposed for the boundary test.
 * Returns SIZE_MAX on any internal overflow (which can never pass a real cap).
 */
size_t TraceHamGatedBytes(long int n, long int nnz_raw, long int k_max,
                          int n_diag);

#endif /* HPHI_EXPEC_TRACE_HAM_H */
