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
 * @file expec_trace_ham.c
 *
 * Phase-3c Task 2: two-pass CSR Hamiltonian collector (see expec_trace_ham.h).
 *
 * The collector installs a static counting sink, then a static fill sink, into
 * the global hamCollectSink hook and runs makeHam() twice in
 * HAM_SINK_TRACE_COLLECT mode:
 *   pass 1 (count): per-row raw-entry counts, prefix-summed in place to rowptr;
 *   pass 2 (fill):  each streamed (irow, jcol, val) written into its row
 *                   segment through a per-row cursor.
 * Each row segment is then stably sorted by column (insertion sort < 16
 * entries, else a bottom-up stable merge sort over the one k_max-sized
 * workspace), adjacent duplicates are summed, and the surviving entries are
 * forward-compacted so colidx/val hold the merged matrix with rowptr[n]==nnz.
 *
 * The sink hook and iHamSinkMode are saved on entry and restored on EVERY exit
 * path (success, gate demotion, injected/real allocation failure), so a
 * caller's dense-mode state is never left corrupted.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <stdint.h>
#include <limits.h>
#include <assert.h>

#include "struct.h"
#include "global.h"
#include "DefCommon.h"
#include "hamstore.h"
#include "makeHam.h"
#include "wrapperMPI.h"
#include "expec_trace_ham.h"

/* ------------------------------------------------------------------ *
 *  File-local collector state (shared with the two static sinks).
 *  TraceHamCollect() runs the two makeHam() passes strictly serially, so
 *  a single set of module globals is sufficient; the collector is not
 *  re-entrant, matching every other makeHam()-driven code path.
 * ------------------------------------------------------------------ */
static long int         g_n;         /* matrix dimension (idim_max)          */
static long int        *g_rowptr;    /* counts -> prefix-sum offsets, n+1     */
static long int        *g_cursor;    /* per-row fill cursors, n+1             */
static long int        *g_colidx;    /* raw then merged column indices        */
static double complex  *g_val;       /* raw then merged values                */
static int              g_overflow;  /* set if a per-row count overflows      */

/* One sort-workspace entry: a (column, value) pair kept together so the
 * stable merge sort permutes both arrays in lockstep. */
typedef struct { long int col; double complex val; } TraceEnt;

/* Counting sink: bump the raw entry count of row `irow` (1-based). */
static void trace_count_sink(long int irow, long int jcol, double complex val) {
  (void)jcol;
  (void)val;
  if (irow < 1 || irow > g_n) { g_overflow = 1; return; }
  if (g_rowptr[irow] == LONG_MAX) { g_overflow = 1; return; }
  g_rowptr[irow]++;
}

/* Fill sink: append (jcol, val) into row `irow`'s segment (0-based column). */
static void trace_fill_sink(long int irow, long int jcol, double complex val) {
  long int pos;
  if (irow < 1 || irow > g_n) return;
  pos = g_cursor[irow]++;
  g_colidx[pos] = jcol - 1;
  g_val[pos] = val;
}

/* ------------------------------------------------------------------ *
 *  Checked byte arithmetic (mirrors phase-3b's guards).
 * ------------------------------------------------------------------ */
static int checked_bytes(uintmax_t a, uintmax_t b, size_t *out) {
  if (a != 0 && b > UINTMAX_MAX / a) return 0;
  uintmax_t r = a * b;
  if (r > SIZE_MAX) return 0;
  *out = (size_t)r;
  return 1;
}

static int add_size(size_t *acc, size_t add) {
  if (*acc > SIZE_MAX - add) return 0;
  *acc += add;
  return 1;
}

size_t TraceHamGatedBytes(long int n, long int nnz_raw, long int k_max,
                          int n_diag) {
  size_t total = 0, t;
  uintmax_t np1 = (uintmax_t)n + 1u;
  uintmax_t nnz = (uintmax_t)nnz_raw;
  uintmax_t nn  = (uintmax_t)n;
  uintmax_t nd  = (uintmax_t)(n_diag > 0 ? n_diag : 0);
  uintmax_t km  = (uintmax_t)k_max;
  uintmax_t ws_elem = (uintmax_t)sizeof(long int) + (uintmax_t)sizeof(double complex);

  /* rowptr */
  if (!checked_bytes(np1, sizeof(long int), &t) || !add_size(&total, t)) return SIZE_MAX;
  /* cursors */
  if (!checked_bytes(np1, sizeof(long int), &t) || !add_size(&total, t)) return SIZE_MAX;
  /* colidx */
  if (!checked_bytes(nnz, sizeof(long int), &t) || !add_size(&total, t)) return SIZE_MAX;
  /* val */
  if (!checked_bytes(nnz, sizeof(double complex), &t) || !add_size(&total, t)) return SIZE_MAX;
  /* y */
  if (!checked_bytes(nn, sizeof(double complex), &t) || !add_size(&total, t)) return SIZE_MAX;
  /* diagonal coefficient arrays: n_diag * n * sizeof(double) */
  {
    size_t nd_n;
    if (!checked_bytes(nd, nn, &nd_n)) return SIZE_MAX;
    if (!checked_bytes((uintmax_t)nd_n, sizeof(double), &t) || !add_size(&total, t)) return SIZE_MAX;
  }
  /* sort workspace: k_max * (sizeof(long int) + sizeof(double complex)) */
  if (!checked_bytes(km, ws_elem, &t) || !add_size(&total, t)) return SIZE_MAX;

  return total;
}

/* ------------------------------------------------------------------ *
 *  Frozen per-model number of diagonal coefficient arrays (design §3b).
 * ------------------------------------------------------------------ */
static int trace_model_n_diag(int iCalcModel) {
  switch (iCalcModel) {
    case Hubbard:
    case HubbardNConserved:
    case Kondo:
    case KondoNConserved:
    case tJ:
    case tJNConserved:
    case HubbardGC:
    case KondoGC:
    case tJGC:
      return 3;      /* D(k), N(k), S(k) */
    case SpinGC:
      return 1;      /* S(k) */
    case Spin:
    default:
      return 0;      /* canonical Spin (and any ineligible model) */
  }
}

/* ------------------------------------------------------------------ *
 *  Gated allocator: allocation number *nalloc (0-based) fails when it
 *  equals fail_at, driving the demotion path from the unit test.
 * ------------------------------------------------------------------ */
static void *galloc(size_t sz, int *nalloc, int fail_at) {
  int idx = (*nalloc)++;
  if (idx == fail_at) return NULL;
  if (sz == 0) sz = 1;
  return malloc(sz);
}

/* Stable insertion sort of a (col,val) segment. */
static void insertion_seg(long int *col, double complex *v, long int len) {
  long int i, j;
  for (i = 1; i < len; i++) {
    long int c = col[i];
    double complex val = v[i];
    j = i - 1;
    while (j >= 0 && col[j] > c) {
      col[j + 1] = col[j];
      v[j + 1] = v[j];
      j--;
    }
    col[j + 1] = c;
    v[j + 1] = val;
  }
}

/* Bottom-up stable merge sort of a (col,val) segment, ws sized >= len. */
static void merge_sort_seg(long int *col, double complex *v, long int len,
                           TraceEnt *ws) {
  long int width;
  for (width = 1; width < len; width *= 2) {
    long int i;
    for (i = 0; i < len; i += 2 * width) {
      long int lo = i;
      long int mid = (i + width < len) ? i + width : len;
      long int hi = (i + 2 * width < len) ? i + 2 * width : len;
      long int a = lo, b = mid, k = lo;
      while (a < mid && b < hi) {
        if (col[a] <= col[b]) { ws[k].col = col[a]; ws[k].val = v[a]; a++; }
        else                  { ws[k].col = col[b]; ws[k].val = v[b]; b++; }
        k++;
      }
      while (a < mid) { ws[k].col = col[a]; ws[k].val = v[a]; a++; k++; }
      while (b < hi)  { ws[k].col = col[b]; ws[k].val = v[b]; b++; k++; }
      for (k = lo; k < hi; k++) { col[k] = ws[k].col; v[k] = ws[k].val; }
    }
  }
}

/* Restore the saved sink hook and mode; used on every exit path. */
static void trace_restore_sink(int prev_mode,
                               void (*prev_hook)(long int, long int,
                                                 double complex)) {
  hamCollectSink = prev_hook;
  iHamSinkMode = prev_mode;
}

int TraceHamCollect(struct BindStruct *X, size_t cap_bytes,
                    int fail_alloc_at, TraceHamCsr *csr) {
  int prev_mode;
  void (*prev_hook)(long int, long int, double complex);
  long int n, nnz_raw, k_max, w;
  int n_diag, nalloc = 0, di;
  size_t want;
  TraceEnt *ws = NULL;
  double complex *y = NULL;
  double *diag[3] = { NULL, NULL, NULL };
  long int r;

  memset(csr, 0, sizeof(*csr));
  n = (long int)X->Check.idim_max;
  n_diag = trace_model_n_diag(X->Def.iCalcModel);

  prev_mode = iHamSinkMode;
  prev_hook = hamCollectSink;

  /* Alloc 0: rowptr (also the per-row counter during pass 1). */
  g_rowptr = (long int *)galloc(sizeof(long int) * (size_t)(n + 1), &nalloc, fail_alloc_at);
  /* Alloc 1: cursors. */
  g_cursor = (long int *)galloc(sizeof(long int) * (size_t)(n + 1), &nalloc, fail_alloc_at);
  if (g_rowptr == NULL || g_cursor == NULL) {
    free(g_rowptr); g_rowptr = NULL;
    free(g_cursor); g_cursor = NULL;
    return 0; /* sink not yet touched; nothing to restore */
  }

  /* -------- Pass 1: count raw entries per row. -------- */
  g_n = n;
  g_overflow = 0;
  for (r = 0; r <= n; r++) g_rowptr[r] = 0;

  iHamSinkMode = HAM_SINK_TRACE_COLLECT;
  hamCollectSink = trace_count_sink;
  if (makeHam(X) != 0) {
    fprintf(stdoutMPI, "ERROR: ExpecMode 2 Hamiltonian re-enumeration failed\n");
    trace_restore_sink(prev_mode, prev_hook);
    exitMPI(-1);
  }

  if (g_overflow) {
    trace_restore_sink(prev_mode, prev_hook);
    free(g_rowptr); g_rowptr = NULL;
    free(g_cursor); g_cursor = NULL;
    return 0;
  }

  /* k_max = largest raw row; prefix-sum counts into rowptr end-offsets.
     The prefix sum is guarded exactly like trace_count_sink's per-increment
     LONG_MAX check: a running total that would overflow signed long int is a
     demotion (free, restore sink, return 0), never signed-overflow UB. */
  k_max = 0;
  for (r = 1; r <= n; r++) if (g_rowptr[r] > k_max) k_max = g_rowptr[r];
  g_rowptr[0] = 0;
  for (r = 1; r <= n; r++) {
    if (g_rowptr[r] > LONG_MAX - g_rowptr[r - 1]) {
      trace_restore_sink(prev_mode, prev_hook);
      free(g_rowptr); g_rowptr = NULL;
      free(g_cursor); g_cursor = NULL;
      return 0;
    }
    g_rowptr[r] += g_rowptr[r - 1];
  }
  nnz_raw = g_rowptr[n];

  /* -------- Memory gate. -------- */
  want = TraceHamGatedBytes(n, nnz_raw, k_max, n_diag);
  if (want > cap_bytes) {
    trace_restore_sink(prev_mode, prev_hook);
    free(g_rowptr); g_rowptr = NULL;
    free(g_cursor); g_cursor = NULL;
    return 0;
  }

  /* Alloc 2: colidx.  Alloc 3: val.  Alloc 4: sort workspace.  Alloc 5: y.
     Alloc 6..: the n_diag coefficient arrays. */
  g_colidx = (long int *)galloc(sizeof(long int) * (size_t)nnz_raw, &nalloc, fail_alloc_at);
  g_val    = (double complex *)galloc(sizeof(double complex) * (size_t)nnz_raw, &nalloc, fail_alloc_at);
  ws       = (TraceEnt *)galloc(sizeof(TraceEnt) * (size_t)k_max, &nalloc, fail_alloc_at);
  y        = (double complex *)galloc(sizeof(double complex) * (size_t)n, &nalloc, fail_alloc_at);
  for (di = 0; di < n_diag; di++) {
    diag[di] = (double *)galloc(sizeof(double) * (size_t)n, &nalloc, fail_alloc_at);
  }
  {
    int alloc_ok = (g_colidx != NULL && g_val != NULL && ws != NULL && y != NULL);
    for (di = 0; di < n_diag; di++) if (diag[di] == NULL) alloc_ok = 0;
    if (!alloc_ok) {
      trace_restore_sink(prev_mode, prev_hook);
      free(g_rowptr); g_rowptr = NULL;
      free(g_cursor); g_cursor = NULL;
      free(g_colidx); g_colidx = NULL;
      free(g_val);    g_val = NULL;
      free(ws);
      free(y);
      for (di = 0; di < 3; di++) free(diag[di]);
      return 0;
    }
  }

  /* Zero the energy-family coefficient arrays (Task 3 fills them). */
  for (di = 0; di < n_diag; di++) {
    memset(diag[di], 0, sizeof(double) * (size_t)n);
  }

  /* -------- Pass 2: fill row segments. -------- */
  for (r = 1; r <= n; r++) g_cursor[r] = g_rowptr[r - 1];
  hamCollectSink = trace_fill_sink;
  if (makeHam(X) != 0) {
    fprintf(stdoutMPI, "ERROR: ExpecMode 2 Hamiltonian re-enumeration failed\n");
    trace_restore_sink(prev_mode, prev_hook);
    exitMPI(-1);
  }

  /* Sinks no longer needed; restore the caller's dense-mode state. */
  trace_restore_sink(prev_mode, prev_hook);

  /* -------- Per-row stable sort + adjacent-duplicate sum + compaction. --------
     g_cursor now holds the raw END of each row (== raw rowptr). Snapshot the
     raw boundaries there while rowptr is overwritten with merged offsets. */
  for (r = 0; r <= n; r++) g_cursor[r] = g_rowptr[r];
  w = 0;
  for (r = 1; r <= n; r++) {
    long int raw_s = g_cursor[r - 1];
    long int raw_e = g_cursor[r];
    long int len = raw_e - raw_s;
    long int *col = g_colidx + raw_s;
    double complex *v = g_val + raw_s;
    long int p;
    if (len < 16) insertion_seg(col, v, len);
    else          merge_sort_seg(col, v, len, ws);
    p = 0;
    while (p < len) {
      long int c = col[p];
      double complex acc = v[p];
      p++;
      while (p < len && col[p] == c) { acc += v[p]; p++; }
      g_colidx[w] = c;
      g_val[w] = acc;
      w++;
    }
    g_rowptr[r] = w;
  }
  g_rowptr[0] = 0;

  free(ws);
  free(g_cursor); g_cursor = NULL;

  /* Debug invariants. */
  assert(g_rowptr[0] == 0);
#ifndef NDEBUG
  for (r = 1; r <= n; r++) assert(g_rowptr[r] >= g_rowptr[r - 1]);
#endif
  assert(g_rowptr[n] == w);

  csr->n = n;
  csr->nnz = w;
  csr->rowptr = g_rowptr;
  csr->colidx = g_colidx;   /* raw capacity retained */
  csr->val = g_val;
  csr->y = y;
  csr->n_diag = n_diag;
  for (di = 0; di < 3; di++) csr->diag[di] = diag[di];

  g_rowptr = NULL;
  g_colidx = NULL;
  g_val = NULL;
  return 1;
}

void TraceHamFree(TraceHamCsr *csr) {
  int di;
  if (csr == NULL) return;
  free(csr->rowptr);
  free(csr->colidx);
  free(csr->val);
  free(csr->y);
  for (di = 0; di < 3; di++) free(csr->diag[di]);
  memset(csr, 0, sizeof(*csr));
}
