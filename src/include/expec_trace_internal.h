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
 * @file expec_trace_internal.h
 *
 * @brief ExpecMode 2 (trace-kernel) KERNEL-INTERNAL API (phase 3b Task 2).
 *
 * This header is deliberately SEPARATE from the frozen orchestration-boundary
 * header src/include/expec_trace.h (Task 1). It exposes:
 *   - the TraceMap operator-mapping type,
 *   - the basis-mapping extraction drivers (TraceMapExtractOneBody/TwoBody,
 *     TraceMapFree) that live in src/expec_trace.c, and
 *   - the private mapping-probe adapters (`<orig>_TraceProbe`) that are
 *     co-located next to their element functions in
 *     src/mltplyHubbardCore.c / src/mltplySpinCore.c.
 *
 * The unit test (test/unit/expec_trace_map_check.c) includes THIS header to
 * call the kernel internals directly, bypassing the production dispatch (whose
 * capability table stays all-FALSE until Task 5). HPhi installs neither header
 * externally; the split marks the orchestration boundary explicitly.
 *
 * ---- Mapping-probe contract (binding, shared by every `*_TraceProbe`) ----
 *
 * A one-body/two-body Green's-function operator O acts on basis state k
 * (1-based vector index j = k+1) by sending it to a single destination basis
 * state k' with an amplitude a, or annihilating it. Each probe reports, for
 * the given input index j:
 *   return 1  => the transition survives; *kprime_out is the 0-based
 *               destination in [0, idim_max-1] and *amp_out is its (complex)
 *               amplitude a. For a DIAGONAL operator the probe always returns
 *               1 with *kprime_out = j-1 and *amp_out = the bit/occupation
 *               result (which may itself be 0).
 *   return 0  => the transition is annihilated; *kprime_out = -1, *amp_out = 0.
 *
 * Vector-index invariant (the reason a single 0-based kprime suffices for both
 * the canonical GetOffComp path and the grand-canonical bare-bit path): the
 * original element function always reads/writes its result vector at index
 * kprime+1. For the canonical off-diagonal path GetOffComp yields a 1-based
 * index off, and kprime = off-1; for the grand-canonical path the bare bit
 * pattern tmp_off is 0-based and kprime = tmp_off. In every case kprime+1 is
 * the vector slot the pre-refactor code touched, so the extracted `*_map`
 * cores drive both the unchanged original AND the probe from one code path.
 *
 * Coupling-constant convention: probes carry the SAME tmp_V the Mode-1 element
 * call would (one-body: implicit 1.0, so one-body probes take no tmp_V; two-
 * body Hubbard: the fixed 1.0; two-body Spin/SpinGC-half: the +-1 factor that
 * Rearray_Interactions folds in). The two-body probes therefore fold tmp_V
 * into *amp_out. The lone exception, faithfully preserved from Mode 1, is the
 * canonical Spin-half EXCHANGE branch: expec_cisajscktalt_SpinHalf multiplies
 * that branch's contribution by the raw child_exchange_spin_element sign only,
 * NOT by tmp_V, so child_exchange_spin_element_TraceProbe reports the bare sign
 * and its driver branch does not reintroduce tmp_V. See inventory doc 2c.
 */
#pragma once
#include <complex.h>
#include "struct.h"

/**
 * @brief Precomputed basis mapping k -> (k', amplitude) for one GF operator.
 *
 * n == X->Check.idim_max (the process-local Hilbert dimension, canonical or
 * grand-canonical). kprime[k] in [-1, n-1]: 0-based destination for source
 * basis state k (0-based; vector slot k+1), or -1 if O annihilates it. amp[k]
 * is meaningful only where kprime[k] >= 0. is_diagonal is a hint (kprime[k]==k
 * wherever it survives).
 *
 * Sentinel: n == 0 with kprime == amp == NULL means "this operator pair is
 * irregular (Rearray_Interactions rejected it): the Mode-1 path writes a 0.0
 * row for it". A regular operator that merely happens to annihilate every
 * state instead has n == idim_max with every kprime == -1 (also streams to 0).
 */
typedef struct {
  long int n;
  long int *kprime;
  double complex *amp;
  int is_diagonal;
} TraceMap;

/**
 * @brief Build the one-body (cisajs) mapping for operator X->Def.CisAjt[ipair].
 * Mirrors the LOCAL (intra-process) dispatch of expec_cisajs.c with
 * X->Large.mode = M_CORR, then restores X->Large so extraction is side-effect
 * free on X. @return 0 on success, -1 on failure (e.g. allocation).
 */
int TraceMapExtractOneBody(struct BindStruct *X, int ipair, TraceMap *map);

/**
 * @brief Build the two-body (cisajscktaltdc) mapping for operator
 * X->Def.CisAjtCkuAlvDC[ipair]. Mirrors the LOCAL dispatch of
 * expec_cisajscktaltdc.c (incl. Rearray_Interactions for the Spin/SpinGC-half
 * models). @return 0 on success, -1 on failure.
 */
int TraceMapExtractTwoBody(struct BindStruct *X, int ipair, TraceMap *map);

/** @brief Release TraceMap buffers (NULL-safe; leaves the struct zeroed). */
void TraceMapFree(TraceMap *map);

/* ------------------------------------------------------------------ */
/* Mapping-probe adapters (defined next to their element functions).   */
/* One-body.                                                           */
/* ------------------------------------------------------------------ */
int CisAjt_TraceProbe(long unsigned int j, struct BindStruct *X,
    long unsigned int is1_spin, long unsigned int is2_spin,
    long unsigned int sum_spin, long unsigned int diff_spin,
    long int *kprime_out, double complex *amp_out);
int GC_CisAjt_TraceProbe(long unsigned int j, struct BindStruct *X,
    long unsigned int is1_spin, long unsigned int is2_spin,
    long unsigned int sum_spin, long unsigned int diff_spin,
    long int *kprime_out, double complex *amp_out);
int GC_CisAis_TraceProbe(long unsigned int j, struct BindStruct *X,
    long unsigned int is1_spin,
    long int *kprime_out, double complex *amp_out);
int child_Spin_CisAis_TraceProbe(long unsigned int j, struct BindStruct *X,
    long unsigned int is1_spin, long unsigned int sigma1,
    long int *kprime_out, double complex *amp_out);
int child_SpinGC_CisAis_TraceProbe(long unsigned int j, struct BindStruct *X,
    long unsigned int is1_spin, long unsigned int sigma1,
    long int *kprime_out, double complex *amp_out);
int child_SpinGC_CisAit_TraceProbe(long unsigned int j, struct BindStruct *X,
    long unsigned int is1_spin, long unsigned int sigma2,
    long int *kprime_out, double complex *amp_out);

/* Two-body: Hubbard canonical (amp folds in tmp_V). */
int CisAisCisAis_element_TraceProbe(long unsigned int j,
    long unsigned int isite1, long unsigned int isite3,
    double complex tmp_V, struct BindStruct *X,
    long int *kprime_out, double complex *amp_out);
int CisAisCjtAku_element_TraceProbe(long unsigned int j,
    long unsigned int isite1, long unsigned int isite3, long unsigned int isite4,
    long unsigned int Bsum, long unsigned int Bdiff,
    double complex tmp_V, struct BindStruct *X,
    long int *kprime_out, double complex *amp_out);
int CisAjtCkuAku_element_TraceProbe(long unsigned int j,
    long unsigned int isite1, long unsigned int isite2, long unsigned int isite3,
    long unsigned int Asum, long unsigned int Adiff,
    double complex tmp_V, struct BindStruct *X,
    long int *kprime_out, double complex *amp_out);
int CisAjtCkuAlv_element_TraceProbe(long unsigned int j,
    long unsigned int isite1, long unsigned int isite2,
    long unsigned int isite3, long unsigned int isite4,
    long unsigned int Asum, long unsigned int Adiff,
    long unsigned int Bsum, long unsigned int Bdiff,
    double complex tmp_V, struct BindStruct *X,
    long int *kprime_out, double complex *amp_out);

/* Two-body: Hubbard grand-canonical. */
int GC_CisAisCisAis_element_TraceProbe(long unsigned int j,
    long unsigned int isite1, long unsigned int isite3,
    double complex tmp_V, struct BindStruct *X,
    long int *kprime_out, double complex *amp_out);
int GC_CisAisCjtAku_element_TraceProbe(long unsigned int j,
    long unsigned int isite1, long unsigned int isite3, long unsigned int isite4,
    long unsigned int Bsum, long unsigned int Bdiff,
    double complex tmp_V, struct BindStruct *X,
    long int *kprime_out, double complex *amp_out);
int GC_CisAjtCkuAku_element_TraceProbe(long unsigned int j,
    long unsigned int isite1, long unsigned int isite2, long unsigned int isite3,
    long unsigned int Asum, long unsigned int Adiff,
    double complex tmp_V, struct BindStruct *X,
    long int *kprime_out, double complex *amp_out);
int GC_CisAjtCkuAlv_element_TraceProbe(long unsigned int j,
    long unsigned int isite1, long unsigned int isite2,
    long unsigned int isite3, long unsigned int isite4,
    long unsigned int Asum, long unsigned int Adiff,
    long unsigned int Bsum, long unsigned int Bdiff,
    double complex tmp_V, struct BindStruct *X,
    long int *kprime_out, double complex *amp_out);

/* Two-body: canonical Spin-half reachable families. */
int CisAisCisAis_spin_element_TraceProbe(long unsigned int j,
    long unsigned int isA_up, long unsigned int isB_up,
    long unsigned int org_sigma2, long unsigned int org_sigma4,
    double complex tmp_V, struct BindStruct *X,
    long int *kprime_out, double complex *amp_out);
int child_exchange_spin_element_TraceProbe(long unsigned int j,
    struct BindStruct *X,
    long unsigned int isA_up, long unsigned int isB_up,
    long unsigned int sigmaA, long unsigned int sigmaB,
    long int *kprime_out, double complex *amp_out);

/* Two-body: SpinGC-half. */
int GC_CisAisCisAis_spin_element_TraceProbe(long unsigned int j,
    long unsigned int isA_up, long unsigned int isB_up,
    long unsigned int org_sigma2, long unsigned int org_sigma4,
    double complex tmp_V, struct BindStruct *X,
    long int *kprime_out, double complex *amp_out);
int GC_CisAisCitAiu_spin_element_TraceProbe(long unsigned int j,
    long unsigned int org_sigma2, long unsigned int org_sigma4,
    long unsigned int isA_up, long unsigned int isB_up,
    double complex tmp_V, struct BindStruct *X,
    long int *kprime_out, double complex *amp_out);
int GC_CisAitCiuAiu_spin_element_TraceProbe(long unsigned int j,
    long unsigned int org_sigma2, long unsigned int org_sigma4,
    long unsigned int isA_up, long unsigned int isB_up,
    double complex tmp_V, struct BindStruct *X,
    long int *kprime_out, double complex *amp_out);
int GC_CisAitCiuAiv_spin_element_TraceProbe(long unsigned int j,
    long unsigned int org_sigma2, long unsigned int org_sigma4,
    long unsigned int isA_up, long unsigned int isB_up,
    double complex tmp_V, struct BindStruct *X,
    long int *kprime_out, double complex *amp_out);
