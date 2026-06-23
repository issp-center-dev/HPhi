/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 The University of Tokyo */

/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

#ifdef MPI
#include <mpi.h>
#endif
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include "bitcalc.h"
#include "nbody_interall.h"
#include "mltplyCommon.h"
#include "wrapperMPI.h"

static int parse_unsigned_token(const char **pp, unsigned int *value)
{
  char *end = NULL;
  unsigned long v;
  const char *p = *pp;
  while (isspace((unsigned char)*p)) p++;
  if (*p == '\0' || *p == '-' || *p == '+') return -1;
  errno = 0;
  v = strtoul(p, &end, 10);
  if (errno != 0 || end == p || v > UINT_MAX) return -1;
  *value = (unsigned int)v;
  *pp = end;
  return 0;
}

static int parse_int_token(const char **pp, int *value)
{
  char *end = NULL;
  long v;
  const char *p = *pp;
  while (isspace((unsigned char)*p)) p++;
  if (*p == '\0') return -1;
  errno = 0;
  v = strtol(p, &end, 10);
  if (errno != 0 || end == p || v < INT_MIN || v > INT_MAX) return -1;
  *value = (int)v;
  *pp = end;
  return 0;
}

static int parse_double_token(const char **pp, double *value)
{
  char *end = NULL;
  double v;
  const char *p = *pp;
  while (isspace((unsigned char)*p)) p++;
  if (*p == '\0') return -1;
  errno = 0;
  v = strtod(p, &end);
  if (errno != 0 || end == p || !isfinite(v)) return -1;
  *value = v;
  *pp = end;
  return 0;
}

int ParseNBodyInterAllLine(
  const char *line,
  unsigned int *N,
  int **factors,
  double *re,
  double *im
) {
  const char *p = line;
  unsigned int n;
  unsigned int k;
  int *buf;
  size_t nints;

  if (parse_unsigned_token(&p, &n) != 0 || n == 0) {
    fprintf(stdoutMPI, "Error: NBodyInterAll line has an invalid factor count.\n");
    return -1;
  }
  if (n > UINT_MAX / 4) {
    fprintf(stdoutMPI, "Error: NBodyInterAll line is too large.\n");
    return -1;
  }
#if SIZE_MAX < UINT_MAX
  if ((size_t)n > SIZE_MAX / 4 / sizeof(int)) {
    fprintf(stdoutMPI, "Error: NBodyInterAll line is too large.\n");
    return -1;
  }
#endif
  nints = (size_t)4 * n;
  buf = (int *)malloc(nints * sizeof(int));
  if (buf == NULL) {
    fprintf(stdoutMPI, "Error: Failed to allocate NBodyInterAll parser buffer.\n");
    return -1;
  }
  for (k = 0; k < 4 * n; k++) {
    if (parse_int_token(&p, &buf[k]) != 0) {
      fprintf(stdoutMPI, "Error: NBodyInterAll line has too few integer fields.\n");
      free(buf);
      return -1;
    }
  }
  if (parse_double_token(&p, re) != 0 || parse_double_token(&p, im) != 0) {
    fprintf(stdoutMPI, "Error: NBodyInterAll line has invalid coefficient fields.\n");
    free(buf);
    return -1;
  }
  while (isspace((unsigned char)*p)) p++;
  if (*p != '\0') {
    fprintf(stdoutMPI, "Error: NBodyInterAll line has extra fields.\n");
    free(buf);
    return -1;
  }

  *N = n;
  *factors = buf;
  return 0;
}

static int nbody_is_supported_spin_model(const struct DefineList *D)
{
  if (D->iCalcModel == SpinGC) return TRUE;
  if (D->iCalcModel == Spin) return TRUE;
  if (D->iCalcModel == HubbardGC) return TRUE;
  if (D->iCalcModel == Hubbard) return TRUE;
  return FALSE;
}

static int nbody_is_hubbard_model(const struct DefineList *D)
{
  return D->iCalcModel == Hubbard ||
         D->iCalcModel == HubbardGC ||
         D->iCalcModel == HubbardNConserved;
}

static int nbody_is_tj_model(const struct DefineList *D)
{
  return D->iCalcModel == tJ || D->iCalcModel == tJGC;
}

static int nbody_is_kondo_model(const struct DefineList *D)
{
  return D->iCalcModel == Kondo ||
         D->iCalcModel == KondoGC ||
         D->iCalcModel == KondoNConserved;
}

static int nbody_is_spinless_model(const struct DefineList *D)
{
  return D->iCalcModel == SpinlessFermion || D->iCalcModel == SpinlessFermionGC;
}

static int nbody_is_spinful_raw_fermion_model(const struct DefineList *D)
{
  return nbody_is_hubbard_model(D) == TRUE ||
         nbody_is_tj_model(D) == TRUE ||
         nbody_is_kondo_model(D) == TRUE;
}

static int nbody_is_raw_fermion_model(const struct DefineList *D)
{
  return nbody_is_spinful_raw_fermion_model(D) == TRUE ||
         nbody_is_spinless_model(D) == TRUE;
}

static int nbody_is_supported_model(const struct DefineList *D)
{
  if (D->iCalcModel == SpinGC) return TRUE;
  if (D->iCalcModel == Spin) return TRUE;
  if (nbody_is_raw_fermion_model(D) == TRUE) return TRUE;
  return FALSE;
}

static int nbody_uses_hubbard_list_path(const struct DefineList *D)
{
  return D->iCalcModel == Hubbard ||
         D->iCalcModel == HubbardNConserved ||
         D->iCalcModel == tJ ||
         D->iCalcModel == tJGC ||
         D->iCalcModel == Kondo ||
         D->iCalcModel == KondoGC ||
         D->iCalcModel == KondoNConserved;
}

static int nbody_requires_spinful_conservation(const struct DefineList *D)
{
  return D->iCalcModel == Hubbard ||
         D->iCalcModel == tJ ||
         D->iCalcModel == Kondo;
}

static int nbody_is_unsupported_nconserved_model(const struct DefineList *D)
{
  return D->iCalcModel == tJNConserved;
}

static int nbody_is_general_spin(const struct DefineList *D)
{
  return D->iFlgGeneralSpin == TRUE &&
         (D->iCalcModel == Spin || D->iCalcModel == SpinGC);
}

int ValidateNBodyInterAllScope(const struct DefineList *D)
{
  unsigned int t, k;
  if (D->NNBodyInterAll == 0) return 0;
  if (nbody_is_supported_model(D) == FALSE) {
    if (nbody_is_unsupported_nconserved_model(D) == TRUE) {
      fprintf(stdoutMPI,
              "Error: NBodyInterAll does not support tJNConserved. "
              "For tJ standard input, define 2Sz to use the Sz-conserved model.\n");
    }
    else {
      fprintf(stdoutMPI,
              "Error: NBodyInterAll is currently supported only for SpinGC, Spin, "
              "HubbardGC, Hubbard, HubbardNConserved, SpinlessFermionGC, SpinlessFermion, "
              "tJGC, tJ, KondoGC, Kondo, and KondoNConserved.\n");
    }
    return -1;
  }
  if (D->iCalcType == TimeEvolution) {
    fprintf(stdoutMPI, "Error: NBodyInterAll is not yet supported in TimeEvolution.\n");
    return -1;
  }
  if (D->iCalcModel == KondoNConserved && D->iFlgCalcSpec != CALCSPEC_NOT) {
    fprintf(stdoutMPI, "Error: NBodyInterAll does not support KondoNConserved with CalcSpec.\n");
    return -1;
  }
  if (D->iCalcModel == HubbardNConserved && D->iFlgCalcSpec != CALCSPEC_NOT) {
    fprintf(stdoutMPI, "Error: NBodyInterAll does not support HubbardNConserved with CalcSpec.\n");
    return -1;
  }
  if (nbody_is_spinless_model(D) == TRUE && D->iCalcType == FullDiag) {
    fprintf(stdoutMPI, "Error: NBodyInterAll is not yet supported in FullDiag for SpinlessFermion / SpinlessFermionGC.\n");
    return -1;
  }
  if (nbody_is_kondo_model(D) == TRUE) {
    unsigned int site;
    for (site = 0; site < D->Nsite; site++) {
      if (D->LocSpn[site] > LOCSPIN) {
        fprintf(stdoutMPI, "Error: NBodyInterAll does not support general-spin Kondo local spins.\n");
        return -1;
      }
    }
  }
  for (t = 0; t < D->NNBodyInterAll; t++) {
    const unsigned int off = D->NBodyInterAll_Offset[t];
    for (k = 0; k < D->NBodyInterAll_N[t]; k++) {
      const int *f = D->NBodyInterAll_Factors[off + k];
      if (f[0] < 0 || f[0] >= (int)D->Nsite || f[2] < 0 || f[2] >= (int)D->Nsite) {
        fprintf(stdoutMPI, "Error: Site index of NBodyInterAll is incorrect.\n");
        return -1;
      }
      if (nbody_is_raw_fermion_model(D) == FALSE && f[0] != f[2]) {
        fprintf(stdoutMPI, "Error: NBodyInterAll currently requires site_out == site_in for every factor.\n");
        return -1;
      }
      if (nbody_is_kondo_model(D) == TRUE &&
          (D->LocSpn[f[0]] != ITINERANT || D->LocSpn[f[2]] != ITINERANT) &&
          f[0] != f[2]) {
        fprintf(stdoutMPI, "Error: Kondo local-spin NBodyInterAll factors require site_out == site_in.\n");
        return -1;
      }
      if (nbody_is_spinless_model(D) == TRUE) {
        if (f[1] != 0 || f[3] != 0) {
          fprintf(stdoutMPI, "Error: Spin index of NBodyInterAll is incorrect.\n");
          return -1;
        }
      }
      else {
        const int max_spin = (nbody_is_spinful_raw_fermion_model(D) == TRUE) ? 1 :
          (nbody_is_general_spin(D) ? D->LocSpn[f[0]] : 1);
        if (max_spin < 1 || f[1] < 0 || f[1] > max_spin || f[3] < 0 || f[3] > max_spin) {
          fprintf(stdoutMPI, "Error: Spin index of NBodyInterAll is incorrect.\n");
          return -1;
        }
      }
    }
  }
  return 0;
}

static int find_site(const int *sites, unsigned int n, int site)
{
  unsigned int i;
  for (i = 0; i < n; i++) {
    if (sites[i] == site) return (int)i;
  }
  return -1;
}

static void sort_canonical(int *sites, int *outs, int *ins, unsigned int n)
{
  unsigned int i, j;
  for (i = 0; i < n; i++) {
    for (j = i + 1; j < n; j++) {
      if (sites[j] < sites[i]) {
        int ts = sites[i], to = outs[i], ti = ins[i];
        sites[i] = sites[j];
        outs[i] = outs[j];
        ins[i] = ins[j];
        sites[j] = ts;
        outs[j] = to;
        ins[j] = ti;
      }
    }
  }
}

int NormalizeNBodyInterAllTerms(struct DefineList *D)
{
  unsigned int t, k;
  unsigned int total = 0;

  D->NBodyInterAll_TotalCanonicalFactors = 0;
  for (t = 0; t < D->NNBodyInterAll; t++) {
    const unsigned int nraw = D->NBodyInterAll_N[t];
    const unsigned int off = D->NBodyInterAll_Offset[t];
    if (nbody_is_raw_fermion_model(D) == TRUE) {
      D->NBodyInterAll_CanonicalOffset[t] = total;
      D->NBodyInterAll_CanonicalN[t] = nraw;
      for (k = 0; k < nraw; k++) {
        const int *f = D->NBodyInterAll_Factors[off + k];
        D->NBodyInterAll_CanonicalFactors[total + k][0] = f[0];
        D->NBodyInterAll_CanonicalFactors[total + k][1] = f[1];
        D->NBodyInterAll_CanonicalFactors[total + k][2] = f[2];
        D->NBodyInterAll_CanonicalFactors[total + k][3] = f[3];
      }
      total += nraw;
      continue;
    }

    int *sites = (int *)malloc(nraw * sizeof(int));
    int *outs = (int *)malloc(nraw * sizeof(int));
    int *ins = (int *)malloc(nraw * sizeof(int));
    unsigned int ncanon = 0;
    if (sites == NULL || outs == NULL || ins == NULL) {
      fprintf(stdoutMPI, "Error: Failed to allocate NBodyInterAll normalization buffer.\n");
      free(sites);
      free(outs);
      free(ins);
      return -1;
    }

    for (k = 0; k < nraw; k++) {
      const int *f = D->NBodyInterAll_Factors[off + k];
      const int site = f[0];
      const int spin_out = f[1];
      const int spin_in = f[3];
      const int pos = find_site(sites, ncanon, site);
      if (pos < 0) {
        sites[ncanon] = site;
        outs[ncanon] = spin_out;
        ins[ncanon] = spin_in;
        ncanon++;
      }
      else {
        if (ins[pos] != spin_out) {
          fprintf(stdoutMPI, "Error: NBodyInterAll contains a zero same-site operator product.\n");
          free(sites);
          free(outs);
          free(ins);
          return -1;
        }
        ins[pos] = spin_in;
      }
    }

    sort_canonical(sites, outs, ins, ncanon);
    D->NBodyInterAll_CanonicalOffset[t] = total;
    D->NBodyInterAll_CanonicalN[t] = ncanon;
    for (k = 0; k < ncanon; k++) {
      D->NBodyInterAll_CanonicalFactors[total + k][0] = sites[k];
      D->NBodyInterAll_CanonicalFactors[total + k][1] = outs[k];
      D->NBodyInterAll_CanonicalFactors[total + k][2] = sites[k];
      D->NBodyInterAll_CanonicalFactors[total + k][3] = ins[k];
    }
    total += ncanon;
    free(sites);
    free(outs);
    free(ins);
  }
  D->NBodyInterAll_TotalCanonicalFactors = total;
  return 0;
}

int CheckNBodyInterAllSpinConservation(const struct DefineList *D)
{
  unsigned int t, k;
  if (D->NNBodyInterAll == 0) return 0;
  if (D->iCalcModel != Spin) return 0;

  for (t = 0; t < D->NNBodyInterAll; t++) {
    const unsigned int n = D->NBodyInterAll_CanonicalN[t];
    const unsigned int off = D->NBodyInterAll_CanonicalOffset[t];
    int delta_nup = 0;
    for (k = 0; k < n; k++) {
      const int *f = D->NBodyInterAll_CanonicalFactors[off + k];
      delta_nup += f[1] - f[3];
    }
    if (delta_nup != 0) {
      fprintf(stdoutMPI,
              "Error: NBodyInterAll term %u does not conserve total Sz: delta2Sz=%d.\n",
              t + 1, 2 * delta_nup);
      return -1;
    }
  }
  return 0;
}

int CheckNBodyInterAllHubbardConservation(const struct DefineList *D)
{
  unsigned int t, k;
  if (D->NNBodyInterAll == 0) return 0;
  if (nbody_requires_spinful_conservation(D) == FALSE) return 0;

  for (t = 0; t < D->NNBodyInterAll; t++) {
    const unsigned int n = D->NBodyInterAll_CanonicalN[t];
    const unsigned int off = D->NBodyInterAll_CanonicalOffset[t];
    int delta_nup = 0;
    int delta_ndown = 0;
    for (k = 0; k < n; k++) {
      const int *f = D->NBodyInterAll_CanonicalFactors[off + k];
      if (f[1] == 0) delta_nup++;
      else delta_ndown++;
      if (f[3] == 0) delta_nup--;
      else delta_ndown--;
    }
    if (delta_nup != 0 || delta_ndown != 0) {
      fprintf(stdoutMPI,
              "Error: NBodyInterAll term %u does not conserve particle numbers: delta_Nup=%d delta_Ndown=%d.\n",
              t + 1, delta_nup, delta_ndown);
      return -1;
    }
  }
  return 0;
}

int ClassifyNBodyInterAllTerms(struct DefineList *D)
{
  unsigned int t, k;
  D->NNBodyInterAll_Diagonal = 0;
  D->NNBodyInterAll_OffDiagonal = 0;
  for (t = 0; t < D->NNBodyInterAll; t++) {
    const unsigned int n = D->NBodyInterAll_CanonicalN[t];
    const unsigned int off = D->NBodyInterAll_CanonicalOffset[t];
    int diagonal = TRUE;
    for (k = 0; k < n; k++) {
      const int *f = D->NBodyInterAll_CanonicalFactors[off + k];
      if (nbody_is_raw_fermion_model(D) == TRUE) {
        if (f[0] != f[2] || f[1] != f[3]) {
          diagonal = FALSE;
          break;
        }
      }
      else if (f[1] != f[3]) {
        diagonal = FALSE;
        break;
      }
    }
    if (diagonal == TRUE) {
      if (fabs(cimag(D->ParaNBodyInterAll[t])) > eps_CheckImag0) {
        fprintf(stdoutMPI, "Error: Diagonal NBodyInterAll term has a finite imaginary part.\n");
        return -1;
      }
      D->NBodyInterAll_DiagonalIndex[D->NNBodyInterAll_Diagonal++] = t;
    }
    else {
      D->NBodyInterAll_OffDiagonalIndex[D->NNBodyInterAll_OffDiagonal++] = t;
    }
  }
  return 0;
}

int CheckNBodyInterAllHermitePairs(const struct DefineList *D)
{
  unsigned int p, k;
  if (D->NNBodyInterAll_OffDiagonal % 2 != 0) {
    fprintf(stdoutMPI, "Error: Off-diagonal NBodyInterAll terms must appear as adjacent Hermite pairs.\n");
    return -1;
  }
  for (p = 0; p < D->NNBodyInterAll_OffDiagonal; p += 2) {
    const unsigned int t0 = D->NBodyInterAll_OffDiagonalIndex[p];
    const unsigned int t1 = D->NBodyInterAll_OffDiagonalIndex[p + 1];
    const unsigned int n0 = D->NBodyInterAll_CanonicalN[t0];
    const unsigned int n1 = D->NBodyInterAll_CanonicalN[t1];
    const unsigned int off0 = D->NBodyInterAll_CanonicalOffset[t0];
    const unsigned int off1 = D->NBodyInterAll_CanonicalOffset[t1];
    if (t1 != t0 + 1 || n0 != n1) {
      fprintf(stdoutMPI, "Error: Off-diagonal NBodyInterAll terms must appear as adjacent Hermite pairs.\n");
      return -1;
    }
    for (k = 0; k < n0; k++) {
      const int *f0 = D->NBodyInterAll_CanonicalFactors[off0 + k];
      const int *f1 = (nbody_is_raw_fermion_model(D) == TRUE) ?
        D->NBodyInterAll_CanonicalFactors[off1 + n0 - 1 - k] :
        D->NBodyInterAll_CanonicalFactors[off1 + k];
      if (nbody_is_raw_fermion_model(D) == TRUE) {
        if (f0[0] != f1[2] || f0[1] != f1[3] || f0[2] != f1[0] || f0[3] != f1[1]) {
          fprintf(stdoutMPI, "Error: Off-diagonal NBodyInterAll Hermite pair has inconsistent factors.\n");
          return -1;
        }
      }
      else {
        if (f0[0] != f1[0] || f0[1] != f1[3] || f0[3] != f1[1]) {
          fprintf(stdoutMPI, "Error: Off-diagonal NBodyInterAll Hermite pair has inconsistent factors.\n");
          return -1;
        }
      }
    }
    if (cabs(D->ParaNBodyInterAll[t1] - conj(D->ParaNBodyInterAll[t0])) > eps_CheckImag0) {
      fprintf(stdoutMPI, "Error: Off-diagonal NBodyInterAll Hermite pair has inconsistent coefficients.\n");
      return -1;
    }
  }
  return 0;
}

static int apply_nbody_interall_bits(
  const struct BindStruct *X,
  unsigned int term_index,
  unsigned long int intra_in,
  int rank_in,
  unsigned long int *intra_out,
  int *rank_out
) {
  const struct DefineList *D = &X->Def;
  unsigned int k;
  unsigned long int lo = intra_in;
  int ro = rank_in;
  const unsigned int n = D->NBodyInterAll_CanonicalN[term_index];
  const unsigned int off = D->NBodyInterAll_CanonicalOffset[term_index];

  for (k = 0; k < n; k++) {
    const int *f = D->NBodyInterAll_CanonicalFactors[off + k];
    const unsigned int site = (unsigned int)f[0];
    const int spin_out = f[1];
    const int spin_in = f[3];
    const unsigned long int mask = D->Tpow[site];
    int bit;

    if (site < D->Nsite) {
      bit = (lo & mask) ? 1 : 0;
      if (bit != spin_in) return 0;
      if (spin_out != spin_in) {
        if (spin_out == 1) lo |= mask;
        else lo &= ~mask;
      }
    }
    else {
      bit = (ro & (int)mask) ? 1 : 0;
      if (bit != spin_in) return 0;
      if (spin_out != spin_in) {
        if (spin_out == 1) ro |= (int)mask;
        else ro &= ~((int)mask);
      }
    }
  }

  *intra_out = lo;
  *rank_out = ro;
  return 1;
}

static int convert_nbody_general_spin_to_list1(
  const struct BindStruct *X,
  unsigned long int local_out,
  unsigned long int *j_out
) {
  return ConvertToList1GeneralSpin(local_out, X->Check.sdim, j_out);
}

static int apply_nbody_interall_general_spin_gc(
  const struct BindStruct *X,
  unsigned int term_index,
  unsigned long int local_in,
  int rank_in,
  unsigned long int *local_out,
  int *rank_out
) {
  const struct DefineList *D = &X->Def;
  unsigned int k;
  unsigned long int lo = local_in;
  unsigned long int ro = (unsigned long int)rank_in;
  const unsigned int n = D->NBodyInterAll_CanonicalN[term_index];
  const unsigned int off = D->NBodyInterAll_CanonicalOffset[term_index];

  for (k = 0; k < n; k++) {
    const int *f = D->NBodyInterAll_CanonicalFactors[off + k];
    const unsigned int site = (unsigned int)f[0];
    const int spin_out = f[1];
    const int spin_in = f[3];
    unsigned long int next = 0;

    if (site < D->Nsite) {
      if (GetOffCompGeneralSpin(lo, (int)site + 1, spin_in, spin_out,
                                &next, D->SiteToBit, D->Tpow) == FALSE) {
        return 0;
      }
      lo = next;
    }
    else {
      if (GetOffCompGeneralSpin(ro, (int)site + 1, spin_in, spin_out,
                                &next, D->SiteToBit, D->Tpow) == FALSE) {
        return 0;
      }
      ro = next;
    }
  }

  *local_out = lo;
  *rank_out = (int)ro;
  return 1;
}

static int apply_hubbardgc_annihilate_mask(
  unsigned long int mask,
  unsigned long int *state,
  int *sign
) {
  int sgn = 1;
  if ((*state & mask) == 0) return 0;
  SgnBit(*state & (mask - 1), &sgn);
  *sign *= sgn;
  *state &= ~mask;
  return 1;
}

static int apply_hubbardgc_create_mask(
  unsigned long int mask,
  unsigned long int *state,
  int *sign
) {
  int sgn = 1;
  if ((*state & mask) != 0) return 0;
  SgnBit(*state & (mask - 1), &sgn);
  *sign *= sgn;
  *state |= mask;
  return 1;
}

static unsigned long int get_hubbardgc_local_block(
  const struct DefineList *D
) {
  if (D->Nsite == 0) return 1;
  return D->OrgTpow[2 * D->Nsite - 1] * 2;
}

static int apply_nbody_interall_hubbardgc_full(
  const struct BindStruct *X,
  unsigned int term_index,
  unsigned long int local_in,
  int rank_in,
  unsigned long int *local_out,
  int *rank_out,
  int *sign
) {
  const struct DefineList *D = &X->Def;
  unsigned int k;
  unsigned long int state;
  const unsigned long int block = get_hubbardgc_local_block(D);
  const unsigned int n = D->NBodyInterAll_CanonicalN[term_index];
  const unsigned int off = D->NBodyInterAll_CanonicalOffset[term_index];

  state = local_in + block * (unsigned long int)rank_in;
  *sign = 1;

  for (k = n; k > 0; k--) {
    const int *f = D->NBodyInterAll_CanonicalFactors[off + k - 1];
    const unsigned int site_out = (unsigned int)f[0];
    const unsigned int spin_out = (unsigned int)f[1];
    const unsigned int site_in = (unsigned int)f[2];
    const unsigned int spin_in = (unsigned int)f[3];
    const unsigned long int mask_in = D->OrgTpow[2 * site_in + spin_in];
    const unsigned long int mask_out = D->OrgTpow[2 * site_out + spin_out];

    if (apply_hubbardgc_annihilate_mask(mask_in, &state, sign) == 0) return 0;
    if (apply_hubbardgc_create_mask(mask_out, &state, sign) == 0) return 0;
  }

  *local_out = state % block;
  *rank_out = (int)(state / block);
  return 1;
}

static int apply_hubbardgc_rank_annihilate(
  const struct DefineList *D,
  unsigned int site,
  unsigned int spin,
  unsigned long int *rank_state
) {
  unsigned long int mask;
  if (site < D->Nsite) return 1;
  mask = D->Tpow[2 * site + spin];
  if ((*rank_state & mask) == 0) return 0;
  *rank_state &= ~mask;
  return 1;
}

static int apply_hubbardgc_rank_create(
  const struct DefineList *D,
  unsigned int site,
  unsigned int spin,
  unsigned long int *rank_state
) {
  unsigned long int mask;
  if (site < D->Nsite) return 1;
  mask = D->Tpow[2 * site + spin];
  if ((*rank_state & mask) != 0) return 0;
  *rank_state |= mask;
  return 1;
}

static int apply_hubbardgc_rank_factor(
  const struct DefineList *D,
  const int *f,
  int dagger,
  unsigned long int *rank_state
) {
  const unsigned int site_out = (unsigned int)f[0];
  const unsigned int spin_out = (unsigned int)f[1];
  const unsigned int site_in = (unsigned int)f[2];
  const unsigned int spin_in = (unsigned int)f[3];

  if (dagger == FALSE) {
    if (apply_hubbardgc_rank_annihilate(D, site_in, spin_in, rank_state) == 0) return 0;
    if (apply_hubbardgc_rank_create(D, site_out, spin_out, rank_state) == 0) return 0;
  }
  else {
    if (apply_hubbardgc_rank_annihilate(D, site_out, spin_out, rank_state) == 0) return 0;
    if (apply_hubbardgc_rank_create(D, site_in, spin_in, rank_state) == 0) return 0;
  }
  return 1;
}

static int apply_hubbardgc_rank_part(
  const struct BindStruct *X,
  unsigned int term,
  int current_rank,
  int dagger,
  int *partner_rank
) {
  const struct DefineList *D = &X->Def;
  const unsigned int n = D->NBodyInterAll_CanonicalN[term];
  const unsigned int off = D->NBodyInterAll_CanonicalOffset[term];
  unsigned int k;
  unsigned long int rank_state = (unsigned long int)current_rank;

  if (dagger == FALSE) {
    for (k = n; k > 0; k--) {
      if (apply_hubbardgc_rank_factor(D, D->NBodyInterAll_CanonicalFactors[off + k - 1],
                                      FALSE, &rank_state) == 0) {
        return 0;
      }
    }
  }
  else {
    for (k = 0; k < n; k++) {
      if (apply_hubbardgc_rank_factor(D, D->NBodyInterAll_CanonicalFactors[off + k],
                                      TRUE, &rank_state) == 0) {
        return 0;
      }
    }
  }

  *partner_rank = (int)rank_state;
  return 1;
}

static int nbody_interall_hubbardgc_partner_rank(
  const struct BindStruct *X,
  unsigned int term,
  int current_rank,
  int *partner_rank,
  int *active
) {
  if (apply_hubbardgc_rank_part(X, term, current_rank, FALSE, partner_rank) == 1) {
    *active = TRUE;
    return 0;
  }
  if (apply_hubbardgc_rank_part(X, term, current_rank, TRUE, partner_rank) == 1) {
    *active = TRUE;
    return 0;
  }
  *active = FALSE;
  *partner_rank = current_rank;
  return 0;
}

static unsigned long int get_spinless_local_block(
  const struct DefineList *D
) {
  if (D->Nsite == 0) return 1;
  return 1UL << D->Nsite;
}

static int apply_nbody_interall_spinless_full(
  const struct BindStruct *X,
  unsigned int term_index,
  unsigned long int local_in,
  int rank_in,
  unsigned long int *local_out,
  int *rank_out,
  int *sign
) {
  const struct DefineList *D = &X->Def;
  unsigned int k;
  unsigned long int state;
  const unsigned long int block = get_spinless_local_block(D);
  const unsigned int n = D->NBodyInterAll_CanonicalN[term_index];
  const unsigned int off = D->NBodyInterAll_CanonicalOffset[term_index];

  state = local_in + block * (unsigned long int)rank_in;
  *sign = 1;

  for (k = n; k > 0; k--) {
    const int *f = D->NBodyInterAll_CanonicalFactors[off + k - 1];
    const unsigned int site_out = (unsigned int)f[0];
    const unsigned int site_in = (unsigned int)f[2];
    const unsigned long int mask_in = 1UL << site_in;
    const unsigned long int mask_out = 1UL << site_out;

    if (apply_hubbardgc_annihilate_mask(mask_in, &state, sign) == 0) return 0;
    if (apply_hubbardgc_create_mask(mask_out, &state, sign) == 0) return 0;
  }

  *local_out = state % block;
  *rank_out = (int)(state / block);
  return 1;
}

static int apply_spinless_rank_annihilate(
  const struct DefineList *D,
  unsigned int site,
  unsigned long int *rank_state
) {
  unsigned long int mask;
  if (site < D->Nsite) return 1;
  mask = D->Tpow[site];
  if ((*rank_state & mask) == 0) return 0;
  *rank_state &= ~mask;
  return 1;
}

static int apply_spinless_rank_create(
  const struct DefineList *D,
  unsigned int site,
  unsigned long int *rank_state
) {
  unsigned long int mask;
  if (site < D->Nsite) return 1;
  mask = D->Tpow[site];
  if ((*rank_state & mask) != 0) return 0;
  *rank_state |= mask;
  return 1;
}

static int apply_spinless_rank_factor(
  const struct DefineList *D,
  const int *f,
  int dagger,
  unsigned long int *rank_state
) {
  const unsigned int site_out = (unsigned int)f[0];
  const unsigned int site_in = (unsigned int)f[2];

  if (dagger == FALSE) {
    if (apply_spinless_rank_annihilate(D, site_in, rank_state) == 0) return 0;
    if (apply_spinless_rank_create(D, site_out, rank_state) == 0) return 0;
  }
  else {
    if (apply_spinless_rank_annihilate(D, site_out, rank_state) == 0) return 0;
    if (apply_spinless_rank_create(D, site_in, rank_state) == 0) return 0;
  }
  return 1;
}

static int apply_spinless_rank_part(
  const struct BindStruct *X,
  unsigned int term,
  int current_rank,
  int dagger,
  int *partner_rank
) {
  const struct DefineList *D = &X->Def;
  const unsigned int n = D->NBodyInterAll_CanonicalN[term];
  const unsigned int off = D->NBodyInterAll_CanonicalOffset[term];
  unsigned int k;
  unsigned long int rank_state = (unsigned long int)current_rank;

  if (dagger == FALSE) {
    for (k = n; k > 0; k--) {
      if (apply_spinless_rank_factor(D, D->NBodyInterAll_CanonicalFactors[off + k - 1],
                                     FALSE, &rank_state) == 0) {
        return 0;
      }
    }
  }
  else {
    for (k = 0; k < n; k++) {
      if (apply_spinless_rank_factor(D, D->NBodyInterAll_CanonicalFactors[off + k],
                                     TRUE, &rank_state) == 0) {
        return 0;
      }
    }
  }

  *partner_rank = (int)rank_state;
  return 1;
}

static int nbody_interall_spinless_partner_rank(
  const struct BindStruct *X,
  unsigned int term,
  int current_rank,
  int *partner_rank,
  int *active
) {
  if (apply_spinless_rank_part(X, term, current_rank, FALSE, partner_rank) == 1) {
    *active = TRUE;
    return 0;
  }
  if (apply_spinless_rank_part(X, term, current_rank, TRUE, partner_rank) == 1) {
    *active = TRUE;
    return 0;
  }
  *active = FALSE;
  *partner_rank = current_rank;
  return 0;
}

int ApplyNBodyInterAllSpinGC(
  const struct BindStruct *X,
  unsigned int term_index,
  unsigned long int local_in,
  int rank_in,
  unsigned long int *local_out,
  int *rank_out,
  double complex *matrix_element
) {
  int ret;
  if (nbody_is_general_spin(&X->Def) == TRUE) {
    ret = apply_nbody_interall_general_spin_gc(X, term_index, local_in, rank_in, local_out, rank_out);
  }
  else {
    ret = apply_nbody_interall_bits(X, term_index, local_in, rank_in, local_out, rank_out);
  }
  if (ret == 1) *matrix_element = X->Def.ParaNBodyInterAll[term_index];
  return ret;
}

int SetDiagonalNBodyInterAllSpinGC(struct BindStruct *X)
{
  unsigned int i;
  if (X->Def.NNBodyInterAll_Diagonal == 0) return 0;
  if (nbody_is_supported_spin_model(&X->Def) == FALSE) return -1;
  if (nbody_is_hubbard_model(&X->Def) == TRUE) return -1;

  for (i = 0; i < X->Def.NNBodyInterAll_Diagonal; i++) {
    const unsigned int term = X->Def.NBodyInterAll_DiagonalIndex[i];
    const double coeff = creal(X->Def.ParaNBodyInterAll[term]);
    const unsigned long int i_max = X->Check.idim_max;
    unsigned long int j;
    if (X->Def.iCalcModel == Spin) {
#pragma omp parallel for default(none) shared(list_Diagonal, list_1, X) firstprivate(i_max, term, coeff, myrank) private(j)
      for (j = 1; j <= i_max; j++) {
        unsigned long int local_out = 0;
        int rank_out = 0;
        int ret = (X->Def.iFlgGeneralSpin == TRUE) ?
          apply_nbody_interall_general_spin_gc(X, term, list_1[j], myrank, &local_out, &rank_out) :
          apply_nbody_interall_bits(X, term, list_1[j], myrank, &local_out, &rank_out);
        if (ret == 1 && local_out == list_1[j] && rank_out == myrank) {
          list_Diagonal[j] += coeff;
        }
      }
      continue;
    }
#pragma omp parallel for default(none) shared(list_Diagonal, X) firstprivate(i_max, term, coeff, myrank) private(j)
    for (j = 1; j <= i_max; j++) {
      unsigned long int local_out = 0;
      int rank_out = 0;
      double complex me = 0.0;
      int ret = ApplyNBodyInterAllSpinGC(X, term, j - 1, myrank, &local_out, &rank_out, &me);
      if (ret == 1 && local_out == j - 1 && rank_out == myrank) {
        list_Diagonal[j] += coeff;
      }
    }
  }
  return 0;
}

int SetDiagonalNBodyInterAllHubbardGC(struct BindStruct *X)
{
  unsigned int i;
  if (X->Def.NNBodyInterAll_Diagonal == 0) return 0;
  if (X->Def.iCalcModel != HubbardGC) return -1;

  for (i = 0; i < X->Def.NNBodyInterAll_Diagonal; i++) {
    const unsigned int term = X->Def.NBodyInterAll_DiagonalIndex[i];
    const double coeff = creal(X->Def.ParaNBodyInterAll[term]);
    const unsigned long int i_max = X->Check.idim_max;
    unsigned long int j;

#pragma omp parallel for default(none) shared(list_Diagonal, X) firstprivate(i_max, term, coeff, myrank) private(j)
    for (j = 1; j <= i_max; j++) {
      unsigned long int local_out = 0;
      int rank_out = 0;
      int sign = 1;
      int ret = apply_nbody_interall_hubbardgc_full(
        X, term, j - 1, myrank, &local_out, &rank_out, &sign);
      if (ret == 1 && local_out == j - 1 && rank_out == myrank) {
        list_Diagonal[j] += coeff * sign;
      }
    }
  }
  return 0;
}

int SetDiagonalNBodyInterAllHubbard(struct BindStruct *X)
{
  unsigned int i;
  if (X->Def.NNBodyInterAll_Diagonal == 0) return 0;
  if (nbody_uses_hubbard_list_path(&X->Def) == FALSE) return -1;

  for (i = 0; i < X->Def.NNBodyInterAll_Diagonal; i++) {
    const unsigned int term = X->Def.NBodyInterAll_DiagonalIndex[i];
    const double coeff = creal(X->Def.ParaNBodyInterAll[term]);
    const unsigned long int i_max = X->Check.idim_max;
    unsigned long int j;

#pragma omp parallel for default(none) shared(list_Diagonal, list_1, X) firstprivate(i_max, term, coeff, myrank) private(j)
    for (j = 1; j <= i_max; j++) {
      unsigned long int local_out = 0;
      int rank_out = 0;
      int sign = 1;
      int ret = apply_nbody_interall_hubbardgc_full(
        X, term, list_1[j], myrank, &local_out, &rank_out, &sign);
      if (ret == 1 && local_out == list_1[j] && rank_out == myrank) {
        list_Diagonal[j] += coeff * sign;
      }
    }
  }
  return 0;
}

int SetDiagonalNBodyInterAllSpinlessGC(struct BindStruct *X)
{
  unsigned int i;
  if (X->Def.NNBodyInterAll_Diagonal == 0) return 0;
  if (X->Def.iCalcModel != SpinlessFermionGC) return -1;

  for (i = 0; i < X->Def.NNBodyInterAll_Diagonal; i++) {
    const unsigned int term = X->Def.NBodyInterAll_DiagonalIndex[i];
    const double coeff = creal(X->Def.ParaNBodyInterAll[term]);
    const unsigned long int i_max = X->Check.idim_max;
    unsigned long int j;

#pragma omp parallel for default(none) shared(list_Diagonal, X) firstprivate(i_max, term, coeff, myrank) private(j)
    for (j = 1; j <= i_max; j++) {
      unsigned long int local_out = 0;
      int rank_out = 0;
      int sign = 1;
      int ret = apply_nbody_interall_spinless_full(
        X, term, j - 1, myrank, &local_out, &rank_out, &sign);
      if (ret == 1 && local_out == j - 1 && rank_out == myrank) {
        list_Diagonal[j] += coeff * sign;
      }
    }
  }
  return 0;
}

int SetDiagonalNBodyInterAllSpinless(struct BindStruct *X)
{
  unsigned int i;
  if (X->Def.NNBodyInterAll_Diagonal == 0) return 0;
  if (X->Def.iCalcModel != SpinlessFermion) return -1;

  for (i = 0; i < X->Def.NNBodyInterAll_Diagonal; i++) {
    const unsigned int term = X->Def.NBodyInterAll_DiagonalIndex[i];
    const double coeff = creal(X->Def.ParaNBodyInterAll[term]);
    const unsigned long int i_max = X->Check.idim_max;
    unsigned long int j;

#pragma omp parallel for default(none) shared(list_Diagonal, list_1, X) firstprivate(i_max, term, coeff, myrank) private(j)
    for (j = 1; j <= i_max; j++) {
      unsigned long int local_out = 0;
      int rank_out = 0;
      int sign = 1;
      int ret = apply_nbody_interall_spinless_full(
        X, term, list_1[j], myrank, &local_out, &rank_out, &sign);
      if (ret == 1 && local_out == list_1[j] && rank_out == myrank) {
        list_Diagonal[j] += coeff * sign;
      }
    }
  }
  return 0;
}

static int nbody_rank_flip_mask(const struct BindStruct *X, unsigned int term, int *mask)
{
  unsigned int k;
  int m = 0;
  const unsigned int n = X->Def.NBodyInterAll_CanonicalN[term];
  const unsigned int off = X->Def.NBodyInterAll_CanonicalOffset[term];
  for (k = 0; k < n; k++) {
    const int *f = X->Def.NBodyInterAll_CanonicalFactors[off + k];
    const unsigned int site = (unsigned int)f[0];
    if (site >= X->Def.Nsite && f[1] != f[3]) {
      m ^= (int)X->Def.Tpow[site];
    }
  }
  *mask = m;
  return 0;
}

static int nbody_interall_general_spin_partner_rank(
  const struct BindStruct *X,
  unsigned int term,
  int current_rank,
  int *partner_rank,
  int *active
) {
  unsigned int k;
  unsigned long int partner = (unsigned long int)current_rank;
  int side = 0;
  const unsigned int n = X->Def.NBodyInterAll_CanonicalN[term];
  const unsigned int off = X->Def.NBodyInterAll_CanonicalOffset[term];

  *active = TRUE;
  for (k = 0; k < n; k++) {
    const int *f = X->Def.NBodyInterAll_CanonicalFactors[off + k];
    const unsigned int site = (unsigned int)f[0];
    const int spin_out = f[1];
    const int spin_in = f[3];
    int digit;
    int this_side;

    if (site < X->Def.Nsite) continue;

    digit = GetBitGeneral(site + 1, (unsigned long int)current_rank,
                          X->Def.SiteToBit, X->Def.Tpow);
    if (spin_out == spin_in) {
      if (digit != spin_out) {
        *active = FALSE;
        *partner_rank = current_rank;
        return 0;
      }
      continue;
    }

    if (digit == spin_out) this_side = 1;
    else if (digit == spin_in) this_side = -1;
    else {
      *active = FALSE;
      *partner_rank = current_rank;
      return 0;
    }

    if (side == 0) side = this_side;
    else if (side != this_side) {
      *active = FALSE;
      *partner_rank = current_rank;
      return 0;
    }

    if (this_side == 1) {
      partner += ((long int)spin_in - (long int)spin_out) * X->Def.Tpow[site];
    }
    else {
      partner += ((long int)spin_out - (long int)spin_in) * X->Def.Tpow[site];
    }
  }

  *partner_rank = (int)partner;
  return 0;
}

static double complex apply_nbody_term_to_rank(
  struct BindStruct *X,
  unsigned int term,
  double complex *tmp_v0,
  const double complex *src_v1,
  double complex *cur_v1,
  int rank_in
) {
  unsigned long int j;
  double complex dam_pr = 0.0;
  const unsigned long int i_max = X->Check.idim_max;
  const int do_update = (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC);

#pragma omp parallel for default(none) reduction(+:dam_pr) \
  shared(X, tmp_v0, src_v1, cur_v1) firstprivate(i_max, term, rank_in, do_update, myrank) private(j)
  for (j = 1; j <= i_max; j++) {
    unsigned long int local_out = 0;
    int rank_out = 0;
    double complex me = 0.0;
    int ret = ApplyNBodyInterAllSpinGC(X, term, j - 1, rank_in, &local_out, &rank_out, &me);
    if (ret == 1 && rank_out == myrank) {
      const double complex dmv = me * src_v1[j];
      if (do_update) tmp_v0[local_out + 1] += dmv;
      dam_pr += conj(cur_v1[local_out + 1]) * dmv;
    }
  }
  return dam_pr;
}

static double complex apply_nbody_term_to_rank_spin(
  struct BindStruct *X,
  unsigned int term,
  double complex *tmp_v0,
  const unsigned long int *src_list_1,
  const double complex *src_v1,
  double complex *cur_v1,
  unsigned long int src_i_max,
  int rank_in
) {
  unsigned long int j;
  double complex dam_pr = 0.0;
  const int do_update = (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC);

#pragma omp parallel for default(none) reduction(+:dam_pr) \
  shared(X, tmp_v0, src_list_1, src_v1, cur_v1, list_2_1, list_2_2) \
  firstprivate(src_i_max, term, rank_in, do_update, myrank) private(j)
  for (j = 1; j <= src_i_max; j++) {
    unsigned long int intra_out = 0;
    unsigned long int j_out = 0;
    int rank_out = 0;
    double complex me = X->Def.ParaNBodyInterAll[term];
    int ret = (X->Def.iFlgGeneralSpin == TRUE) ?
      apply_nbody_interall_general_spin_gc(X, term, src_list_1[j], rank_in, &intra_out, &rank_out) :
      apply_nbody_interall_bits(X, term, src_list_1[j], rank_in, &intra_out, &rank_out);
    int in_sector = FALSE;
    if (ret == 1 && rank_out == myrank) {
      in_sector = (X->Def.iFlgGeneralSpin == TRUE) ?
        convert_nbody_general_spin_to_list1(X, intra_out, &j_out) :
        GetOffComp(list_2_1, list_2_2, intra_out,
                   X->Large.irght, X->Large.ilft, X->Large.ihfbit, &j_out);
    }
    if (in_sector == TRUE) {
      const double complex dmv = me * src_v1[j];
      if (do_update) tmp_v0[j_out] += dmv;
      dam_pr += conj(cur_v1[j_out]) * dmv;
    }
  }
  return dam_pr;
}

static double complex multiply_nbody_pair_spin_general_spin(
  struct BindStruct *X,
  unsigned int offdiag_pair_pos,
  double complex *tmp_v0,
  double complex *tmp_v1
) {
  const unsigned int term0 = X->Def.NBodyInterAll_OffDiagonalIndex[offdiag_pair_pos];
  const unsigned int term1 = X->Def.NBodyInterAll_OffDiagonalIndex[offdiag_pair_pos + 1];
  int origin = myrank;
  int active = TRUE;
  double complex dam_pr = 0.0;

  if (nbody_interall_general_spin_partner_rank(X, term0, myrank, &origin, &active) != 0) {
    return 0.0;
  }
  if (active == FALSE || origin == myrank) {
    dam_pr += apply_nbody_term_to_rank_spin(
      X, term0, tmp_v0, list_1, tmp_v1, tmp_v1, X->Check.idim_max, myrank);
    dam_pr += apply_nbody_term_to_rank_spin(
      X, term1, tmp_v0, list_1, tmp_v1, tmp_v1, X->Check.idim_max, myrank);
    return dam_pr;
  }

#ifdef MPI
  {
    MPI_Status statusMPI;
    unsigned long int idim_max_buf = 0;
    int ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                            &idim_max_buf,      1, MPI_UNSIGNED_LONG, origin, 0,
                            MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
                        list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf,  idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    dam_pr += apply_nbody_term_to_rank_spin(
      X, term0, tmp_v0, list_1buf, v1buf, tmp_v1, idim_max_buf, origin);
    dam_pr += apply_nbody_term_to_rank_spin(
      X, term1, tmp_v0, list_1buf, v1buf, tmp_v1, idim_max_buf, origin);
  }
#else
  fprintf(stdoutMPI, "Error: NBodyInterAll reached an MPI-only rank flip path without MPI.\n");
  return 0.0;
#endif
  return dam_pr;
}

static double complex multiply_nbody_pair(
  struct BindStruct *X,
  unsigned int offdiag_pair_pos,
  double complex *tmp_v0,
  double complex *tmp_v1
) {
  const unsigned int term0 = X->Def.NBodyInterAll_OffDiagonalIndex[offdiag_pair_pos];
  const unsigned int term1 = X->Def.NBodyInterAll_OffDiagonalIndex[offdiag_pair_pos + 1];
  int mask = 0;
  int origin;
  double complex dam_pr = 0.0;

  nbody_rank_flip_mask(X, term0, &mask);
  origin = myrank ^ mask;
  if (origin == myrank) {
    dam_pr += apply_nbody_term_to_rank(X, term0, tmp_v0, tmp_v1, tmp_v1, myrank);
    dam_pr += apply_nbody_term_to_rank(X, term1, tmp_v0, tmp_v1, tmp_v1, myrank);
    return dam_pr;
  }

#ifdef MPI
  {
    MPI_Status statusMPI;
    int ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                            v1buf,  X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                            MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    dam_pr += apply_nbody_term_to_rank(X, term0, tmp_v0, v1buf, tmp_v1, origin);
    dam_pr += apply_nbody_term_to_rank(X, term1, tmp_v0, v1buf, tmp_v1, origin);
  }
#else
  fprintf(stdoutMPI, "Error: NBodyInterAll reached an MPI-only rank flip path without MPI.\n");
  return 0.0;
#endif
  return dam_pr;
}

static double complex multiply_nbody_pair_general_spin_gc(
  struct BindStruct *X,
  unsigned int offdiag_pair_pos,
  double complex *tmp_v0,
  double complex *tmp_v1
) {
  const unsigned int term0 = X->Def.NBodyInterAll_OffDiagonalIndex[offdiag_pair_pos];
  const unsigned int term1 = X->Def.NBodyInterAll_OffDiagonalIndex[offdiag_pair_pos + 1];
  int origin = myrank;
  int active = TRUE;
  double complex dam_pr = 0.0;

  if (nbody_interall_general_spin_partner_rank(X, term0, myrank, &origin, &active) != 0) {
    return 0.0;
  }
  if (active == FALSE || origin == myrank) {
    dam_pr += apply_nbody_term_to_rank(X, term0, tmp_v0, tmp_v1, tmp_v1, myrank);
    dam_pr += apply_nbody_term_to_rank(X, term1, tmp_v0, tmp_v1, tmp_v1, myrank);
    return dam_pr;
  }

#ifdef MPI
  {
    MPI_Status statusMPI;
    int ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                            v1buf,  X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                            MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    dam_pr += apply_nbody_term_to_rank(X, term0, tmp_v0, v1buf, tmp_v1, origin);
    dam_pr += apply_nbody_term_to_rank(X, term1, tmp_v0, v1buf, tmp_v1, origin);
  }
#else
  fprintf(stdoutMPI, "Error: NBodyInterAll reached an MPI-only rank flip path without MPI.\n");
  return 0.0;
#endif
  return dam_pr;
}

static double complex multiply_nbody_pair_spin(
  struct BindStruct *X,
  unsigned int offdiag_pair_pos,
  double complex *tmp_v0,
  double complex *tmp_v1
) {
  const unsigned int term0 = X->Def.NBodyInterAll_OffDiagonalIndex[offdiag_pair_pos];
  const unsigned int term1 = X->Def.NBodyInterAll_OffDiagonalIndex[offdiag_pair_pos + 1];
  int mask = 0;
  int origin;
  double complex dam_pr = 0.0;

  nbody_rank_flip_mask(X, term0, &mask);
  origin = myrank ^ mask;
  if (origin == myrank) {
    dam_pr += apply_nbody_term_to_rank_spin(
      X, term0, tmp_v0, list_1, tmp_v1, tmp_v1, X->Check.idim_max, myrank);
    dam_pr += apply_nbody_term_to_rank_spin(
      X, term1, tmp_v0, list_1, tmp_v1, tmp_v1, X->Check.idim_max, myrank);
    return dam_pr;
  }

#ifdef MPI
  {
    MPI_Status statusMPI;
    unsigned long int idim_max_buf = 0;
    int ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                            &idim_max_buf,      1, MPI_UNSIGNED_LONG, origin, 0,
                            MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
                        list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf,  idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    dam_pr += apply_nbody_term_to_rank_spin(
      X, term0, tmp_v0, list_1buf, v1buf, tmp_v1, idim_max_buf, origin);
    dam_pr += apply_nbody_term_to_rank_spin(
      X, term1, tmp_v0, list_1buf, v1buf, tmp_v1, idim_max_buf, origin);
  }
#else
  fprintf(stdoutMPI, "Error: NBodyInterAll reached an MPI-only rank flip path without MPI.\n");
  return 0.0;
#endif
  return dam_pr;
}

static double complex apply_nbody_hubbardgc_term_to_rank(
  struct BindStruct *X,
  unsigned int term,
  double complex *tmp_v0,
  const double complex *src_v1,
  double complex *cur_v1,
  unsigned long int src_i_max,
  int rank_in
) {
  unsigned long int j;
  double complex dam_pr = 0.0;
  const int do_update = (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC);

#pragma omp parallel for default(none) reduction(+:dam_pr) \
  shared(X, tmp_v0, src_v1, cur_v1) firstprivate(src_i_max, term, rank_in, do_update, myrank) private(j)
  for (j = 1; j <= src_i_max; j++) {
    unsigned long int local_out = 0;
    int rank_out = 0;
    int sign = 1;
    int ret = apply_nbody_interall_hubbardgc_full(
      X, term, j - 1, rank_in, &local_out, &rank_out, &sign);
    if (ret == 1 && rank_out == myrank) {
      const double complex dmv = X->Def.ParaNBodyInterAll[term] * sign * src_v1[j];
      if (do_update) tmp_v0[local_out + 1] += dmv;
      dam_pr += conj(cur_v1[local_out + 1]) * dmv;
    }
  }
  return dam_pr;
}

static double complex apply_nbody_hubbard_term_to_rank(
  struct BindStruct *X,
  unsigned int term,
  double complex *tmp_v0,
  const unsigned long int *src_list_1,
  const double complex *src_v1,
  double complex *cur_v1,
  unsigned long int src_i_max,
  int rank_in
) {
  unsigned long int j;
  double complex dam_pr = 0.0;
  const int do_update = (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC);

#pragma omp parallel for default(none) reduction(+:dam_pr) \
  shared(X, tmp_v0, src_list_1, src_v1, cur_v1, list_2_1, list_2_2) \
  firstprivate(src_i_max, term, rank_in, do_update, myrank) private(j)
  for (j = 1; j <= src_i_max; j++) {
    unsigned long int local_out = 0;
    unsigned long int j_out = 0;
    int rank_out = 0;
    int sign = 1;
    int ret = apply_nbody_interall_hubbardgc_full(
      X, term, src_list_1[j], rank_in, &local_out, &rank_out, &sign);
    int in_sector = FALSE;
    if (ret == 1 && rank_out == myrank) {
      in_sector = GetOffComp(list_2_1, list_2_2, local_out,
                             X->Large.irght, X->Large.ilft, X->Large.ihfbit, &j_out);
    }
    if (in_sector == TRUE) {
      const double complex dmv = X->Def.ParaNBodyInterAll[term] * sign * src_v1[j];
      if (do_update) tmp_v0[j_out] += dmv;
      dam_pr += conj(cur_v1[j_out]) * dmv;
    }
  }
  return dam_pr;
}

static double complex apply_nbody_spinlessgc_term_to_rank(
  struct BindStruct *X,
  unsigned int term,
  double complex *tmp_v0,
  const double complex *src_v1,
  double complex *cur_v1,
  unsigned long int src_i_max,
  int rank_in
) {
  unsigned long int j;
  double complex dam_pr = 0.0;
  const int do_update = (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC);

#pragma omp parallel for default(none) reduction(+:dam_pr) \
  shared(X, tmp_v0, src_v1, cur_v1) firstprivate(src_i_max, term, rank_in, do_update, myrank) private(j)
  for (j = 1; j <= src_i_max; j++) {
    unsigned long int local_out = 0;
    int rank_out = 0;
    int sign = 1;
    int ret = apply_nbody_interall_spinless_full(
      X, term, j - 1, rank_in, &local_out, &rank_out, &sign);
    if (ret == 1 && rank_out == myrank) {
      const double complex dmv = X->Def.ParaNBodyInterAll[term] * sign * src_v1[j];
      if (do_update) tmp_v0[local_out + 1] += dmv;
      dam_pr += conj(cur_v1[local_out + 1]) * dmv;
    }
  }
  return dam_pr;
}

static double complex apply_nbody_spinless_term_to_rank(
  struct BindStruct *X,
  unsigned int term,
  double complex *tmp_v0,
  const unsigned long int *src_list_1,
  const double complex *src_v1,
  double complex *cur_v1,
  unsigned long int src_i_max,
  int rank_in
) {
  unsigned long int j;
  double complex dam_pr = 0.0;
  const int do_update = (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC);

#pragma omp parallel for default(none) reduction(+:dam_pr) \
  shared(X, tmp_v0, src_list_1, src_v1, cur_v1, list_2_1, list_2_2) \
  firstprivate(src_i_max, term, rank_in, do_update, myrank) private(j)
  for (j = 1; j <= src_i_max; j++) {
    unsigned long int local_out = 0;
    unsigned long int j_out = 0;
    int rank_out = 0;
    int sign = 1;
    int ret = apply_nbody_interall_spinless_full(
      X, term, src_list_1[j], rank_in, &local_out, &rank_out, &sign);
    int in_sector = FALSE;
    if (ret == 1 && rank_out == myrank) {
      in_sector = GetOffComp(list_2_1, list_2_2, local_out,
                             X->Large.irght, X->Large.ilft, X->Large.ihfbit, &j_out);
    }
    if (in_sector == TRUE) {
      const double complex dmv = X->Def.ParaNBodyInterAll[term] * sign * src_v1[j];
      if (do_update) tmp_v0[j_out] += dmv;
      dam_pr += conj(cur_v1[j_out]) * dmv;
    }
  }
  return dam_pr;
}

static double complex multiply_nbody_hubbardgc_term(
  struct BindStruct *X,
  unsigned int term,
  double complex *tmp_v0,
  double complex *tmp_v1
) {
  int origin = myrank;
  int active = TRUE;
  double complex dam_pr = 0.0;

  if (nbody_interall_hubbardgc_partner_rank(X, term, myrank, &origin, &active) != 0) {
    return 0.0;
  }
  if (active == FALSE) return 0.0;
  if (origin == myrank) {
    return apply_nbody_hubbardgc_term_to_rank(
      X, term, tmp_v0, tmp_v1, tmp_v1, X->Check.idim_max, myrank);
  }

#ifdef MPI
  {
    MPI_Status statusMPI;
    unsigned long int idim_max_buf = 0;
    int ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                            &idim_max_buf,      1, MPI_UNSIGNED_LONG, origin, 0,
                            MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf,  idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    dam_pr = apply_nbody_hubbardgc_term_to_rank(
      X, term, tmp_v0, v1buf, tmp_v1, idim_max_buf, origin);
  }
#else
  fprintf(stdoutMPI, "Error: NBodyInterAll reached an MPI-only rank flip path without MPI.\n");
  return 0.0;
#endif
  return dam_pr;
}

static double complex multiply_nbody_hubbard_term(
  struct BindStruct *X,
  unsigned int term,
  double complex *tmp_v0,
  double complex *tmp_v1
) {
  int origin = myrank;
  int active = TRUE;
  double complex dam_pr = 0.0;

  if (nbody_interall_hubbardgc_partner_rank(X, term, myrank, &origin, &active) != 0) {
    return 0.0;
  }
  if (active == FALSE) return 0.0;
  if (origin == myrank) {
    return apply_nbody_hubbard_term_to_rank(
      X, term, tmp_v0, list_1, tmp_v1, tmp_v1, X->Check.idim_max, myrank);
  }

#ifdef MPI
  {
    MPI_Status statusMPI;
    unsigned long int idim_max_buf = 0;
    int ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                            &idim_max_buf,      1, MPI_UNSIGNED_LONG, origin, 0,
                            MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
                        list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf,  idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    dam_pr = apply_nbody_hubbard_term_to_rank(
      X, term, tmp_v0, list_1buf, v1buf, tmp_v1, idim_max_buf, origin);
  }
#else
  fprintf(stdoutMPI, "Error: NBodyInterAll reached an MPI-only rank flip path without MPI.\n");
  return 0.0;
#endif
  return dam_pr;
}

static double complex multiply_nbody_spinlessgc_term(
  struct BindStruct *X,
  unsigned int term,
  double complex *tmp_v0,
  double complex *tmp_v1
) {
  int origin = myrank;
  int active = TRUE;
  double complex dam_pr = 0.0;

  if (nbody_interall_spinless_partner_rank(X, term, myrank, &origin, &active) != 0) {
    return 0.0;
  }
  if (active == FALSE) return 0.0;
  if (origin == myrank) {
    return apply_nbody_spinlessgc_term_to_rank(
      X, term, tmp_v0, tmp_v1, tmp_v1, X->Check.idim_max, myrank);
  }

#ifdef MPI
  {
    MPI_Status statusMPI;
    unsigned long int idim_max_buf = 0;
    int ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                            &idim_max_buf,      1, MPI_UNSIGNED_LONG, origin, 0,
                            MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf,  idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    dam_pr = apply_nbody_spinlessgc_term_to_rank(
      X, term, tmp_v0, v1buf, tmp_v1, idim_max_buf, origin);
  }
#else
  fprintf(stdoutMPI, "Error: NBodyInterAll reached an MPI-only rank flip path without MPI.\n");
  return 0.0;
#endif
  return dam_pr;
}

static double complex multiply_nbody_spinless_term(
  struct BindStruct *X,
  unsigned int term,
  double complex *tmp_v0,
  double complex *tmp_v1
) {
  int origin = myrank;
  int active = TRUE;
  double complex dam_pr = 0.0;

  if (nbody_interall_spinless_partner_rank(X, term, myrank, &origin, &active) != 0) {
    return 0.0;
  }
  if (active == FALSE) return 0.0;
  if (origin == myrank) {
    return apply_nbody_spinless_term_to_rank(
      X, term, tmp_v0, list_1, tmp_v1, tmp_v1, X->Check.idim_max, myrank);
  }

#ifdef MPI
  {
    MPI_Status statusMPI;
    unsigned long int idim_max_buf = 0;
    int ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                            &idim_max_buf,      1, MPI_UNSIGNED_LONG, origin, 0,
                            MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(list_1, X->Check.idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
                        list_1buf, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf,  idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    dam_pr = apply_nbody_spinless_term_to_rank(
      X, term, tmp_v0, list_1buf, v1buf, tmp_v1, idim_max_buf, origin);
  }
#else
  fprintf(stdoutMPI, "Error: NBodyInterAll reached an MPI-only rank flip path without MPI.\n");
  return 0.0;
#endif
  return dam_pr;
}

int MultiplyNBodyInterAllHubbardGC(
  struct BindStruct *X,
  double complex *tmp_v0,
  double complex *tmp_v1
) {
  unsigned int p;
  if (X->Def.NNBodyInterAll_OffDiagonal == 0) return 0;
  if (X->Def.iCalcModel != HubbardGC) return -1;

  for (p = 0; p < X->Def.NNBodyInterAll_OffDiagonal; p++) {
    const unsigned int term = X->Def.NBodyInterAll_OffDiagonalIndex[p];
    X->Large.prdct += multiply_nbody_hubbardgc_term(X, term, tmp_v0, tmp_v1);
  }
  return 0;
}

int MultiplyNBodyInterAllHubbard(
  struct BindStruct *X,
  double complex *tmp_v0,
  double complex *tmp_v1
) {
  unsigned int p;
  if (X->Def.NNBodyInterAll_OffDiagonal == 0) return 0;
  if (nbody_uses_hubbard_list_path(&X->Def) == FALSE) return -1;

  for (p = 0; p < X->Def.NNBodyInterAll_OffDiagonal; p++) {
    const unsigned int term = X->Def.NBodyInterAll_OffDiagonalIndex[p];
    X->Large.prdct += multiply_nbody_hubbard_term(X, term, tmp_v0, tmp_v1);
  }
  return 0;
}

int MultiplyNBodyInterAllSpinlessGC(
  struct BindStruct *X,
  double complex *tmp_v0,
  double complex *tmp_v1
) {
  unsigned int p;
  if (X->Def.NNBodyInterAll_OffDiagonal == 0) return 0;
  if (X->Def.iCalcModel != SpinlessFermionGC) return -1;

  for (p = 0; p < X->Def.NNBodyInterAll_OffDiagonal; p++) {
    const unsigned int term = X->Def.NBodyInterAll_OffDiagonalIndex[p];
    X->Large.prdct += multiply_nbody_spinlessgc_term(X, term, tmp_v0, tmp_v1);
  }
  return 0;
}

int MultiplyNBodyInterAllSpinless(
  struct BindStruct *X,
  double complex *tmp_v0,
  double complex *tmp_v1
) {
  unsigned int p;
  if (X->Def.NNBodyInterAll_OffDiagonal == 0) return 0;
  if (X->Def.iCalcModel != SpinlessFermion) return -1;

  for (p = 0; p < X->Def.NNBodyInterAll_OffDiagonal; p++) {
    const unsigned int term = X->Def.NBodyInterAll_OffDiagonalIndex[p];
    X->Large.prdct += multiply_nbody_spinless_term(X, term, tmp_v0, tmp_v1);
  }
  return 0;
}

int MultiplyNBodyInterAllSpinGC(
  struct BindStruct *X,
  double complex *tmp_v0,
  double complex *tmp_v1
) {
  unsigned int p;
  if (X->Def.NNBodyInterAll_OffDiagonal == 0) return 0;
  if (nbody_is_supported_spin_model(&X->Def) == FALSE) return -1;
  if (nbody_is_hubbard_model(&X->Def) == TRUE) return -1;

  for (p = 0; p < X->Def.NNBodyInterAll_OffDiagonal; p += 2) {
    if (X->Def.iCalcModel == Spin) {
      if (X->Def.iFlgGeneralSpin == TRUE) {
        X->Large.prdct += multiply_nbody_pair_spin_general_spin(X, p, tmp_v0, tmp_v1);
      }
      else {
        X->Large.prdct += multiply_nbody_pair_spin(X, p, tmp_v0, tmp_v1);
      }
    }
    else if (nbody_is_general_spin(&X->Def) == TRUE) {
      X->Large.prdct += multiply_nbody_pair_general_spin_gc(X, p, tmp_v0, tmp_v1);
    }
    else {
      X->Large.prdct += multiply_nbody_pair(X, p, tmp_v0, tmp_v1);
    }
  }
  return 0;
}

int AddNBodyInterAllToHamHubbard(struct BindStruct *X)
{
  unsigned int p;
  unsigned long int j;
  if (X->Def.NNBodyInterAll_OffDiagonal == 0) return 0;
  if (nbody_uses_hubbard_list_path(&X->Def) == FALSE) return -1;

  for (p = 0; p < X->Def.NNBodyInterAll_OffDiagonal; p++) {
    const unsigned int term = X->Def.NBodyInterAll_OffDiagonalIndex[p];
    for (j = 1; j <= X->Check.idim_max; j++) {
      unsigned long int local_out = 0;
      unsigned long int j_out = 0;
      int rank_out = 0;
      int sign = 1;
      int ret = apply_nbody_interall_hubbardgc_full(
        X, term, list_1[j], myrank, &local_out, &rank_out, &sign);
      if (ret == 1) {
        int in_sector;
        if (rank_out != myrank) {
          fprintf(stdoutMPI, "Error: FullDiag NBodyInterAll cannot handle inter-process output.\n");
          return -1;
        }
        in_sector = GetOffComp(list_2_1, list_2_2, local_out,
                               X->Large.irght, X->Large.ilft, X->Large.ihfbit, &j_out);
        if (in_sector == TRUE) {
          Ham[j_out][j] += X->Def.ParaNBodyInterAll[term] * sign;
        }
      }
    }
  }
  return 0;
}

int AddNBodyInterAllToHamHubbardGC(struct BindStruct *X)
{
  unsigned int p;
  unsigned long int j;
  if (X->Def.NNBodyInterAll_OffDiagonal == 0) return 0;
  if (X->Def.iCalcModel != HubbardGC) return -1;

  for (p = 0; p < X->Def.NNBodyInterAll_OffDiagonal; p++) {
    const unsigned int term = X->Def.NBodyInterAll_OffDiagonalIndex[p];
    for (j = 1; j <= X->Check.idim_max; j++) {
      unsigned long int local_out = 0;
      int rank_out = 0;
      int sign = 1;
      int ret = apply_nbody_interall_hubbardgc_full(
        X, term, j - 1, myrank, &local_out, &rank_out, &sign);
      if (ret == 1) {
        if (rank_out != myrank) {
          fprintf(stdoutMPI, "Error: FullDiag NBodyInterAll cannot handle inter-process output.\n");
          return -1;
        }
        Ham[local_out + 1][j] += X->Def.ParaNBodyInterAll[term] * sign;
      }
    }
  }
  return 0;
}

int AddNBodyInterAllToHamSpinGC(struct BindStruct *X)
{
  unsigned int p;
  unsigned long int j;
  if (X->Def.NNBodyInterAll_OffDiagonal == 0) return 0;
  if (nbody_is_supported_spin_model(&X->Def) == FALSE) return -1;
  if (nbody_is_hubbard_model(&X->Def) == TRUE) return -1;

  for (p = 0; p < X->Def.NNBodyInterAll_OffDiagonal; p++) {
    const unsigned int term = X->Def.NBodyInterAll_OffDiagonalIndex[p];
    for (j = 1; j <= X->Check.idim_max; j++) {
      unsigned long int local_out = 0;
      int rank_out = 0;
      double complex me = 0.0;
      const unsigned long int local_in = (X->Def.iCalcModel == Spin) ? list_1[j] : j - 1;
      int ret = ApplyNBodyInterAllSpinGC(X, term, local_in, myrank, &local_out, &rank_out, &me);
      if (ret == 1) {
        if (rank_out != myrank) {
          fprintf(stdoutMPI, "Error: FullDiag NBodyInterAll cannot handle inter-process output.\n");
          return -1;
        }
        if (X->Def.iCalcModel == Spin) {
          unsigned long int j_out = 0;
          const int in_sector = (X->Def.iFlgGeneralSpin == TRUE) ?
            convert_nbody_general_spin_to_list1(X, local_out, &j_out) :
            GetOffComp(list_2_1, list_2_2, local_out,
                       X->Large.irght, X->Large.ilft, X->Large.ihfbit, &j_out);
          if (in_sector == TRUE) {
            Ham[j_out][j] += me;
          }
        }
        else {
          Ham[local_out + 1][j] += me;
        }
      }
    }
  }
  return 0;
}
