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
#include <stdint.h>
#include <stdlib.h>
#include "bitcalc.h"
#include "nbody_correlation.h"
#include "FileIO.h"
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

int ParseNBodyGLine(
  const char *line,
  unsigned int *N,
  int **factors
) {
  const char *p = line;
  unsigned int n;
  unsigned int k;
  int *buf;
  size_t nints;

  if (parse_unsigned_token(&p, &n) != 0 || n == 0) {
    fprintf(stdoutMPI, "Error: NBodyG line has an invalid factor count.\n");
    return -1;
  }
  if (n > UINT_MAX / 4) {
    fprintf(stdoutMPI, "Error: NBodyG line is too large.\n");
    return -1;
  }
#if SIZE_MAX < UINT_MAX
  if ((size_t)n > SIZE_MAX / 4 / sizeof(int)) {
    fprintf(stdoutMPI, "Error: NBodyG line is too large.\n");
    return -1;
  }
#endif
  nints = (size_t)4 * n;
  buf = (int *)malloc(nints * sizeof(int));
  if (buf == NULL) {
    fprintf(stdoutMPI, "Error: Failed to allocate NBodyG parser buffer.\n");
    return -1;
  }
  for (k = 0; k < 4 * n; k++) {
    if (parse_int_token(&p, &buf[k]) != 0) {
      fprintf(stdoutMPI, "Error: NBodyG line has too few integer fields.\n");
      free(buf);
      return -1;
    }
  }
  while (isspace((unsigned char)*p)) p++;
  if (*p != '\0') {
    fprintf(stdoutMPI, "Error: NBodyG line has extra fields.\n");
    free(buf);
    return -1;
  }

  *N = n;
  *factors = buf;
  return 0;
}

int ValidateNBodyGScope(const struct DefineList *D)
{
  unsigned int t, k;
  if (D->NNBodyG == 0) return 0;
  if (D->iCalcModel != SpinGC && D->iCalcModel != Spin) {
    fprintf(stdoutMPI, "Error: NBodyG is currently supported only for spin-1/2 SpinGC/Spin.\n");
    return -1;
  }
  if (D->iFlgGeneralSpin != FALSE) {
    fprintf(stdoutMPI, "Error: NBodyG is currently supported only for spin-1/2 SpinGC/Spin.\n");
    return -1;
  }
  for (t = 0; t < D->NNBodyG; t++) {
    const unsigned int off = D->NBodyG_Offset[t];
    for (k = 0; k < D->NBodyG_N[t]; k++) {
      const int *f = D->NBodyG_Factors[off + k];
      if (f[0] < 0 || f[0] >= (int)D->Nsite || f[2] < 0 || f[2] >= (int)D->Nsite) {
        fprintf(stdoutMPI, "Error: Site index of NBodyG is incorrect.\n");
        return -1;
      }
      if (f[0] != f[2]) {
        fprintf(stdoutMPI, "Error: NBodyG currently requires site_out == site_in for every factor.\n");
        return -1;
      }
      if (f[1] < 0 || f[1] > 1 || f[3] < 0 || f[3] > 1) {
        fprintf(stdoutMPI, "Error: Spin index of NBodyG is incorrect.\n");
        return -1;
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

int NormalizeNBodyGTerms(struct DefineList *D)
{
  unsigned int t, k;
  unsigned int total = 0;

  D->NBodyG_TotalCanonicalFactors = 0;
  for (t = 0; t < D->NNBodyG; t++) {
    const unsigned int nraw = D->NBodyG_N[t];
    const unsigned int off = D->NBodyG_Offset[t];
    int *sites = (int *)malloc(nraw * sizeof(int));
    int *outs = (int *)malloc(nraw * sizeof(int));
    int *ins = (int *)malloc(nraw * sizeof(int));
    unsigned int ncanon = 0;
    int is_zero = FALSE;
    if (sites == NULL || outs == NULL || ins == NULL) {
      fprintf(stdoutMPI, "Error: Failed to allocate NBodyG normalization buffer.\n");
      free(sites);
      free(outs);
      free(ins);
      return -1;
    }

    for (k = 0; k < nraw; k++) {
      const int *f = D->NBodyG_Factors[off + k];
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
          is_zero = TRUE;
          break;
        }
        ins[pos] = spin_in;
      }
    }

    D->NBodyG_IsZero[t] = is_zero;
    D->NBodyG_CanonicalOffset[t] = total;
    if (is_zero == TRUE) {
      D->NBodyG_CanonicalN[t] = 0;
      free(sites);
      free(outs);
      free(ins);
      continue;
    }

    sort_canonical(sites, outs, ins, ncanon);
    D->NBodyG_CanonicalN[t] = ncanon;
    for (k = 0; k < ncanon; k++) {
      D->NBodyG_CanonicalFactors[total + k][0] = sites[k];
      D->NBodyG_CanonicalFactors[total + k][1] = outs[k];
      D->NBodyG_CanonicalFactors[total + k][2] = sites[k];
      D->NBodyG_CanonicalFactors[total + k][3] = ins[k];
    }
    total += ncanon;
    free(sites);
    free(outs);
    free(ins);
  }
  D->NBodyG_TotalCanonicalFactors = total;
  return 0;
}

int CheckNBodyGSpinConservation(const struct DefineList *D)
{
  unsigned int t, k;
  if (D->NNBodyG == 0) return 0;
  if (D->iCalcModel != Spin || D->iFlgGeneralSpin != FALSE) return 0;

  for (t = 0; t < D->NNBodyG; t++) {
    const unsigned int n = D->NBodyG_CanonicalN[t];
    const unsigned int off = D->NBodyG_CanonicalOffset[t];
    int delta_nup = 0;
    if (D->NBodyG_IsZero[t] == TRUE) continue;
    for (k = 0; k < n; k++) {
      const int *f = D->NBodyG_CanonicalFactors[off + k];
      delta_nup += f[1] - f[3];
    }
    if (delta_nup != 0) {
      fprintf(stdoutMPI,
              "Error: NBodyG term %u does not conserve total Sz: delta2Sz=%d.\n",
              t + 1, 2 * delta_nup);
      return -1;
    }
  }
  return 0;
}

static int apply_nbodyg_spingc(
  const struct BindStruct *X,
  unsigned int term,
  unsigned long int local_in,
  int rank_in,
  unsigned long int *local_out,
  int *rank_out
) {
  const struct DefineList *D = &X->Def;
  unsigned int k;
  unsigned long int lo = local_in;
  int ro = rank_in;
  const unsigned int n = D->NBodyG_CanonicalN[term];
  const unsigned int off = D->NBodyG_CanonicalOffset[term];

  for (k = 0; k < n; k++) {
    const int *f = D->NBodyG_CanonicalFactors[off + k];
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

  *local_out = lo;
  *rank_out = ro;
  return 1;
}

static int nbodyg_rank_flip_mask(const struct BindStruct *X, unsigned int term, int *mask)
{
  unsigned int k;
  int m = 0;
  const unsigned int n = X->Def.NBodyG_CanonicalN[term];
  const unsigned int off = X->Def.NBodyG_CanonicalOffset[term];
  for (k = 0; k < n; k++) {
    const int *f = X->Def.NBodyG_CanonicalFactors[off + k];
    const unsigned int site = (unsigned int)f[0];
    if (site >= X->Def.Nsite && f[1] != f[3]) {
      m ^= (int)X->Def.Tpow[site];
    }
  }
  *mask = m;
  return 0;
}

static double complex expec_nbodyg_term_to_rank(
  struct BindStruct *X,
  unsigned int term,
  const double complex *src_vec,
  const double complex *bra_vec,
  int rank_in
) {
  unsigned long int j;
  double complex dam_pr = 0.0;
  const unsigned long int i_max = X->Check.idim_max;

#pragma omp parallel for default(none) reduction(+:dam_pr) \
  shared(X, src_vec, bra_vec) firstprivate(i_max, term, rank_in, myrank) private(j)
  for (j = 1; j <= i_max; j++) {
    unsigned long int local_out = 0;
    int rank_out = 0;
    int ret = apply_nbodyg_spingc(X, term, j - 1, rank_in, &local_out, &rank_out);
    if (ret == 1 && rank_out == myrank) {
      dam_pr += conj(bra_vec[local_out + 1]) * src_vec[j];
    }
  }
  return dam_pr;
}

static double complex expec_nbodyg_term_to_rank_spin(
  struct BindStruct *X,
  unsigned int term,
  const unsigned long int *src_list_1,
  const double complex *src_vec,
  const double complex *bra_vec,
  unsigned long int src_i_max,
  int rank_in
) {
  unsigned long int j;
  double complex dam_pr = 0.0;

#pragma omp parallel for default(none) reduction(+:dam_pr) \
  shared(X, src_list_1, src_vec, bra_vec, list_2_1, list_2_2) \
  firstprivate(src_i_max, term, rank_in, myrank) private(j)
  for (j = 1; j <= src_i_max; j++) {
    unsigned long int local_out = 0;
    unsigned long int j_out = 0;
    int rank_out = 0;
    int ret = apply_nbodyg_spingc(X, term, src_list_1[j], rank_in, &local_out, &rank_out);
    if (ret == 1 && rank_out == myrank &&
        GetOffComp(list_2_1, list_2_2, local_out,
                   X->Large.irght, X->Large.ilft, X->Large.ihfbit, &j_out) == TRUE) {
      dam_pr += conj(bra_vec[j_out]) * src_vec[j];
    }
  }
  return dam_pr;
}

static double complex calc_nbodyg_term(struct BindStruct *X, unsigned int term, double complex *vec)
{
  int mask = 0;
  int origin;
  double complex dam_pr = 0.0;

  nbodyg_rank_flip_mask(X, term, &mask);
  origin = myrank ^ mask;
  if (origin == myrank) {
    dam_pr = expec_nbodyg_term_to_rank(X, term, vec, vec, myrank);
    return SumMPI_dc(dam_pr);
  }

#ifdef MPI
  {
    MPI_Status statusMPI;
    int ierr = MPI_Sendrecv(vec, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                            v1buf, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                            MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    dam_pr = expec_nbodyg_term_to_rank(X, term, v1buf, vec, origin);
  }
#else
  fprintf(stdoutMPI, "Error: NBodyG reached an MPI-only rank flip path without MPI.\n");
  return 0.0;
#endif
  return SumMPI_dc(dam_pr);
}

static double complex calc_nbodyg_term_spin(struct BindStruct *X, unsigned int term, double complex *vec)
{
  int mask = 0;
  int origin;
  double complex dam_pr = 0.0;

  nbodyg_rank_flip_mask(X, term, &mask);
  origin = myrank ^ mask;
  if (origin == myrank) {
    dam_pr = expec_nbodyg_term_to_rank_spin(
      X, term, list_1, vec, vec, X->Check.idim_max, myrank);
    return SumMPI_dc(dam_pr);
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
    ierr = MPI_Sendrecv(vec, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    dam_pr = expec_nbodyg_term_to_rank_spin(
      X, term, list_1buf, v1buf, vec, idim_max_buf, origin);
  }
#else
  fprintf(stdoutMPI, "Error: NBodyG reached an MPI-only rank flip path without MPI.\n");
  return 0.0;
#endif
  return SumMPI_dc(dam_pr);
}

static int write_nbodyg_line(FILE *fp, const struct DefineList *D, unsigned int term, double complex value)
{
  unsigned int k;
  const unsigned int n = D->NBodyG_N[term];
  const unsigned int off = D->NBodyG_Offset[term];
  fprintf(fp, "%u", n);
  for (k = 0; k < n; k++) {
    const int *f = D->NBodyG_Factors[off + k];
    fprintf(fp, " %4d %4d %4d %4d", f[0], f[1], f[2], f[3]);
  }
  fprintf(fp, " %.10lf %.10lf\n", creal(value), cimag(value));
  return 0;
}

static int get_nbodyg_filename(struct BindStruct *X, char *sdt)
{
  switch (X->Def.iCalcType) {
  case Lanczos:
    sprintf(sdt, cFileNameNBodyG_Lanczos, X->Def.CDataFileHead);
    break;
  case TPQCalc:
  case cTPQ:
    sprintf(sdt, cFileNameNBodyG_TPQ, X->Def.CDataFileHead, X->Def.irand, X->Def.istep);
    break;
  case TimeEvolution:
    sprintf(sdt, cFileNameNBodyG_TE, X->Def.CDataFileHead, X->Def.istep);
    break;
  case FullDiag:
  case CG:
    sprintf(sdt, cFileNameNBodyG_FullDiag, X->Def.CDataFileHead, X->Phys.eigen_num);
    break;
  default:
    fprintf(stdoutMPI, "Error: NBodyG does not support this calculation type.\n");
    return -1;
  }
  return 0;
}

int expec_nbodyg(struct BindStruct *X, double complex *vec)
{
  FILE *fp;
  char sdt[D_FileNameMax];
  unsigned int t;

  if (X->Def.NNBodyG < 1) return 0;
  if ((X->Def.iCalcModel != SpinGC && X->Def.iCalcModel != Spin) ||
      X->Def.iFlgGeneralSpin != FALSE) {
    fprintf(stdoutMPI, "Error: NBodyG is currently supported only for spin-1/2 SpinGC/Spin.\n");
    return -1;
  }
  if (get_nbodyg_filename(X, sdt) != 0) return -1;
  if (childfopenMPI(sdt, "w", &fp) != 0) return -1;

  for (t = 0; t < X->Def.NNBodyG; t++) {
    double complex value = 0.0;
    if (X->Def.NBodyG_IsZero[t] == FALSE) {
      if (X->Def.iCalcModel == Spin) value = calc_nbodyg_term_spin(X, t, vec);
      else value = calc_nbodyg_term(X, t, vec);
    }
    write_nbodyg_line(fp, &X->Def, t, value);
  }

  fclose(fp);
  return 0;
}
