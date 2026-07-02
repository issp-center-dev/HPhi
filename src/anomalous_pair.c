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
#include <stdlib.h>
#include "anomalous_pair.h"
#include "bitcalc.h"
#include "FileIO.h"
#include "mltplyCommon.h"
#include "wrapperMPI.h"

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

static int require_line_end(const char *p)
{
  while (isspace((unsigned char)*p)) p++;
  return (*p == '\0') ? 0 : -1;
}

int ParseAnomalousTermLine(
  const char *line,
  int term[5],
  double *re,
  double *im
) {
  const char *p = line;
  int i;
  for (i = 0; i < 5; i++) {
    if (parse_int_token(&p, &term[i]) != 0) {
      fprintf(stdoutMPI, "Error: AnomalousTerm line has too few integer fields.\n");
      return -1;
    }
  }
  if (parse_double_token(&p, re) != 0 || parse_double_token(&p, im) != 0) {
    fprintf(stdoutMPI, "Error: AnomalousTerm line has invalid coefficient fields.\n");
    return -1;
  }
  if (require_line_end(p) != 0) {
    fprintf(stdoutMPI, "Error: AnomalousTerm line has extra fields.\n");
    return -1;
  }
  return 0;
}

int ParseAnomalousGLine(
  const char *line,
  int term[5]
) {
  const char *p = line;
  int i;
  for (i = 0; i < 5; i++) {
    if (parse_int_token(&p, &term[i]) != 0) {
      fprintf(stdoutMPI, "Error: AnomalousG line has too few integer fields.\n");
      return -1;
    }
  }
  if (require_line_end(p) != 0) {
    fprintf(stdoutMPI, "Error: AnomalousG line has extra fields.\n");
    return -1;
  }
  return 0;
}

static int validate_anomalous_pair(const struct DefineList *D, const int term[5], const char *name)
{
  const int type = term[0];
  const int site1 = term[1];
  const int spin1 = term[2];
  const int site2 = term[3];
  const int spin2 = term[4];

  if (type != 0 && type != 1) {
    fprintf(stdoutMPI, "Error: %s type must be 0 or 1.\n", name);
    return -1;
  }
  if (site1 < 0 || site1 >= (int)D->Nsite ||
      site2 < 0 || site2 >= (int)D->Nsite) {
    fprintf(stdoutMPI, "Error: Site index of %s is incorrect.\n", name);
    return -1;
  }
  if (spin1 < 0 || spin1 > 1 || spin2 < 0 || spin2 > 1) {
    fprintf(stdoutMPI, "Error: Spin index of %s is incorrect.\n", name);
    return -1;
  }
  if (site1 == site2 && spin1 == spin2) {
    fprintf(stdoutMPI, "Error: %s cannot use the same fermion operator twice in one pair.\n", name);
    return -1;
  }
  return 0;
}

static int validate_anomalous_scope_common(const struct DefineList *D, unsigned int n, const char *name)
{
  unsigned int t;
  int **terms = (strcmp(name, "AnomalousTerm") == 0) ? D->AnomalousTerm : D->AnomalousG;
  if (n == 0) return 0;
  if (D->iCalcModel != HubbardGC) {
    fprintf(stdoutMPI, "Error: %s is currently supported only for HubbardGC.\n", name);
    return -1;
  }
  if (D->iFlgCalcSpec != CALCSPEC_NOT) {
    fprintf(stdoutMPI, "Error: %s does not support CalcSpec.\n", name);
    return -1;
  }
  if (nproc > 1 && D->iCalcType == FullDiag) {
    fprintf(stdoutMPI, "Error: %s does not support MPI FullDiag.\n", name);
    return -1;
  }
  for (t = 0; t < n; t++) {
    if (validate_anomalous_pair(D, terms[t], name) != 0) return -1;
  }
  return 0;
}

int ValidateAnomalousTermScope(const struct DefineList *D)
{
  if (D->NAnomalousTerm == 0) return 0;
  if (D->iCalcType == TimeEvolution) {
    fprintf(stdoutMPI, "Error: AnomalousTerm is not supported in TimeEvolution.\n");
    return -1;
  }
  return validate_anomalous_scope_common(D, D->NAnomalousTerm, "AnomalousTerm");
}

int ValidateAnomalousGScope(const struct DefineList *D)
{
  return validate_anomalous_scope_common(D, D->NAnomalousG, "AnomalousG");
}

int CheckAnomalousTermHermitePairs(const struct DefineList *D)
{
  unsigned int t;
  if (D->NAnomalousTerm == 0) return 0;
  if (D->NAnomalousTerm % 2 != 0) {
    fprintf(stdoutMPI, "Error: AnomalousTerm terms must appear as adjacent Hermite pairs.\n");
    return -1;
  }
  for (t = 0; t < D->NAnomalousTerm; t += 2) {
    const int *a = D->AnomalousTerm[t];
    const int *b = D->AnomalousTerm[t + 1];
    if (b[0] != 1 - a[0] ||
        b[1] != a[3] || b[2] != a[4] ||
        b[3] != a[1] || b[4] != a[2]) {
      fprintf(stdoutMPI, "Error: AnomalousTerm Hermite pair has inconsistent operators.\n");
      return -1;
    }
    if (cabs(D->ParaAnomalousTerm[t + 1] - conj(D->ParaAnomalousTerm[t])) > eps_CheckImag0) {
      fprintf(stdoutMPI, "Error: AnomalousTerm Hermite pair has inconsistent coefficients.\n");
      return -1;
    }
  }
  return 0;
}

static int apply_annihilate_mask(unsigned long int mask, unsigned long int *state, int *sign)
{
  int sgn = 1;
  if ((*state & mask) == 0) return 0;
  SgnBit(*state & (mask - 1), &sgn);
  *sign *= sgn;
  *state &= ~mask;
  return 1;
}

static int apply_create_mask(unsigned long int mask, unsigned long int *state, int *sign)
{
  int sgn = 1;
  if ((*state & mask) != 0) return 0;
  SgnBit(*state & (mask - 1), &sgn);
  *sign *= sgn;
  *state |= mask;
  return 1;
}

static unsigned long int get_hubbardgc_local_block(const struct DefineList *D)
{
  if (D->Nsite == 0) return 1;
  return D->OrgTpow[2 * D->Nsite - 1] * 2;
}

static int apply_anomalous_rank_annihilate(
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

static int apply_anomalous_rank_create(
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

static int apply_anomalous_rank_term(
  const struct DefineList *D,
  const int term[5],
  unsigned long int *rank_state
) {
  const int type = term[0];
  const unsigned int site1 = (unsigned int)term[1];
  const unsigned int spin1 = (unsigned int)term[2];
  const unsigned int site2 = (unsigned int)term[3];
  const unsigned int spin2 = (unsigned int)term[4];

  if (type == 1) {
    if (apply_anomalous_rank_create(D, site2, spin2, rank_state) == 0) return 0;
    if (apply_anomalous_rank_create(D, site1, spin1, rank_state) == 0) return 0;
  }
  else {
    if (apply_anomalous_rank_annihilate(D, site2, spin2, rank_state) == 0) return 0;
    if (apply_anomalous_rank_annihilate(D, site1, spin1, rank_state) == 0) return 0;
  }
  return 1;
}

static int anomalous_hubbardgc_partner_rank(
  const struct DefineList *D,
  const int term[5],
  int current_rank,
  int *partner_rank,
  int *active
) {
  int dagger_term[5];
  unsigned long int rank_state = (unsigned long int)current_rank;

  if (apply_anomalous_rank_term(D, term, &rank_state) == 1) {
    *active = TRUE;
    *partner_rank = (int)rank_state;
    return 0;
  }

  dagger_term[0] = 1 - term[0];
  dagger_term[1] = term[3];
  dagger_term[2] = term[4];
  dagger_term[3] = term[1];
  dagger_term[4] = term[2];
  rank_state = (unsigned long int)current_rank;
  if (apply_anomalous_rank_term(D, dagger_term, &rank_state) == 1) {
    *active = TRUE;
    *partner_rank = (int)rank_state;
    return 0;
  }

  *active = FALSE;
  *partner_rank = current_rank;
  return 0;
}

int ApplyAnomalousPairHubbardGC(
  const struct DefineList *D,
  const int term[5],
  unsigned long int local_in,
  int rank_in,
  unsigned long int *local_out,
  int *rank_out,
  int *sign
) {
  const int type = term[0];
  const unsigned int site1 = (unsigned int)term[1];
  const unsigned int spin1 = (unsigned int)term[2];
  const unsigned int site2 = (unsigned int)term[3];
  const unsigned int spin2 = (unsigned int)term[4];
  const unsigned long int mask1 = D->OrgTpow[2 * site1 + spin1];
  const unsigned long int mask2 = D->OrgTpow[2 * site2 + spin2];
  const unsigned long int block = get_hubbardgc_local_block(D);
  unsigned long int state = local_in + block * (unsigned long int)rank_in;

  *sign = 1;
  if (type == 1) {
    if (apply_create_mask(mask2, &state, sign) == 0) return 0;
    if (apply_create_mask(mask1, &state, sign) == 0) return 0;
  }
  else {
    if (apply_annihilate_mask(mask2, &state, sign) == 0) return 0;
    if (apply_annihilate_mask(mask1, &state, sign) == 0) return 0;
  }

  *local_out = state % block;
  *rank_out = (int)(state / block);
  return 1;
}

static double complex apply_anomalous_term_to_rank(
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
  shared(X, tmp_v0, src_v1, cur_v1) \
  firstprivate(src_i_max, term, rank_in, do_update, myrank) private(j)
  for (j = 1; j <= src_i_max; j++) {
    unsigned long int local_out = 0;
    int rank_out = 0;
    int sign = 1;
    int ret = ApplyAnomalousPairHubbardGC(
      &X->Def, X->Def.AnomalousTerm[term], j - 1, rank_in, &local_out, &rank_out, &sign);
    if (ret == 1 && rank_out == myrank) {
      const double complex dmv = X->Def.ParaAnomalousTerm[term] * sign * src_v1[j];
      if (do_update) tmp_v0[local_out + 1] += dmv;
      dam_pr += conj(cur_v1[local_out + 1]) * dmv;
    }
  }
  return dam_pr;
}

static double complex multiply_anomalous_term(
  struct BindStruct *X,
  unsigned int term,
  double complex *tmp_v0,
  double complex *tmp_v1
) {
  int partner = myrank;
  int active = TRUE;
  double complex dam_pr = 0.0;

  if (anomalous_hubbardgc_partner_rank(&X->Def, X->Def.AnomalousTerm[term],
                                       myrank, &partner, &active) != 0) {
    return 0.0;
  }
  if (active == FALSE) return 0.0;
  if (partner == myrank) {
    return apply_anomalous_term_to_rank(
      X, term, tmp_v0, tmp_v1, tmp_v1, X->Check.idim_max, myrank);
  }

#ifdef MPI
  {
    MPI_Status statusMPI;
    unsigned long int idim_max_buf = 0;
    int ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, partner, 0,
                            &idim_max_buf,      1, MPI_UNSIGNED_LONG, partner, 0,
                            MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, partner, 0,
                        v1buf,  idim_max_buf + 1, MPI_DOUBLE_COMPLEX, partner, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    dam_pr = apply_anomalous_term_to_rank(
      X, term, tmp_v0, v1buf, tmp_v1, idim_max_buf, partner);
  }
#else
  fprintf(stdoutMPI, "Error: AnomalousTerm reached an MPI-only rank flip path without MPI.\n");
  return 0.0;
#endif
  return dam_pr;
}

int MultiplyAnomalousTermHubbardGC(
  struct BindStruct *X,
  double complex *tmp_v0,
  double complex *tmp_v1
) {
  unsigned int t;
  if (X->Def.NAnomalousTerm == 0) return 0;
  if (X->Def.iCalcModel != HubbardGC) return -1;

  for (t = 0; t < X->Def.NAnomalousTerm; t++) {
    X->Large.prdct += multiply_anomalous_term(X, t, tmp_v0, tmp_v1);
  }
  return 0;
}

int AddAnomalousTermToHamHubbardGC(struct BindStruct *X)
{
  unsigned int t;
  unsigned long int j;
  if (X->Def.NAnomalousTerm == 0) return 0;
  if (X->Def.iCalcModel != HubbardGC) return -1;

  for (t = 0; t < X->Def.NAnomalousTerm; t++) {
    for (j = 1; j <= X->Check.idim_max; j++) {
      unsigned long int local_out = 0;
      int rank_out = 0;
      int sign = 1;
      int ret = ApplyAnomalousPairHubbardGC(
        &X->Def, X->Def.AnomalousTerm[t], j - 1, myrank, &local_out, &rank_out, &sign);
      if (ret == 1) {
        if (rank_out != myrank) {
          fprintf(stdoutMPI, "Error: FullDiag AnomalousTerm cannot handle inter-process output.\n");
          return -1;
        }
        Ham[local_out + 1][j] += X->Def.ParaAnomalousTerm[t] * sign;
      }
    }
  }
  return 0;
}

static double complex calc_anomalousg_term_hubbardgc_rank(
  struct BindStruct *X,
  unsigned int term,
  const double complex *src_vec,
  const double complex *cur_vec,
  unsigned long int src_i_max,
  int rank_in
) {
  unsigned long int j;
  double complex value = 0.0;

#pragma omp parallel for default(none) reduction(+:value) \
  shared(X, src_vec, cur_vec) firstprivate(src_i_max, term, rank_in, myrank) private(j)
  for (j = 1; j <= src_i_max; j++) {
    unsigned long int local_out = 0;
    int rank_out = 0;
    int sign = 1;
    int ret = ApplyAnomalousPairHubbardGC(
      &X->Def, X->Def.AnomalousG[term], j - 1, rank_in, &local_out, &rank_out, &sign);
    if (ret == 1 && rank_out == myrank) {
      value += conj(cur_vec[local_out + 1]) * sign * src_vec[j];
    }
  }
  return value;
}

static double complex calc_anomalousg_term_hubbardgc(
  struct BindStruct *X,
  unsigned int term,
  double complex *vec
) {
  int partner = myrank;
  int active = TRUE;
  double complex value = 0.0;

  if (anomalous_hubbardgc_partner_rank(&X->Def, X->Def.AnomalousG[term],
                                       myrank, &partner, &active) != 0) {
    return SumMPI_dc(0.0);
  }
  if (active == FALSE) {
    return SumMPI_dc(0.0);
  }
  if (partner == myrank) {
    value = calc_anomalousg_term_hubbardgc_rank(
      X, term, vec, vec, X->Check.idim_max, myrank);
    return SumMPI_dc(value);
  }

#ifdef MPI
  {
    MPI_Status statusMPI;
    unsigned long int idim_max_buf = 0;
    int ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, partner, 0,
                            &idim_max_buf,      1, MPI_UNSIGNED_LONG, partner, 0,
                            MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    ierr = MPI_Sendrecv(vec, X->Check.idim_max + 1, MPI_DOUBLE_COMPLEX, partner, 0,
                        v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, partner, 0,
                        MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);
    value = calc_anomalousg_term_hubbardgc_rank(
      X, term, v1buf, vec, idim_max_buf, partner);
  }
#else
  fprintf(stdoutMPI, "Error: AnomalousG reached an MPI-only rank flip path without MPI.\n");
  return SumMPI_dc(0.0);
#endif
  return SumMPI_dc(value);
}

static int write_anomalousg_line(FILE *fp, const struct DefineList *D, unsigned int term, double complex value)
{
  const int *op = D->AnomalousG[term];
  fprintf(fp, "%4d %4d %4d %4d %4d %.10lf %.10lf\n",
          op[0], op[1], op[2], op[3], op[4], creal(value), cimag(value));
  return 0;
}

static int get_anomalousg_filename(struct BindStruct *X, char *sdt)
{
  switch (X->Def.iCalcType) {
  case Lanczos:
  case CG:
    sprintf(sdt, cFileNameAnomalousG_Lanczos, X->Def.CDataFileHead);
    break;
  case TPQCalc:
  case cTPQ:
    sprintf(sdt, cFileNameAnomalousG_TPQ, X->Def.CDataFileHead, X->Def.irand, X->Def.istep);
    break;
  case TimeEvolution:
    sprintf(sdt, cFileNameAnomalousG_TE, X->Def.CDataFileHead, X->Def.istep);
    break;
  case FullDiag:
    sprintf(sdt, cFileNameAnomalousG_FullDiag, X->Def.CDataFileHead, X->Phys.eigen_num);
    break;
  default:
    fprintf(stdoutMPI, "Error: AnomalousG does not support this calculation type.\n");
    return -1;
  }
  return 0;
}

int expec_anomalousg(struct BindStruct *X, double complex *vec)
{
  FILE *fp;
  char sdt[D_FileNameMax];
  unsigned int t;

  if (X->Def.NAnomalousG < 1) return 0;
  if (X->Def.iCalcModel != HubbardGC) {
    fprintf(stdoutMPI, "Error: AnomalousG is currently supported only for HubbardGC.\n");
    return -1;
  }
  if (get_anomalousg_filename(X, sdt) != 0) return -1;
  if (childfopenMPI(sdt, "w", &fp) != 0) return -1;

  for (t = 0; t < X->Def.NAnomalousG; t++) {
    const double complex value = calc_anomalousg_term_hubbardgc(X, t, vec);
    write_anomalousg_line(fp, &X->Def, t, value);
  }

  fclose(fp);
  return 0;
}
