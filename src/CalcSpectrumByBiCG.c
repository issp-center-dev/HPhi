/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 Takahiro Misawa, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Youhei Yamaji, Synge Todo, Naoki Kawashima */

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
/**@file
@author Mitsuaki Kawamura (The University of Tokyo)
@brief  File for givinvg functions of calculating spectrum by Lanczos
*/
#include "Common.h"
#include "CalcSpectrumByLanczos.h"
#include "Lanczos_EigenValue.h"
#include "FileIO.h"
#include "wrapperMPI.h"
#include "common/setmemory.h"
#include "komega/komega.h"
#include "mltply.h"
#include <math.h>
#ifdef MPI
#include <mpi.h>
#endif

#define BICG_STATUS_NONFINITE 5
#define BICG_DIAG_SAMPLES 3
#define BICG_SPIKE_RATIO_WARN 1.0e3

typedef struct {
  double sum2;
  unsigned long nbad;
  int nsample;
  unsigned long sample_index[BICG_DIAG_SAMPLES];
  double sample_real[BICG_DIAG_SAMPLES];
  double sample_imag[BICG_DIAG_SAMPLES];
} BiCGVectorDiag;

static int IsFiniteComplex(double complex z) {
  return isfinite(creal(z)) && isfinite(cimag(z));
}

static void AnalyzeBiCGVector(const double complex *v, unsigned long n, BiCGVectorDiag *diag) {
  unsigned long i;
  diag->sum2 = 0.0;
  diag->nbad = 0;
  diag->nsample = 0;
  for (i = 0; i < n; i++) {
    const double real = creal(v[i]);
    const double imag = cimag(v[i]);
    if (IsFiniteComplex(v[i])) {
      diag->sum2 += real * real + imag * imag;
    } else {
      if (diag->nsample < BICG_DIAG_SAMPLES) {
        diag->sample_index[diag->nsample] = i + 1;
        diag->sample_real[diag->nsample] = real;
        diag->sample_imag[diag->nsample] = imag;
        diag->nsample++;
      }
      diag->nbad++;
    }
  }
}

static void PrintBiCGVectorDiag(const char *name, const BiCGVectorDiag *diag) {
  int i;
  fprintf(stderr, "%s_sum2=%25.17e %s_nbad=%lu", name, diag->sum2, name, diag->nbad);
  for (i = 0; i < diag->nsample; i++) {
    fprintf(stderr, " %s_bad[%d]=(idx=%lu,value=%25.17e,%25.17e)",
            name, i, diag->sample_index[i], diag->sample_real[i], diag->sample_imag[i]);
  }
}

static const char *BiCGStatusReason(int status) {
  switch (status) {
    case 0: return "converged";
    case 1: return "not converged";
    case 2: return "alpha breakdown";
    case 3: return "pi breakdown";
    case 4: return "rho breakdown";
    case BICG_STATUS_NONFINITE: return "non-finite detected";
    default: return "unknown";
  }
}

static void PrintBiCGIteration1Diag(
  const BiCGVectorDiag *v2_diag,
  const BiCGVectorDiag *v12_diag,
  int status0,
  int status1,
  int status2,
  double residual
) {
  fprintf(stderr,
          "BiCG iteration-1 diagnostic: rank=%d status=(%d,%d,%d) residual=%25.17e ",
          myrank, status0, status1, status2, residual);
  PrintBiCGVectorDiag("v2", v2_diag);
  fprintf(stderr, " ");
  PrintBiCGVectorDiag("Hv2", v12_diag);
  fprintf(stderr, "\n");
  fflush(stderr);
}

static double BiCGResidualRatio(double residual, double initial_residual) {
  if (isfinite(residual) == FALSE || isfinite(initial_residual) == FALSE ||
      initial_residual <= 0.0) {
    return NAN;
  }
  return residual / initial_residual;
}

static void PrintBiCGStatusHeader(
  FILE *fp,
  const struct EDMainCalStruct *X,
  int nBra,
  int Nomega,
  const double complex *dcomega,
  double initial_residual
) {
  fprintf(fp, "# HPhi BiCG residual diagnostics\n");
  fprintf(fp, "# nproc=%d nthreads=%d mpi_batching=%s nBra=%d Nomega=%d idim_max=%lu Lanczos_max=%u threshold=%25.17e\n",
          nproc, nthreads, iFlgMPIBatch ? "ON" : "OFF", nBra, Nomega,
          X->Bind.Check.idim_max, X->Bind.Def.Lanczos_max, eps_Lanczos);
  fprintf(fp, "# omega_org=(%25.17e,%25.17e) omega_first=(%25.17e,%25.17e) omega_last=(%25.17e,%25.17e)\n",
          creal(X->Bind.Def.dcOmegaOrg), cimag(X->Bind.Def.dcOmegaOrg),
          creal(dcomega[0]), cimag(dcomega[0]),
          creal(dcomega[Nomega - 1]), cimag(dcomega[Nomega - 1]));
  fprintf(fp, "# initial_residual_2_norm=%25.17e spike_warning_ratio=%25.17e\n",
          initial_residual, BICG_SPIKE_RATIO_WARN);
  fprintf(fp, "# columns: iter status seed unshifted_resnorm initial_resnorm max_unshifted_resnorm ratio_to_initial first_spike_iter shifted_residual_max shifted_residual_min shifted_below_threshold shifted_nonfinite\n");
}

static void PrintBiCGStatusTrace(
  FILE *fp,
  int stp,
  const int status[3],
  double unshifted_residual,
  double initial_residual,
  double max_unshifted_residual,
  int first_spike_iter,
  int Nomega,
  const double *shifted_residual
) {
  int iomega;
  double ratio = BiCGResidualRatio(max_unshifted_residual, initial_residual);
  double shifted_max = NAN;
  double shifted_min = NAN;
  int n_below_threshold = 0;
  int n_nonfinite = 0;
  for (iomega = 0; iomega < Nomega; iomega++) {
    const double res = shifted_residual == NULL ? NAN : shifted_residual[iomega];
    if (isfinite(res) == FALSE) {
      n_nonfinite++;
      continue;
    }
    if (isfinite(shifted_max) == FALSE || res > shifted_max) shifted_max = res;
    if (isfinite(shifted_min) == FALSE || res < shifted_min) shifted_min = res;
    if (res < eps_Lanczos) n_below_threshold++;
  }
  fprintf(fp,
          "%7d %7d %7d %25.17e %25.17e %25.17e %25.17e %7d "
          "%25.17e %25.17e %7d %7d\n",
          stp, status[1], status[2], unshifted_residual,
          initial_residual, max_unshifted_residual, ratio, first_spike_iter,
          shifted_max, shifted_min, n_below_threshold, n_nonfinite);
}

static void PrintBiCGResidualSummary(
  FILE *stream,
  const char *prefix,
  int last_iter,
  const int status[3],
  double initial_residual,
  double max_unshifted_residual,
  int first_spike_iter
) {
  const double ratio = BiCGResidualRatio(max_unshifted_residual, initial_residual);
  fprintf(stream,
          "%sBiCG residual summary: last_iter=%d status=%d [%s] seed=%d "
          "initial_resnorm=%25.17e max_resnorm=%25.17e max_over_initial=%25.17e "
          "first_spike_iter=%d spike_warning_ratio=%25.17e\n",
          prefix, last_iter, status[1], BiCGStatusReason(status[1]), status[2],
          initial_residual, max_unshifted_residual, ratio, first_spike_iter,
          BICG_SPIKE_RATIO_WARN);
}
/**@brief
Read @f$\alpha, \beta@f$, projected residual for restart
*/
void ReadTMComponents_BiCG(
  struct EDMainCalStruct *X,//!<[inout]
  double complex *v2,//!<[inout] [CheckList::idim_max] Residual vector
  double complex *v4,//!<[inout] [CheckList::idim_max] Shadow esidual vector
  double complex *v12,//!<[inout] [CheckList::idim_max] Old residual vector
  double complex *v14,//!<[inout] [CheckList::idim_max] Old shadow residual vector
  int nBra,//!<[in] Number of left (bra) vectors (Komega nl). Restart paths only ever run with nBra==1.
  int Nomega,//!<[in] Number of frequencies
  double complex *dcSpectrum,//!<[inout] [nBra*Nomega] Projected result vector, spectrum
  double complex *dcomega//!<[in] [Nomega] Frequency
) {
  char sdt[D_FileNameMax];
  char ctmp[256];

  int status[3], idim_max2int, max_step, iter_old;
  unsigned long int idx;
  double complex *alphaCG, *betaCG, *res_save, z_seed;
  double z_seed_r, z_seed_i, alpha_r, alpha_i, beta_r, beta_i, res_r, res_i;
  FILE *fp;
  int comm;

#if defined(MPI)
  comm = MPI_Comm_c2f(MPI_COMM_WORLD);
#else
  comm = 0;
#endif
  idim_max2int = (int)X->Bind.Check.idim_max;

  if (X->Bind.Def.iFlgCalcSpec == RECALC_FROM_TMComponents ||
      X->Bind.Def.iFlgCalcSpec == RECALC_FROM_TMComponents_VEC ||
      X->Bind.Def.iFlgCalcSpec == RECALC_INOUT_TMComponents_VEC)  {
    sprintf(sdt, cFileNameTridiagonalMatrixComponents, X->Bind.Def.CDataFileHead);
    if (childfopenALL(sdt, "rb", &fp) != 0) {
      fprintf(stdoutMPI, "INFO: File for the restart is not found.\n");
      fprintf(stdoutMPI, "      Start from SCRATCH.\n");
      max_step = (int)X->Bind.Def.Lanczos_max;
      komega_bicg_init(&idim_max2int, &nBra, &Nomega, dcSpectrum, dcomega, &max_step, &eps_Lanczos, &comm);
    }
    else {
      fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);
      sscanf(ctmp, "%d", &iter_old);
      if (X->Bind.Def.iFlgCalcSpec > RECALC_FROM_TMComponents) {
        alphaCG = (double complex*)malloc((iter_old + X->Bind.Def.Lanczos_max) * sizeof(double complex));
        betaCG = (double complex*)malloc((iter_old + X->Bind.Def.Lanczos_max) * sizeof(double complex));
        res_save = (double complex*)malloc((iter_old + X->Bind.Def.Lanczos_max) * sizeof(double complex));
      }
      else {
        alphaCG = (double complex*)malloc(iter_old * sizeof(double complex));
        betaCG = (double complex*)malloc(iter_old * sizeof(double complex));
        res_save = (double complex*)malloc(iter_old * sizeof(double complex));
      }
      fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);
      sscanf(ctmp, "%lf %lf\n", &z_seed_r, &z_seed_i);
      z_seed = z_seed_r + I * z_seed_i;

      idx = 0;
      while (fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp) != NULL) {
        sscanf(ctmp, "%lf %lf %lf %lf %lf %lf\n",
          &alpha_r, &alpha_i, &beta_r, &beta_i, &res_r, &res_i);
        alphaCG[idx] = alpha_r + I * alpha_i;
        betaCG[idx] = beta_r + I * beta_i;
        res_save[idx] = res_r + I * res_i;
        idx += 1;
      }
      fclose(fp);

      if (X->Bind.Def.iFlgCalcSpec == RECALC_FROM_TMComponents) X->Bind.Def.Lanczos_max = 0;
      max_step = (int)(iter_old + X->Bind.Def.Lanczos_max);

      komega_bicg_restart(&idim_max2int, &nBra, &Nomega, dcSpectrum, dcomega, &max_step, &eps_Lanczos, status,
        &iter_old, &v2[1], &v12[1], &v4[1], &v14[1], alphaCG, betaCG, &z_seed, res_save, &comm);
      free(alphaCG);
      free(betaCG);
      free(res_save);
    }/*if (childfopenALL(sdt, "rb", &fp) == 0)*/
  }/*if (X->Bind.Def.iFlgCalcSpec > RECALC_NOT)*/
  else {
    max_step = (int)X->Bind.Def.Lanczos_max;
    komega_bicg_init(&idim_max2int, &nBra, &Nomega, dcSpectrum, dcomega, &max_step, &eps_Lanczos, &comm);
  }

}/*int ReadTMComponents_BiCG*/
/**@brief
write @f$\alpha, \beta@f$, projected residual for restart
*/
int OutputTMComponents_BiCG(
  struct EDMainCalStruct *X,//!<[inout]
  int liLanczosStp//!<[in] the BiCG step
)
{
  char sdt[D_FileNameMax];
  unsigned long int stp;
  FILE *fp;
  double complex *alphaCG, *betaCG, *res_save, z_seed;

  alphaCG = (double complex*)malloc(liLanczosStp * sizeof(double complex));
  betaCG = (double complex*)malloc(liLanczosStp * sizeof(double complex));
  res_save = (double complex*)malloc(liLanczosStp * sizeof(double complex));

  komega_bicg_getcoef(alphaCG, betaCG, &z_seed, res_save);

  sprintf(sdt, cFileNameTridiagonalMatrixComponents, X->Bind.Def.CDataFileHead);
  childfopenMPI(sdt, "w", &fp);
  fprintf(fp, "%d \n", liLanczosStp);
  fprintf(fp, "%.10lf %.10lf\n", creal(z_seed), cimag(z_seed));
  for (stp = 0; stp < liLanczosStp; stp++) {
    fprintf(fp, "%25.16le %25.16le %25.16le %25.16le %25.16le %25.16le\n",
      creal(alphaCG[stp]), cimag(alphaCG[stp]),
      creal(betaCG[stp]), cimag(betaCG[stp]),
      creal(res_save[stp]), cimag(res_save[stp]));
  }
  fclose(fp);
  free(alphaCG);
  free(betaCG);
  free(res_save);

  return TRUE;
}/*int OutputTMComponents_BiCG*/
/**@brief
Initialize Shadow Residual as a random vector (Experimental)
*/
void InitShadowRes(
  struct BindStruct *X,//!<[inout]
  double complex *v4//!<[out] [CheckList::idim_max] shadow residual vector
)
{

  long int iv;
  long unsigned int idim;
  int mythread;
  double dnorm;
  /*
  For DSFMT
  */
  long unsigned int u_long_i;
  dsfmt_t dsfmt;

  iv = X->Def.initial_iv;
#pragma omp parallel default(none) private(idim, u_long_i, mythread, dsfmt) \
              shared(v4, iv, X, nthreads, myrank)
  {
    /*
     Initialise MT
    */
#ifdef _OPENMP
    mythread = omp_get_thread_num();
#else
    mythread = 0;
#endif
    u_long_i = 123432 + labs(iv) + mythread + nthreads * myrank;
    dsfmt_init_gen_rand(&dsfmt, u_long_i);

#pragma omp for
    for (idim = 1; idim <= X->Check.idim_max; idim++)
      v4[idim] = 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5)
               + 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5)*I;
  }/*#pragma omp parallel*/

  dnorm = sqrt(creal(VecProdMPI(X->Check.idim_max, v4, v4)));
#pragma omp parallel for default(none) shared(X,v4,dnorm) private(idim)
  for (idim = 1; idim <= X->Check.idim_max; idim++) v4[idim] /= dnorm;

}/*void InitShadowRes*/
/** 
 * @brief A main function to calculate spectrum by BiCG method
 * In this function, the @f$K\omega@f$ library is used.
 * The detailed procedure is written in the document of @f$K\omega@f$.
 * https://issp-center-dev.github.io/Komega/library/en/_build/html/komega_workflow_en.html#the-schematic-workflow-of-shifted-bicg-library
 * 
 * @retval 0 normally finished
 * @retval -1 error
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 * 
 */
int CalcSpectrumByBiCG(
  struct EDMainCalStruct *X,//!<[inout]
  double complex *vrhs,//!<[in] [CheckList::idim_max] Right hand side vector, excited (ket) state A|phi>.
  double complex **vlhs_Bra,//!<[in] [nBra][CheckList::idim_max] Left (bra) states B_b|phi>, projection <B_b phi|r>. Pass {vrhs} for the diagonal G_AA.
  int nBra,//!<[in] Number of left (bra) vectors. One BiCG solve yields the spectrum for every bra (Komega nl projections).
  double complex *v2,//!<[inout] [CheckList::idim_max] Work space for residual vector @f${\bf r}@f$
  double complex *v4,//!<[inout] [CheckList::idim_max] Work space for shadow residual vector @f${\bf {\tilde r}}@f$
  int Nomega,//!<[in] Number of Frequencies
  double complex *dcSpectrum,//!<[out] [nBra*Nomega] Spectrum, Fortran layout x(nBra, Nomega): dcSpectrum[iomega*nBra + ibra]
  double complex *dcomega//!<[in] [Nomega] Frequency
)
{
  char sdt[D_FileNameMax];
  unsigned long int idim, i_max;
  FILE *fp, *fp_status;
  size_t byte_size;
  int iret, ibra;
  unsigned long int liLanczosStp_vec = 0;
  double complex *v12, *v14, *res_proj;
  int stp, status[3], iomega;
  double *resz;
  double initial_residual_sq;
  int ran_bicg_loop = FALSE;
  int bicg_failed = FALSE;
  double final_shifted_max_residual = 0.0;
  double initial_residual = 0.0;
  double max_unshifted_residual = 0.0;
  int first_spike_iter = 0;
  BiCGVectorDiag iter1_v2_diag, iter1_v12_diag;
  int have_iter1_diag = FALSE;

  fprintf(stdoutMPI, "#####  Spectrum calculation with BiCG  #####\n\n");
  status[0] = 0;
  status[1] = 0;
  status[2] = 0;
  /* Defense in depth (independent of the top-level SpectrumNumBra validation): the
     tridiagonal-component restart format stores ONE projected residual stream per BiCG step,
     so multi-bra (nBra>1) is only valid for CalcSpec=Normal. Fail hard before touching any
     restart buffer or file rather than under-allocating res_save by a factor of nBra. */
  if (nBra > 1 && X->Bind.Def.iFlgCalcSpec != RECALC_NOT) {
    fprintf(stderr, "Error: multi-bra BiCG (SpectrumNumBra>1) supports only CalcSpec=\"Normal\" "
                    "(no restart/recalc); the restart format carries one residual stream per step.\n");
    return FALSE;
  }
  /**
  <ul>
  <li>Malloc vector for old residual vector (@f${\bf r}_{\rm old}@f$)
  and old shadow residual vector (@f${\bf {\tilde r}}_{\rm old}@f$).</li>
  */
  v12 = (double complex*)malloc((X->Bind.Check.idim_max + 1) * sizeof(double complex));
  v14 = (double complex*)malloc((X->Bind.Check.idim_max + 1) * sizeof(double complex));
  resz = (double*)malloc(Nomega * sizeof(double));
  /* One projected residual per bra; Komega advances all nBra spectra from a single solve. */
  res_proj = (double complex*)malloc(nBra * sizeof(double complex));
  if (v12 == NULL || v14 == NULL || resz == NULL || res_proj == NULL) {
    fprintf(stderr, "Error: out of memory in CalcSpectrumByBiCG (nBra=%d, Nomega=%d).\n", nBra, Nomega);
    free(v12); free(v14); free(resz); free(res_proj);
    return FALSE;
  }
  /**
  <li>Set initial result vector(+shadow result vector)
  Read residual vectors if restart</li>
  */
  if (X->Bind.Def.iFlgCalcSpec == RECALC_FROM_TMComponents_VEC ||
      X->Bind.Def.iFlgCalcSpec == RECALC_INOUT_TMComponents_VEC) {
    fprintf(stdoutMPI, "  Start: Input vectors for recalculation.\n");
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_InputSpectrumRecalcvecStart, "a");

    sprintf(sdt, cFileNameOutputRestartVec, X->Bind.Def.CDataFileHead, myrank);
    if (childfopenALL(sdt, "rb", &fp) != 0) {
      fprintf(stdoutMPI, "INFO: File for the restart is not found.\n");
      fprintf(stdoutMPI, "      Start from SCRATCH.\n");
#pragma omp parallel for default(none) shared(v2,v4,vrhs,X) private(idim)
      for (idim = 1; idim <= X->Bind.Check.idim_max; idim++) {
        v2[idim] = vrhs[idim];
        v4[idim] = vrhs[idim];
      }
      //InitShadowRes(&(X->Bind), v4);
    }
    else {
      byte_size = fread(&liLanczosStp_vec, sizeof(int), 1, fp);
      byte_size = fread(&i_max, sizeof(i_max), 1, fp);
      if (i_max != X->Bind.Check.idim_max) {
        fprintf(stderr, "Error: The size of the input vector is incorrect.\n");
        printf("%s %ld %ld %ld\n", sdt, i_max, X->Bind.Check.idim_max, liLanczosStp_vec);
        exitMPI(-1);
      }
      byte_size = fread(v2, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
      byte_size = fread(v12, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
      byte_size = fread(v4, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
      byte_size = fread(v14, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
      fclose(fp);
      fprintf(stdoutMPI, "  End:   Input vectors for recalculation.\n");
      TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_InputSpectrumRecalcvecEnd, "a");
      if (byte_size == 0) printf("byte_size : %d\n", (int)byte_size);
    }/*if (childfopenALL(sdt, "rb", &fp) == 0)*/
  }/*if (X->Bind.Def.iFlgCalcSpec > RECALC_FROM_TMComponents)*/
  else {
#pragma omp parallel for default(none) shared(v2,v4,vrhs,X) private(idim)
    for (idim = 1; idim <= X->Bind.Check.idim_max; idim++) {
      v2[idim] = vrhs[idim];
      v4[idim] = vrhs[idim];
    }
    //InitShadowRes(&(X->Bind), v4);
  }
  /**
  <li>Input @f$\alpha, \beta@f$, projected residual, or start from scratch</li>
  */
  ReadTMComponents_BiCG(X, v2, v4, v12, v14, nBra, Nomega, dcSpectrum, dcomega);
  initial_residual_sq = creal(VecProdMPI(X->Bind.Check.idim_max, v2, v2));
  if (isfinite(initial_residual_sq) == TRUE && initial_residual_sq >= 0.0) {
    initial_residual = sqrt(initial_residual_sq);
  } else {
    initial_residual = NAN;
  }
  max_unshifted_residual = initial_residual;
  /**
  <li>@b DO BiCG loop</li>
  <ul>
  */
  fprintf(stdoutMPI, "    Start: Calculate tridiagonal matrix components.\n");
  TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_GetTridiagonalStart, "a");
  fprintf(stdoutMPI, "\n  Iteration     Status     Seed     Residual-2-Norm\n");
  childfopenMPI("residual.dat", "w", &fp);
  childfopenMPI("bicg_status.dat", "w", &fp_status);
  PrintBiCGStatusHeader(fp_status, X, nBra, Nomega, dcomega, initial_residual);
  fflush(fp_status);

  for (stp = 1; stp <= X->Bind.Def.Lanczos_max; stp++) {
    ran_bicg_loop = TRUE;
    /**
    <li>@f${\bf v}_{2}={\hat H}{\bf v}_{12}, {\bf v}_{4}={\hat H}{\bf v}_{14}@f$,
    where @f${\bf v}_{12}, {\bf v}_{14}@f$ are old (shadow) residual vector.</li>
    */
#pragma omp parallel for default(none) shared(X,v12,v14) private(idim)
    for (idim = 1; idim <= X->Bind.Check.idim_max; idim++) {
      v12[idim] = 0.0;
      v14[idim] = 0.0;
    }
    iret = mltply(&X->Bind, v12, v2);
    if (iret == -1){
      return FALSE;
    }
    iret = mltply(&X->Bind, v14, v4);
    if (iret == -1) return FALSE;

    if (stp == 1) {
      AnalyzeBiCGVector(&v2[1], X->Bind.Check.idim_max, &iter1_v2_diag);
      AnalyzeBiCGVector(&v12[1], X->Bind.Check.idim_max, &iter1_v12_diag);
      have_iter1_diag = TRUE;
    }

    for (ibra = 0; ibra < nBra; ibra++)
      res_proj[ibra] = VecProdMPI(X->Bind.Check.idim_max, vlhs_Bra[ibra], v2);
    /**
    <li>Update projected result vector dcSpectrum (all nBra spectra at once).</li>
    */

    komega_bicg_update(&v12[1], &v2[1], &v14[1], &v4[1], dcSpectrum, res_proj, status);

    if (stp == 1 && have_iter1_diag == TRUE &&
        (status[1] == BICG_STATUS_NONFINITE || IsFiniteComplex(v12[1]) == FALSE)) {
      if (status[0] >= 0) status[0] = -stp;
      if (status[1] == 0) status[1] = BICG_STATUS_NONFINITE;
      PrintBiCGIteration1Diag(&iter1_v2_diag, &iter1_v12_diag,
        status[0], status[1], status[2], creal(v12[1]));
    }

    if (status[1] < 2) {
      komega_bicg_getresidual(resz);
    } else {
      for (iomega = 0; iomega < Nomega; iomega++) resz[iomega] = NAN;
    }

    if (isfinite(creal(v12[1])) == TRUE) {
      if (isfinite(max_unshifted_residual) == FALSE ||
          creal(v12[1]) > max_unshifted_residual) {
        max_unshifted_residual = creal(v12[1]);
      }
      if (first_spike_iter == 0 &&
          BiCGResidualRatio(creal(v12[1]), initial_residual) > BICG_SPIKE_RATIO_WARN) {
        first_spike_iter = stp;
      }
    }
    PrintBiCGStatusTrace(fp_status, stp, status, creal(v12[1]), initial_residual,
                         max_unshifted_residual, first_spike_iter, Nomega, resz);
    fflush(fp_status);

    /**
    <li>Output residuals at each frequency for some analysis.  Keep the
    historical 10-step cadence, and add the first spike/final iteration so
    rare failures carry the decisive per-frequency values without making
    large production spectrum runs write Lanczos_max*Nomega lines.</li>
    */
    if (stp % 10 == 0 || status[0] < 0 || first_spike_iter == stp) {
      for (iomega = 0; iomega < Nomega; iomega++) {
        /* Report the first bra's spectrum as a representative trace; resz is per-frequency. */
        fprintf(fp, "%7i %20.10e %20.10e %20.10e %20.10e\n",
          stp, creal(dcomega[iomega]),
          creal(dcSpectrum[iomega*nBra]), cimag(dcSpectrum[iomega*nBra]),
          resz[iomega]);
      }
      fprintf(fp, "\n");
      fflush(fp);
    }

    fprintf(stdoutMPI, "  %9d  %9d %8d %25.15e\n", abs(status[0]), status[1], status[2], creal(v12[1]));
    if (status[0] < 0) break;
  }/*for (stp = 0; stp <= X->Bind.Def.Lanczos_max; stp++)*/
  if (ran_bicg_loop == TRUE) {
    PrintBiCGResidualSummary(fp_status, "# ", abs(status[0]), status, initial_residual,
                             max_unshifted_residual, first_spike_iter);
  }
  fclose(fp);
  fclose(fp_status);
  if (ran_bicg_loop == TRUE) {
    PrintBiCGResidualSummary(stdoutMPI, "  ", abs(status[0]), status, initial_residual,
                             max_unshifted_residual, first_spike_iter);
    if (first_spike_iter > 0) {
      fprintf(stdoutMPI,
              "  WARNING: BiCG residual spike exceeded %.1e times the initial residual "
              "at iteration %d. See bicg_status.dat and residual.dat for diagnostics.\n",
              BICG_SPIKE_RATIO_WARN, first_spike_iter);
    }
  }
  if (ran_bicg_loop == TRUE && (status[0] >= 0 || status[1] != 0)) {
    bicg_failed = TRUE;
    if (status[1] >= 2) {
      final_shifted_max_residual = NAN;
    } else {
      final_shifted_max_residual = 0.0;
      komega_bicg_getresidual(resz);
      for (iomega = 0; iomega < Nomega; iomega++) {
        if (isfinite(resz[iomega]) == FALSE) {
          final_shifted_max_residual = resz[iomega];
          break;
        }
        if (resz[iomega] > final_shifted_max_residual) final_shifted_max_residual = resz[iomega];
      }
    }
  }
  /**
  </ul>
  <li>@b END @b DO BiCG loop</li>
  */
  fprintf(stdoutMPI, "    End:   Calculate tridiagonal matrix components.\n\n");
  TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_GetTridiagonalEnd, "a");

  if (bicg_failed == TRUE && status[1] >= 2) {
    fprintf(stderr,
      "Error: BiCG spectrum did not finish successfully within Lanczos_max=%u "
      "(last iteration=%d, status=%d [%s], seed=%d, max residual=%25.15e).\n",
      X->Bind.Def.Lanczos_max, abs(status[0]), status[1],
      BiCGStatusReason(status[1]), status[2], final_shifted_max_residual);
    komega_bicg_finalize();
    free(resz);
    free(res_proj);
    free(v12);
    free(v14);
    return FALSE;
  }

  /**
  <li>Save @f$\alpha, \beta@f$, projected residual</li>
  */
  /* The tridiagonal-component restart file stores ONE projected residual per step
     (komega_bicg_getcoef copies r_l_save(nl,step)); its writer is hard-coded to nl=1.
     Multi-bra (nBra>1) is Normal-only and does not support restart/recalc, so skip it
     rather than overrun the nl=1-sized buffer. */
  if (X->Bind.Def.iFlgCalcSpec != RECALC_FROM_TMComponents && nBra == 1)
    OutputTMComponents_BiCG(X, abs(status[0]));
  /**
  <li>output vectors for recalculation</li>
  </ul>
  */
  if (X->Bind.Def.iFlgCalcSpec == RECALC_OUTPUT_TMComponents_VEC ||
      X->Bind.Def.iFlgCalcSpec == RECALC_INOUT_TMComponents_VEC) {
    fprintf(stdoutMPI, "    Start: Output vectors for recalculation.\n");
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_OutputSpectrumRecalcvecStart, "a");

    komega_bicg_getvec(&v12[1], &v14[1]);

    sprintf(sdt, cFileNameOutputRestartVec, X->Bind.Def.CDataFileHead, myrank);
    if (childfopenALL(sdt, "wb", &fp) != 0) {
      exitMPI(-1);
    }
    byte_size = fwrite(&status[0], sizeof(status[0]), 1, fp);
    byte_size = fwrite(&X->Bind.Check.idim_max, sizeof(X->Bind.Check.idim_max), 1, fp);
    byte_size = fwrite(v2, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
    byte_size = fwrite(v12, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
    byte_size = fwrite(v4, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
    byte_size = fwrite(v14, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
    fclose(fp);

    fprintf(stdoutMPI, "    End:   Output vectors for recalculation.\n");
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_OutputSpectrumRecalcvecEnd, "a");
  }/*if (X->Bind.Def.iFlgCalcSpec > RECALC_FROM_TMComponents)*/

  if (bicg_failed == TRUE) {
    fprintf(stderr,
      "Error: BiCG spectrum did not finish successfully within Lanczos_max=%u "
      "(last iteration=%d, status=%d [%s], seed=%d, max residual=%25.15e).\n",
      X->Bind.Def.Lanczos_max, abs(status[0]), status[1],
      BiCGStatusReason(status[1]), status[2], final_shifted_max_residual);
    komega_bicg_finalize();
    free(resz);
    free(res_proj);
    free(v12);
    free(v14);
    return FALSE;
  }

  komega_bicg_finalize();

  free(resz);
  free(res_proj);
  free(v12);
  free(v14);
  return TRUE;
}/*int CalcSpectrumByBiCG*/
