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
#ifdef MPI
#include <mpi.h>
#endif
/**@brief
Read @f$\alpha, \beta@f$, projected residual for restart
*/
void ReadTMComponents_BiCG(
  struct EDMainCalStruct *X,//!<[inout]
  double complex *v2,//!<[inout] [CheckList::idim_max] Residual vector
  double complex *v4,//!<[inout] [CheckList::idim_max] Shadow esidual vector
  double complex *v12,//!<[inout] [CheckList::idim_max] Old residual vector
  double complex *v14,//!<[inout] [CheckList::idim_max] Old shadow residual vector
  int Nomega,//!<[in] Number of frequencies
  double complex *dcSpectrum,//!<[inout] [Nomega] Projected result vector, spectrum
  double complex *dcomega//!<[in] [Nomega] Frequency
) {
  char sdt[D_FileNameMax];
  char ctmp[256];

  int one = 1, status[3], idim_max2int, max_step, iter_old;
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
      komega_bicg_init(&idim_max2int, &one, &Nomega, dcSpectrum, dcomega, &max_step, &eps_Lanczos, &comm);
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

      komega_bicg_restart(&idim_max2int, &one, &Nomega, dcSpectrum, dcomega, &max_step, &eps_Lanczos, status,
        &iter_old, &v2[1], &v12[1], &v4[1], &v14[1], alphaCG, betaCG, &z_seed, res_save, &comm);
      free(alphaCG);
      free(betaCG);
      free(res_save);
    }/*if (childfopenALL(sdt, "rb", &fp) == 0)*/
  }/*if (X->Bind.Def.iFlgCalcSpec > RECALC_NOT)*/
  else {
    max_step = (int)X->Bind.Def.Lanczos_max;
    komega_bicg_init(&idim_max2int, &one, &Nomega, dcSpectrum, dcomega, &max_step, &eps_Lanczos, &comm);
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
  double complex *vrhs,//!<[in] [CheckList::idim_max] Right hand side vector, excited state.
  double complex *v2,//!<[inout] [CheckList::idim_max] Work space for residual vector @f${\bf r}@f$
  double complex *v4,//!<[inout] [CheckList::idim_max] Work space for shadow residual vector @f${\bf {\tilde r}}@f$
  int Nomega,//!<[in] Number of Frequencies
  double complex *dcSpectrum,//!<[out] [Nomega] Spectrum
  double complex *dcomega//!<[in] [Nomega] Frequency
)
{
  char sdt[D_FileNameMax];
  unsigned long int idim, i_max;
  FILE *fp;
  size_t byte_size;
  int iret;
  unsigned long int liLanczosStp_vec = 0;
  double complex *v12, *v14, res_proj;
  int stp, status[3], iomega;
  double *resz;

  fprintf(stdoutMPI, "#####  Spectrum calculation with BiCG  #####\n\n");
  /**
  <ul>
  <li>Malloc vector for old residual vector (@f${\bf r}_{\rm old}@f$)
  and old shadow residual vector (@f${\bf {\tilde r}}_{\rm old}@f$).</li>
  */
  v12 = (double complex*)malloc((X->Bind.Check.idim_max + 1) * sizeof(double complex));
  v14 = (double complex*)malloc((X->Bind.Check.idim_max + 1) * sizeof(double complex));
  resz = (double*)malloc(Nomega * sizeof(double));
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
  ReadTMComponents_BiCG(X, v2, v4, v12, v14, Nomega, dcSpectrum, dcomega);
  /**
  <li>@b DO BiCG loop</li>
  <ul>
  */
  fprintf(stdoutMPI, "    Start: Calculate tridiagonal matrix components.\n");
  TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_GetTridiagonalStart, "a");
  fprintf(stdoutMPI, "\n  Iteration     Status     Seed     Residual-2-Norm\n");
  childfopenMPI("residual.dat", "w", &fp);

  for (stp = 1; stp <= X->Bind.Def.Lanczos_max; stp++) {
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

    res_proj = VecProdMPI(X->Bind.Check.idim_max, vrhs, v2);
    /**
    <li>Update projected result vector dcSpectrum.</li>
    */

    komega_bicg_update(&v12[1], &v2[1], &v14[1], &v4[1], dcSpectrum, &res_proj, status);

    /**
    <li>Output residuals at each frequency for some analysis</li>
    */
    if (stp % 10 == 0) {
      komega_bicg_getresidual(resz);

      for (iomega = 0; iomega < Nomega; iomega++) {
        fprintf(fp, "%7i %20.10e %20.10e %20.10e %20.10e\n", 
          stp, creal(dcomega[iomega]), 
          creal(dcSpectrum[iomega]), cimag(dcSpectrum[iomega]),
          resz[iomega]);
      }
      fprintf(fp, "\n");
    }

    fprintf(stdoutMPI, "  %9d  %9d %8d %25.15e\n", abs(status[0]), status[1], status[2], creal(v12[1]));
    if (status[0] < 0) break;
  }/*for (stp = 0; stp <= X->Bind.Def.Lanczos_max; stp++)*/
  fclose(fp);
  /**
  </ul>
  <li>@b END @b DO BiCG loop</li>
  */
  fprintf(stdoutMPI, "    End:   Calculate tridiagonal matrix components.\n\n");
  TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_GetTridiagonalEnd, "a");
  /**
  <li>Save @f$\alpha, \beta@f$, projected residual</li>
  */
  if (X->Bind.Def.iFlgCalcSpec != RECALC_FROM_TMComponents)
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

  komega_bicg_finalize();

  free(resz);
  free(v12);
  free(v14);
  return TRUE;
}/*int CalcSpectrumByBiCG*/
