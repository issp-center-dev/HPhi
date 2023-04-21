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
#include "FileIO.h"
#include "wrapperMPI.h"
#include "common/setmemory.h"
#include "mltply.h"
#include "CalcSpectrum.h"
#include "mltplyCommon.h"
#ifdef MPI
#include <mpi.h>
#endif
/*
*@brief Solve Shifted equation
*/
void ShiftedEq(
  int iter,
  int Nomega,
  int NdcSpectrum,
  int* lz_conv,
  double complex* alpha,
  double complex* beta,
  double complex* dcomega,
  double complex z_seed,
  double complex** pBiCG,
  double complex** res_proj,
  double complex** pi,
  double complex** dcSpectrum
) {
  int iomega, idcSpectrum;
  double complex pi_2;

  for (iomega = 0; iomega < Nomega; iomega++) {

    if (lz_conv[iomega] == 1) {
      pi[iter][iomega] = pi[iter - 1][iomega];
      continue;
    }

    if (iter == 1)
      pi_2 = 1.0;
    else
      pi_2 = pi[iter - 2][iomega];

    pi[iter][iomega] = (1.0 + alpha[iter] * (dcomega[iomega] - z_seed)) * pi[iter - 1][iomega]
      - alpha[iter] * beta[iter] / alpha[iter - 1] * (pi_2 - pi[iter - 1][iomega]);
    for (idcSpectrum = 0; idcSpectrum < NdcSpectrum; idcSpectrum++) {
      pBiCG[iomega][idcSpectrum] = res_proj[iter][idcSpectrum] / pi[iter - 1][iomega]
        + (pi_2 / pi[iter - 1][iomega]) * (pi_2 / pi[iter - 1][iomega]) * beta[iter] * pBiCG[iomega][idcSpectrum];
      dcSpectrum[iomega][idcSpectrum] +=
        pi[iter - 1][iomega] / pi[iter][iomega] * alpha[iter] * pBiCG[iomega][idcSpectrum];
    }
  }/*for (iomega = 0; iomega < Nomega; iomega++)*/
}
/**
@brief Perform Seed Switch
*/
void SeedSwitch(
  int istate,
  int iter,
  int Nomega,
  int NdcSpectrum,
  int* lz_conv,
  int* iz_seed,
  double complex* z_seed,
  double complex* rho,
  double complex* dcomega,
  long int ndim,
  double complex** v2,
  double complex** v3,
  double complex** v4,
  double complex** v5,
  double complex** pi,
  double complex* alpha,
  double complex* beta,
  double complex** res_proj
) {
  double pi_min;
  double complex pi_seed;
  int iz_seed0, iomega, jter, idcSpectrum;
  long int idim;
  //
  // Initialize for min
  //
  iz_seed0 = -1;
  for (iomega = 0; iomega < Nomega; iomega++)
    if (lz_conv[iomega] == 0) {
      iz_seed0 = iomega;
      pi_min = cabs(pi[iter][iz_seed0]);
    }
  if (iz_seed0 == -1) return;
  //
  // Search min.
  //
  for (iomega = 0; iomega < Nomega; iomega++) {
    if (lz_conv[iomega] == 0)
      if (cabs(pi[iter][iomega]) < pi_min) {
        iz_seed0 = iomega;
        pi_min = cabs(pi[iter][iomega]);
      }
  }/*for (iomega = 0; iomega < Nomega; iomega++)*/

  if (cabs(pi[iter][iz_seed0]) < 1.0e-50) {
    printf("Error : pi at seed (%d) is 0.", iz_seed0);
    exitMPI(-1);
  }

  if (iz_seed0 != *iz_seed) {

    *iz_seed = iz_seed0;
    *z_seed = dcomega[iz_seed0];

    *rho /= (pi[iter - 1][iz_seed0] * pi[iter - 1][iz_seed0]);

    for (idim = 1; idim <= ndim; idim++) {
      v2[idim][istate] /= pi[iter][iz_seed0];
      v4[idim][istate] /= conj(pi[iter][iz_seed0]);
      v3[idim][istate] /= pi[iter - 1][iz_seed0];
      v5[idim][istate] /= conj(pi[iter - 1][iz_seed0]);
    }
    /*
    For restarting
    */
    for (jter = 1; jter <= iter; jter++) {
      alpha[jter] *= pi[jter - 1][iz_seed0] / pi[jter][iz_seed0];
      if (jter != 1)
        beta[jter] *= (pi[jter - 2][iz_seed0] / pi[jter - 1][iz_seed0])* (pi[jter - 2][iz_seed0] / pi[jter - 1][iz_seed0]);
      for (idcSpectrum = 0; idcSpectrum < NdcSpectrum; idcSpectrum++) {
        res_proj[jter][idcSpectrum] /= pi[jter - 1][iz_seed0];
      }
    }

    for (jter = 1; jter <= iter; jter++) {
      pi_seed = pi[jter][iz_seed0];
      for (iomega = 0; iomega < Nomega; iomega++)
        pi[jter][iomega] /= pi_seed;
    }
  }

}/*void SeedSwitch*/
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
  double complex **v2,//!<[in] [CheckList::idim_max] Right hand side vector, excited state.
  double complex **v4,//!<[inout] [CheckList::idim_max] Work space for residual vector @f${\bf r}@f$
  int Nomega,//!<[in] Number of Frequencies
  int NdcSpectrum,
  double complex ***dcSpectrum,//!<[out] [Nomega] Spectrum
  double complex **dcomega,//!<[in] [Nomega] Frequency
  double complex **v1Org
)
{
  char sdt[D_FileNameMax], ctmp[256];
  unsigned long int idim, i_max;
  FILE *fp;
  size_t byte_size;
  int idcSpectrum, iter_old;
  unsigned long int liLanczosStp_vec = 0;
  double complex **vL, **v12, **v14, **v3, **v5, *z_seed, *rho, *rho_old,
    ** alpha, ** beta, *** res_proj, *** pi, *** pBiCG, alpha_denom;
  int stp, iomega, istate, *iz_seed, ** lz_conv, *lz_conv_state, lz_conv_all;
  double resz, *resnorm, dtmp[4];

  fprintf(stdoutMPI, "#####  Spectrum calculation with BiCG  #####\n\n");
  /**
  <ul>
  <li>Malloc vector for old residual vector (@f${\bf r}_{\rm old}@f$)
  and old shadow residual vector (@f${\bf {\tilde r}}_{\rm old}@f$).</li>
  */
  z_seed = cd_1d_allocate(X->Bind.Def.k_exct);
  iz_seed = i_1d_allocate(X->Bind.Def.k_exct);
  rho = cd_1d_allocate(X->Bind.Def.k_exct);
  rho_old = cd_1d_allocate(X->Bind.Def.k_exct);
  resnorm = d_1d_allocate(X->Bind.Def.k_exct);
  for (istate = 0; istate < X->Bind.Def.k_exct; istate++) {
    iz_seed[istate] = 0;
    z_seed[istate] = dcomega[istate][iz_seed[istate]];
    rho[istate] = 1.0;
  }
  pBiCG = cd_3d_allocate(X->Bind.Def.k_exct, Nomega, NdcSpectrum);
  v3 =  cd_2d_allocate(X->Bind.Check.idim_max + 1, X->Bind.Def.k_exct);
  v5 =  cd_2d_allocate(X->Bind.Check.idim_max + 1, X->Bind.Def.k_exct);
  v12 = cd_2d_allocate(X->Bind.Check.idim_max + 1, X->Bind.Def.k_exct);
  v14 = cd_2d_allocate(X->Bind.Check.idim_max + 1, X->Bind.Def.k_exct);
  vL = cd_2d_allocate(X->Bind.Check.idim_max + 1, X->Bind.Def.k_exct);
  lz_conv = i_2d_allocate(X->Bind.Def.k_exct, Nomega);
  lz_conv_state = i_1d_allocate(X->Bind.Def.k_exct);
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
      zclear(X->Bind.Check.idim_max, &v2[1][0]);
      GetExcitedState(&(X->Bind), 1, v2, v1Org, 0);
#pragma omp parallel for default(none) shared(v2,v4,v1Org,X) private(idim)
      for (idim = 1; idim <= X->Bind.Check.idim_max; idim++) 
        v4[idim][0] = v2[idim][0];
    }
    else {
      byte_size = fread(&liLanczosStp_vec, sizeof(int), 1, fp);
      byte_size = fread(&i_max, sizeof(i_max), 1, fp);
      if (i_max != X->Bind.Check.idim_max) {
        fprintf(stderr, "Error: The size of the input vector is incorrect.\n");
        printf("%s %ld %ld %ld\n", sdt, i_max, X->Bind.Check.idim_max, liLanczosStp_vec);
        exitMPI(-1);
      }
      byte_size = fread(&v2[0][0], sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
      byte_size = fread(&v12[0][0], sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
      byte_size = fread(&v4[0][0], sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
      byte_size = fread(&v14[0][0], sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
      fclose(fp);
      fprintf(stdoutMPI, "  End:   Input vectors for recalculation.\n");
      TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_InputSpectrumRecalcvecEnd, "a");
      if (byte_size == 0) printf("byte_size : %d\n", (int)byte_size);
    }/*if (childfopenALL(sdt, "rb", &fp) == 0)*/
  }/*if (X->Bind.Def.iFlgCalcSpec > RECALC_FROM_TMComponents)*/
  else {
    zclear(X->Bind.Check.idim_max, &v2[1][0]);
    GetExcitedState(&(X->Bind), X->Bind.Def.k_exct, v2, v1Org, 0);
#pragma omp parallel for default(none) shared(v2,v4,v1Org,X) private(idim,istate)
    for (idim = 1; idim <= X->Bind.Check.idim_max; idim++)
      for (istate = 0;istate < X->Bind.Def.k_exct;istate++)
        v4[idim][istate] = v2[idim][istate];
  }
  /**
  <li>Input @f$\alpha, \beta@f$, projected residual, or start from scratch</li>
  */
  if (X->Bind.Def.iFlgCalcSpec == RECALC_FROM_TMComponents ||
    X->Bind.Def.iFlgCalcSpec == RECALC_FROM_TMComponents_VEC ||
    X->Bind.Def.iFlgCalcSpec == RECALC_INOUT_TMComponents_VEC) {
    sprintf(sdt, cFileNameTridiagonalMatrixComponents, X->Bind.Def.CDataFileHead);
    if (childfopenALL(sdt, "rb", &fp) != 0) {
      fprintf(stdoutMPI, "INFO: File for the restart is not found.\n");
      fprintf(stdoutMPI, "      Start from SCRATCH.\n");
    }
    else {
      fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);
      sscanf(ctmp, "%d", &iter_old);
      if (X->Bind.Def.iFlgCalcSpec > RECALC_FROM_TMComponents) {
        alpha = cd_2d_allocate(X->Bind.Def.k_exct, iter_old + X->Bind.Def.Lanczos_max);
        beta = cd_2d_allocate(X->Bind.Def.k_exct, iter_old + X->Bind.Def.Lanczos_max);
        res_proj = cd_3d_allocate(X->Bind.Def.k_exct, iter_old + X->Bind.Def.Lanczos_max, NdcSpectrum);
        pi = cd_3d_allocate(X->Bind.Def.k_exct, iter_old + X->Bind.Def.Lanczos_max, Nomega);
      }
      else {
        alpha = cd_2d_allocate(X->Bind.Def.k_exct, iter_old);
        beta = cd_2d_allocate(X->Bind.Def.k_exct, iter_old);
        res_proj = cd_3d_allocate(X->Bind.Def.k_exct, iter_old, NdcSpectrum);
        pi = cd_3d_allocate(X->Bind.Def.k_exct, iter_old, Nomega);
      }
      fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);
      for (istate = 0; istate < X->Bind.Def.k_exct; istate++) {
        sscanf(ctmp, "%lf %lf\n", &dtmp[0], &dtmp[1]);
        z_seed[istate] = dtmp[0] + I * dtmp[1];
      }

      for (istate = 0; istate < X->Bind.Def.k_exct; istate++) {
        for (stp = 0; stp < iter_old; stp++) {
          fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);
          sscanf(ctmp, "%lf %lf %lf %lf\n",
            &dtmp[0], &dtmp[1], &dtmp[2], &dtmp[3]);
          alpha[istate][stp] = dtmp[0] + I * dtmp[1];
          beta[istate][stp] = dtmp[2] + I * dtmp[3];
          for (idcSpectrum = 0; idcSpectrum < NdcSpectrum; idcSpectrum++) {
            fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);
            sscanf(ctmp, "%lf %lf\n", &dtmp[0], &dtmp[1]);
            res_proj[istate][stp][idcSpectrum] = dtmp[0] + I * dtmp[1];
          }
        }
      }
      fclose(fp);

      for (stp = 1; stp <= iter_old; stp++)
        for (istate = 0; istate < X->Bind.Def.k_exct; istate++)
          ShiftedEq(stp, Nomega, NdcSpectrum, lz_conv[istate], alpha[istate], beta[istate], dcomega[istate],
            z_seed[istate], pBiCG[istate], res_proj[istate], pi[istate], dcSpectrum[istate]);

      MultiVecProdMPI(X->Bind.Check.idim_max, X->Bind.Def.k_exct, v5, v3, rho);

      for (istate = 0; istate < X->Bind.Def.k_exct; istate++) {
        SeedSwitch(istate, stp, Nomega, NdcSpectrum, lz_conv[istate], &iz_seed[istate],
          &z_seed[istate], &rho[istate], dcomega[istate], X->Bind.Def.k_exct, v2, v3, v4, v5,
          pi[istate], alpha[istate], beta[istate], res_proj[istate]);
      }

      resnorm = d_1d_allocate(X->Bind.Def.k_exct);
      NormMPI_dv(X->Bind.Check.idim_max, X->Bind.Def.k_exct, v2, resnorm);

      for (istate = 0; istate < X->Bind.Def.k_exct; istate++)
        for (iomega = 0; iomega < Nomega; iomega++)
          if (fabs(resnorm[istate] / pi[istate][stp][iomega]) < eps_Lanczos)
            lz_conv[istate][idcSpectrum] = 1;
      free_d_1d_allocate(resnorm);

      if (X->Bind.Def.iFlgCalcSpec == RECALC_FROM_TMComponents) X->Bind.Def.Lanczos_max = 0;

    }/*if (childfopenALL(sdt, "rb", &fp) == 0)*/
  }/*if (X->Bind.Def.iFlgCalcSpec > RECALC_NOT)*/
  else {
    iter_old = 0;
    alpha = cd_2d_allocate(X->Bind.Def.k_exct, X->Bind.Def.Lanczos_max);
    beta = cd_2d_allocate(X->Bind.Def.k_exct, X->Bind.Def.Lanczos_max);
    res_proj = cd_3d_allocate(X->Bind.Def.k_exct, X->Bind.Def.Lanczos_max, NdcSpectrum);
    pi = cd_3d_allocate(X->Bind.Def.k_exct, X->Bind.Def.Lanczos_max, Nomega);
    for (istate = 0; istate < X->Bind.Def.k_exct; istate++) {
      alpha[istate][0] = 1.0;
      beta[istate][0] = 0.0;
      for (iomega = 0; iomega < Nomega; iomega++) {
        pi[istate][0][iomega] = 1.0;
        for (idcSpectrum = 0; idcSpectrum < NdcSpectrum; idcSpectrum++) {
          pBiCG[istate][iomega][idcSpectrum] = 0.0;
          dcSpectrum[istate][iomega][idcSpectrum] = 0.0;
        }
      }/*for (iomega = 0; iomega < Nomega; iomega++)*/
    }/*for (istate = 0; istate < nstate; istate++)*/
  }
  /**
  <li>@b DO BiCG loop</li>
  <ul>
  */
  fprintf(stdoutMPI, "    Start: Calculate tridiagonal matrix components.\n");
  TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_GetTridiagonalStart, "a");
  fprintf(stdoutMPI, "\n  Iteration     Status     Seed     Residual-2-Norm\n");
  childfopenMPI("residual.dat", "w", &fp);

  for (stp = iter_old + 1; stp <= iter_old + X->Bind.Def.Lanczos_max; stp++) {
    /**
    <li>@f${\bf v}_{2}={\hat H}{\bf v}_{12}, {\bf v}_{4}={\hat H}{\bf v}_{14}@f$,
    where @f${\bf v}_{12}, {\bf v}_{14}@f$ are old (shadow) residual vector.</li>
    */
    zclear(X->Bind.Check.idim_max * X->Bind.Def.k_exct, &v12[1][0]);
    zclear(X->Bind.Check.idim_max * X->Bind.Def.k_exct, &v14[1][0]);
    mltply(&X->Bind, X->Bind.Def.k_exct, v12, v2);
    mltply(&X->Bind, X->Bind.Def.k_exct, v14, v4);
    
    for (idcSpectrum = 0; idcSpectrum < NdcSpectrum; idcSpectrum++) {
      zclear(X->Bind.Check.idim_max * X->Bind.Def.k_exct, &vL[1][0]);
      GetExcitedState(&(X->Bind), X->Bind.Def.k_exct, vL, v1Org, idcSpectrum + 1);
      MultiVecProdMPI(X->Bind.Check.idim_max, X->Bind.Def.k_exct, vL, v2, rho_old);
      for (istate = 0; istate < X->Bind.Def.k_exct; istate++)res_proj[istate][stp][idcSpectrum] = rho_old[istate];
    }
    /**
    <li>Update projected result vector dcSpectrum.</li>
    */
    for (istate = 0; istate < X->Bind.Def.k_exct; istate++)rho_old[istate] = rho[istate];
    MultiVecProdMPI(X->Bind.Check.idim_max, X->Bind.Def.k_exct, v4, v2, rho);
    for (istate = 0; istate < X->Bind.Def.k_exct; istate++) {
      lz_conv_all *= lz_conv[istate][iomega];

      if (stp == 1)
        beta[istate][stp] = 0.0;
      else
        beta[istate][stp] = rho[istate] / rho_old[istate];

      for (idim = 1; idim <= X->Bind.Check.idim_max; idim++) {
        v12[idim][istate] = z_seed[istate] * v2[idim][istate] - v12[idim][istate];
        v14[idim][istate] = conj(z_seed[istate]) * v4[idim][istate] - v14[idim][istate];
      }
    }/*for (istate = 0; istate < nstate; istate++)*/

    MultiVecProdMPI(X->Bind.Check.idim_max, X->Bind.Def.k_exct, v4, v12, rho_old);

    for (istate = 0; istate < X->Bind.Def.k_exct; istate++) {
      if (lz_conv_state[istate] == 1) continue;

      alpha_denom = rho_old[istate] - beta[istate][stp] * rho[istate] / alpha[istate][stp - 1];

      if (cabs(alpha_denom) < 1.0e-50) {
        printf("Error : The denominator of alpha is zero.\n");
        exitMPI(-1);
      }
      else if (cabs(rho[istate]) < 1.0e-50) {
        printf("Error : rho is zero.\n");
        exitMPI(-1);
      }
      alpha[istate][stp] = rho[istate] / alpha_denom;
      /*
      Shifted equation
      */
      ShiftedEq(stp, Nomega, NdcSpectrum, lz_conv[istate], alpha[istate], beta[istate], dcomega[istate],
        z_seed[istate], pBiCG[istate], res_proj[istate], pi[istate], dcSpectrum[istate]);
      /*
      Update residual
      */
      for (idim = 1; idim <= X->Bind.Check.idim_max; idim++) {
        v12[idim][istate] = (1.0 + alpha[istate][stp] * beta[istate][stp] / alpha[istate][stp - 1]) * v2[idim][istate]
          - alpha[istate][stp] * v12[idim][istate]
          - alpha[istate][stp] * beta[istate][stp] / alpha[istate][stp - 1] * v3[idim][istate];
        v3[idim][istate] = v2[idim][istate];
        v2[idim][istate] = v12[idim][istate];
        v14[idim][istate] = (1.0 + conj(alpha[istate][stp] * beta[istate][stp] / alpha[istate][stp - 1])) * v4[idim][istate]
          - conj(alpha[istate][stp]) * v14[idim][istate]
          - conj(alpha[istate][stp] * beta[istate][stp] / alpha[istate][stp - 1]) * v5[idim][istate];
        v5[idim][istate] = v4[idim][istate];
        v4[idim][istate] = v14[idim][istate];
      }/*for (idim = 0; idim < Check::idim_maxs; idim++)*/
      /*
      Seed Switching
      */
      SeedSwitch(istate, stp, Nomega, NdcSpectrum, lz_conv[istate], &iz_seed[istate],
        &z_seed[istate], &rho[istate], dcomega[istate], X->Bind.Check.idim_max, v2, v3, v4, v5,
        pi[istate], alpha[istate], beta[istate], res_proj[istate]);
    }/*for (istate = 0; istaet < nstate; istate++)*/
    /*
    Convergence check
    */
    NormMPI_dv(X->Bind.Check.idim_max, X->Bind.Def.k_exct, v2, resnorm);

    lz_conv_all = 1;
    fprintf(stdoutMPI, "  %9d  ", stp);
    for (istate = 0; istate < X->Bind.Def.k_exct; istate++) {
      lz_conv_state[istate] = 1;
      for (iomega = 0; iomega < Nomega; iomega++) {
        if (lz_conv[istate][iomega] == 0)
          if(fabs(resnorm[istate] / pi[istate][stp][iomega]) < eps_Lanczos)
            lz_conv[istate][iomega] = 1;
        lz_conv_state[istate] *= lz_conv[istate][iomega];
      }
      lz_conv_all *= lz_conv_state[istate];

      fprintf(stdoutMPI, "%9d %25.15e", iz_seed[istate], resnorm[istate]);
    }/*for (istate = 0; istate < nstate; istate++)*/
    fprintf(stdoutMPI, "\n");

    if (lz_conv_all == 1) break;

    /**
    <li>Output residuals at each frequency for some analysis</li>
    */
    if (stp % 10 == 0) {

      for (iomega = 0; iomega < Nomega; iomega++) {
        fprintf(fp, "%7i ", stp);
        for (istate = 0; istate < X->Bind.Def.k_exct; istate++) {
          resz = resnorm[istate] / cabs(pi[istate][stp][iomega]);//FIXME

          fprintf(fp, "%20.10e %20.10e %20.10e %20.10e ",
            creal(dcomega[istate][iomega]),
            creal(dcSpectrum[istate][iomega][0]), cimag(dcSpectrum[istate][iomega][0]),resz);
        }
        fprintf(fp, "\n");
      }
      fprintf(fp, "\n");
    }
    
  }/*for (stp = 0; stp <= X->Bind.Def.Lanczos_max; stp++)*/
  fclose(fp);

  if (stp >= iter_old + X->Bind.Def.Lanczos_max)
    fprintf(stdoutMPI, "Remark : Not converged in iteration %d.", stp);
  iter_old = stp;
  /**
  </ul>
  <li>@b END @b DO BiCG loop</li>
  */
  fprintf(stdoutMPI, "    End:   Calculate tridiagonal matrix components.\n\n");
  TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_GetTridiagonalEnd, "a");
  /**
  <li>Save @f$\alpha, \beta@f$, projected residual</li>
  */
  if (X->Bind.Def.iFlgCalcSpec != RECALC_FROM_TMComponents) {
    sprintf(sdt, cFileNameTridiagonalMatrixComponents, X->Bind.Def.CDataFileHead);
    childfopenMPI(sdt, "w", &fp);
    fprintf(fp, "%d \n", iter_old);
    for (istate = 0; istate < X->Bind.Def.k_exct; istate++)
      fprintf(fp, "%.10lf %.10lf\n", creal(z_seed[istate]), cimag(z_seed[istate]));
    for (istate = 0; istate < X->Bind.Def.k_exct; istate++) {
      for (stp = 0; stp < iter_old; stp++) {
        fprintf(fp, "%25.16le %25.16le %25.16le %25.16le\n",
          creal(alpha[istate][stp]), cimag(alpha[istate][stp]),
          creal(beta[istate][stp]), cimag(beta[istate][stp]));
        for (idcSpectrum = 0; idcSpectrum < NdcSpectrum; idcSpectrum++) {
          fprintf(fp, "%25.16le %25.16le\n",
            creal(res_proj[istate][stp][idcSpectrum]), cimag(res_proj[istate][stp][idcSpectrum]));
        }
      }
    }
    fclose(fp);
  }
  /**
  <li>output vectors for recalculation</li>
  </ul>
  */
  if (X->Bind.Def.iFlgCalcSpec == RECALC_OUTPUT_TMComponents_VEC ||
      X->Bind.Def.iFlgCalcSpec == RECALC_INOUT_TMComponents_VEC) {
    fprintf(stdoutMPI, "    Start: Output vectors for recalculation.\n");
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_OutputSpectrumRecalcvecStart, "a");

    sprintf(sdt, cFileNameOutputRestartVec, X->Bind.Def.CDataFileHead, myrank);
    if (childfopenALL(sdt, "wb", &fp) != 0) {
      exitMPI(-1);
    }
    byte_size = fwrite(&iter_old, sizeof(iter_old), 1, fp);
    byte_size = fwrite(&X->Bind.Check.idim_max, sizeof(X->Bind.Check.idim_max), 1, fp);
    byte_size = fwrite(&v2[0][0], sizeof(complex double), (X->Bind.Check.idim_max + 1)* X->Bind.Def.k_exct, fp);
    byte_size = fwrite(&v12[0][0], sizeof(complex double), (X->Bind.Check.idim_max + 1) * X->Bind.Def.k_exct, fp);
    byte_size = fwrite(&v4[0][0], sizeof(complex double), (X->Bind.Check.idim_max + 1) * X->Bind.Def.k_exct, fp);
    byte_size = fwrite(&v14[0][0], sizeof(complex double), (X->Bind.Check.idim_max + 1) * X->Bind.Def.k_exct, fp);
    fclose(fp);

    fprintf(stdoutMPI, "    End:   Output vectors for recalculation.\n");
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_OutputSpectrumRecalcvecEnd, "a");
  }/*if (X->Bind.Def.iFlgCalcSpec > RECALC_FROM_TMComponents)*/

  free_cd_1d_allocate(z_seed);
  free_i_1d_allocate(iz_seed);
  free_cd_1d_allocate(rho);
  free_cd_1d_allocate(rho_old);
  free_d_1d_allocate(resnorm);
  free_cd_3d_allocate(pBiCG);
  free_cd_2d_allocate(v3);
  free_cd_2d_allocate(v5);
  free_cd_2d_allocate(v12);
  free_cd_2d_allocate(v14);
  free_cd_2d_allocate(vL);
  free_i_2d_allocate(lz_conv);
  free_i_1d_allocate(lz_conv_state);
  free_cd_2d_allocate(alpha);
  free_cd_2d_allocate(beta);
  free_cd_3d_allocate(pi);
  free_cd_3d_allocate(res_proj);
  return TRUE;
}/*int CalcSpectrumByBiCG*/
