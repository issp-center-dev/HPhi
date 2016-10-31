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
#include "CalcSpectrumByLanczos.h"
#include "Lanczos_EigenValue.h"
#include "FileIO.h"
#include "wrapperMPI.h"
#include "mfmemory.h"
#ifdef MPI
#include "komega/pkomega_bicg.h"
#else
#include "komega/komega_bicg.h"
#endif
#include "mltply.h"

/**
 * @file   CalcSpectrumByBiCG.c
 * @version 1.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 * @brief  File for givinvg functions of calculating spectrum by Lanczos
 * 
 * 
 */


int ReadTMComponents_BiCG(
  struct EDMainCalStruct *X,
  double complex *v2,
  double complex *v4,
  double complex *v12,
  double complex *v14,
  int Nomega,
  double complex *dcSpectrum,
  double complex *dcomega
) {
  char sdt[D_FileNameMax];
  char ctmp[256];

  int one = 1, status[3];
  unsigned long int idx;
  unsigned long int liLanczosStp, liLanczosStp2;
  double complex *alphaCG, *betaCG, *res_save, z_seed;
  double z_seed_r, z_seed_i, alpha_r, alpha_i, beta_r, beta_i, res_r, res_i;
  FILE *fp;
#ifdef MPI
  MPI_Comm comm = MPI_COMM_WORLD;
#endif

  idx = 0;
  sprintf(sdt, cFileNameTridiagonalMatrixComponents, X->Bind.Def.CDataFileHead);
  childfopenMPI(sdt, "r", &fp);

  fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);
  sscanf(ctmp, "%ld \n", &liLanczosStp);
  if (X->Bind.Def.iFlgCalcSpec == RECALC_INOUT_TMComponents_VEC ||
    X->Bind.Def.iFlgCalcSpec == RECALC_FROM_TMComponents_VEC) {
    alphaCG = (double complex*)malloc((liLanczosStp + X->Bind.Def.Lanczos_max) * sizeof(double complex));
    betaCG = (double complex*)malloc((liLanczosStp + X->Bind.Def.Lanczos_max) * sizeof(double complex));
    res_save = (double complex*)malloc((liLanczosStp + X->Bind.Def.Lanczos_max) * sizeof(double complex));
  }
  else if (X->Bind.Def.iFlgCalcSpec == RECALC_FROM_TMComponents) {
    alphaCG = (double complex*)malloc(liLanczosStp * sizeof(double complex));
    betaCG = (double complex*)malloc(liLanczosStp * sizeof(double complex));
    res_save = (double complex*)malloc(liLanczosStp * sizeof(double complex));
  }
  fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);
  sscanf(ctmp, "%lf %lf\n", &z_seed_r, &z_seed_i);
  z_seed = z_seed_r + I * z_seed_i;
  while (fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp) != NULL) {
    sscanf(ctmp, "%lf %lf %lf %lf %lf %lf\n",
      &alpha_r, &alpha_i, &beta_r, &beta_i, &res_r, &res_i);
    alphaCG[idx] = alpha_r + I * alpha_i;
    betaCG[idx] = beta_r + I * beta_i;
    res_save[idx] = res_r + I * res_i;
    idx += 1;
  }
  fclose(fp);

  if (X->Bind.Def.iFlgCalcSpec == RECALC_FROM_TMComponents) {
    X->Bind.Def.Lanczos_restart = 0;
  }
  else if (X->Bind.Def.iFlgCalcSpec == RECALC_INOUT_TMComponents_VEC ||
    X->Bind.Def.iFlgCalcSpec == RECALC_FROM_TMComponents_VEC) {
    X->Bind.Def.Lanczos_restart = liLanczosStp;
    liLanczosStp2 = liLanczosStp + X->Bind.Def.Lanczos_max;
  }
#ifdef MPI
  pkomega_bicg_restart(&X->Bind.Check.idim_max, &one, &Nomega, dcSpectrum, dcomega, &liLanczosStp2, &eps_Lanczos, &comm, status,
    &liLanczosStp, &v2[1], &v12[1], &v4[1], &v14[1], alphaCG, betaCG, &z_seed, res_save);
#else
  komega_bicg_restart(&X->Bind.Check.idim_max, &one, &Nomega, dcSpectrum, dcomega, &liLanczosStp2, &eps_Lanczos, status,
    &liLanczosStp, &v2[1], &v12[1], &v4[1], &v14[1], alphaCG, betaCG, &z_seed, res_save);
#endif
  free(alphaCG);
  free(betaCG);
  free(res_save);

  return TRUE;
}/*int ReadTMComponents_BiCG*/


int OutputTMComponents_BiCG(
  struct EDMainCalStruct *X,
  int liLanczosStp
)
{
  char sdt[D_FileNameMax];
  unsigned long int i;
  FILE *fp;
  double complex *alphaCG, *betaCG, *res_save, z_seed;

  printf("debug %d\n", liLanczosStp);
  alphaCG = (double complex*)malloc(liLanczosStp * sizeof(double complex));
  betaCG = (double complex*)malloc(liLanczosStp * sizeof(double complex));
  res_save = (double complex*)malloc(liLanczosStp * sizeof(double complex));
#ifdef MPI
  pkomega_bicg_getcoef(alphaCG, betaCG, &z_seed, res_save);
#else
  komega_bicg_getcoef(alphaCG, betaCG, &z_seed, res_save);
#endif

  sprintf(sdt, cFileNameTridiagonalMatrixComponents, X->Bind.Def.CDataFileHead);
  childfopenMPI(sdt, "w", &fp);
  fprintf(fp, "%d \n", liLanczosStp);
  fprintf(fp, "%.10lf %.10lf\n", creal(z_seed), cimag(z_seed));
  for (i = 0; i < liLanczosStp; i++) {
    fprintf(fp, "%25.16le %25.16le %25.16le %25.16le %25.16le %25.16le\n",
      creal(alphaCG[i]), cimag(alphaCG[i]),
      creal(betaCG[i]), cimag(betaCG[i]),
      creal(res_save[i]), cimag(res_save[i]));
  }
  fclose(fp);
  free(alphaCG);
  free(betaCG);
  free(res_save);

  return TRUE;
}/*int OutputTMComponents_BiCG*/

/** 
 * @brief A main function to calculate spectrum by BiCG method
 * 
 * @param[in,out] X CalcStruct list for getting and pushing calculation information 
 * @retval 0 normally finished
 * @retval -1 error
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 * 
 */
int CalcSpectrumByBiCG(		 
  struct EDMainCalStruct *X,
  double complex *vrhs,
  double complex *v2,
  double complex *v4,
  int Nomega,
  double complex *dcSpectrum,
  double complex *dcomega
)
{
  char sdt[D_FileNameMax];
  unsigned long int i, i_max;
  FILE *fp;
  int iret;
  unsigned long int liLanczosStp_vec = 0;
  double complex *v12, *v14, res_proj;
  int stp, one = 1, status[3];
#ifdef MPI
  MPI_Comm comm = MPI_COMM_WORLD;
#endif

  fprintf(stdoutMPI, "#####  Spectrum calculation with BiCG  #####\n\n");
  
  v12 = (double complex*)malloc((X->Bind.Check.idim_max + 1) * sizeof(double complex));
  v14 = (double complex*)malloc((X->Bind.Check.idim_max + 1) * sizeof(double complex));

  /*
    Read residual vectors
  */
  if (X->Bind.Def.iFlgCalcSpec == RECALC_FROM_TMComponents_VEC ||
      X->Bind.Def.iFlgCalcSpec == RECALC_INOUT_TMComponents_VEC) {
    fprintf(stdoutMPI, "  Start: Input vectors for recalculation.\n");
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_InputSpectrumRecalcvecStart, "a");

    sprintf(sdt, cFileNameOutputRestartVec, X->Bind.Def.CDataFileHead, myrank);
    if (childfopenALL(sdt, "rb", &fp) != 0) {
      exitMPI(-1);
    }
    fread(&liLanczosStp_vec, sizeof(liLanczosStp_vec), 1, fp);
    fread(&i_max, sizeof(long int), 1, fp);
    if (i_max != X->Bind.Check.idim_max) {
      fprintf(stderr, "Error: The size of the input vector is incorrect.\n");
      exitMPI(-1);
    }
    fread(v2, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
    fread(v12, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
    fread(v4, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
    fread(v14, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
    fclose(fp);
    fprintf(stdoutMPI, "  End:   Input vectors for recalculation.\n");
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_InputSpectrumRecalcvecEnd, "a");
  }
  else {
    i_max = X->Bind.Check.idim_max;
    for (i = 1; i <= i_max; i++) {
      v2[i] = vrhs[i];
      v4[i] = vrhs[i];
    }
  }
  /*
    Input alpha, beta, projected residual, or start from scratch
  */
  if (X->Bind.Def.iFlgCalcSpec == RECALC_FROM_TMComponents ||
      X->Bind.Def.iFlgCalcSpec == RECALC_FROM_TMComponents_VEC ||
      X->Bind.Def.iFlgCalcSpec == RECALC_INOUT_TMComponents_VEC)
  {
    iret = ReadTMComponents_BiCG(X,v2,v4,v12,v14,Nomega,dcSpectrum,dcomega);
    if (!iret == TRUE) {
      fprintf(stdoutMPI, "  Error: Fail to read TMcomponents\n");
      return FALSE;
    }
  }
  else {
#ifdef MPI
    pkomega_bicg_init(&i_max, &one, &Nomega, dcSpectrum, dcomega, &X->Bind.Def.Lanczos_max, &eps_Lanczos, &comm);
#else
    komega_bicg_init(&i_max, &one, &Nomega, dcSpectrum, dcomega, &X->Bind.Def.Lanczos_max, &eps_Lanczos);
#endif
  }

  /*
  BiCG loop
  */
  fprintf(stdoutMPI, "    Start: Calculate tridiagonal matrix components.\n");
  TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_GetTridiagonalStart, "a");

  for (stp = 0; stp <= X->Bind.Def.Lanczos_max; stp++) {

    for (i = 1; i <= i_max; i++) v12[i] = 0.0;
    iret = mltply(&X->Bind, v12, v2);
    for (i = 1; i <= i_max; i++) v14[i] = 0.0;
    iret = mltply(&X->Bind, v14, v4);

    res_proj = 0.0;
    for (i = 1; i <= i_max; i++) res_proj += conj(vrhs[i]) * v2[i];
    res_proj = SumMPI_dc(res_proj);

    printf("debug %15.5e %15.5e\n", creal(res_proj), cimag(res_proj));
#ifdef MPI
    pkomega_bicg_update(&v12[1], &v2[1], &v14[1], &v4[1], dcSpectrum, &res_proj, status);
#else
    komega_bicg_update(&v12[1], &v2[1], &v14[1], &v4[1], dcSpectrum, &res_proj, status);
#endif

    fprintf(stdoutMPI, "%d %d %d %25.15e\n", status[0], status[1], status[2], creal(v12[1]));
    if (status[0] < 0) break;
  }/*for (stp = 0; stp <= X->Bind.Def.Lanczos_max; stp++)*/

  fprintf(stdoutMPI, "    End:   Calculate tridiagonal matrix components.\n\n");
  TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_GetTridiagonalEnd, "a");
  /*
   Save alpha, beta, projected residual
  */
  if (X->Bind.Def.iFlgCalcSpec != RECALC_FROM_TMComponents)
    OutputTMComponents_BiCG(X, abs(status[0]));
  /*
    output vectors for recalculation
  */
  if (X->Bind.Def.iFlgCalcSpec == RECALC_OUTPUT_TMComponents_VEC ||
      X->Bind.Def.iFlgCalcSpec == RECALC_INOUT_TMComponents_VEC) {

    fprintf(stdoutMPI, "    Start: Output vectors for recalculation.\n");
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_OutputSpectrumRecalcvecStart, "a");

#ifdef MPI
    pkomega_bicg_getvec(&v12[1], &v14[1]);
#else
    komega_bicg_getvec(&v12[1], &v14[1]);
#endif

    sprintf(sdt, cFileNameOutputRestartVec, X->Bind.Def.CDataFileHead, myrank);
    if (childfopenALL(sdt, "wb", &fp) != 0) {
      exitMPI(-1);
    }
    fwrite(&status[0], sizeof(status[0]), 1, fp);
    fwrite(&X->Bind.Check.idim_max, sizeof(X->Bind.Check.idim_max), 1, fp);
    fwrite(v2, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
    fwrite(v12, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
    fwrite(v4, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
    fwrite(v14, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
    fclose(fp);

    fprintf(stdoutMPI, "    End:   Output vectors for recalculation.\n");
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_OutputSpectrumRecalcvecEnd, "a");
  }

#ifdef MPI
  pkomega_bicg_finalize();
#else
  komega_bicg_finalize();
#endif

  free(v12);
  free(v14);

  return TRUE;
}/*int CalcSpectrumByBiCG*/
