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
/**@file
@brief Functions to perform calculations with the
localy optimal block (preconditioned) conjugate gradient method.
*/
#include "Common.h"
#include "xsetmem.h"
#include "mltply.h"
#include "FileIO.h"
#include "wrapperMPI.h"
#include "expec_cisajs.h"
#include "expec_cisajscktaltdc.h"
#include "expec_totalspin.h"
#include "expec_energy_flct.h"
#include "phys.h"
#include <math.h>
#include "./common/setmemory.h"

void zheevd_(char *jobz, char *uplo, int *n, double complex *a, int *lda, double *w, double complex *work, int *lwork, double *rwork, int * lrwork, int *iwork, int *liwork, int *info);
void zgemm_(char *transa, char *transb, int *m, int *n, int *k, double complex *alpha, double complex *a, int *lda, double complex *b, int *ldb, double complex *beta, double complex *c, int *ldc);
/**
@brief Solve the generalized eigenvalue problem
@f[
{\hat H} |\phi\rangle = \varepsilon {\hat O} |\phi\rangle
@f]
with the Lowdin's orthogonalization
@return the truncated dimension, nsub2
*/
static int diag_ovrp(
  int nsub,//!<[in] Original dimension of subspace
  double complex *hsub,//!<[inout] (nsub*nsub) subspace hamiltonian -> eigenvector
  double complex *ovlp,//!<[inout] (nsub*nsub) Overrap matrix -> @f${\hat O}^{1/2}@f$
  double *eig//!<[out] (nsub) Eigenvalue
)
{
  int *iwork, info, isub, jsub, nsub2;
  char jobz = 'V', uplo = 'U', transa = 'N', transb = 'N';
  double *rwork;
  double complex *work, *mat;
  int liwork, lwork, lrwork;
  double complex one = 1.0, zero = 0.0;

  liwork = 5 * nsub + 3;
  lwork = nsub*nsub + 2 * nsub;
  lrwork = 3 * nsub*nsub + (4 + (int)log2(nsub) + 1) * nsub + 1;

  iwork = (int*)malloc(liwork * sizeof(int));
  rwork = (double*)malloc(lrwork * sizeof(double));
  work = (double complex*)malloc(lwork * sizeof(double complex));
  mat = (double complex*)malloc(nsub*nsub * sizeof(double complex));
  for (isub = 0; isub < nsub*nsub; isub++)mat[isub] = 0.0;
  /**@brief
  (1) Compute @f${\hat O}^{-1/2}@f$ with diagonalizing overrap matrix
  */
  zheevd_(&jobz, &uplo, &nsub, ovlp, &nsub, eig, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
  /**@brief
  @f[
   {\hat O}^{-1/2} = \left(\frac{|O_1\rangle}{\sqrt{o_1}}, \frac{|O_2\rangle}{\sqrt{o_2}},
   ...\right)
   \\
   {\hat O} |O_i\rangle = o_i |O_i\rangle
  @f]
  if @f$o_i@f$ is very small that dimension is ignored. Therefore @f${\hat O}^{-1/2}@f$
  is nsub*nsub2 matrix.
  */
  nsub2 = 0;
  for (isub = 0; isub < nsub; isub++) {
    if (eig[isub] > 1.0e-10) { /*to be changed default 1.0e-14*/
      for (jsub = 0; jsub < nsub; jsub++)
        ovlp[jsub + nsub*nsub2] = ovlp[jsub + nsub*isub] / sqrt(eig[isub]);
      nsub2 += 1;
    }
  }
  for (isub = nsub2; isub < nsub; isub++) 
    for (jsub = 0; jsub < nsub; jsub++)
      ovlp[jsub + nsub*isub] = 0.0;
  /**
  (2) Transform @f${\hat H}'\equiv {\hat O}^{-1/2 \dagger}{\hat H}{\hat O}^{-1/2}@f$.
  @f${\hat H}'@f$ is nsub2*nsub2 matrix.
  */
  transa = 'N';
  zgemm_(&transa, &transb, &nsub, &nsub, &nsub, &one, hsub, &nsub, ovlp, &nsub, &zero, mat, &nsub);
  transa = 'C';
  zgemm_(&transa, &transb, &nsub, &nsub, &nsub, &one, ovlp, &nsub, mat, &nsub, &zero, hsub, &nsub);
  /**
  (3) Diagonalize @f${\hat H}'@f$. It is the standard eigenvalue problem.
  @f[
  {\hat H}' |\phi'_i\rangle = \varepsilon_i |\phi'_i\rangle
  @f]
  */
  zheevd_(&jobz, &uplo, &nsub2, hsub, &nsub, eig, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
  /**
  (4) Transform eigenvector into the original nsub space as
  @f[
  |\phi_i\rangle = {\hat O}^{-1/2} |\phi'_i\rangle
  @f]
  */
  transa = 'N';
  zgemm_(&transa, &transb, &nsub, &nsub, &nsub, &one, ovlp, &nsub, hsub, &nsub, &zero, mat, &nsub);
 // printf("%d %d %15.5f %15.5f %15.5f\n", info, nsub2, eig[0], eig[1], eig[2]);
  for (isub = 0; isub < nsub*nsub; isub++)hsub[isub] = mat[isub];

  free(mat);
  free(work);
  free(rwork);
  free(iwork);

  return(nsub2);
}/*void diag_ovrp*/
/**@brief
Compute adaptively shifted preconditionar written in
S. Yamada, et al., Transactions of JSCES, Paper No. 20060027 (2006).
@return adaptive shift for preconditioning
*/
static double calc_preshift(
  double eig,//!<[in] Eigenvalue in this step
  double res,//!<[in] Residual 2-norm in this step
  double eps_LOBPCG//!<[in] Convergence threshold
)
{
  double k, i;
  double preshift;

  if (res < 1.0) {
    k = trunc(log10(fabs(eig)));
    if (eps_LOBPCG > res) i = ceil(log10(eps_LOBPCG));
    else i = ceil(log10(res));

    preshift = trunc(eig / pow(10.0, k + i))*pow(10.0, k + i);
  }
  else preshift = 0.0;

  return(preshift);
}/*void calc_preshift*/
/*
@brief Compute initial guess for LOBPCG.
If this is resuterting run, read from files.
*/
static void Initialize_wave(
  struct BindStruct *X,//!<[inout]
  double complex **wave//!<[out] [exct][CheckList::idim_max] initial eigenvector
) 
{
  FILE *fp;
  char sdt[D_FileNameMax];
  size_t byte_size = 0;

  int iproc, ie, ierr;
  long int idim, iv, i_max;
  unsigned long int i_max_tmp, sum_i_max;
  int mythread;
  double dnorm;
  /*
  For DSFMT
  */
  long unsigned int u_long_i;
  dsfmt_t dsfmt;
  /**@brief
  (A) For restart: Read saved eigenvector files (as binary files) from each processor
  */
  if (X->Def.iReStart == RESTART_INOUT || X->Def.iReStart == RESTART_IN) {
    //StartTimer(3600);
    //TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecStart, "a", rand_i, step_i);
    fprintf(stdoutMPI, "%s", cLogInputVecStart);

    ierr = 0;
    for (ie = 0; ie < X->Def.k_exct; ie++) {

      sprintf(sdt, cFileNameInputVector, ie, myrank);
      childfopenALL(sdt, "rb", &fp);
      if (fp == NULL) {
        fprintf(stdout, "Restart file is not found.\n");
        fprintf(stdout, "Start from scratch.\n");
        ierr = 1;
        break;
      }
      else {
        byte_size = fread(&iproc, sizeof(int), 1, fp);
        byte_size = fread(&i_max, sizeof(unsigned long int), 1, fp);
        //fprintf(stdoutMPI, "Debug: i_max=%ld, step_i=%d\n", i_max, step_i);
        if (i_max != X->Check.idim_max) {
          fprintf(stderr, "Error: Invalid restart file.\n");
          exitMPI(-1);
        }
        byte_size = fread(wave[ie], sizeof(complex double), X->Check.idim_max + 1, fp);
        fclose(fp);
      }
    }/*for (ie = 0; ie < X->Def.k_exct; ie++)*/

    if (ierr == 0) {
      //TimeKeeperWithRandAndStep(X, cFileNameTPQStep, cOutputVecFinish, "a", rand_i, step_i);
      fprintf(stdoutMPI, "%s", cLogInputVecFinish);
      //StopTimer(3600);
      if (byte_size == 0) printf("byte_size: %d \n", (int)byte_size);
      return;
    }/*if (ierr == 0)*/

  }/*X->Def.iReStart == RESTART_INOUT || X->Def.iReStart == RESTART_IN*/

  /**@brief
  (B) For scratch (including the case that restart files are not found):
  initialize eigenvectors in the same way as TPQ and Lanczos.
  */
  i_max = X->Check.idim_max;

  if (initial_mode == 0) {

    for (ie = 0; ie < X->Def.k_exct; ie++) {

      sum_i_max = SumMPI_li(X->Check.idim_max);
      X->Large.iv = (sum_i_max / 2 + X->Def.initial_iv + ie) % sum_i_max + 1;
      iv = X->Large.iv;
      fprintf(stdoutMPI, "  initial_mode=%d normal: iv = %ld i_max=%ld k_exct =%d\n\n", 
        initial_mode, iv, i_max, X->Def.k_exct);
#pragma omp parallel for default(none) private(idim) shared(wave,i_max,ie)
      for (idim = 1; idim <= i_max; idim++) wave[ie][idim] = 0.0;
      
      sum_i_max = 0;
      for (iproc = 0; iproc < nproc; iproc++) {

        i_max_tmp = BcastMPI_li(iproc, i_max);
        if (sum_i_max <= iv && iv < sum_i_max + i_max_tmp) {

          if (myrank == iproc) {
            wave[ie][iv - sum_i_max + 1] = 1.0;
            if (X->Def.iInitialVecType == 0) {
              wave[ie][iv - sum_i_max + 1] += 1.0*I;
              wave[ie][iv - sum_i_max + 1] /= sqrt(2.0);
            }
          }/*if (myrank == iproc)*/
        }/*if (sum_i_max <= iv && iv < sum_i_max + i_max_tmp)*/

        sum_i_max += i_max_tmp;

      }/*for (iproc = 0; iproc < nproc; iproc++)*/
    }/*for (ie = 0; ie < X->Def.k_exct; ie++)*/
  }/*if(initial_mode == 0)*/
  else if (initial_mode == 1) {
    iv = X->Def.initial_iv;
    fprintf(stdoutMPI, "  initial_mode=%d (random): iv = %ld i_max=%ld k_exct =%d\n\n",
      initial_mode, iv, i_max, X->Def.k_exct);
#pragma omp parallel default(none) private(idim, u_long_i, mythread, dsfmt, ie) \
              shared(wave, iv, X, nthreads, myrank, i_max)
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

      for (ie = 0; ie < X->Def.k_exct; ie++) {
        if (X->Def.iInitialVecType == 0) {
#pragma omp for
          for (idim = 1; idim <= i_max; idim++)
            wave[ie][idim] = 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5) + 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5)*I;
        }
        else {
#pragma omp for
          for (idim = 1; idim <= i_max; idim++)
            wave[ie][idim] = 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5);
        }
      }/*for (ie = 0; ie < X->Def.k_exct; ie++)*/

    }/*#pragma omp parallel*/

    for (ie = 0; ie < X->Def.k_exct; ie++) {
      dnorm = sqrt(creal(VecProdMPI(i_max, wave[ie], wave[ie])));
#pragma omp parallel for default(none) shared(i_max,wave,dnorm,ie) private(idim)
      for (idim = 1; idim <= i_max; idim++) wave[ie][idim] /= dnorm;
    }
  }/*else if(initial_mode==1)*/
}/*static void Initialize_wave*/
/**
@brief Output eigenvectors for restart LOBPCG method
*/
static void Output_restart(
  struct BindStruct *X,//!<[inout]
  double complex **wave//!<[in] [exct][CheckList::idim_max] initial eigenvector
)
{
  FILE *fp;
  size_t byte_size = 0;
  char sdt[D_FileNameMax];
  int ie;

  //TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecStart, "a", rand_i, step_i);
  fprintf(stdoutMPI, "%s", cLogOutputVecStart);
  
  for (ie = 0; ie < X->Def.k_exct; ie++) {
    sprintf(sdt, cFileNameOutputVector, ie, myrank);
    if (childfopenALL(sdt, "wb", &fp) != 0) exitMPI(-1);
    byte_size = fwrite(&X->Large.itr, sizeof(X->Large.itr), 1, fp);
    byte_size = fwrite(&X->Check.idim_max, sizeof(X->Check.idim_max), 1, fp);
    byte_size = fwrite(wave[ie], sizeof(complex double), X->Check.idim_max + 1, fp);
    fclose(fp);
  }/*for (ie = 0; ie < X->Def.k_exct; ie++)*/
  //TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecFinish, "a", rand_i, step_i);
  fprintf(stdoutMPI, "%s", cLogOutputVecFinish);
  if(byte_size == 0) printf("byte_size : %d\n", (int)byte_size);
}/*static void Output_restart*/
/**@brief
Core routine for the LOBPCG method
This method is introduced in 
-# S. Yamada, et al., Transactions of JSCES, Paper No. 20060027 (2006).
   https://www.jstage.jst.go.jp/article/jsces/2006/0/2006_0_20060027/_pdf
-# A. V. Knyazev, SIAM J. Sci.  Compute. 23, 517 (2001).
   http://epubs.siam.org/doi/pdf/10.1137/S1064827500366124
*/
int LOBPCG_Main(
  struct BindStruct *X//!<[inout]
)
{
  char sdt[D_FileNameMax], sdt_2[D_FileNameMax];
  FILE *fp;
  int iconv = -1;
  long int idim, i_max;
  int ii, jj, ie, je, nsub, stp, mythread, nsub_cut;
  double complex ***wxp/*[0] w, [1] x, [2] p of Ref.1*/, 
    ***hwxp/*[0] h*w, [1] h*x, [2] h*p of Ref.1*/,
    *hsub, *ovlp /*Subspace Hamiltonian and Overlap*/,
    **work;
  double *eig, dnorm, eps_LOBPCG, eigabs_max, preshift, precon, dnormmax, *eigsub, eig_pos_shift;

  nsub = 3 * X->Def.k_exct;
  eig_pos_shift = LargeValue * X->Def.NsiteMPI;

  eig = d_1d_allocate(X->Def.k_exct);
  eigsub = d_1d_allocate(nsub);
  hsub = cd_1d_allocate(nsub*nsub);
  ovlp = cd_1d_allocate(nsub*nsub);
  work = cd_2d_allocate(nthreads, nsub);

  i_max = X->Check.idim_max;

  free(v0);
  free(v1);
  free(vg);
  wxp = cd_3d_allocate(3, X->Def.k_exct, X->Check.idim_max + 1);
  hwxp = cd_3d_allocate(3, X->Def.k_exct, X->Check.idim_max + 1);
  /**@brief
  <ul>
  <li>Set initial guess of wavefunction: 
  @f${\bf x}=@f$initial guess</li>
  */
  Initialize_wave(X, wxp[1]);

  TimeKeeper(X, cFileNameTimeKeep, cLanczos_EigenValueStart, "a");

  for (ie = 0; ie < X->Def.k_exct; ie++)
    for (idim = 1; idim <= i_max; idim++) hwxp[1][ie][idim] = 0.0;
  for (ie = 0; ie < X->Def.k_exct; ie++) 
    mltply(X, hwxp[1][ie], wxp[1][ie]);
  stp = 1;
  TimeKeeperWithStep(X, cFileNameTimeKeep, cLanczos_EigenValueStep, "a", 0);

  for (ie = 0; ie < X->Def.k_exct; ie++){
    for (idim = 1; idim <= i_max; idim++) wxp[2][ie][idim] = 0.0;
    for (idim = 1; idim <= i_max; idim++) hwxp[2][ie][idim] = 0.0;

    eig[ie] = creal(VecProdMPI(i_max, wxp[1][ie], hwxp[1][ie]));
  }/*for (ie = 0; ie < X->Def.k_exct; ie++)*/

  sprintf(sdt_2, cFileNameLanczosStep, X->Def.CDataFileHead);
  childfopenMPI(sdt_2, "w", &fp);
  fprintf(stdoutMPI, "    Step   Residual-2-norm     Threshold      Energy\n");
  fprintf(fp, "    Step   Residual-2-norm     Threshold      Energy\n");
  fclose(fp);

  nsub_cut = nsub;
  /**@brief
  <li>@b DO LOBPCG loop</li>
  <ul>
  */
  for (stp = 1; stp <= X->Def.Lanczos_max; stp++) {
    /**@brief
    <li>Scale convergence threshold with the absolute value of eigenvalue
    for numerical stability</li>
    */
    eigabs_max = 0.0;
    for (ie = 0; ie < X->Def.k_exct; ie++)
      if (fabs(eig[ie]) > eigabs_max) eigabs_max = fabs(eig[ie]);
    eps_LOBPCG = pow(10, -0.5 *X->Def.LanczosEps);
    if (eigabs_max > 1.0) eps_LOBPCG *= eigabs_max;
    /**@brief
    <li>@b DO each eigenvector</li>
    <ul>
    */
    dnormmax = 0.0;
    for (ie = 0; ie < X->Def.k_exct; ie++) {
      /**@brief
      <li>Compute residual vectors: @f${\bf w}={\bf X}-\mu {\bf x}@f$</li>
      */
#pragma omp parallel for default(none) shared(i_max,wxp,hwxp,eig,ie) private(idim)
      for (idim = 1; idim <= i_max; idim++)
        wxp[0][ie][idim] = hwxp[1][ie][idim] - eig[ie] * wxp[1][ie][idim];
      dnorm = sqrt(creal(VecProdMPI(i_max, wxp[0][ie], wxp[0][ie])));
      if (dnorm > dnormmax) dnormmax = dnorm;

      if (stp != 1) {
        /**@brief
        <li>Preconditioning (Point Jacobi): @f${\bf w}={\hat T}^{-1} {\bf w}@f$</li>
        */
        if (X->Def.PreCG == 1) {
          preshift = calc_preshift(eig[ie]+ eig_pos_shift, dnorm, eps_LOBPCG) - eig_pos_shift;
#pragma omp parallel for default(none) shared(wxp,ie,list_Diagonal,preshift,i_max,eps_LOBPCG) private(idim,precon)
          for (idim = 1; idim <= i_max; idim++) {
            precon = list_Diagonal[idim] - preshift;
            if(fabs(precon) > eps_LOBPCG) wxp[0][ie][idim] /= precon;
          }
        }/*if(X->Def.PreCG == 1)*/
        /**@brief
        <li>Normalize residual vector: @f${\bf w}={\bf w}/|w|@f$
        */
        dnorm = sqrt(creal(VecProdMPI(i_max, wxp[0][ie], wxp[0][ie])));
#pragma omp parallel for default(none) shared(i_max,wxp,dnorm,ie) private(idim)
        for (idim = 1; idim <= i_max; idim++) wxp[0][ie][idim] /= dnorm;
      }/*if (stp /= 1)*/

    }/*for (ie = 0; ie < X->Def.k_exct; ie++)*/
    /**@brief
    </ul>
    <li>@b END @b DO each eigenvector</li>
    <li>Convergence check</li>
    */
    childfopenMPI(sdt_2, "a", &fp);
    fprintf(stdoutMPI, "%9d %15.5e %15.5e      ", stp, dnormmax, eps_LOBPCG);
    fprintf(fp, "%9d %15.5e %15.5e      ", stp, dnormmax, eps_LOBPCG);
    for (ie = 0; ie < X->Def.k_exct; ie++) {
      fprintf(stdoutMPI, " %15.5e", eig[ie]);
      fprintf(fp, " %15.5e", eig[ie]);
    }
    if(nsub_cut == 0) printf("nsub_cut : %d", nsub_cut);
    fprintf(stdoutMPI, "\n");
    fprintf(fp, "\n");
    fclose(fp);

    if (dnormmax < eps_LOBPCG) {
      iconv = 0;
      break;
    }
    /**@brief
    <li>@f${\bf W}={\hat H}{\bf w}@f$</li>
    */
#pragma omp parallel default(none) shared(hwxp,i_max,X) private(idim,ie)
    {
#pragma omp for nowait
      for (ie = 0; ie < X->Def.k_exct; ie++)
        for (idim = 1; idim <= i_max; idim++)  hwxp[0][ie][idim] = 0.0;
    }
    for (ie = 0; ie < X->Def.k_exct; ie++)
      mltply(X, hwxp[0][ie], wxp[0][ie]);

    TimeKeeperWithStep(X, cFileNameTimeKeep, cLanczos_EigenValueStep, "a", stp);
    /**@brief
    <li>Compute subspace Hamiltonian and overrap matrix:
    @f${\hat H}_{\rm sub}=\{{\bf w},{\bf x},{\bf p}\}^\dagger \{{\bf W},{\bf X},{\bf P}\}@f$, 
    @f${\hat O}=\{{\bf w},{\bf x},{\bf p}\}^\dagger \{{\bf w},{\bf x},{\bf p}\}@f$,
    </li>
    */
    for (ii = 0; ii < 3; ii++) {
      for (ie = 0; ie < X->Def.k_exct; ie++){
        for (jj = 0; jj < 3; jj++) {
          for (je = 0; je < X->Def.k_exct; je++){
            hsub[je + jj*X->Def.k_exct + ie * nsub + ii * nsub*X->Def.k_exct]
              = VecProdMPI(i_max, wxp[jj][je], hwxp[ii][ie]);
            ovlp[je + jj*X->Def.k_exct + ie * nsub + ii * nsub*X->Def.k_exct]
              = VecProdMPI(i_max, wxp[jj][je], wxp[ii][ie]);
          }/*for (je = 0; je < X->Def.k_exct; je++)*/
        }/*for (jj = 0; jj < 3; jj++)*/
      }/*for (ie = 0; ie < X->Def.k_exct; ie++)*/
    }/*for (ii = 0; ii < 3; ii++)*/
    for (ie = 0; ie < X->Def.k_exct; ie++)
      eig[ie] =
      creal(hsub[ie + 1 * X->Def.k_exct + ie * nsub + 1 * nsub*X->Def.k_exct]);
    /**@brief
    <li>Subspace diagonalization with the Lowdin's orthogonalization for
        generalized eigenvalue problem: @f${\hat H}_{\rm sub}{\bf v}={\hat O}\mu_{\rm sub}{\bf v}@f$,
        @f${\bf v}=(\alpha, \beta, \gamma)@f$</li>
    */
    nsub_cut = diag_ovrp(nsub, hsub, ovlp, eigsub);
    /**@brief
    <li>Update @f$\mu=(\mu+\mu_{\rm sub})/2@f$</li>
    */
    for (ie = 0; ie < X->Def.k_exct; ie++)
      eig[ie] = 0.5 * (eig[ie] + eigsub[ie]);

#pragma omp parallel default(none) shared(i_max,X,wxp,hwxp,hsub,nsub,work) private(idim,ie,je,jj,mythread)
    {
#if defined(_OPENMP)
      mythread = omp_get_thread_num();
#else
      mythread = 0;
#endif

#pragma omp for
      for (idim = 1; idim <= i_max; idim++) {
        /**@brief
        <li>@f${\bf x}=\alpha {\bf w}+\beta {\bf x}+\gamma {\bf p}@f$,
        Normalize @f${\bf x}@f$</li>
        */
        for (ie = 0; ie < X->Def.k_exct; ie++) {
          work[mythread][ie] = 0.0;
          for (jj = 0; jj < 3; jj++)
            for (je = 0; je < X->Def.k_exct; je++)
              work[mythread][ie] += wxp[jj][je][idim] * hsub[je + jj *X->Def.k_exct + nsub*ie];
        }
        for (ie = 0; ie < X->Def.k_exct; ie++) wxp[1][ie][idim] = work[mythread][ie];
        /**@brief
        <li>@f${\bf X}=\alpha {\bf W}+\beta {\bf X}+\gamma {\bf P}@f$,
        Normalize @f${\bf X}@f$</li>
        */
        for (ie = 0; ie < X->Def.k_exct; ie++) {
          work[mythread][ie] = 0.0;
          for (jj = 0; jj < 3; jj++)
            for (je = 0; je < X->Def.k_exct; je++)
              work[mythread][ie] += hwxp[jj][je][idim] * hsub[je + jj *X->Def.k_exct + nsub*ie];
        }
        for (ie = 0; ie < X->Def.k_exct; ie++) hwxp[1][ie][idim] = work[mythread][ie];
        /**@brief
        <li>@f${\bf p}=\alpha {\bf w}+\gamma {\bf p}@f$,
        Normalize @f${\bf p}@f$</li>
        */
        for (ie = 0; ie < X->Def.k_exct; ie++) {
          work[mythread][ie] = 0.0;
          for (jj = 0; jj < 3; jj += 2) {
            for (je = 0; je < X->Def.k_exct; je++)
              work[mythread][ie] += wxp[jj][je][idim] * hsub[je + jj *X->Def.k_exct + nsub*ie];
          }
        }
        for (ie = 0; ie < X->Def.k_exct; ie++) wxp[2][ie][idim] = work[mythread][ie];
        /**@brief
        <li>@f${\bf P}=\alpha {\bf W}+\gamma {\bf P}@f$,
        Normalize @f${\bf P}@f$</li>
        */
        for (ie = 0; ie < X->Def.k_exct; ie++) {
          work[mythread][ie] = 0.0;
          for (jj = 0; jj < 3; jj += 2)
            for (je = 0; je < X->Def.k_exct; je++)
              work[mythread][ie] += hwxp[jj][je][idim] * hsub[je + jj *X->Def.k_exct + nsub*ie];
        }
        for (ie = 0; ie < X->Def.k_exct; ie++) hwxp[2][ie][idim] = work[mythread][ie];

      }/*for (idim = 1; idim <= i_max; idim++)*/
    }/*pragma omp parallel*/
    /**@brief
    <li>Normalize @f${\bf w}@f$ and @f${\bf W}@f$</li>
    */
    for (ii = 1; ii < 3; ii++) {
      for (ie = 0; ie < X->Def.k_exct; ie++) {
        dnorm = sqrt(creal(VecProdMPI(i_max, wxp[ii][ie], wxp[ii][ie])));
#pragma omp parallel for default(none) shared(i_max,wxp,hwxp,dnorm,ie,ii) private(idim)
        for (idim = 1; idim <= i_max; idim++) {
          wxp[ii][ie][idim] /= dnorm;
          hwxp[ii][ie][idim] /= dnorm;
        }
      }/* for (ie = 0; ie < X->Def.k_exct; ie++)*/
    }/*for (ii = 1; ii < 3; ii++)*/

  }/*for (stp = 1; stp <= X->Def.Lanczos_max; stp++)*/
  /**@brief
  </ul>
  <li>@b END @b DO LOBPCG iteration
  */
  //fclose(fp);

  X->Large.itr = stp;
  sprintf(sdt, cFileNameTimeKeep, X->Def.CDataFileHead);

  TimeKeeper(X, cFileNameTimeKeep, cLanczos_EigenValueFinish, "a");
  fprintf(stdoutMPI, "%s", cLogLanczos_EigenValueEnd);

  free_d_1d_allocate(eig);
  free_d_1d_allocate(eigsub);
  free_cd_1d_allocate(hsub);
  free_cd_1d_allocate(ovlp);
  free_cd_2d_allocate(work);
  free_cd_3d_allocate(hwxp);
  /**@brief
  <li>Output resulting vectors for restart</li>
  */
  if (X->Def.iReStart == RESTART_OUT || X->Def.iReStart == RESTART_INOUT){
      Output_restart(X, wxp[1]);
      if(iconv != 0) {
          sprintf(sdt, "%s", cLogLanczos_EigenValueNotConverged);
          return 1;
      }
  }
  /**@brief
  <li>Just Move wxp[1] into ::L_vec. The latter must be start from 0-index (the same as FullDiag)</li>
  </ul>
  */
  L_vec = cd_2d_allocate(X->Def.k_exct, X->Check.idim_max + 1);
#pragma omp parallel default(none) shared(i_max,wxp,L_vec,X) private(idim,ie)
  {
    for (ie = 0; ie < X->Def.k_exct; ie++) {
#pragma omp for nowait
      for (idim = 0; idim < i_max; idim++)
        L_vec[ie][idim] = wxp[1][ie][idim + 1];
    }/*for (ie = 0; ie < X->Def.k_exct; ie++)*/
  }/*X->Def.k_exct, X->Check.idim_max + 1);*/
  free_cd_3d_allocate(wxp);

  v0 = cd_1d_allocate(X->Check.idim_max + 1);
  v1 = cd_1d_allocate(X->Check.idim_max + 1);
  vg = cd_1d_allocate(X->Check.idim_max + 1);
  if (iconv != 0) {
    sprintf(sdt, "%s", cLogLanczos_EigenValueNotConverged);
    return -1;
  }
  else {
    return 0;
  }
}/*int LOBPCG_Main*/
/**
@brief Driver routine for LOB(P)CG method.
*/
int CalcByLOBPCG(
  struct EDMainCalStruct *X//![inout]
)
{
  char sdt[D_FileNameMax];
  size_t byte_size = 0;
  long int i_max = 0, ie, idim;
  FILE *fp;

  fprintf(stdoutMPI, "######  Eigenvalue with LOBPCG  #######\n\n");

  if (X->Bind.Def.iInputEigenVec == FALSE) {

    // this part will be modified
    switch (X->Bind.Def.iCalcModel) {
    case HubbardGC:
    case SpinGC:
    case KondoGC:
    case SpinlessFermionGC:
      initial_mode = 1; // 1 -> random initial vector
      break;
    case Hubbard:
    case Kondo:
    case Spin:
    case SpinlessFermion:

      if (X->Bind.Def.iFlgGeneralSpin == TRUE) {
        initial_mode = 1;
      }
      else {
        if (X->Bind.Def.initial_iv>0) {
          initial_mode = 0; // 0 -> only v[iv] = 1
        }
        else {
          initial_mode = 1; // 1 -> random initial vector
        }
      }
      break;
    default:
      //fclose(fp);
      exitMPI(-1);
    }

    int iret = LOBPCG_Main(&(X->Bind));
    if (iret != 0) {
      if(iret ==1) return (TRUE);
      else{
          fprintf(stdoutMPI, "  LOBPCG is not converged in this process.\n");
          return(FALSE);
      }
    }
  }/*if (X->Bind.Def.iInputEigenVec == FALSE)*/
  else {// X->Bind.Def.iInputEigenVec=true :input v1:
    /**@brief
    If this run is for spectrum calculation, eigenvectors are not computed
    and read from files.
    */
    fprintf(stdoutMPI, "An Eigenvector is inputted.\n");
    L_vec = cd_2d_allocate(X->Bind.Def.k_exct, X->Bind.Check.idim_max + 1);
    for (ie = 0; ie < X->Bind.Def.k_exct; ie++) {
      TimeKeeper(&(X->Bind), cFileNameTimeKeep, cReadEigenVecStart, "a");
      sprintf(sdt, cFileNameInputEigen, X->Bind.Def.CDataFileHead, ie, myrank);
      childfopenALL(sdt, "rb", &fp);
      if (fp == NULL) {
        fprintf(stderr, "Error: Inputvector file is not found.\n");
        exitMPI(-1);
      }
      byte_size = fread(&step_i, sizeof(int), 1, fp);
      byte_size = fread(&i_max, sizeof(long int), 1, fp);
      if (i_max != X->Bind.Check.idim_max) {
        fprintf(stderr, "Error: Invalid Inputvector file.\n");
        exitMPI(-1);
      }
      byte_size = fread(v1, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
#pragma omp parallel for default(none) shared(L_vec, v1) firstprivate(i_max, ie), private(idim)
      for (idim = 0; idim < i_max; idim++) {
        L_vec[ie][idim] = v1[idim + 1];
      }
    }/*for (ie = 0; ie < X->Def.k_exct; ie++)*/
    fclose(fp);
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, cReadEigenVecFinish, "a");

    if(byte_size == 0) printf("byte_size : %d\n", (int)byte_size);
  }/*X->Bind.Def.iInputEigenVec == TRUE*/

  fprintf(stdoutMPI, "%s", cLogLanczos_EigenVecEnd);
  /**@brief
    Compute & Output physical variables to a file
    the same function as FullDiag [phys()] is used.
  */
  phys(&(X->Bind), X->Bind.Def.k_exct);

  X->Bind.Def.St=1;
  if (X->Bind.Def.St == 0) {
    sprintf(sdt, cFileNameEnergy_Lanczos, X->Bind.Def.CDataFileHead);
  }
  else if (X->Bind.Def.St == 1) {
    sprintf(sdt, cFileNameEnergy_CG, X->Bind.Def.CDataFileHead);
  }

  if (childfopenMPI(sdt, "w", &fp) != 0) {
    exitMPI(-1);
  }
  for (ie = 0; ie < X->Bind.Def.k_exct; ie++) {
    //phys(&(X->Bind), ie);
    fprintf(fp, "State %ld\n", ie);
    fprintf(fp, "  Energy  %.16lf \n", X->Bind.Phys.all_energy[ie]);
    fprintf(fp, "  Doublon  %.16lf \n", X->Bind.Phys.all_doublon[ie]);
    fprintf(fp, "  Sz  %.16lf \n", X->Bind.Phys.all_sz[ie]);
    //fprintf(fp, "  S^2  %.16lf \n", X->Bind.Phys.all_s2[ie]);
    //fprintf(fp, "  N_up  %.16lf \n", X->Bind.Phys.all_num_up[ie]);
    //fprintf(fp, "  N_down  %.16lf \n", X->Bind.Phys.all_num_down[ie]);
    fprintf(fp, "\n");
  }
  fclose(fp);
  /*
   Output Eigenvector to a file
  */
  if (X->Bind.Def.iOutputEigenVec == TRUE) {
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, cOutputEigenVecStart, "a");

    for (ie = 0; ie < X->Bind.Def.k_exct; ie++) {

#pragma omp parallel for default(none) shared(X,v1,L_vec,ie) private(idim)
      for (idim = 0; idim < X->Bind.Check.idim_max; idim++)
        v1[idim + 1] = L_vec[ie][idim];
      
      sprintf(sdt, cFileNameOutputEigen, X->Bind.Def.CDataFileHead, ie, myrank);
      if (childfopenALL(sdt, "wb", &fp) != 0) exitMPI(-1);
      byte_size = fwrite(&X->Bind.Large.itr, sizeof(X->Bind.Large.itr), 1, fp);
      byte_size = fwrite(&X->Bind.Check.idim_max, sizeof(X->Bind.Check.idim_max), 1, fp);
      byte_size = fwrite(v1, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
      fclose(fp);
    }/*for (ie = 0; ie < X->Bind.Def.k_exct; ie++)*/

    TimeKeeper(&(X->Bind), cFileNameTimeKeep, cOutputEigenVecStart, "a");
  }/*if (X->Bind.Def.iOutputEigenVec == TRUE)*/

  return TRUE;

}/*int CalcByLOBPCG*/
