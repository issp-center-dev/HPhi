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
#include "mltplyCommon.h"
#include "./common/setmemory.h"

void debug_print(int num, double complex *var){
  int i;
  for (i=0;i<num;i++){
    printf("debug %d %f %f\n", i, creal(var[i]), cimag(var[i]));
  }
}

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

  iwork = i_1d_allocate(liwork);
  rwork = d_1d_allocate(lrwork);
  mat = cd_1d_allocate(nsub*nsub);
#ifdef FUJITSU
  void *vptr;
  posix_memalign(&vptr, 256, lwork * sizeof(double complex));
  work = (double complex*)vptr;
#else
  work = cd_1d_allocate(lwork);
#endif
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

  free_cd_1d_allocate(mat);
  free_cd_1d_allocate(work);
  free_d_1d_allocate(rwork);
  free_i_1d_allocate(iwork);

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
  double complex **wave//!<[out] [CheckList::idim_max][exct] initial eigenvector
) 
{
  FILE *fp;
  char sdt[D_FileNameMax];
  size_t byte_size = 0;
  double complex *vin;
  int iproc, ie, ierr;
  long int idim, iv, i_max;
  unsigned long int i_max_tmp, sum_i_max;
  int mythread;
  double *dnorm;
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
    vin = cd_1d_allocate(X->Check.idim_max + 1);
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
        byte_size = fread(vin, sizeof(complex double), X->Check.idim_max + 1, fp);
        for (idim = 1; idim <= i_max; idim++) wave[idim][ie] = vin[idim];
        fclose(fp);
      }
    }/*for (ie = 0; ie < X->Def.k_exct; ie++)*/
    free_cd_1d_allocate(vin);

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
      for (idim = 1; idim <= i_max; idim++) wave[idim][ie] = 0.0;
      
      sum_i_max = 0;
      for (iproc = 0; iproc < nproc; iproc++) {

        i_max_tmp = BcastMPI_li(iproc, i_max);
        if (sum_i_max <= iv && iv < sum_i_max + i_max_tmp) {

          if (myrank == iproc) {
            wave[iv - sum_i_max + 1][ie] = 1.0;
            if (X->Def.iInitialVecType == 0) {
              wave[iv - sum_i_max + 1][ie] += 1.0*I;
              wave[iv - sum_i_max + 1][ie] /= sqrt(2.0);
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
            wave[idim][ie] = 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5) + 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5)*I;
        }
        else {
#pragma omp for
          for (idim = 1; idim <= i_max; idim++)
            wave[idim][ie] = 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5);
        }
      }/*for (ie = 0; ie < X->Def.k_exct; ie++)*/

    }/*#pragma omp parallel*/

    dnorm = d_1d_allocate(X->Def.k_exct);
    NormMPI_dv(i_max, X->Def.k_exct, wave, dnorm);
#pragma omp parallel for default(none) shared(i_max,wave,dnorm,X) private(idim,ie)
    for (idim = 1; idim <= i_max; idim++) 
      for (ie = 0; ie < X->Def.k_exct; ie++) wave[idim][ie] /= dnorm[ie];
    free_d_1d_allocate(dnorm);
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
  long unsigned int idim;
  double complex *vout;

  //TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecStart, "a", rand_i, step_i);
  fprintf(stdoutMPI, "%s", cLogOutputVecStart);
  
  vout = cd_1d_allocate(X->Check.idim_max + 1);
  for (ie = 0; ie < X->Def.k_exct; ie++) {
    sprintf(sdt, cFileNameOutputVector, ie, myrank);
    if (childfopenALL(sdt, "wb", &fp) != 0) exitMPI(-1);
    byte_size = fwrite(&X->Large.itr, sizeof(X->Large.itr), 1, fp);
    byte_size = fwrite(&X->Check.idim_max, sizeof(X->Check.idim_max), 1, fp);
    for (idim = 1; idim <= X->Check.idim_max; idim++) vout[idim] = wave[idim][ie];
    byte_size = fwrite(vout, sizeof(complex double), X->Check.idim_max + 1, fp);
    fclose(fp);
  }/*for (ie = 0; ie < X->Def.k_exct; ie++)*/
  free_cd_1d_allocate(vout);

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
  int iconv = -1, i4_max;
  long int idim, i_max;
  int ii, jj, ie, nsub, stp, nsub_cut, nstate;
  double complex ***wxp/*[0] w, [1] x, [2] p of Ref.1*/,
    ***hwxp/*[0] h*w, [1] h*x, [2] h*p of Ref.1*/,
    ****hsub, ****ovlp; /*Subspace Hamiltonian and Overlap*/
  double *eig, *dnorm, eps_LOBPCG, eigabs_max, preshift, precon, dnormmax, *eigsub, eig_pos_shift;
  char tN = 'N', tC = 'C';
  double complex one = 1.0, zero = 0.0;

  nsub = 3 * X->Def.k_exct;
  nstate = X->Def.k_exct;
  eig_pos_shift = LargeValue * X->Def.NsiteMPI;

  eig = d_1d_allocate(X->Def.k_exct);
  dnorm = d_1d_allocate(X->Def.k_exct);
  eigsub = d_1d_allocate(nsub);
  hsub = cd_4d_allocate(3, X->Def.k_exct, 3, X->Def.k_exct);
  ovlp = cd_4d_allocate(3, X->Def.k_exct, 3, X->Def.k_exct);

  i_max = X->Check.idim_max;
  i4_max = (int)i_max;

  free_cd_2d_allocate(v0);
  free_cd_2d_allocate(v1);
  wxp = cd_3d_allocate(3, X->Check.idim_max + 1, X->Def.k_exct);
  hwxp = cd_3d_allocate(3, X->Check.idim_max + 1, X->Def.k_exct);
  /**@brief
  <ul>
  <li>Set initial guess of wavefunction: 
  @f${\bf x}=@f$initial guess</li>
  */
  Initialize_wave(X, wxp[1]);

  TimeKeeper(X, cFileNameTimeKeep, cLanczos_EigenValueStart, "a");

  zclear(i_max*X->Def.k_exct, &hwxp[1][1][0]);
  mltply(X, X->Def.k_exct, hwxp[1], wxp[1]);
  stp = 1;
  TimeKeeperWithStep(X, cFileNameTimeKeep, cLanczos_EigenValueStep, "a", 0);

  zclear(i_max*X->Def.k_exct, &wxp[2][1][0]);
  zclear(i_max*X->Def.k_exct, &hwxp[2][1][0]);
  for (ie = 0; ie < X->Def.k_exct; ie++) eig[ie] = 0.0;
  for (idim = 1; idim <= i_max; idim++) {
    for (ie = 0; ie < X->Def.k_exct; ie++) {
      wxp[2][idim][ie] = 0.0;
      hwxp[2][idim][ie] = 0.0;
      eig[ie] += conj(wxp[1][idim][ie]) * hwxp[1][idim][ie];
    }
  }
  SumMPI_dv(X->Def.k_exct, eig);

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
    /**@brief
     <li>Compute residual vectors: @f${\bf w}={\bf X}-\mu {\bf x}@f$</li>
    */
#pragma omp parallel for default(none) shared(i_max,wxp,hwxp,eig,X) private(idim,ie) 
    for (idim = 1; idim <= i_max; idim++) {
      for (ie = 0; ie < X->Def.k_exct; ie++) {
        wxp[0][idim][ie] = hwxp[1][idim][ie] - eig[ie] * wxp[1][idim][ie];
      }
    }        
    NormMPI_dv(i_max, X->Def.k_exct, wxp[0], dnorm);

    dnormmax = 0.0;
    for (ie = 0; ie < X->Def.k_exct; ie++) 
      if (dnorm[ie] > dnormmax) dnormmax = dnorm[ie];
    /**@brief
    <li>Preconditioning (Point Jacobi): @f${\bf w}={\hat T}^{-1} {\bf w}@f$</li>
    */
    if (stp != 1) {
      if (X->Def.PreCG == 1) {
        for (ie = 0; ie < X->Def.k_exct; ie++) 
          preshift = calc_preshift(eig[ie] + eig_pos_shift, dnorm[ie], eps_LOBPCG) - eig_pos_shift;
#pragma omp parallel for default(none) shared(wxp,list_Diagonal,preshift,i_max,eps_LOBPCG,X) \
private(idim,precon,ie)
        for (idim = 1; idim <= i_max; idim++) {
          for (ie = 0; ie < X->Def.k_exct; ie++){
            precon = list_Diagonal[idim] - preshift;
            if (fabs(precon) > eps_LOBPCG) wxp[0][idim][ie] /= precon;
          }
        }
      }/*if(X->Def.PreCG == 1)*/
      /**@brief
        <li>Normalize residual vector: @f${\bf w}={\bf w}/|w|@f$
      */
      NormMPI_dv(i_max, X->Def.k_exct, wxp[0], dnorm);
#pragma omp parallel for default(none) shared(i_max,wxp,dnorm,X) private(idim,ie)
      for (idim = 1; idim <= i_max; idim++)
        for (ie = 0; ie < X->Def.k_exct; ie++)
          wxp[0][idim][ie] /= dnorm[ie];
    }/*if (stp /= 1)*/
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
    zclear(i_max*X->Def.k_exct, &hwxp[0][1][0]);
    mltply(X, X->Def.k_exct, hwxp[0], wxp[0]);

    TimeKeeperWithStep(X, cFileNameTimeKeep, cLanczos_EigenValueStep, "a", stp);
    /**@brief
    <li>Compute subspace Hamiltonian and overrap matrix:
    @f${\hat H}_{\rm sub}=\{{\bf w},{\bf x},{\bf p}\}^\dagger \{{\bf W},{\bf X},{\bf P}\}@f$, 
    @f${\hat O}=\{{\bf w},{\bf x},{\bf p}\}^\dagger \{{\bf w},{\bf x},{\bf p}\}@f$,
    </li>
    */
    for (ii = 0; ii < 3; ii++) {
      for (jj = 0; jj < 3; jj++) {
        zgemm_(&tN, &tC, &nstate, &nstate, &i4_max, &one,
          &wxp[ii][1][0], &nstate, &wxp[jj][1][0], &nstate, &zero, &ovlp[jj][0][ii][0], &nsub);
        zgemm_(&tN, &tC, &nstate, &nstate, &i4_max, &one,
          &wxp[ii][1][0], &nstate, &hwxp[jj][1][0], &nstate, &zero, &hsub[jj][0][ii][0], &nsub);
      }
    }
    SumMPI_cv(nsub*nsub, &ovlp[0][0][0][0]);
    SumMPI_cv(nsub*nsub, &hsub[0][0][0][0]);

    for (ie = 0; ie < X->Def.k_exct; ie++)
      eig[ie] = creal(hsub[1][ie][1][ie]);
    /**@brief
    <li>Subspace diagonalization with the Lowdin's orthogonalization for
        generalized eigenvalue problem: @f${\hat H}_{\rm sub}{\bf v}={\hat O}\mu_{\rm sub}{\bf v}@f$,
        @f${\bf v}=(\alpha, \beta, \gamma)@f$</li>
    */
    nsub_cut = diag_ovrp(nsub, &hsub[0][0][0][0], &ovlp[0][0][0][0], eigsub);
    /**@brief
    <li>Update @f$\mu=(\mu+\mu_{\rm sub})/2@f$</li>
    */
    for (ie = 0; ie < X->Def.k_exct; ie++)
      eig[ie] = 0.5 * (eig[ie] + eigsub[ie]);
    /**@brief
      <li>@f${\bf x}=\alpha {\bf w}+\beta {\bf x}+\gamma {\bf p}@f$,
      Normalize @f${\bf x}@f$</li>
    */
    zclear(i_max*X->Def.k_exct, &v1buf[1][0]);
    for (ii = 0; ii < 3; ii++) {
      zgemm_(&tC, &tN, &nstate, &i4_max, &nstate, &one,
        &hsub[0][0][ii][0], &nsub, &wxp[ii][1][0], &nstate, &one, &v1buf[1][0], &nstate);
    }
    for (idim = 1; idim <= i_max; idim++) for (ie = 0; ie < X->Def.k_exct; ie++)
      wxp[1][idim][ie] = v1buf[idim][ie];
    /**@brief
    <li>@f${\bf X}=\alpha {\bf W}+\beta {\bf X}+\gamma {\bf P}@f$,
    Normalize @f${\bf X}@f$</li>
    */
    zclear(i_max*X->Def.k_exct, &v1buf[1][0]);
    for (ii = 0; ii < 3; ii++) {
      zgemm_(&tC, &tN, &nstate, &i4_max, &nstate, &one,
        &hsub[0][0][ii][0], &nsub, &hwxp[ii][1][0], &nstate, &one, &v1buf[1][0], &nstate);
    }
    for (idim = 1; idim <= i_max; idim++) for (ie = 0; ie < X->Def.k_exct; ie++)
      hwxp[1][idim][ie] = v1buf[idim][ie];
    /**@brief
    <li>@f${\bf p}=\alpha {\bf w}+\gamma {\bf p}@f$,
    Normalize @f${\bf p}@f$</li>
    */
    zclear(i_max*X->Def.k_exct, &v1buf[1][0]);
    for (ii = 0; ii < 3; ii += 2) {
      zgemm_(&tC, &tN, &nstate, &i4_max, &nstate, &one,
        &hsub[0][0][ii][0], &nsub, &wxp[ii][1][0], &nstate, &one, &v1buf[1][0], &nstate);
    }
    for (idim = 1; idim <= i_max; idim++) for (ie = 0; ie < X->Def.k_exct; ie++)
      wxp[2][idim][ie] = v1buf[idim][ie];
    /**@brief
    <li>@f${\bf P}=\alpha {\bf W}+\gamma {\bf P}@f$,
    Normalize @f${\bf P}@f$</li>
    */
    zclear(i_max*X->Def.k_exct, &v1buf[1][0]);
    for (ii = 0; ii < 3; ii += 2) {
      zgemm_(&tC, &tN, &nstate, &i4_max, &nstate, &one,
        &hsub[0][0][ii][0], &nsub, &hwxp[ii][1][0], &nstate, &one, &v1buf[1][0], &nstate);
    }
    for (idim = 1; idim <= i_max; idim++) for (ie = 0; ie < X->Def.k_exct; ie++)
      hwxp[2][idim][ie] = v1buf[idim][ie];
    /**@brief
    <li>Normalize @f${\bf w}@f$ and @f${\bf W}@f$</li>
    */
    for (ii = 1; ii < 3; ii++) {
      NormMPI_dv(i_max, X->Def.k_exct, wxp[ii], dnorm);
#pragma omp parallel for default(none) shared(i_max,wxp,hwxp,dnorm,ii,X) private(idim,ie)
      for (idim = 1; idim <= i_max; idim++) {
        for (ie = 0; ie < X->Def.k_exct; ie++) {
          wxp[ii][idim][ie] /= dnorm[ie];
          hwxp[ii][idim][ie] /= dnorm[ie];
        }/* for (ie = 0; ie < X->Def.k_exct; ie++)*/
      }
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
  free_d_1d_allocate(dnorm);
  free_d_1d_allocate(eigsub);
  free_cd_4d_allocate(hsub);
  free_cd_4d_allocate(ovlp);
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
  <li>Just Move wxp[1] into ::v1. The latter must be start from 0-index (the same as FullDiag)</li>
  </ul>
  */
  v0 = cd_2d_allocate(X->Check.idim_max + 1, X->Def.k_exct);
#pragma omp parallel for default(none) shared(i_max,wxp,v0,X) private(idim,ie)
  for (idim = 1; idim <= i_max; idim++)
    for (ie = 0; ie < X->Def.k_exct; ie++) 
      v0[idim][ie] = wxp[1][idim][ie];
  free_cd_3d_allocate(wxp);
  v1 = cd_2d_allocate(X->Check.idim_max + 1, X->Def.k_exct);

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
  double complex *vin;

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
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, cReadEigenVecStart, "a");
    vin = cd_1d_allocate(X->Bind.Check.idim_max + 1);
    for (ie = 0; ie < X->Bind.Def.k_exct; ie++) {
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
      byte_size = fread(vin, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
#pragma omp parallel for default(none) shared(v1,vin, i_max, ie), private(idim)
      for (idim = 1; idim <= i_max; idim++) {
        v1[idim][ie] = vin[idim];
      }
      fclose(fp);
    }/*for (ie = 0; ie < X->Def.k_exct; ie++)*/
    free_cd_1d_allocate(vin);
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
    fprintf(fp, "  Energy  %.16lf \n", X->Bind.Phys.energy[ie]);
    fprintf(fp, "  Doublon  %.16lf \n", X->Bind.Phys.doublon[ie]);
    fprintf(fp, "  Sz  %.16lf \n", X->Bind.Phys.Sz[ie]);
    //fprintf(fp, "  S^2  %.16lf \n", X->Bind.Phys.s2[ie]);
    //fprintf(fp, "  N_up  %.16lf \n", X->Bind.Phys.num_up[ie]);
    //fprintf(fp, "  N_down  %.16lf \n", X->Bind.Phys.num_down[ie]);
    fprintf(fp, "\n");
  }
  fclose(fp);
  /*
   Output Eigenvector to a file
  */
  if (X->Bind.Def.iOutputEigenVec == TRUE) {
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, cOutputEigenVecStart, "a");

    vin = cd_1d_allocate(X->Bind.Check.idim_max + 1);
    for (ie = 0; ie < X->Bind.Def.k_exct; ie++) {

#pragma omp parallel for default(none) shared(X,v1,ie,vin) private(idim)
      for (idim = 1; idim <= X->Bind.Check.idim_max; idim++)
        vin[idim] = v1[idim][ie];
      
      sprintf(sdt, cFileNameOutputEigen, X->Bind.Def.CDataFileHead, ie, myrank);
      if (childfopenALL(sdt, "wb", &fp) != 0) exitMPI(-1);
      byte_size = fwrite(&X->Bind.Large.itr, sizeof(X->Bind.Large.itr), 1, fp);
      byte_size = fwrite(&X->Bind.Check.idim_max, sizeof(X->Bind.Check.idim_max), 1, fp);
      byte_size = fwrite(vin, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
      fclose(fp);
    }/*for (ie = 0; ie < X->Bind.Def.k_exct; ie++)*/
    free_cd_1d_allocate(vin);

    TimeKeeper(&(X->Bind), cFileNameTimeKeep, cOutputEigenVecStart, "a");
  }/*if (X->Bind.Def.iOutputEigenVec == TRUE)*/

  return TRUE;

}/*int CalcByLOBPCG*/
