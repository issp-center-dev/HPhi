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

#include "Common.h"
#include "mfmemory.h"
#include "xsetmem.h"
#include "mltply.h"
#include "FileIO.h"
#include "wrapperMPI.h"
#include "expec_cisajs.h"
#include "expec_cisajscktaltdc.h"
#include "expec_totalspin.h"
#include "expec_energy.h"
#include <math.h>

void zhegv_(int *itype, char *jobz, char *uplo, int *n, double complex *a, int *lda, double complex *b, int *ldb, double *w, double complex *work, int *lwork, double *rwork, int *info);
void zheevd_(char *jobz, char *uplo, int *n, double complex *a, int *lda, double *w, double complex *work, int *lwork, double *rwork, int * lrwork, int *iwork, int *liwork, int *info);
void zgemm_(char *transa, char *transb, int *m, int *n, int *k, double complex *alpha, double complex *a, int *lda, double complex *b, int *ldb, double complex *beta, double complex *c, int *ldc);

/*
Conjgate vector product
*/
double complex vec_prod(
  long unsigned int ndim,
  double complex *v1,
  double complex *v2)
{
  long unsigned int idim;
  double complex prod;

  prod = 0.0;
#pragma omp parallel for default(none) shared(v1,v2,ndim) private(idim) reduction(+: prod)
  for (idim = 1; idim <= ndim; idim++) prod += conj(v1[idim]) * v2[idim];
  prod = SumMPI_dc(prod);

  return(prod);
}/*double complex vec_prod*/

void diag_ovrp(
  int nsub,
  double complex *hsub,
  double complex *ovrp,
  double *eig
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

  zheevd_(&jobz, &uplo, &nsub, ovrp, &nsub, eig, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
  /*
   Compute O^{-1/2}
  */
  nsub2 = 0;
  for (isub = 0; isub < nsub; isub++) {
    if (eig[isub] > sqrt(eps_Lanczos)) {
      for (jsub = 0; jsub < nsub; jsub++)
        ovrp[jsub + nsub*nsub2] = ovrp[jsub + nsub*isub] / sqrt(eig[isub]);
      nsub2 += 1;
    }
  }
  for (isub = nsub2; isub < nsub; isub++) 
    for (jsub = 0; jsub < nsub; jsub++)
      ovrp[jsub + nsub*isub] = 0.0;
  /*
  Compute orthogonalized hsub
  */
  transa = 'N';
  zgemm_(&transa, &transb, &nsub, &nsub, &nsub, &one, hsub, &nsub, ovrp, &nsub, &zero, mat, &nsub);
  transa = 'C';
  zgemm_(&transa, &transb, &nsub, &nsub, &nsub, &one, ovrp, &nsub, mat, &nsub, &zero, hsub, &nsub);
  /*
   Diagonalize orthogonalized one
  */
  zheevd_(&jobz, &uplo, &nsub2, hsub, &nsub, eig, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
  transa = 'N';
  zgemm_(&transa, &transb, &nsub, &nsub, &nsub, &one, ovrp, &nsub, hsub, &nsub, &zero, mat, &nsub);
 // printf("%d %d %15.5f %15.5f %15.5f\n", info, nsub2, eig[0], eig[1], eig[2]);
  for (isub = 0; isub < nsub*nsub; isub++)hsub[isub] = mat[isub];

  free(mat);
  free(work);
  free(rwork);
  free(iwork);
}/*void diag_ovrp*/

double calc_preshift(
  double eig,
  double res,
  double eps_LOBPCG
)
{
  double k, i;
  double preshift;

  if (fabs(eig) > 10.0) k = trunc(log10(fabs(eig)));
  else k = 1.0;

  if (res < 1.0) {
    if (eps_LOBPCG > res) i = ceil(log10(eps_LOBPCG));
    else i = ceil(log10(res));

    preshift = trunc(eig / pow(10.0, k + i))*pow(10.0, k + i);
  }
  else preshift = 0.0;

  return(preshift);

}/*void calc_preshift*/

/*
Core routine for the LOBPCG method
*/
int LOBPCG_Main(
  struct BindStruct *X
)
{
  FILE *fp;
  char sdt[D_FileNameMax], sdt_2[D_FileNameMax];
  int stp, iproc;
  long int idim, iv, i_max;
  unsigned long int i_max_tmp, sum_i_max;
  int iconv = -1;
  int mythread;

  int ii, jj, ie, je, nsub;
  double complex ***wxp/*[0] w, [1] x, [2] p */, 
    ***hwxp/*[0] h*w, [1] h*x, [2] h*p */,
    *hsub, *ovrp;
  double *eig, dnorm, eps_LOBPCG, eigabs_max, preshift, precon, dnormmax;
  /*
  Variables for zhegv
  */
  int lwork, info, lda, ldb;
  int itype = 1, nsub2;
  double *rwork, *eigsub;
  double complex *work;
  char jobz = 'V', uplo = 'U';
  /*
  For DSFMT
  */
  long unsigned int u_long_i;
  dsfmt_t dsfmt;

  sprintf(sdt_2, cFileNameLanczosStep, X->Def.CDataFileHead);
  nsub = 3 * X->Def.k_exct;
  lda = nsub;
  ldb = nsub;
  lwork = 2 * nsub - 1;
  c_malloc1(work, lwork);
  d_malloc1(rwork, 3 * nsub - 2);
  d_malloc1(eigsub, nsub);

  d_malloc1(eig, X->Def.k_exct);
  c_malloc1(hsub, nsub*nsub);
  c_malloc1(ovrp, nsub*nsub);

  i_max = X->Check.idim_max;

  free(v0);
  free(v1);
  free(vg);
  c_malloc3(wxp, 3, X->Def.k_exct, X->Check.idim_max + 1);
  c_malloc3(hwxp, 3, X->Def.k_exct, X->Check.idim_max + 1);
  /*
   Set initial guess of wavefunction
  */
  if (initial_mode == 0) {

    for (ie = 0; ie < X->Def.k_exct; ie++) {

      sum_i_max = SumMPI_li(X->Check.idim_max);
      X->Large.iv = (sum_i_max / 2 + X->Def.initial_iv + ie) % sum_i_max + 1;
      iv = X->Large.iv;
      fprintf(stdoutMPI, "  initial_mode=%d normal: iv = %ld i_max=%ld k_exct =%d \n\n", initial_mode, iv, i_max, X->Def.k_exct);
#pragma omp parallel for default(none) private(idim) shared(wxp,i_max,ie)
      for (idim = 1; idim <= i_max; idim++) {
        wxp[1][ie][idim] = 0.0;
      }

      sum_i_max = 0;
      for (iproc = 0; iproc < nproc; iproc++) {

        i_max_tmp = BcastMPI_li(iproc, i_max);
        if (sum_i_max <= iv && iv < sum_i_max + i_max_tmp) {

          if (myrank == iproc) {
            wxp[1][ie][iv - sum_i_max + 1] = 1.0;
            if (X->Def.iInitialVecType == 0) {
              wxp[1][ie][iv - sum_i_max + 1] += 1.0*I;
              wxp[1][ie][iv - sum_i_max + 1] /= sqrt(2.0);
            }
          }/*if (myrank == iproc)*/
        }/*if (sum_i_max <= iv && iv < sum_i_max + i_max_tmp)*/

        sum_i_max += i_max_tmp;

      }/*for (iproc = 0; iproc < nproc; iproc++)*/
    }/*for (ie = 0; ie < X->Def.k_exct; ie++)*/
  }/*if(initial_mode == 0)*/
  else if (initial_mode == 1) {
    iv = X->Def.initial_iv;
    fprintf(stdoutMPI, "  initial_mode=%d (random): iv = %ld i_max=%ld k_exct =%d \n\n", initial_mode, iv, i_max, X->Def.k_exct);
#pragma omp parallel default(none) private(idim, u_long_i, mythread, dsfmt, ie) \
            shared(wxp, iv, X, nthreads, myrank, i_max)
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
            wxp[1][ie][idim] = 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5) + 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5)*I;
        }
        else {
#pragma omp for
          for (idim = 1; idim <= i_max; idim++)
            wxp[1][ie][idim] = 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5);
        }
      }/*for (ie = 0; ie < X->Def.k_exct; ie++)*/

    }/*#pragma omp parallel*/

    for (ie = 0; ie < X->Def.k_exct; ie++) {
      dnorm = sqrt(creal(vec_prod(i_max, wxp[1][ie], wxp[1][ie])));
#pragma omp parallel for default(none) private(idim) shared(wxp,i_max, dnorm, ie)
      for (idim = 1; idim <= i_max; idim++)  wxp[1][ie][idim] = wxp[1][ie][idim] / dnorm;
    }
  }/*else if(initial_mode==1)*/

   //Eigenvalues by Lanczos method
  TimeKeeper(X, cFileNameTimeKeep, cLanczos_EigenValueStart, "a");

  for (ie = 0; ie < X->Def.k_exct; ie++)
    for (idim = 1; idim <= i_max; idim++)  hwxp[1][ie][idim] = 0.0;

  for (ie = 0; ie < X->Def.k_exct; ie++) 
    mltply(X, hwxp[1][ie], wxp[1][ie]);
  stp = 1;
  TimeKeeperWithStep(X, cFileNameTimeKeep, cLanczos_EigenValueStep, "a", stp);

  for (ie = 0; ie < X->Def.k_exct; ie++){
    for (idim = 1; idim <= i_max; idim++) wxp[2][ie][idim] = 0.0;
    for (idim = 1; idim <= i_max; idim++) hwxp[2][ie][idim] = 0.0;

    eig[ie] = creal(vec_prod(i_max, wxp[1][ie], hwxp[1][ie]));
  }/*for (ie = 0; ie < X->Def.k_exct; ie++)*/

  fprintf(stdoutMPI, "    Step   Residual-2-norm     Threshold   zhegv-info      Energy\n");

  info = 0;
  for (stp = 1; stp <= X->Def.Lanczos_max; stp++) {
    /*
     Convergence threshold
    */
    eigabs_max = 0.0;
    for (ie = 0; ie < X->Def.k_exct; ie++)
      if (fabs(eig[ie]) > eigabs_max) eigabs_max = fabs(eig[ie]);
    eps_LOBPCG = pow(10, -0.5 *X->Def.LanczosEps);
    if (eigabs_max > 1.0) eps_LOBPCG *= eigabs_max;
    /*
     Compute residual vectors
    */
    dnormmax = 0.0;
    for (ie = 0; ie < X->Def.k_exct; ie++) {
      for (idim = 1; idim <= i_max; idim++)
        wxp[0][ie][idim] = hwxp[1][ie][idim] - eig[ie] * wxp[1][ie][idim];
      dnorm = sqrt(creal(vec_prod(i_max, wxp[0][ie], wxp[0][ie])));
      if (dnorm > dnormmax) dnormmax = dnorm;

      if (stp /= 1) {
        /*
         Preconditioning (Point Jacobi)
        */
        preshift = calc_preshift(eig[ie], dnorm, eps_LOBPCG);
        for (idim = 1; idim <= i_max; idim++) {
          precon = list_Diagonal[idim] - preshift;
          if (fabs(precon) > eps_LOBPCG) wxp[0][ie][idim] /= precon;
        }
        /*
         Normalize residual vector
        */
        dnorm = sqrt(creal(vec_prod(i_max, wxp[0][ie], wxp[0][ie])));
        for (idim = 1; idim <= i_max; idim++) wxp[0][ie][idim] = wxp[0][ie][idim] / dnorm;
      }/*if (stp /= 1)*/

    }/*for (ie = 0; ie < X->Def.k_exct; ie++)*/
    /*
    Convergence check
    */

    fprintf(stdoutMPI, "%9d %15.5e %15.5e      %d   ",
      stp, dnormmax, eps_LOBPCG, info);
    for (ie = 0; ie < X->Def.k_exct; ie++)
      fprintf(stdoutMPI, " %15.5e", eig[ie]);
    fprintf(stdoutMPI, "\n");

    if (dnormmax < eps_LOBPCG) {
      iconv = 0;
      break;
    }

    for (ie = 0; ie < X->Def.k_exct; ie++)
      for (idim = 1; idim <= i_max; idim++)  hwxp[0][ie][idim] = 0.0;
    for (ie = 0; ie < X->Def.k_exct; ie++)
      mltply(X, hwxp[0][ie], wxp[0][ie]);

    for (ii = 0; ii < 3; ii++) {
      for (ie = 0; ie < X->Def.k_exct; ie++){
        for (jj = 0; jj < 3; jj++) {
          for (je = 0; je < X->Def.k_exct; je++){
            hsub[je + jj*X->Def.k_exct + ie * nsub + ii * nsub*X->Def.k_exct]
              = vec_prod(i_max, wxp[jj][je], hwxp[ii][ie]);
            ovrp[je + jj*X->Def.k_exct + ie * nsub + ii * nsub*X->Def.k_exct]
              = vec_prod(i_max, wxp[jj][je], wxp[ii][ie]);
          }
        }
      }
    }
    for (ie = 0; ie < X->Def.k_exct; ie++)
      eig[ie] =
      creal(hsub[ie + 1 * X->Def.k_exct + ie * nsub + 1 * nsub*X->Def.k_exct]);

    /*
    if (stp == 1) {
      nsub2 = 2 * X->Def.k_exct;
      zhegv_(&itype, &jobz, &uplo, &nsub2, hsub, &lda, ovrp, &ldb, eigsub, work, &lwork, rwork, &info);
    }
    else {
      nsub2 = nsub;
      zhegv_(&itype, &jobz, &uplo, &nsub2, hsub, &lda, ovrp, &ldb, eigsub, work, &lwork, rwork, &info);
    }
    */
    
    diag_ovrp(nsub, hsub, ovrp, eigsub);
   
    for (ie = 0; ie < X->Def.k_exct; ie++)
      eig[ie] = 0.5 * (eig[ie] + eigsub[ie]);

    for (idim = 1; idim <= i_max; idim++) {

      for (ie = 0; ie < X->Def.k_exct; ie++) {
        work[ie] = 0.0;
        for (jj = 0; jj < 3; jj++)
          for (je = 0; je < X->Def.k_exct; je++)
            work[ie] += wxp[jj][je][idim] * hsub[je + jj *X->Def.k_exct + nsub*ie];
      }
      for (ie = 0; ie < X->Def.k_exct; ie++) wxp[1][ie][idim] = work[ie];

      for (ie = 0; ie < X->Def.k_exct; ie++) {
        work[ie] = 0.0;
        for (jj = 0; jj < 3; jj++)
          for (je = 0; je < X->Def.k_exct; je++)
            work[ie] += hwxp[jj][je][idim] * hsub[je + jj *X->Def.k_exct + nsub*ie];
      }
      for (ie = 0; ie < X->Def.k_exct; ie++) hwxp[1][ie][idim] = work[ie];

      for (ie = 0; ie < X->Def.k_exct; ie++) {
        work[ie] = 0.0;
        for (jj = 0; jj < 3; jj += 2) {
          for (je = 0; je < X->Def.k_exct; je++)
            work[ie] += wxp[jj][je][idim] * hsub[je + jj *X->Def.k_exct + nsub*ie];
        }
      }
      for (ie = 0; ie < X->Def.k_exct; ie++) wxp[2][ie][idim] = work[ie];

      for (ie = 0; ie < X->Def.k_exct; ie++) {
        work[ie] = 0.0;
        for (jj = 0; jj < 3; jj += 2)
          for (je = 0; je < X->Def.k_exct; je++)
            work[ie] += hwxp[jj][je][idim] * hsub[je + jj *X->Def.k_exct + nsub*ie];
      }
      for (ie = 0; ie < X->Def.k_exct; ie++) hwxp[2][ie][idim] = work[ie];

    }/*for (idim = 1; idim <= i_max; idim++)*/

    for (ii = 1; ii < 3; ii++) {
      for (ie = 0; ie < X->Def.k_exct; ie++) {
        dnorm = sqrt(creal(vec_prod(i_max, wxp[ii][ie], wxp[ii][ie])));
        for (idim = 1; idim <= i_max; idim++) {
          wxp[ii][ie][idim] = wxp[ii][ie][idim] / dnorm;
          hwxp[ii][ie][idim] = hwxp[ii][ie][idim] / dnorm;
        }/*for (idim = 1; idim <= i_max; idim++)*/
      }/* for (ie = 0; ie < X->Def.k_exct; ie++)*/
    }/*for (ii = 1; ii < 3; ii++)*/

  }/*for (stp = 1; stp <= X->Def.Lanczos_max; stp++)*/

  sprintf(sdt, cFileNameTimeKeep, X->Def.CDataFileHead);
  if (iconv != 0) {
    sprintf(sdt, "%s", cLogLanczos_EigenValueNotConverged);
    return -1;
  }

  TimeKeeper(X, cFileNameTimeKeep, cLanczos_EigenValueFinish, "a");
  fprintf(stdoutMPI, "%s", cLogLanczos_EigenValueEnd);

  c_free1(work, lwork);
  d_free1(rwork, 3 * nsub - 2);
  d_free1(eigsub, nsub);

  d_free1(eig, X->Def.CDataFileHead);
  c_free1(hsub, nsub*nsub);
  c_free1(ovrp, nsub*nsub);

  c_free3(hwxp, 3, X->Def.k_exct, X->Check.idim_max + 1);
  c_malloc1(v0, X->Check.idim_max + 1);
  for (idim = 1; idim <= i_max; idim++) v0[idim] = wxp[1][X->Def.k_exct-1][idim];
  c_free3(wxp, 3, X->Def.k_exct, X->Check.idim_max + 1);
  c_malloc1(v1, X->Check.idim_max + 1);
  c_malloc1(vg, X->Check.idim_max + 1);

  return 0;
}/*int LOBPCG_Main*/
 
int CalcByLOBPCG(
  struct EDMainCalStruct *X
)
{
  char sdt[D_FileNameMax];
  double var;
  long int i_max = 0;
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

    if (LOBPCG_Main(&(X->Bind)) != 0) {
      fprintf(stderr, "  LOBPCG is not converged in this process.\n");
      return(FALSE);
    }

    expec_energy(&(X->Bind));
    var = fabs(X->Bind.Phys.var - X->Bind.Phys.energy*X->Bind.Phys.energy) / fabs(X->Bind.Phys.var);
    fprintf(stdoutMPI, "\n");
    fprintf(stdoutMPI, "  LOBPCG Accuracy check !!!\n");
    fprintf(stdoutMPI, "    variance = %.14e \n ", var);
    fprintf(stdoutMPI, "\n");

  }
  else {// X->Bind.Def.iInputEigenVec=true :input v1:
    fprintf(stdoutMPI, "An Eigenvector is inputted.\n");
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, cReadEigenVecStart, "a");
    sprintf(sdt, cFileNameInputEigen, X->Bind.Def.CDataFileHead, X->Bind.Def.k_exct - 1, myrank);
    childfopenALL(sdt, "rb", &fp);
    if (fp == NULL) {
      fprintf(stderr, "Error: A file of Inputvector does not exist.\n");
      exitMPI(-1);
    }
    fread(&step_i, sizeof(long int), 1, fp);
    fread(&i_max, sizeof(long int), 1, fp);
    if (i_max != X->Bind.Check.idim_max) {
      fprintf(stderr, "Error: A file of Inputvector is incorrect.\n");
      exitMPI(-1);
    }
    fread(v1, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
    fclose(fp);
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, cReadEigenVecFinish, "a");
  }

  fprintf(stdoutMPI, "%s", cLogLanczos_EigenVecEnd);
  // v1 is eigen vector

  if (!expec_cisajs(&(X->Bind), v1) == 0) {
    fprintf(stderr, "Error: calc OneBodyG.\n");
    exitMPI(-1);
  }

  if (!expec_cisajscktaltdc(&(X->Bind), v1) == 0) {
    fprintf(stderr, "Error: calc TwoBodyG.\n");
    exitMPI(-1);
  }

  if (!expec_totalSz(&(X->Bind), v1) == 0) {
    fprintf(stderr, "Error: calc TotalSz.\n");
    exitMPI(-1);
  }
  /*
   Output physical variables to a file
  */
  if (X->Bind.Def.St == 0) {
    sprintf(sdt, cFileNameEnergy_Lanczos, X->Bind.Def.CDataFileHead);
  }
  else if (X->Bind.Def.St == 1) {
    sprintf(sdt, cFileNameEnergy_CG, X->Bind.Def.CDataFileHead);
  }

  if (childfopenMPI(sdt, "w", &fp) != 0) {
    exitMPI(-1);
  }
  fprintf(fp, "Energy  %.16lf \n", X->Bind.Phys.energy);
  fprintf(fp, "Doublon  %.16lf \n", X->Bind.Phys.doublon);
  fprintf(fp, "Sz  %.16lf \n", X->Bind.Phys.sz);
  fclose(fp);
  /*
   Output Eigenvector to a file
  */
  if (X->Bind.Def.iOutputEigenVec == TRUE) {
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, cOutputEigenVecStart, "a");
    sprintf(sdt, cFileNameOutputEigen, X->Bind.Def.CDataFileHead, X->Bind.Def.k_exct - 1, myrank);
    if (childfopenALL(sdt, "wb", &fp) != 0) {
      exitMPI(-1);
    }
    fwrite(&X->Bind.Large.itr, sizeof(X->Bind.Large.itr), 1, fp);
    fwrite(&X->Bind.Check.idim_max, sizeof(X->Bind.Check.idim_max), 1, fp);
    fwrite(v1, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
    fclose(fp);
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, cOutputEigenVecStart, "a");
  }/*if (X->Bind.Def.iOutputEigenVec == TRUE)*/

  return TRUE;

}/*int CalcByLOBPCG*/
