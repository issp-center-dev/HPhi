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

void zhegv_(int *itype, char *jobz, char *uplo, int *n, double complex *a, int *lda, double complex *b, int *ldb, double *w, double complex *work, int *lwork, double *rwork, int *info);

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
  long int i, iv, i_max;
  unsigned long int i_max_tmp, sum_i_max;
  int iconv = -1;
  int mythread;

  int ii, jj, jtarget, k_exct;
  double complex **wxp, **hwxp, hsub[9], ovrp[9];
  double eig, dnorm;
  /*
  Variables for zhegv
  */
  int lwork = 5, info, lda = 3, ldb = 3;
  int itype = 1, nsub;
  double rwork[7], eigsub[3];
  double complex work[5];
  char jobz = 'V', uplo = 'U';
  /*
  For DSFMT
  */
  long unsigned int u_long_i;
  dsfmt_t dsfmt;

  sprintf(sdt_2, cFileNameLanczosStep, X->Def.CDataFileHead);

  i_max = X->Check.idim_max;

  c_malloc2(wxp, 3, X->Check.idim_max + 1);
  c_malloc2(hwxp, 3, X->Check.idim_max + 1);

  if (initial_mode == 0) {

    sum_i_max = SumMPI_li(X->Check.idim_max);
    X->Large.iv = (sum_i_max / 2 + X->Def.initial_iv) % sum_i_max + 1;
    iv = X->Large.iv;
    fprintf(stdoutMPI, "  initial_mode=%d normal: iv = %ld i_max=%ld k_exct =%d \n\n", initial_mode, iv, i_max, k_exct);
#pragma omp parallel for default(none) private(i) shared(wxp,i_max)
    for (i = 1; i <= i_max; i++) {
      wxp[1][i] = 0.0;
    }

    sum_i_max = 0;
    for (iproc = 0; iproc < nproc; iproc++) {

      i_max_tmp = BcastMPI_li(iproc, i_max);
      if (sum_i_max <= iv && iv < sum_i_max + i_max_tmp) {

        if (myrank == iproc) {
          wxp[1][iv - sum_i_max + 1] = 1.0;
          if (X->Def.iInitialVecType == 0) {
            wxp[1][iv - sum_i_max + 1] += 1.0*I;
            wxp[1][iv - sum_i_max + 1] /= sqrt(2.0);
          }
        }/*if (myrank == iproc)*/
      }/*if (sum_i_max <= iv && iv < sum_i_max + i_max_tmp)*/

      sum_i_max += i_max_tmp;

    }/*for (iproc = 0; iproc < nproc; iproc++)*/
  }/*if(initial_mode == 0)*/
  else if (initial_mode == 1) {
    iv = X->Def.initial_iv;
    fprintf(stdoutMPI, "  initial_mode=%d (random): iv = %ld i_max=%ld k_exct =%d \n\n", initial_mode, iv, i_max, k_exct);
#pragma omp parallel default(none) private(i, u_long_i, mythread, dsfmt) \
            shared(wxp, iv, X, nthreads, myrank) firstprivate(i_max)
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

      if (X->Def.iInitialVecType == 0) {
#pragma omp for
        for (i = 1; i <= i_max; i++)
          wxp[1][i] = 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5) + 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5)*I;
      }
      else {
#pragma omp for
        for (i = 1; i <= i_max; i++)
          wxp[1][i] = 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5);
      }

    }/*#pragma omp parallel*/

    dnorm = sqrt(creal(vec_prod(i_max, wxp[1], wxp[1])));
#pragma omp parallel for default(none) private(i) shared(wxp,i_max, dnorm)
    for (i = 1; i <= i_max; i++)  wxp[1][i] = wxp[1][i] / dnorm;
  }/*else if(initial_mode==1)*/

   //Eigenvalues by Lanczos method
  TimeKeeper(X, cFileNameTimeKeep, cLanczos_EigenValueStart, "a");

  for (i = 1; i <= i_max; i++)  hwxp[1][i] = 0.0;
  mltply(X, hwxp[1], wxp[1]);
  stp = 1;
  TimeKeeperWithStep(X, cFileNameTimeKeep, cLanczos_EigenValueStep, "a", stp);

  for (i = 1; i <= i_max; i++) wxp[2][i] = 0.0;
  for (i = 1; i <= i_max; i++) hwxp[2][i] = 0.0;

  eig = creal(vec_prod(i_max, wxp[1], hwxp[1]));

  for (stp = 1; stp <= X->Def.Lanczos_max; stp++) {

    for (i = 1; i <= i_max; i++) wxp[0][i] = hwxp[1][i] - eig * wxp[1][i];
    dnorm = creal(vec_prod(i_max, wxp[0], wxp[0]));

    /*
    Convergence check
    */
    fprintf(stdoutMPI, "%8d %15.5e %15.5e\n", stp, dnorm, eig);
    if (dnorm < eps_Lanczos) break;

    if (stp /= 1)
      for (i = 1; i <= i_max; i++) wxp[0][i] = wxp[0][i] / dnorm;

    for (i = 1; i <= i_max; i++)  hwxp[0][i] = 0.0;
    mltply(X, hwxp[0], wxp[0]);

    for (ii = 0; ii < 3; ii++) {
      for (jj = 0; jj < 3; jj++) {
        hsub[jj + ii * 3] = vec_prod(i_max, wxp[jj], hwxp[ii]);
        ovrp[jj + ii * 3] = vec_prod(i_max, wxp[jj], wxp[ii]);
      }
    }
    eig = creal(hsub[2 + 2 * 3]);

    if (stp == 1) {
      nsub = 2;
      zhegv_(&itype, &jobz, &uplo, &nsub, hsub, &lda, ovrp, &ldb, eigsub, work, &lwork, rwork, &info);
      jtarget = 0;
    }
    else {
      nsub = 3;
      zhegv_(&itype, &jobz, &uplo, &nsub, hsub, &lda, ovrp, &ldb, eigsub, work, &lwork, rwork, &info);
      jtarget = 0;
    }

    eig = 0.5 * (eig + eigsub[jtarget]);

    for (i = 1; i <= i_max; i++) {
      wxp[1][i] = hsub[0 + jtarget * 3] * wxp[0][i]
                + hsub[1 + jtarget * 3] * wxp[1][i]
                + hsub[2 + jtarget * 3] * wxp[2][i];
      hwxp[1][i] = hsub[0 + jtarget * 3] * hwxp[0][i]
                 + hsub[1 + jtarget * 3] * hwxp[1][i]
                 + hsub[2 + jtarget * 3] * hwxp[2][i];

      wxp[2][i] = hsub[0 + jtarget * 3] * wxp[0][i]
                + hsub[2 + jtarget * 3] * wxp[2][i];
      hwxp[2][i] = hsub[0 + jtarget * 3] * hwxp[0][i]
                 + hsub[2 + jtarget * 3] * hwxp[2][i];
    }

    for (ii = 1; ii < 3; ii++) {
      dnorm = sqrt(creal(vec_prod(i_max, wxp[ii], wxp[ii])));
      for (i = 1; i <= i_max; i++) {
        wxp[ii][i] = wxp[ii][i] / dnorm;
        hwxp[ii][i] = hwxp[ii][i] / dnorm;
      }/*for (i = 1; i <= i_max; i++)*/
    }/*for (ii = 1; ii < 3; ii++)*/

  }/*for (stp = 1; stp <= X->Def.Lanczos_max; stp++)*/

  sprintf(sdt, cFileNameTimeKeep, X->Def.CDataFileHead);
  if (iconv != 0) {
    sprintf(sdt, "%s", cLogLanczos_EigenValueNotConverged);
    return -1;
  }

  TimeKeeper(X, cFileNameTimeKeep, cLanczos_EigenValueFinish, "a");
  fprintf(stdoutMPI, "%s", cLogLanczos_EigenValueEnd);

  for (i = 1; i <= i_max; i++) v0[i] = wxp[1][i];
  c_free2(wxp, 3, X->Check.idim_max + 1);
  c_free2(hwxp, 3, X->Check.idim_max + 1);

  return 0;
}/*int LOBPCG_Main*/
 
int CalcByLOBPCG(
  struct EDMainCalStruct *X
)
{
  char sdt[D_FileNameMax];
  double diff_ene, var;
  long int i_max = 0;
  FILE *fp;

  if (X->Bind.Def.iInputEigenVec == FALSE) {

    if (LOBPCG_Main(&(X->Bind)) != 0) {
      fprintf(stderr, "  LOBPCG is not converged in this process.\n");
      return(FALSE);
    }

    expec_energy(&(X->Bind));
    var = fabs(X->Bind.Phys.var - X->Bind.Phys.energy*X->Bind.Phys.energy) / fabs(X->Bind.Phys.var);
    diff_ene = fabs(X->Bind.Phys.Target_energy - X->Bind.Phys.energy) / fabs(X->Bind.Phys.Target_energy);
    fprintf(stdoutMPI, "\n");
    fprintf(stdoutMPI, "  CG Accuracy check !!!\n");
    fprintf(stdoutMPI, "  LanczosEnergy = %.14e\n  EnergyByVec   = %.14e\n  diff_ene      = %.14e\n  var           = %.14e \n ", X->Bind.Phys.Target_energy, X->Bind.Phys.energy, diff_ene, var);
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
