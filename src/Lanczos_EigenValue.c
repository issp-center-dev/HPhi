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
#include "mltply.h"
#include "vec12.h"
#include "bisec.h"
#include "FileIO.h"
#include "matrixlapack.h"
#include "Lanczos_EigenValue.h"
#include "wrapperMPI.h"
#include "CalcTime.h"
/** 
 * 
 * 
 * @param X 
 * 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @return 
 */
int Lanczos_EigenValue(struct BindStruct *X) {

  fprintf(stdoutMPI, "%s", cLogLanczos_EigenValueStart);
  FILE *fp;
  char sdt[D_FileNameMax], sdt_2[D_FileNameMax];
  int stp;
  long int i, i_max;
  unsigned long int i_max_tmp;
  int k_exct, Target;
  int iconv = -1;
  double beta1, alpha1; //beta,alpha1 should be real
  double complex temp1, temp2;
  double complex cbeta1;
  double E[5], ebefor, E_target;

// for GC
  double dnorm=0.0;

  double **tmp_mat;
  double *tmp_E;
  int int_i, int_j, mfint[7];
  int iret=0;
  sprintf(sdt_2, cFileNameLanczosStep, X->Def.CDataFileHead);

  i_max = X->Check.idim_max;
  k_exct = X->Def.k_exct;
  unsigned long int liLanczosStp = X->Def.Lanczos_max;
  unsigned long int liLanczosStp_vec=0;

    if (X->Def.iReStart == RESTART_INOUT || X->Def.iReStart == RESTART_IN){
      if(ReadInitialVector( X, v0, v1, &liLanczosStp_vec)!=0){
        fprintf(stdoutMPI, "  Error: Fail to read initial vectors\n");
        return -2;
      }
      iret=ReadTMComponents(X, &dnorm, &liLanczosStp, 0);
      if(iret !=TRUE){
        fprintf(stdoutMPI, "  Error: Fail to read TMcomponents\n");
        return -2;
      }
      if(liLanczosStp_vec !=liLanczosStp){
        fprintf(stdoutMPI, "  Error: Input files for vector and TMcomponents are incoorect.\n");
        fprintf(stdoutMPI, "  Error: Input vector %ld th stps, TMcomponents  %ld th stps.\n", liLanczosStp_vec, liLanczosStp);
        return -2;
      }
      X->Def.Lanczos_restart=liLanczosStp;
      liLanczosStp = liLanczosStp+X->Def.Lanczos_max;
      alpha1=alpha[X->Def.Lanczos_restart];
      beta1=beta[X->Def.Lanczos_restart];

    }/*X->Def.iReStart == RESTART_INOUT || X->Def.iReStart == RESTART_IN*/
    else {
      SetInitialVector(X, v0, v1);
      //Eigenvalues by Lanczos method
      TimeKeeper(X, cFileNameTimeKeep, cLanczos_EigenValueStart, "a");
      StartTimer(4101);
      mltply(X, v0, v1);
      StopTimer(4101);
      stp = 1;
      TimeKeeperWithStep(X, cFileNameTimeKeep, cLanczos_EigenValueStep, "a", stp);

      alpha1 = creal(X->Large.prdct);// alpha = v^{\dag}*H*v

      alpha[1] = alpha1;
      cbeta1 = 0.0;

#pragma omp parallel for reduction(+:cbeta1) default(none) private(i) shared(v0, v1) firstprivate(i_max, alpha1)
      for (i = 1; i <= i_max; i++) {
        cbeta1 += conj(v0[i] - alpha1 * v1[i]) * (v0[i] - alpha1 * v1[i]);
      }
      cbeta1 = SumMPI_dc(cbeta1);
      beta1 = creal(cbeta1);
      beta1 = sqrt(beta1);
      beta[1] = beta1;
      ebefor = 0;
      liLanczosStp = X->Def.Lanczos_max;
      X->Def.Lanczos_restart =1;
    }/*else restart*/

/*
      Set Maximum number of loop to the dimention of the Wavefunction
    */
  i_max_tmp = SumMPI_li(i_max);
  if (i_max_tmp < liLanczosStp) {
    liLanczosStp = i_max_tmp;
  }
  if (i_max_tmp < X->Def.LanczosTarget) {
    liLanczosStp = i_max_tmp;
  }
  if (i_max_tmp == 1) {
    E[1] = alpha[1];
    StartTimer(4102);
    vec12(alpha, beta, stp, E, X);
    StopTimer(4102);
    X->Large.itr = stp;
    X->Phys.Target_energy = E[k_exct];
    iconv = 0;
    fprintf(stdoutMPI, "  LanczosStep  E[1] \n");
    fprintf(stdoutMPI, "  stp=%d %.10lf \n", stp, E[1]);
  }

  fprintf(stdoutMPI, "  LanczosStep  E[1] E[2] E[3] E[4] Target:E[%d] E_Max/Nsite\n", X->Def.LanczosTarget + 1);
  for (stp = X->Def.Lanczos_restart+1; stp <= liLanczosStp; stp++) {
#pragma omp parallel for default(none) private(i,temp1, temp2) shared(v0, v1) firstprivate(i_max, alpha1, beta1)
    for (i = 1; i <= i_max; i++) {
      temp1 = v1[i];
      temp2 = (v0[i] - alpha1 * v1[i]) / beta1;
      v0[i] = -beta1 * temp1;
      v1[i] = temp2;
    }

    StartTimer(4101);
    mltply(X, v0, v1);
    StopTimer(4101);
    TimeKeeperWithStep(X, cFileNameTimeKeep, cLanczos_EigenValueStep, "a", stp);
    alpha1 = creal(X->Large.prdct);
    alpha[stp] = alpha1;
    cbeta1 = 0.0;

#pragma omp parallel for reduction(+:cbeta1) default(none) private(i) shared(v0, v1) firstprivate(i_max, alpha1)
    for (i = 1; i <= i_max; i++) {
      cbeta1 += conj(v0[i] - alpha1 * v1[i]) * (v0[i] - alpha1 * v1[i]);
    }
    cbeta1 = SumMPI_dc(cbeta1);
    beta1 = creal(cbeta1);
    beta1 = sqrt(beta1);
    beta[stp] = beta1;

    Target = X->Def.LanczosTarget;

    if (stp == 2) {
      d_malloc2(tmp_mat, stp, stp);
      d_malloc1(tmp_E, stp + 1);

      for (int_i = 0; int_i < stp; int_i++) {
        for (int_j = 0; int_j < stp; int_j++) {
          tmp_mat[int_i][int_j] = 0.0;
        }
      }
      tmp_mat[0][0] = alpha[1];
      tmp_mat[0][1] = beta[1];
      tmp_mat[1][0] = beta[1];
      tmp_mat[1][1] = alpha[2];
      DSEVvalue(stp, tmp_mat, tmp_E);
      E[1] = tmp_E[0];
      E[2] = tmp_E[1];
      if (Target < 2) {
        E_target = tmp_E[Target];
        ebefor = E_target;
      }
      d_free1(tmp_E, stp + 1);
      d_free2(tmp_mat, stp, stp);

      childfopenMPI(sdt_2, "w", &fp);

      fprintf(fp, "LanczosStep  E[1] E[2] E[3] E[4] Target:E[%d] E_Max/Nsite\n", Target + 1);
      if (Target < 2) {
        fprintf(stdoutMPI, "  stp = %d %.10lf %.10lf xxxxxxxxxx xxxxxxxxx  %.10lf xxxxxxxxx \n", stp, E[1], E[2],
                E_target);
        fprintf(fp, "stp = %d %.10lf %.10lf xxxxxxxxxx xxxxxxxxx %.10lf xxxxxxxxx \n", stp, E[1], E[2], E_target);
      } else {
        fprintf(stdoutMPI, "  stp = %d %.10lf %.10lf xxxxxxxxxx xxxxxxxxx xxxxxxxxx xxxxxxxxx \n", stp, E[1], E[2]);
        fprintf(fp, "stp = %d %.10lf %.10lf xxxxxxxxxx xxxxxxxxx xxxxxxxxx xxxxxxxxx \n", stp, E[1], E[2]);
      }

      fclose(fp);
    }

    if (stp > 2 && stp % 2 == 0) {

      childfopenMPI(sdt_2, "a", &fp);

      d_malloc2(tmp_mat, stp, stp);
      d_malloc1(tmp_E, stp + 1);

      for (int_i = 0; int_i < stp; int_i++) {
        for (int_j = 0; int_j < stp; int_j++) {
          tmp_mat[int_i][int_j] = 0.0;
        }
      }
      tmp_mat[0][0] = alpha[1];
      tmp_mat[0][1] = beta[1];
      for (int_i = 1; int_i < stp - 1; int_i++) {
        tmp_mat[int_i][int_i] = alpha[int_i + 1];
        tmp_mat[int_i][int_i + 1] = beta[int_i + 1];
        tmp_mat[int_i][int_i - 1] = beta[int_i];
      }
      tmp_mat[int_i][int_i] = alpha[int_i + 1];
      tmp_mat[int_i][int_i - 1] = beta[int_i];
      StartTimer(4103);
      DSEVvalue(stp, tmp_mat, tmp_E);
      StopTimer(4103);
      E[1] = tmp_E[0];
      E[2] = tmp_E[1];
      E[3] = tmp_E[2];
      E[4] = tmp_E[3];
      E[0] = tmp_E[stp - 1];
      if (stp > Target) {
        E_target = tmp_E[Target];
      }
      d_free1(tmp_E, stp + 1);
      d_free2(tmp_mat, stp, stp);
      if (stp > Target) {
        fprintf(stdoutMPI, "  stp = %d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n", stp, E[1], E[2], E[3], E[4],
                E_target, E[0] / (double) X->Def.NsiteMPI);
        fprintf(fp, "stp=%d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n", stp, E[1], E[2], E[3], E[4], E_target,
                E[0] / (double) X->Def.NsiteMPI);
      } else {
        fprintf(stdoutMPI, "  stp = %d %.10lf %.10lf %.10lf %.10lf xxxxxxxxx %.10lf\n", stp, E[1], E[2], E[3], E[4],
                E[0] / (double) X->Def.NsiteMPI);
        fprintf(fp, "stp=%d %.10lf %.10lf %.10lf %.10lf xxxxxxxxx %.10lf\n", stp, E[1], E[2], E[3], E[4],
                E[0] / (double) X->Def.NsiteMPI);
      }
      fclose(fp);
      if (stp > Target) {
        if (fabs((E_target - ebefor) / E_target) < eps_Lanczos || fabs(beta[stp]) < pow(10.0, -14)) {
          if(X->Def.iReStart == RESTART_INOUT ||X->Def.iReStart == RESTART_OUT){
            break;
          }
          d_malloc1(tmp_E, stp + 1);
          StartTimer(4102);
          vec12(alpha, beta, stp, tmp_E, X);
          StopTimer(4102);
          X->Large.itr = stp;
          X->Phys.Target_energy = E_target;
          X->Phys.Target_CG_energy = E[k_exct]; //for CG
          iconv = 0;
          d_free1(tmp_E, stp + 1);
          break;
        }
        ebefor = E_target;
      }

    }
  }
  if (X->Def.iReStart == RESTART_INOUT ||X->Def.iReStart == RESTART_OUT ){
    OutputTMComponents(X, alpha,beta, dnorm, stp-1);
    OutputLanczosVector(X, v0, v1, stp-1);
    return 0;
  }

  sprintf(sdt, cFileNameTimeKeep, X->Def.CDataFileHead);
  if (iconv != 0) {
    sprintf(sdt, "%s", cLogLanczos_EigenValueNotConverged);
    fprintf(stderr, "  Lanczos Eigenvalue is not converged in this process.\n");
    return -1;
  }

  TimeKeeper(X, cFileNameTimeKeep, cLanczos_EigenValueFinish, "a");
  fprintf(stdoutMPI, "%s", cLogLanczos_EigenValueEnd);

  return 0;
}

/** 
 * @brief Calculate tridiagonal matrix components by Lanczos method
 * 
 * @param _alpha 
 * @param _beta 
 * @param _v1 
 * @param Lanczos_step 
 * 
 * @return 0
 */
int Lanczos_GetTridiagonalMatrixComponents(
        struct BindStruct *X,
        double *_alpha,
        double *_beta,
        double complex *tmp_v1,
        unsigned long int *liLanczos_step
 )
{

 char sdt[D_FileNameMax];
 int stp;
 long int i,i_max;
 i_max=X->Check.idim_max;
 
 unsigned long int i_max_tmp;
 double beta1,alpha1; //beta,alpha1 should be real
 double  complex temp1,temp2;
 double complex cbeta1;

 sprintf(sdt, cFileNameLanczosStep, X->Def.CDataFileHead);  
  
  /*
    Set Maximum number of loop to the dimension of the Wavefunction
  */
 i_max_tmp = SumMPI_li(i_max);
 if(i_max_tmp < *liLanczos_step){
    *liLanczos_step = i_max_tmp;
  }
  if(i_max_tmp < X->Def.LanczosTarget){
    *liLanczos_step = i_max_tmp;
  }

  if(X->Def.Lanczos_restart==0) {
#pragma omp parallel for default(none) private(i) shared(v0, v1, tmp_v1) firstprivate(i_max)
        for (i = 1; i <= i_max; i++) {
            v0[i] = 0.0;
            v1[i] = tmp_v1[i];
        }
        stp = 0;
        mltply(X, v0, tmp_v1);
        TimeKeeperWithStep(X, cFileNameTimeKeep, c_Lanczos_SpectrumStep, "a", stp);
        alpha1 = creal(X->Large.prdct);// alpha = v^{\dag}*H*v
        _alpha[1] = alpha1;
        cbeta1 = 0.0;

#pragma omp parallel for reduction(+:cbeta1) default(none) private(i) shared(v0, v1) firstprivate(i_max, alpha1)
        for (i = 1; i <= i_max; i++) {
            cbeta1 += conj(v0[i] - alpha1 * v1[i]) * (v0[i] - alpha1 * v1[i]);
        }
        cbeta1 = SumMPI_dc(cbeta1);
        beta1 = creal(cbeta1);
        beta1 = sqrt(beta1);
        _beta[1] = beta1;
        X->Def.Lanczos_restart =1;
    }
else{
        alpha1=alpha[X->Def.Lanczos_restart];
        beta1=beta[X->Def.Lanczos_restart];
    }
//  for(stp = 2; stp <= *liLanczos_step; stp++){
    for(stp = X->Def.Lanczos_restart+1; stp <= *liLanczos_step; stp++){
      if(fabs(beta[stp-1])<pow(10.0, -14)){
          *liLanczos_step=stp-1;
          break;
      }

#pragma omp parallel for default(none) private(i,temp1, temp2) shared(v0, v1) firstprivate(i_max, alpha1, beta1)
      for(i=1;i<=i_max;i++){
	temp1 = v1[i];
	temp2 = (v0[i]-alpha1*v1[i])/beta1;
	v0[i] = -beta1*temp1;
	v1[i] =  temp2;
      }

      mltply(X, v0, v1);
      TimeKeeperWithStep(X, cFileNameTimeKeep, c_Lanczos_SpectrumStep, "a", stp);
      alpha1=creal(X->Large.prdct);
      _alpha[stp]=alpha1;
      cbeta1=0.0;
      
#pragma omp parallel for reduction(+:cbeta1) default(none) private(i) shared(v0, v1) firstprivate(i_max, alpha1)
      for(i=1;i<=i_max;i++){
	    cbeta1+=conj(v0[i]-alpha1*v1[i])*(v0[i]-alpha1*v1[i]);
      }
      cbeta1 = SumMPI_dc(cbeta1);
      beta1=creal(cbeta1);
      beta1=sqrt(beta1);
      _beta[stp]=beta1;
  }

  return TRUE;
}

int ReadInitialVector(struct BindStruct *X, double complex* _v0, double complex *_v1, unsigned long int *liLanczosStp_vec)
{
  size_t byte_size;
  char sdt[D_FileNameMax];
  FILE *fp;
  unsigned long int i_max;

  fprintf(stdoutMPI, "  Start: Input vectors for recalculation.\n");
  TimeKeeper(X, cFileNameTimeKeep, c_InputSpectrumRecalcvecStart, "a");
  sprintf(sdt, cFileNameOutputRestartVec, X->Def.CDataFileHead, myrank);
  if (childfopenALL(sdt, "rb", &fp) != 0) {
    return -1;
  }
  byte_size = fread(liLanczosStp_vec, sizeof(*liLanczosStp_vec),1,fp);
  byte_size = fread(&i_max, sizeof(long int), 1, fp);
  if(i_max != X->Check.idim_max){
    fprintf(stderr, "Error: A size of Inputvector is incorrect.\n");
    return -1;
  }
  byte_size = fread(_v0, sizeof(complex double), X->Check.idim_max + 1, fp);
  byte_size = fread(_v1, sizeof(complex double), X->Check.idim_max + 1, fp);
  fclose(fp);
  fprintf(stdoutMPI, "  End:   Input vectors for recalculation.\n");
  TimeKeeper(X, cFileNameTimeKeep, c_InputSpectrumRecalcvecEnd, "a");
  if (byte_size == 0) printf("byte_size: %d \n", (int)byte_size);
  return 0;
}

int OutputLanczosVector(struct BindStruct *X,
                        double complex* tmp_v0,
                        double complex *tmp_v1,
                        unsigned long int liLanczosStp_vec){
  char sdt[D_FileNameMax];
  FILE *fp;

  fprintf(stdoutMPI, "    Start: Output vectors for recalculation.\n");
  TimeKeeper(X, cFileNameTimeKeep, c_OutputSpectrumRecalcvecStart, "a");

  sprintf(sdt, cFileNameOutputRestartVec, X->Def.CDataFileHead, myrank);
  if(childfopenALL(sdt, "wb", &fp)!=0){
    return -1;
  }
  fwrite(&liLanczosStp_vec, sizeof(liLanczosStp_vec),1,fp);
  fwrite(&X->Check.idim_max, sizeof(X->Check.idim_max),1,fp);
  fwrite(tmp_v0, sizeof(complex double),X->Check.idim_max+1, fp);
  fwrite(tmp_v1, sizeof(complex double),X->Check.idim_max+1, fp);
  fclose(fp);

  fprintf(stdoutMPI, "    End:   Output vectors for recalculation.\n");
  TimeKeeper(X, cFileNameTimeKeep, c_OutputSpectrumRecalcvecEnd, "a");
  return 0;
}


int SetInitialVector(struct BindStruct *X, double complex* tmp_v0, double complex *tmp_v1) {
  int iproc;
  long int i, iv, i_max;
  unsigned long int i_max_tmp, sum_i_max;
  int mythread;

// for GC
  double dnorm;
  double complex cdnorm;
  long unsigned int u_long_i;
  dsfmt_t dsfmt;

  i_max = X->Check.idim_max;
  if (initial_mode == 0) {
    if(X->Def.iFlgMPI==0) {
      sum_i_max = SumMPI_li(X->Check.idim_max);
    }
    else{
      sum_i_max =X->Check.idim_max;
    }
    X->Large.iv = (sum_i_max / 2 + X->Def.initial_iv) % sum_i_max + 1;
    iv = X->Large.iv;
    fprintf(stdoutMPI, "  initial_mode=%d normal: iv = %ld i_max=%ld k_exct =%d \n\n", initial_mode, iv, i_max,
            X->Def.k_exct);
#pragma omp parallel for default(none) private(i) shared(tmp_v0, tmp_v1) firstprivate(i_max)
    for (i = 1; i <= i_max; i++) {
      tmp_v0[i] = 0.0;
      tmp_v1[i] = 0.0;
    }

    sum_i_max = 0;
    if(X->Def.iFlgMPI==0) {
      for (iproc = 0; iproc < nproc; iproc++) {
        i_max_tmp = BcastMPI_li(iproc, i_max);
        if (sum_i_max <= iv && iv < sum_i_max + i_max_tmp) {
          if (myrank == iproc) {
            tmp_v1[iv - sum_i_max + 1] = 1.0;
            if (X->Def.iInitialVecType == 0) {
              tmp_v1[iv - sum_i_max + 1] += 1.0 * I;
              tmp_v1[iv - sum_i_max + 1] /= sqrt(2.0);
            }
          }/*if (myrank == iproc)*/
        }/*if (sum_i_max <= iv && iv < sum_i_max + i_max_tmp)*/

        sum_i_max += i_max_tmp;

      }/*for (iproc = 0; iproc < nproc; iproc++)*/
    }
    else {
      tmp_v1[iv + 1] = 1.0;
      if (X->Def.iInitialVecType == 0) {
        tmp_v1[iv + 1] += 1.0 * I;
        tmp_v1[iv + 1] /= sqrt(2.0);
      }
    }
  }/*if(initial_mode == 0)*/
  else if (initial_mode == 1) {
    iv = X->Def.initial_iv;
    fprintf(stdoutMPI, "  initial_mode=%d (random): iv = %ld i_max=%ld k_exct =%d \n\n", initial_mode, iv, i_max,
            X->Def.k_exct);
#pragma omp parallel default(none) private(i, u_long_i, mythread, dsfmt) \
            shared(tmp_v0, tmp_v1, iv, X, nthreads, myrank) firstprivate(i_max)
    {

#pragma omp for
      for (i = 1; i <= i_max; i++) {
        tmp_v0[i] = 0.0;
      }
      /*
       Initialise MT
      */
#ifdef _OPENMP
      mythread = omp_get_thread_num();
#else
      mythread = 0;
#endif
      if(X->Def.iFlgMPI==0) {
        u_long_i = 123432 + labs(iv) + mythread + nthreads * myrank;
      }
      else{
        u_long_i = 123432 + labs(iv)+ mythread ;
      }
      dsfmt_init_gen_rand(&dsfmt, u_long_i);

      if (X->Def.iInitialVecType == 0) {
#pragma omp for
        for (i = 1; i <= i_max; i++)
          tmp_v1[i] = 2.0 * (dsfmt_genrand_close_open(&dsfmt) - 0.5) +
                      2.0 * (dsfmt_genrand_close_open(&dsfmt) - 0.5) * I;
      } else {
#pragma omp for
        for (i = 1; i <= i_max; i++)
          tmp_v1[i] = 2.0 * (dsfmt_genrand_close_open(&dsfmt) - 0.5);
      }

    }/*#pragma omp parallel*/

    cdnorm = 0.0;
#pragma omp parallel for default(none) private(i) shared(tmp_v1, i_max) reduction(+: cdnorm)
    for (i = 1; i <= i_max; i++) {
      cdnorm += conj(tmp_v1[i]) * tmp_v1[i];
    }
    if(X->Def.iFlgMPI==0) {
      cdnorm = SumMPI_dc(cdnorm);
    }
    dnorm = creal(cdnorm);
    dnorm = sqrt(dnorm);
#pragma omp parallel for default(none) private(i) shared(tmp_v1) firstprivate(i_max, dnorm)
    for (i = 1; i <= i_max; i++) {
      tmp_v1[i] = tmp_v1[i] / dnorm;
    }
  }/*else if(initial_mode==1)*/
  return 0;
}

int ReadTMComponents(
        struct BindStruct *X,
        double *_dnorm,
        unsigned long int *_i_max,
        int iFlg
){
  char sdt[D_FileNameMax];
  char ctmp[256];

  unsigned long int idx;
  unsigned long int i_max;
  double dnorm;
  FILE *fp;
  idx=1;
  sprintf(sdt, cFileNameTridiagonalMatrixComponents, X->Def.CDataFileHead);
  childfopenMPI(sdt,"r", &fp);

  fgetsMPI(ctmp, sizeof(ctmp)/sizeof(char), fp);
  sscanf(ctmp,"%ld \n", &i_max);
  if(iFlg==0) {
    alpha = (double *) realloc(alpha, sizeof(double) * (i_max + X->Def.Lanczos_max + 1));
    beta = (double *) realloc(beta, sizeof(double) * (i_max + X->Def.Lanczos_max + 1));
  }
  else if(iFlg==1){
    alpha=(double*)realloc(alpha, sizeof(double)*(i_max + 1));
    beta=(double*)realloc(beta, sizeof(double)*(i_max + 1));
  }
  else{
    fclose(fp);
    return -1;
  }
  fgetsMPI(ctmp, sizeof(ctmp)/sizeof(char), fp);
  sscanf(ctmp,"%lf \n", &dnorm);
  while(fgetsMPI(ctmp, sizeof(ctmp)/sizeof(char), fp) != NULL){
    sscanf(ctmp,"%lf %lf \n", &alpha[idx], &beta[idx]);
    idx++;
  }
  fclose(fp);
  *_dnorm=dnorm;
  *_i_max=i_max;
  return TRUE;
}


int OutputTMComponents(
        struct BindStruct *X,
        double *_alpha,
        double *_beta,
        double _dnorm,
        int liLanczosStp
)
{
  char sdt[D_FileNameMax];
  unsigned long int i;
  FILE *fp;

  sprintf(sdt, cFileNameTridiagonalMatrixComponents, X->Def.CDataFileHead);
  childfopenMPI(sdt,"w", &fp);
  fprintf(fp, "%d \n",liLanczosStp);
  fprintf(fp, "%.10lf \n",_dnorm);
  for( i = 1 ; i <= liLanczosStp; i++){
//    fprintf(stdoutMPI, "Debug: %d, %.10lf %.10lf \n",i, _alpha[i], _beta[i]);
    fprintf(fp,"%.10lf %.10lf \n", _alpha[i], _beta[i]);
  }
  fclose(fp);
  return TRUE;
}