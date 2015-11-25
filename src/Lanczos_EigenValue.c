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

#include "Common.h"
#include "mfmemory.h"
#include "Lanczos_EigenValue.h"
#include "wrapperMPI.h"

/** 
 * 
 * 
 * @param X 
 * 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @return 
 */
int Lanczos_EigenValue(struct BindStruct *X)
{

  fprintf(stdoutMPI, "%s", cLogLanczos_EigenValueStart);
  FILE *fp;
  char sdt[D_FileNameMax],sdt_2[D_FileNameMax];
  int stp;
  long int i,iv,i_max;      
  int k_exct,Target;
  int iconv=-1;
  double beta1,alpha1; //beta,alpha1 should be real
  double  complex temp1,temp2;
  double complex cbeta1;
  double E[5],ebefor;

// for GC
  double dnorm;
  double complex cdnorm;
  long unsigned int u_long_i;
  dsfmt_t dsfmt;

#ifdef lapack
  double **tmp_mat;
  double *tmp_E;
  int    int_i,int_j,mfint[7];
#endif
      
  sprintf(sdt_2, cFileNameLanczosStep, X->Def.CDataFileHead);

  i_max=X->Check.idim_max;      
  k_exct = X->Def.k_exct;

  if(initial_mode == 0){
    X->Large.iv=(X->Check.idim_max/3+X->Def.initial_iv)%X->Check.idim_max+1;
    if(X->Def.iCalcModel==Spin || X->Def.iCalcModel==Kondo){
      X->Large.iv=(X->Check.idim_max/2+X->Def.initial_iv)%X->Check.idim_max+1;
    }
    iv=X->Large.iv;
    fprintf(stdoutMPI, "initial_mode=%d normal: iv = %ld i_max=%ld k_exct =%d \n",initial_mode,iv,i_max,k_exct);       
#pragma omp parallel for default(none) private(i) shared(v0, v1) firstprivate(i_max)
    for(i = 1; i <= i_max; i++){
      v0[i]=0.0;
      v1[i]=0.0;
    }
    v1[iv]=1.0;
    if(X->Def.iInitialVecType==0){
      v1[iv]+=1.0*I;
      v1[iv]/=sqrt(2.0);
    }
  }else if(initial_mode==1){
    iv = X->Def.initial_iv;
    fprintf(stdoutMPI, "initial_mode=%d (random): iv = %ld i_max=%ld k_exct =%d \n",initial_mode,iv,i_max,k_exct);       
    #pragma omp parallel for default(none) private(i) shared(v0, v1) firstprivate(i_max)
    for(i = 1; i <= i_max; i++){
      v0[i]=0.0;
    }
    u_long_i = 123432 + abs(iv);
    dsfmt_init_gen_rand(&dsfmt, u_long_i);
    if(X->Def.iInitialVecType==0){
      for(i = 1; i <= i_max; i++){
	v1[i]=2.0*(dsfmt_genrand_close_open(&dsfmt)-0.5)+2.0*(dsfmt_genrand_close_open(&dsfmt)-0.5)*I;
      }
    }
    else{
       for(i = 1; i <= i_max; i++){
	 v1[i]=2.0*(dsfmt_genrand_close_open(&dsfmt)-0.5);
      }
    }
    cdnorm=0.0;
#pragma omp parallel for default(none) private(i) shared(v1, i_max) reduction(+: cdnorm) 
    for(i=1;i<=i_max;i++){
     cdnorm += conj(v1[i])*v1[i];
    }
    cdnorm = SumMPI_dc(cdnorm);
    dnorm=creal(cdnorm);
    dnorm=sqrt(dnorm);
    #pragma omp parallel for default(none) private(i) shared(v1) firstprivate(i_max, dnorm)
    for(i=1;i<=i_max;i++){
      v1[i] = v1[i]/dnorm;
    }
  }
  //Eigenvalues by Lanczos method

  TimeKeeper(X, cFileNameTimeKeep, cLanczos_EigenValueStart, "a");
  mltply(X, v0, v1);
  stp=1;

  TimeKeeperWithStep(X, cFileNameTimeKeep, cLanczos_EigenValueStep, "a", stp);          
  alpha1=creal(X->Large.prdct) ;// alpha = v^{\dag}*H*v
  alpha[1]=alpha1;
  cbeta1=0.0;
        
#pragma omp parallel for reduction(+:cbeta1) default(none) private(i) shared(v0, v1) firstprivate(i_max, alpha1)
  for(i = 1; i <= i_max; i++){
    cbeta1+=conj(v0[i]-alpha1*v1[i])*(v0[i]-alpha1*v1[i]);
  }
  cbeta1 = SumMPI_dc(cbeta1);
  beta1=creal(cbeta1);
  beta1=sqrt(beta1);
  beta[1]=beta1;
  ebefor=0;
  //  fprintf(stdoutMPI, "alpha[%d]=%lf, beta[%d]=%lf\n", 1, alpha1, 1, beta1);
  
  for(stp = 2; stp <= X->Def.Lanczos_max; stp++){
#pragma omp parallel for default(none) private(i,temp1, temp2) shared(v0, v1) firstprivate(i_max, alpha1, beta1)
    for(i=1;i<=i_max;i++){
      temp1 = v1[i];
      temp2 = (v0[i]-alpha1*v1[i])/beta1;
      v0[i] = -beta1*temp1;
      v1[i] =  temp2;
    }

    mltply(X, v0, v1);

    TimeKeeperWithStep(X, cFileNameTimeKeep, cLanczos_EigenValueStep, "a", stp);
   
    alpha1=creal(X->Large.prdct);
    alpha[stp]=alpha1;
    cbeta1=0.0;

#pragma omp parallel for reduction(+:cbeta1) default(none) private(i) shared(v0, v1) firstprivate(i_max, alpha1)
    for(i=1;i<=i_max;i++){
      cbeta1+=conj(v0[i]-alpha1*v1[i])*(v0[i]-alpha1*v1[i]);
    }
    cbeta1 = SumMPI_dc(cbeta1);
    beta1=creal(cbeta1);
    beta1=sqrt(beta1);
    beta[stp]=beta1;
    
    Target  = X->Def.LanczosTarget;

    //    fprintf(stdoutMPI, "alpha[%d]=%lf, beta[%d]=%lf\n", stp, alpha1, stp, beta1);
    
    if(stp==2){      
     #ifdef lapack
      d_malloc2(tmp_mat,stp,stp);
      d_malloc1(tmp_E,stp+1);

       for(int_i=0;int_i<stp;int_i++){
         for(int_j=0;int_j<stp;int_j++){
           tmp_mat[int_i][int_j] = 0.0;
         }
       } 
       tmp_mat[0][0]   = alpha[1]; 
       tmp_mat[0][1]   = beta[1]; 
       tmp_mat[1][0]   = beta[1]; 
       tmp_mat[1][1]   = alpha[2]; 
       DSEVvalue(stp,tmp_mat,tmp_E);
       E[1] = tmp_E[0];
       E[2] = tmp_E[1];
       E[3] = tmp_E[2];
       E[4] = tmp_E[3];
       d_free1(tmp_E,stp+1);
       d_free2(tmp_mat,stp,stp);
     #else
       bisec(alpha,beta,stp,E,4,eps_Bisec);
     #endif
       ebefor=E[Target];
    }
            
    if(stp>2 && stp%2==0){
      
#ifdef lapack
      d_malloc2(tmp_mat,stp,stp);
      d_malloc1(tmp_E,stp+1);

       for(int_i=0;int_i<stp;int_i++){
         for(int_j=0;int_j<stp;int_j++){
           tmp_mat[int_i][int_j] = 0.0;
         }
       } 
       tmp_mat[0][0]   = alpha[1]; 
       tmp_mat[0][1]   = beta[1]; 
       for(int_i=1;int_i<stp-1;int_i++){
         tmp_mat[int_i][int_i]     = alpha[int_i+1]; 
         tmp_mat[int_i][int_i+1]   = beta[int_i+1]; 
         tmp_mat[int_i][int_i-1]   = beta[int_i]; 
       }
       tmp_mat[int_i][int_i]       = alpha[int_i+1]; 
       tmp_mat[int_i][int_i-1]     = beta[int_i]; 
       DSEVvalue(stp,tmp_mat,tmp_E);
       E[1] = tmp_E[0];
       E[2] = tmp_E[1];
       E[3] = tmp_E[2];
       E[4] = tmp_E[3];
       d_free1(tmp_E,stp+1);
       d_free2(tmp_mat,stp,stp);
#else
       bisec(alpha,beta,stp,E,4,eps_Bisec);
#endif
      

       fprintf(stdoutMPI, "stp=%d %.10lf %.10lf %.10lf %.10lf \n",stp,E[1],E[2],E[3],E[4]);
       if(stp==4){
	 childfopenMPI(sdt_2,"w", &fp);
       }
       else{
	 childfopenMPI(sdt_2,"a", &fp);
       }
       fprintf(fp,"stp=%d %.10lf %.10lf %.10lf %.10lf\n",stp,E[1],E[2],E[3],E[4]);
       fclose(fp);

      if(fabs((E[Target]-ebefor)/E[Target])<eps_Lanczos){
        vec12(alpha,beta,stp,E,X);		
        X->Large.itr=stp;       
        X->Phys.Target_energy=E[k_exct];
	iconv=0;
	break;
      }

      if(abs(beta[stp])<pow(10.0, -14)){
	beta[stp]=pow(10.0, -14)*dShiftBeta;
      }

      ebefor=E[Target];            
    }
  }        


  sprintf(sdt,cFileNameTimeKeep,X->Def.CDataFileHead);
  if(iconv!=0){
    sprintf(sdt,  cLogLanczos_EigenValueNotConverged);
    return -1;
  }

  TimeKeeper(X, cFileNameTimeKeep, cLanczos_EigenValueFinish, "a");
  fprintf(stdoutMPI, "%s", cLogLanczos_EigenValueEnd);
  
  return 0;  
}
