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
@brief Inversed power method with CG
*/
#include "CG_EigenVector.h"
#include "FileIO.h"
#include "mltply.h"
#include "wrapperMPI.h"
#include "CalcTime.h"
/**
@brief inversed power method with CG
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
@return -1 if file can not be opened, 0 for other.
*/
int CG_EigenVector(struct BindStruct *X/**<[inout]*/){

  fprintf(stdoutMPI, "%s", cLogCG_EigenVecStart);
  TimeKeeper(X, cFileNameTimeKeep, cCG_EigenVecStart, "a");

  time_t start,mid;
  FILE *fp_0;
  char sdt_1[D_FileNameMax];
  dsfmt_t dsfmt;
  long unsigned int u_long_i;
  int mythread;

  long int i,j;
  double Eig;
  int i_itr,itr,iv,itr_max;
  int t_itr;
  double bnorm,xnorm,rnorm,rnorm2;
  double complex alpha,beta,xb,rp,yp,gosa1,tmp_r,gosa2;
  double complex *y,*b;
  long int L_size;
  long int i_max;

  i_max=X->Check.idim_max;
  Eig=X->Phys.Target_CG_energy;
    
  strcpy(sdt_1, cFileNameTimeEV_CG);
  if(childfopenMPI(sdt_1, "w", &fp_0) !=0){
    return -1;
  }
    
  L_size=sizeof(double complex)*(i_max+1);
    
  b=(double complex *)malloc(L_size);
  y=(double complex *)malloc(L_size);
    
  if(y==NULL){
    fprintf(fp_0,"BAD in CG_EigenVector  \n");
  }else{
    fprintf(fp_0,"allocate succeed !!! \n");
  }
  fclose(fp_0);
        
  start=time(NULL);  
  /*
    add random components
  */
  iv = X->Def.initial_iv;
  bnorm = 0.0;
#pragma omp parallel default(none) private(i, u_long_i, mythread, dsfmt) \
  shared(v0, v1, iv, X, nthreads, myrank, b, bnorm) firstprivate(i_max)
  {
    /*
      Initialise MT
    */
#ifdef _OPENMP
    mythread = omp_get_thread_num();
#else
    mythread = 0;
#endif
    u_long_i = 123432 + abs(iv) + mythread + nthreads * myrank;
    dsfmt_init_gen_rand(&dsfmt, u_long_i);

#pragma omp for
    for (i = 1; i <= i_max; i++) {
      v0[i] = v1[i];
      b[i] = v0[i];
    }

#pragma omp for reduction(+:bnorm)
    for (i = 1; i <= i_max; i++) {
      b[i] += 2.0*(dsfmt_genrand_close_open(&dsfmt) - 0.5)*0.001;
      bnorm += conj(b[i])*b[i];
    }
  }/*#pragma omp*/
  /*
    Normalize b
  */
  bnorm = SumMPI_d(bnorm);
  bnorm=sqrt(bnorm);
  
#pragma omp parallel for default(none) private(i) shared(b) firstprivate(i_max,bnorm)
  for(i=1;i<=i_max;i++){
    b[i]=b[i]/bnorm;
  }

  t_itr=0;
  for(i_itr=0;i_itr<=50;i_itr++){
    //CG start!!
    bnorm=0.0;
    //initialization
#pragma omp parallel for reduction(+:bnorm) default(none) private(i) shared(b, v1, vg, v0) firstprivate(i_max)
    for(i=1;i<=i_max;i++){
      bnorm+=conj(b[i])*b[i];
      v1[i]=b[i];
      vg[i]=b[i];
      v0[i]=0.0;
    }
    bnorm = SumMPI_d(bnorm);
    if(iv >= 0){
      childfopenMPI(sdt_1,"a",&fp_0);    
      fprintf(fp_0,"b[%d]=%lf bnorm== %lf \n ",iv,creal(b[iv]),bnorm);
      fclose(fp_0);           
    }
    //iteration
    if(i_itr==0){
      itr_max=500;
    }else{
      itr_max=500;
    }
  
    for(itr=1;itr<=itr_max;itr++){
#pragma omp parallel for default(none) private(j) shared(y, vg) firstprivate(i_max, Eig,eps_CG)
      for(j=1;j<=i_max;j++){  
        y[j]=(-Eig+eps_CG)*vg[j];   //y = -E*p
      }
      StartTimer(4401);
      mltply(X,y,vg);      // y += H*p
      StopTimer(4401);
      // (H-E)p=y finish!
      rp=0.0;
      yp=0.0;
#pragma omp parallel for reduction(+:rp, yp) default(none) private(i) shared(v1, vg, y) firstprivate(i_max) 
      for(i=1;i<=i_max;i++){
        rp+=v1[i]*conj(v1[i]);
        yp+=y[i]*conj(vg[i]);
      }
      rp = SumMPI_dc(rp);
      yp = SumMPI_dc(yp);
      alpha=rp/yp;
      rnorm=0.0;
#pragma omp parallel for reduction(+:rnorm) default(none) private(i) shared(v0, v1, vg)firstprivate(i_max, alpha) 
      for(i=1;i<=i_max;i++){
        v0[i]+=alpha*vg[i];
        rnorm+=conj(v1[i])*v1[i];
      }
      rnorm = SumMPI_d(rnorm);
      rnorm2=0.0;
      gosa1=0.0;
#pragma omp parallel for reduction(+:rnorm2, gosa1) default(none) private(i) shared(v1 , y) firstprivate(i_max, alpha) private(tmp_r) 
      for(i=1;i<=i_max;i++){
        tmp_r=v1[i]-alpha*y[i];
        gosa1+=conj(tmp_r)*v1[i];// old r and new r should be orthogonal -> gosa1=0
        v1[i]=tmp_r; 
        rnorm2+=conj(v1[i])*v1[i];
      }
      gosa1 = SumMPI_dc(gosa1);
      rnorm2 = SumMPI_d(rnorm2);
            
      gosa2=0.0;
#pragma omp parallel for reduction(+:gosa2) default(none) private(i) shared(v1, vg) firstprivate(i_max) 
      for(i=1;i<=i_max;i++){
        gosa2+=v1[i]*conj(vg[i]); // new r and old p should be orthogonal
      }
      gosa2 = SumMPI_dc(gosa2);
            
      beta=rnorm2/rnorm;
#pragma omp parallel for default(none) shared(v1, vg) firstprivate(i_max, beta)
      for(i=1;i<=i_max;i++){
        vg[i]=v1[i]+beta*vg[i];
      }
      // if(itr%5==0){
      childfopenMPI(sdt_1,"a", &fp_0);
      fprintf(fp_0,"i_itr=%d itr=%d %.10lf %.10lf \n ",
              i_itr,itr,sqrt(rnorm2),pow(10,-5)*sqrt(bnorm));
      fclose(fp_0);                
      if(sqrt(rnorm2)<eps*sqrt(bnorm)){
        t_itr+=itr;
        childfopenMPI(sdt_1,"a", &fp_0);
        fprintf(fp_0,"CG OK:   t_itr=%d \n ",t_itr);
        fclose(fp_0); 
        break;
      }
      //}
    }
    //CG finish!!
    xnorm=0.0;
#pragma omp parallel for reduction(+:xnorm) default(none) private(i) shared(v0) firstprivate(i_max)
    for(i=1;i<=i_max;i++){
      xnorm+=conj(v0[i])*v0[i];
    }
    xnorm = SumMPI_d(xnorm);
    xnorm=sqrt(xnorm);

#pragma omp parallel for default(none) shared(v0) firstprivate(i_max, xnorm) 
    for(i=1;i<=i_max;i++){
      v0[i]=v0[i]/xnorm;
    }
    xb=0.0;

#pragma omp parallel for default(none) reduction(+:xb) private(i) shared(v0, b) firstprivate(i_max)
    for(i=1;i<=i_max;i++){
      xb+=conj(v0[i])*b[i];
    }
    xb = SumMPI_dc(xb);
       
    mid=time(NULL);
       
    childfopenMPI(sdt_1,"a", &fp_0);
    fprintf(fp_0,"i_itr=%d itr=%d time=%lf  fabs(fabs(xb)-1.0)=%.16lf\n"
            ,i_itr,itr,difftime(mid,start),fabs(cabs(xb)-1.0));
    fclose(fp_0);
        
    if(fabs(cabs(xb)-1.0)<eps){
      childfopenMPI(sdt_1,"a", &fp_0);
      fprintf(fp_0,"number of iterations in inv1:i_itr=%d itr=%d t_itr=%d %lf\n ",
              i_itr,itr,t_itr,fabs(cabs(xb)-1.0));
      fclose(fp_0);
      
      break;
    }else{
#pragma omp parallel for default(none) private(i) shared(b, v0) firstprivate(i_max)
      for(i=1;i<=i_max;i++){
        b[i]=v0[i];
      }
    }       
    //inv1 routine finish!!
  }
    
  free(b);
  free(y);

  TimeKeeper(X, cFileNameTimeKeep, cCG_EigenVecFinish, "a");
  fprintf(stdoutMPI, "%s", cLogCG_EigenVecEnd);
  
  return 0;
}/*int CG_EigenVector*/
