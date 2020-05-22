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
#include "PowerLanczos.h"
#include "mltply.h"
#include "wrapperMPI.h"

int PowerLanczos(struct BindStruct *X){
  
  long int i,j;

  double dnorm, dnorm_inv;
  double complex dam_pr1,dam_pr2a,dam_pr2b,dam_pr3,dam_pr4;
  double E1,E2a,E2b,E3,E4;
  double alpha_p,alpha_m;
  double Lz_Var;
  double Lz_Ene_p,Lz_Ene_m;
  double Lz_Var_p,Lz_Var_m;
  long int i_max,i_Lz;
    
  i_max=X->Check.idim_max;    

// v1 is eigenvector
// v0 = H*v1
// this subroutine 
  for(i_Lz=0;i_Lz<50;i_Lz++){
    //v1 -> eigen_vec
    //v0 -> v0=H*v1
//  if(i_Lz>1){
    #pragma omp parallel for default(none) private(i) shared(v0) firstprivate(i_max)
    for(i = 1; i <= i_max; i++){
      v0[i]=0.0+0.0*I;
    }
    mltply(X, v0, v1); // v0+=H*v1

    dam_pr1=0.0;
    dam_pr2a=0.0;
    #pragma omp parallel for default(none) reduction(+:dam_pr1,dam_pr2a) private(j) shared(v0, v1)firstprivate(i_max) 
    for(j=1;j<=i_max;j++){
      dam_pr1   += conj(v1[j])*v0[j]; // E   = <v1|H|v1>=<v1|v0>
      dam_pr2a  += conj(v0[j])*v0[j]; // E^2 = <v1|H*H|v1>=<v0|v0>
     //v0[j]=v1[j]; v1-> orginal v0=H*v1
    }  
    dam_pr1 = SumMPI_dc(dam_pr1);
    dam_pr2a = SumMPI_dc(dam_pr2a);
    E1    = creal(dam_pr1); // E
    E2a   = creal(dam_pr2a);// E^2

    #pragma omp parallel for default(none) private(i) shared(vg) firstprivate(i_max)
    for(i = 1; i <= i_max; i++){
      vg[i]=0.0;
    }
  
    mltply(X, vg, v0); // vg=H*v0=H*H*v1
    dam_pr2b = 0.0;
    dam_pr3  = 0.0;
    dam_pr4  = 0.0;
    #pragma omp parallel for default(none) reduction(+:dam_pr2b,dam_pr3, dam_pr4) private(j) shared(v0,v1, vg)firstprivate(i_max) 
    for(j=1;j<=i_max;j++){
      dam_pr2b   += conj(vg[j])*v1[j]; // E^2   = <v1|H^2|v1>=<vg|v1>
      dam_pr3    += conj(vg[j])*v0[j]; // E^3   = <v1|H^3|v1>=<vg|v0>
      dam_pr4    += conj(vg[j])*vg[j]; // E^4   = <v1|H^4|v1>=<vg|vg>
    } 
    dam_pr2b = SumMPI_dc(dam_pr2b);
    dam_pr3 = SumMPI_dc(dam_pr3);
    dam_pr4 = SumMPI_dc(dam_pr4);
    //E1    = X->Phys.energy;// E^1
    //E2a   = X->Phys.var ;// E^2 = <v1|H*H|v1>
    E2b   = creal(dam_pr2b) ;// E^2 = (<v1|H^2)|v1>
    E3    = creal(dam_pr3) ;// E^3
    E4    = creal(dam_pr4) ;// E^4
  
    if(solve_2ndPolinomial(X,&alpha_p,&alpha_m,E1,E2a,E2b,E3,E4)!=TRUE){
      fprintf(stdoutMPI,"Power Lanczos break \n");
      return 0;
    }
    //printf("E1=%.16lf E2a=%.16lf E2b=%.16lf E3=%.16lf E4=%.16lf \n",E1,E2a,E2b,E3,E4);
    
    Lz(X,alpha_p,&Lz_Ene_p,&Lz_Var_p,E1,E2a,E3,E4);
    Lz(X,alpha_m,&Lz_Ene_m,&Lz_Var_m,E1,E2a,E3,E4);
  
    if(Lz_Ene_p < Lz_Ene_m){
      Lz_Var=Lz_Var_p;
      fprintf(stdoutMPI,"Power Lanczos (P): %.16lf %.16lf \n",Lz_Ene_p,Lz_Var_p);
      #pragma omp parallel for default(none)  private(j) shared(v0, v1) firstprivate(i_max,alpha_p) 
      for(j=1;j<=i_max;j++){
        v1[j]   = v1[j]+alpha_p*v0[j];   // (1+alpha*H)v1=v1+alpha*v0
      }
    }else{
      Lz_Var=Lz_Var_m;
      fprintf(stdoutMPI,"Power Lanczos (M): %.16lf %.16lf \n",Lz_Ene_m,Lz_Var_m);
      #pragma omp parallel for default(none)  private(j) shared(v0, v1)firstprivate(i_max,alpha_m) 
      for(j=1;j<=i_max;j++){
        v1[j]   = v1[j]+alpha_m*v0[j]; // (1+alpha*H)v1=v1+alpha*v0
      }
    } 
    //normalization
    dnorm=0.0;
    #pragma omp parallel for default(none) reduction(+:dnorm) private(j) shared(v1) firstprivate(i_max)
    for(j=1;j<=i_max;j++){
      dnorm += conj(v1[j])*v1[j];
    }
    dnorm = SumMPI_d(dnorm);
    dnorm=sqrt(dnorm);
    dnorm_inv=1.0/dnorm;
#pragma omp parallel for default(none) private(j) shared(v1) firstprivate(i_max, dnorm_inv)
    for(j=1;j<=i_max;j++){
      v1[j] = v1[j]*dnorm_inv;
    }
    if(Lz_Var < eps_Energy){
      fprintf(stdoutMPI,"Power Lanczos break \n");
      return 1;
      //break;
    }
  }
  return 0;
}

int  solve_2ndPolinomial(struct BindStruct *X,double *alpha_p,double *alpha_m,double E1,double E2a,double E2b,double E3,double E4){
  double a,b,c,d;
  double tmp_AA,tmp_BB,tmp_CC;
 
  //not solving 2nd Polinomial
  //approximate linear equation is solved


  a = E1;
  b = E2a;
  c = E2a;
  d = E3;

  tmp_AA  = b*(b+c)-2*a*d;
  tmp_BB  = -a*b+d;
  tmp_CC  = b*((b+c)*(b+c))-(a*a)*b*(b+2*c)+4*(a*a*a)*d-2*a*(2*b+c)*d+d*d;
  if(tmp_AA< pow(b, 2)*pow(10.0, -15)){
    return FALSE;
  }
  //printf("XXX: %.16lf %.16lf %.16lf %.16lf \n",a, b, c, d);
  //printf("XXX: %.16lf %.16lf %.16lf \n",tmp_AA,tmp_BB,tmp_CC);
  if(tmp_CC>=0){
  *alpha_p =  (tmp_BB+sqrt((tmp_CC)))/tmp_AA; 
  *alpha_m =  (tmp_BB-sqrt((tmp_CC)))/tmp_AA; 
  //printf("YYY: %.16lf %.16lf  \n",*alpha_p,*alpha_m);
  }
  else{
  //*alpha_m =  -E2a/E3*(1-E1*E1/E2a)/(1-E1*E2a/E3);
    //*alpha_p =  0.0;
  //*alpha_m =  -E3/E4*(1+2*E1*E1*E1/E3-3*E1*E2a/E3)/(1-5*E2a*E2a/E4+4*E1*E1*E2a/E4);
    *alpha_p =  cabs((tmp_BB+csqrt((tmp_CC))))/tmp_AA; 
    *alpha_m =  cabs((tmp_BB-csqrt((tmp_CC))))/tmp_AA; 
  }
  return TRUE;
}

void  Lz(struct BindStruct *X,double alpha,double *Lz_Ene,double *Lz_Var,double E1,double E2,double E3,double E4){

  double tmp_ene,tmp_var;

  tmp_ene        = (E1+2*alpha*E2+alpha*alpha*E3)/(1+2*alpha*E1+alpha*alpha*E2);
  *Lz_Ene        = tmp_ene;
  X->Phys.energy = tmp_ene;
  tmp_var        = (E2+2*alpha*E3+alpha*alpha*E4)/(1+2*alpha*E1+alpha*alpha*E2);
  X->Phys.var    = tmp_var;
  *Lz_Var        = fabs(tmp_var-tmp_ene*tmp_ene)/tmp_var;

}

