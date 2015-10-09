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

/** 
 * 
 * 
 * @param alpha 
 * @param beta 
 * @param ndim 
 * @param E 
 * @param X 
 * @author Takahiro Misawa (The University of Tokyo)
 */
void vec12(double alpha[],double beta[],int ndim,
	   double E[],struct BindStruct *X){
  
  int j,k,l,km,nvec;
  int i_deg,j_deg;
  double tmp_deg;
  double *di,*bl,*bu,*bv;
  double *cm,*lex;
  double s,dnorm;
    
  nvec = X->Def.nvec;
  di=(double *)malloc((X->Def.Lanczos_max+1)*sizeof(double));
  bl=(double *)malloc((X->Def.Lanczos_max+1)*sizeof(double));
  bu=(double *)malloc((X->Def.Lanczos_max+1)*sizeof(double));
  bv=(double *)malloc((X->Def.Lanczos_max+1)*sizeof(double));
  cm=(double *)malloc((X->Def.Lanczos_max+1)*sizeof(double));
  lex=(double *)malloc((X->Def.Lanczos_max+1)*sizeof(double));

  for(k=1;k<=nvec;k++){
#pragma omp parallel for default(none) firstprivate(ndim, k, E, alpha, beta) private(j) shared(di, bl, bu)
    for(j=1;j<=ndim;j++){
      di[j]=E[k]-alpha[j];
      bl[j]=-beta[j];
      bu[j]=-beta[j];
    }

    //LU decomposition
    for(j=1;j<=ndim-1;j++){         
      if(fabs(di[j]>fabs(bl[j]))){
	//non pivoting    
	lex[j]=0;
	if(fabs(di[j]<eps_vec12)){
	  di[j]=eps_vec12;
	}   
	cm[j+1]=bl[j]/di[j];
	di[j+1]=di[j+1]-cm[j+1]*bu[j];
	bv[j]=0.0;
      }else{
	//pivoting
	lex[j]=1.0;
	cm[j+1]=di[j]/bl[j];
	di[j]=bl[j];
	s=bu[j];
	bu[j]=di[j+1];
	bv[j]=bu[j+1];
	di[j+1]=s-cm[j+1]*bu[j];
	bu[j+1]=-cm[j+1]*bv[j];
      }
    }
    
    if(fabs(di[ndim])<eps_vec12){
      di[ndim]=eps_vec12;
    }

    //initial vector
#pragma omp parallel for default(none) firstprivate(k, ndim) private(j) shared(vec)
    for(j=1;j<=ndim;j++){
      vec[k][j]=1.0/(1.0*j*5.0);
    }
    
    //degeneracy check up
    if(k==1){
      km=k;   
    }else if(fabs(E[k]-E[km]) > eps_Energy ){
      km=k;
    }else{
      printf("DEGENERACY IN vec12.c \n");
      for(i_deg=km;i_deg<=k-1;i_deg++){
	tmp_deg = 0.0;
		
#pragma omp parallel for default(none) reduction(+:tmp_deg) firstprivate(ndim, k, i_deg) private(j_deg) shared(vec)
	for(j_deg=1;j_deg<=ndim;j_deg++){
	  tmp_deg += vec[i_deg][j_deg]*vec[k][i_deg]; 
	}
#pragma omp parallel for default(none) firstprivate(ndim, k, i_deg, tmp_deg) private(j_deg) shared(vec)
	for(j_deg=1;j_deg<=ndim;j_deg++){
	  vec[k][j_deg] += -tmp_deg*vec[i_deg][j_deg];
	}	
      }			
    }		
    //inverse iteration
    for(l=1;l<=k-km+3;l++){
      if(l!=1 || k!=km){
	// forward substition
	for(j=1;j<=ndim-1;j++){
	  if(lex[j]==0){
	    vec[k][j+1]+=-vec[k][j]*cm[j+1];
	  }else{
	    s=vec[k][j];
	    vec[k][j]=vec[k][j+1];
	    vec[k][j+1]=s-vec[k][j]*cm[j+1];
	  }
	}
      }
        
      //backward substition
      for(j=ndim;j>=1;j--){
	s=vec[k][j];
	if(j<ndim-1){
	  s+=-vec[k][j+1]*bu[j];
	}
	if(j<ndim-2){
	  s+=-vec[k][j+2]*bv[j];
	}
	vec[k][j]=s/di[j];
      }
        
      //normalization
      dnorm=0.0;
#pragma omp parallel for default(none) reduction(+:dnorm)  firstprivate(ndim, k) private(j) shared(vec)
      for(j=1;j<=ndim;j++){
	dnorm+=pow(vec[k][j],2);
      }

      if(dnorm> eps_vec12){
	dnorm=1.0/sqrt(dnorm);
      }

#pragma omp parallel for default(none) firstprivate(ndim, k, dnorm) private(j) shared(vec)
      for(j=1;j<=ndim;j++){
	vec[k][j]=vec[k][j]*dnorm;
      }
    }
  }

  free(di);
  free(bl);
  free(bu);
  free(bv);
  free(cm);
  free(lex);

}   
