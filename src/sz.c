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

#include <bitcalc.h>
#include "mfmemory.h"
#include "FileIO.h"
#include "sz.h"
#include "wrapperMPI.h"

/**
 * @file   sz.c
 * 
 * @brief  File List 
 * 
 * @version 0.2
 * @details 
 *
 * @version 0.1
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 */

/** 
 * 
 * 
 * @param X 
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int sz
(
 struct BindStruct *X
 )
{
  FILE *fp,*fp_err;
  char sdt[D_FileNameMax],sdt_err[D_FileNameMax];
    
  long unsigned int i,icnt; 
  long unsigned int ib,jb;
    
  long unsigned int j;
  long unsigned int div;
  long unsigned int num_up,num_down;
  long unsigned int irght,ilft,ihfbit;

  //*[s] for omp parall
  int  all_up,all_down,tmp_res,num_threads;
  long unsigned int tmp_1,tmp_2,tmp_3;
  int mfint[7];
  long int **comb;
  //*[e] for omp parall

  // [s] for Kondo
  int N_all_up, N_all_down;
  int all_loc;
  long unsigned int num_loc, div_down;
  int num_loc_up;
  int icheck_loc;
  int ihfSpinDown=0;
  // [e] for Kondo
    
  long unsigned int i_max;
  double idim=0.0;
  long unsigned int div_up;

  int iSpnup, iMinup,iAllup;

  int N2=0;
  int N=0;
  fprintf(stdoutMPI, "%s", cProStartCalcSz);
  TimeKeeper(X, cFileNameSzTimeKeep, cInitalSz, "w");

  if(X->Check.idim_max!=0){
 
  switch(X->Def.iCalcModel){
  case HubbardGC:
  case HubbardNConserved:
  case Hubbard:
    N2=2*X->Def.Nsite;
    idim = pow(2.0,N2);
    break;
  case KondoGC:
  case Kondo:
    N2  = 2*X->Def.Nsite;
    N  =  X->Def.Nsite;
    idim = pow(2.0,N2);
    for(j=0;j<N;j++){
      fprintf(stdoutMPI, cStateLocSpin,j,X->Def.LocSpn[j]);
    }
    break;
  case SpinGC:
  case Spin:
    N=X->Def.Nsite;
    if(X->Def.iFlgGeneralSpin==FALSE){
      idim = pow(2.0, N);
    }
    else{
      idim=1;
      for(j=0; j<N; j++){
	idim *= X->Def.SiteToBit[j];
      }
    }
    break;
  }
  li_malloc2(comb, X->Def.Nsite+1,X->Def.Nsite+1);
  i_max=X->Check.idim_max;
  
  switch(X->Def.iCalcModel){
  case HubbardNConserved:
  case Hubbard:
  case KondoGC:
  case Kondo:
  case Spin:
    if(X->Def.iFlgGeneralSpin==FALSE){
      if(GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit)!=0){
	exitMPI(-1);
      }
      //fprintf(stdoutMPI, "idim=%lf irght=%ld ilft=%ld ihfbit=%ld \n",idim,irght,ilft,ihfbit);
    }
    else{
      ihfbit=X->Check.sdim;
      //fprintf(stdoutMPI, "idim=%lf ihfbit=%ld \n",idim, ihfbit);
    }
  break;
 default:
   break;
}   
  
  icnt=1;
  jb=0;

  if(X->Def.READ==1){
    if(!Read_sz(X, irght, ilft, ihfbit, &i_max)==0){
      exitMPI(-1);
    }
  }
  else{ 
    sprintf(sdt, cFileNameSzTimeKeep, X->Def.CDataFileHead);
    #ifdef _OPENMP
    num_threads  = omp_get_max_threads();
    #else
    num_threads=1;
    #endif
    childfopenMPI(sdt,"a", &fp);
    fprintf(fp, "num_threads==%d\n",num_threads);
    fclose(fp);
    
    //*[s] omp parallel

    TimeKeeper(X, cFileNameSzTimeKeep, cOMPSzStart, "a");
    
    switch(X->Def.iCalcModel){
    case HubbardGC:
      icnt = X->Def.Tpow[2*X->Def.Nsite-1]*2+0;/*Tpow[2*X->Def.Nsit]=1*/
      break;
      
    case SpinGC:
      if(X->Def.iFlgGeneralSpin==FALSE){
	icnt = X->Def.Tpow[X->Def.Nsite-1]*2+0;/*Tpow[X->Def.Nsit]=1*/
      }
      else{
	icnt = X->Def.Tpow[X->Def.Nsite-1]*X->Def.SiteToBit[X->Def.Nsite-1];
      }
      break;
      
    case KondoGC:
   // this part can not be parallelized
      jb = 0;
      num_loc=0;
      for(j=X->Def.Nsite/2; j< X->Def.Nsite ;j++){
	if(X->Def.LocSpn[j] != ITINERANT){
	  num_loc += 1;
	}
      }

      for(ib=0;ib<X->Check.sdim;ib++){
	list_jb[ib]=jb;
	i=ib*ihfbit;
	icheck_loc=1;
	
	for(j=(X->Def.Nsite+1)/2; j< X->Def.Nsite ;j++){
	  div_up    = i & X->Def.Tpow[2*j];
	  div_up    = div_up/X->Def.Tpow[2*j];
	  div_down  = i & X->Def.Tpow[2*j+1];
	  div_down  = div_down/X->Def.Tpow[2*j+1];

	  if(X->Def.LocSpn[j] != ITINERANT){
	     if(X->Def.Nsite%2==1 && j==(X->Def.Nsite/2)){
	      icheck_loc= icheck_loc;
	    }
	     else{
	       icheck_loc   = icheck_loc*(div_up^div_down);// exclude doubllly ocupited site
	     }
	  }
	}
	if(icheck_loc == 1){
	  if(X->Def.Nsite%2==1 && X->Def.LocSpn[X->Def.Nsite/2] != ITINERANT){
	    jb +=X->Def.Tpow[X->Def.Nsite-1-(X->Def.NLocSpn-num_loc)];
	  }else{
	    jb +=X->Def.Tpow[X->Def.Nsite-(X->Def.NLocSpn-num_loc)];
	  }
	}
      }

      icnt = 0; 
#pragma omp parallel for default(none) reduction(+:icnt) private(ib) firstprivate(ihfbit, N2, X) shared(list_1)
      for(ib=0;ib<X->Check.sdim;ib++){
	icnt+=child_omp_sz_KondoGC(ib,ihfbit,N2,X);
      }      
    break;
      
 case Hubbard:
      // this part can not be parallelized
      jb = 0;

      for(ib=0;ib<X->Check.sdim;ib++){
	list_jb[ib]=jb;

	i=ib*ihfbit;
	num_up=0;
	for(j=0;j<=N2-2;j+=2){
	  div=i & X->Def.Tpow[j];
	  div=div/X->Def.Tpow[j];
	  num_up+=div;
	}
	num_down=0;
	for(j=1;j<=N2-1;j+=2){
	  div=i & X->Def.Tpow[j];
	  div=div/X->Def.Tpow[j];
	  num_down+=div;
	}
	
	tmp_res  = X->Def.Nsite%2; // even Ns-> 0, odd Ns -> 1
	all_up   = (X->Def.Nsite+tmp_res)/2;
	all_down = (X->Def.Nsite-tmp_res)/2;

	tmp_1 = Binomial(all_up,X->Def.Nup-num_up,comb,all_up);
	tmp_2 = Binomial(all_down,X->Def.Ndown-num_down,comb,all_down);
	jb   += tmp_1*tmp_2;
      }

      //#pragma omp barrier
      TimeKeeper(X, cFileNameSzTimeKeep, cOMPSzFinish, "a");
 
      icnt = 0;
      #pragma omp parallel for default(none) reduction(+:icnt) private(ib) firstprivate(ihfbit, N2, X) 
      for(ib=0;ib<X->Check.sdim;ib++){
	icnt+=child_omp_sz(ib,ihfbit,N2,X);
      }
      break;

    case HubbardNConserved:
      // this part can not be parallelized
      jb = 0;
      iSpnup=0;
      iMinup=0;
      iAllup=X->Def.Ne;
      if(X->Def.Ne > X->Def.Nsite){
	iMinup = X->Def.Ne-X->Def.Nsite;
	iAllup = X->Def.Nsite;
      }
      
      for(ib=0;ib<X->Check.sdim;ib++){
	list_jb[ib]=jb;

	i=ib*ihfbit;
	num_up=0;
	for(j=0;j<=N2-2;j+=2){
	  div=i & X->Def.Tpow[j];
	  div=div/X->Def.Tpow[j];
	  num_up+=div;
	}
	num_down=0;
	for(j=1;j<=N2-1;j+=2){
	  div=i & X->Def.Tpow[j];
	  div=div/X->Def.Tpow[j];
	  num_down+=div;
	}
	
	tmp_res  = X->Def.Nsite%2; // even Ns-> 0, odd Ns -> 1
	all_up   = (X->Def.Nsite+tmp_res)/2;
	all_down = (X->Def.Nsite-tmp_res)/2;

	for(iSpnup=iMinup; iSpnup<= iAllup; iSpnup++){
	  tmp_1 = Binomial(all_up, iSpnup-num_up,comb,all_up);
	  tmp_2 = Binomial(all_down, X->Def.Ne-iSpnup-num_down,comb,all_down);
	  jb   += tmp_1*tmp_2;
	}
      }
      //#pragma omp barrier
      TimeKeeper(X, cFileNameSzTimeKeep, cOMPSzFinish, "a");
	
      icnt = 0;
#pragma omp parallel for default(none) reduction(+:icnt) private(ib) firstprivate(ihfbit, N2, X) 
      for(ib=0;ib<X->Check.sdim;ib++){
	icnt+=child_omp_sz(ib,ihfbit,N2,X);
      }
      break;
            
    case Kondo:
      // this part can not be parallelized
      N_all_up   = X->Def.Nup;
      N_all_down = X->Def.Ndown;
      fprintf(stdoutMPI, cStateNupNdown, N_all_up,N_all_down);

      jb = 0;
      num_loc=0;
      for(j=X->Def.Nsite/2; j< X->Def.Nsite ;j++){
	if(X->Def.LocSpn[j] != ITINERANT){
	  num_loc += 1;
	}
      }

      for(ib=0;ib<X->Check.sdim;ib++){
	list_jb[ib]=jb;
	i=ib*ihfbit;
	num_up=0;
	num_down=0;	
	icheck_loc=1;

	for(j=X->Def.Nsite/2; j< X->Def.Nsite ;j++){
	  div_up    = i & X->Def.Tpow[2*j];
	  div_up    = div_up/X->Def.Tpow[2*j];
	  div_down  = i & X->Def.Tpow[2*j+1];
	  div_down  = div_down/X->Def.Tpow[2*j+1];
	  if(X->Def.LocSpn[j] == ITINERANT){
	    num_up   += div_up;        
	    num_down += div_down;  
	  }else{    
	    num_up   += div_up;     
	    num_down += div_down;
	    if(X->Def.Nsite%2==1 && j==(X->Def.Nsite/2)){
	      icheck_loc= icheck_loc;
	      ihfSpinDown=div_down;
	      if(div_down ==0){
		num_up += 1;
	      }
	    }
	    else{
	      icheck_loc   = icheck_loc*(div_up^div_down);// exclude doubllly ocupited site
	    }
	  }
	}

	if(icheck_loc == 1){
	  tmp_res  = X->Def.Nsite%2; // even Ns-> 0, odd Ns -> 1
	  all_loc =  X->Def.NLocSpn-num_loc;
	  all_up   = (X->Def.Nsite+tmp_res)/2-all_loc;
	  all_down = (X->Def.Nsite-tmp_res)/2-all_loc;
	  
	  if(X->Def.Nsite%2==1 && X->Def.LocSpn[X->Def.Nsite/2] != ITINERANT){
	    all_up   = (X->Def.Nsite)/2-all_loc;
	    all_down = (X->Def.Nsite)/2-all_loc;
	  }
	  
	  for(num_loc_up=0; num_loc_up <= all_loc; num_loc_up++){
	    tmp_1 = Binomial(all_loc, num_loc_up, comb, all_loc);
	    
	    if( X->Def.Nsite%2==1 && X->Def.LocSpn[X->Def.Nsite/2] != ITINERANT){
	      if(ihfSpinDown !=0){
		tmp_2 = Binomial(all_up, X->Def.Nup-num_up-num_loc_up, comb, all_up);
		tmp_3 = Binomial(all_down, X->Def.Ndown-num_down-(all_loc-num_loc_up), comb, all_down);
		
	      }
	      else{
		tmp_2 = Binomial(all_up, X->Def.Nup-num_up-num_loc_up, comb, all_up);
		tmp_3 = Binomial(all_down, X->Def.Ndown-num_down-(all_loc-num_loc_up), comb, all_down);
	      }
	    }
	    else{
	      tmp_2 = Binomial(all_up, X->Def.Nup-num_up-num_loc_up, comb, all_up);
	      tmp_3 = Binomial(all_down, X->Def.Ndown-num_down-(all_loc-num_loc_up), comb, all_down);
	    }
	    jb   += tmp_1*tmp_2*tmp_3;  
	  }
	}

      }
      //#pragma omp barrier
      
      TimeKeeper(X, cFileNameSzTimeKeep, cOMPSzMid, "a");
 
      icnt = 0;
#pragma omp parallel for default(none) reduction(+:icnt) private(ib) firstprivate(ihfbit, N2, X) 
      for(ib=0;ib<X->Check.sdim;ib++){
	icnt+=child_omp_sz_Kondo(ib,ihfbit,N2,X);
      }
      break;

    case Spin:
      // this part can not be parallelized
      if(X->Def.iFlgGeneralSpin==FALSE){
	jb = 0;
	//fprintf(stdoutMPI, "Check.sdim=%ld, ihfbit=%ld\n", X->Check.sdim, ihfbit);
	for(ib=0;ib<X->Check.sdim;ib++){
	  list_jb[ib]=jb;
	  i=ib*ihfbit;
	  num_up=0;
	  for(j=0;j<N; j++){
	    div_up = i & X->Def.Tpow[j];
	    div_up = div_up/X->Def.Tpow[j];
	    num_up +=div_up;
	  }
	  all_up   = (X->Def.Nsite+1)/2;
	  tmp_1 = Binomial(all_up,X->Def.Ne-num_up,comb,all_up);
	  jb   += tmp_1;
	}
	//#pragma omp barrier

	TimeKeeper(X, cFileNameSzTimeKeep, cOMPSzMid, "a");
 
	icnt = 0;
#pragma omp parallel for default(none) reduction(+:icnt) private(ib) firstprivate(ihfbit, N, X)
	for(ib=0;ib<X->Check.sdim;ib++){
	  icnt+=child_omp_sz_spin(ib,ihfbit,N,X);
	}
      }
      else{
	int Max2Sz=0;
	int irghtsite=1;
	long int itmpSize=1;
	int i2Sz=0;
	for(j=0; j<X->Def.Nsite; j++){
	  itmpSize *= X->Def.SiteToBit[j];
	  if(itmpSize==ihfbit){
	    break;
	  }
	  irghtsite++;
	}
        for(j=0; j<X->Def.Nsite; j++){
          Max2Sz += X->Def.LocSpn[j];
        }
	
	lui_malloc1(HilbertNumToSz, 2*Max2Sz+1);
	for(ib=0; ib<2*Max2Sz+1; ib++){
	  HilbertNumToSz[ib]=0;
	}
	
	for(ib =0; ib<ihfbit; ib++){
	  i2Sz=0;
	  for(j=1; j<= irghtsite; j++){
	    i2Sz += GetLocal2Sz(j,ib, X->Def.SiteToBit, X->Def.Tpow);
	  }
	  list_2_1_Sz[ib]=i2Sz;
	  HilbertNumToSz[i2Sz+Max2Sz]++;
	}
        jb = 0;
	long int ilftdim=(X->Def.Tpow[X->Def.Nsite-1]*X->Def.SiteToBit[X->Def.Nsite-1])/ihfbit;
	for(ib=0;ib<ilftdim;ib++){
	  list_jb[ib]=jb;
	  i2Sz=0;
	  for(j=1;j<=(N-irghtsite); j++){
	    i2Sz += GetLocal2Sz(j,ib, X->Def.SiteToBit, X->Def.Tpow);
	  }
	  list_2_2_Sz[ib]=i2Sz;
	  if(X->Def.Total2Sz- i2Sz +Max2Sz>=0 && X->Def.Total2Sz- i2Sz <= Max2Sz){
	    jb += HilbertNumToSz[X->Def.Total2Sz- i2Sz +Max2Sz];
	  }
	}
	
	TimeKeeper(X, cFileNameSzTimeKeep, cOMPSzMid, "a");

	icnt = 0;
#pragma omp parallel for default(none) reduction(+:icnt) private(ib) firstprivate(ilftdim, ihfbit, N, X) 
	for(ib=0;ib<ilftdim; ib++){
	  icnt+=child_omp_sz_GeneralSpin(ib,ihfbit,N,X);
	}
		
	i_free1(HilbertNumToSz, 2*Max2Sz+1);	
      }
      
      break;
    default:
      exitMPI(-1);
       
    }    
    i_max=icnt;
    //fprintf(stdoutMPI, "Xicnt=%ld \n",icnt);
    TimeKeeper(X, cFileNameSzTimeKeep, cOMPSzFinish, "a");
  }

  if(X->Def.iCalcModel==HubbardNConserved){
    X->Def.iCalcModel=Hubbard;
  }
  
  //Error message
  //i_max=i_max+1;
  if(i_max!=X->Check.idim_max){
    fprintf(stderr, "%s", cErrSz);
    fprintf(stderr, cErrSz_ShowDim, i_max, X->Check.idim_max);
    strcpy(sdt_err,cFileNameErrorSz);
    if(childfopenMPI(sdt_err,"a",&fp_err)!=0){
      exitMPI(-1);
    }
    fprintf(fp_err,cErrSz_OutFile);
    fclose(fp_err);
    exitMPI(-1);
  }
  
  i_free2(comb, X->Def.Nsite+1,X->Def.Nsite+1);
  }
  fprintf(stdoutMPI, "%s", cProEndCalcSz);
  return 0;    
}

/** 
 * 
 * 
 * @param n 
 * @param k 
 * @param comb 
 * @param Nsite 
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
long int Binomial(int n,int k,long int **comb,int Nsite){
  // nCk, Nsite=max(n)
  int tmp_i,tmp_j;

  if(n==0 && k==0){
    return 1;
  } 
  else if(n<0 || k<0 || n<k){
    return 0;
  }
  
  for(tmp_i=0;tmp_i<=Nsite;tmp_i++){
    for(tmp_j=0;tmp_j<=Nsite;tmp_j++){
      comb[tmp_i][tmp_j] = 0;
    }
  }

  comb[0][0] = 1;
  comb[1][0] = 1;
  comb[1][1] = 1;
  for(tmp_i=2;tmp_i<=n;tmp_i++){
    for(tmp_j=0;tmp_j<=tmp_i;tmp_j++){
      if(tmp_j==0){
        comb[tmp_i][tmp_j] = 1;
      }else if(tmp_j==tmp_i){
        comb[tmp_i][tmp_j] = 1;
      }else{
        comb[tmp_i][tmp_j] = comb[tmp_i-1][tmp_j-1]+comb[tmp_i-1][tmp_j];
      }
    }
  }
  return comb[n][k];
}

/** 
 * 
 * 
 * @param ib 
 * @param ihfbit 
 * @param N2 
 * @param X 
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int child_omp_sz(long unsigned int ib, long unsigned int ihfbit,int N2,struct BindStruct *X){

  long unsigned int i,j; 
  long unsigned int ia,ja,jb;
  long unsigned int div_down, div_up;
  long unsigned int num_up,num_down;
  long unsigned int tmp_num_up,tmp_num_down;
    
  jb = list_jb[ib];
  i  = ib*ihfbit;
    
  num_up   = 0;
  num_down = 0;
  for(j=0;j< X->Def.Nsite ;j++){
    div_up    = i & X->Def.Tpow[2*j];
    div_up    = div_up/X->Def.Tpow[2*j];
    div_down  = i & X->Def.Tpow[2*j+1];
    div_down  = div_down/X->Def.Tpow[2*j+1];
    num_up += div_up;
    num_down += div_down;
  }
  
  ja=1;
  tmp_num_up   = num_up;
  tmp_num_down = num_down;

  if(X->Def.iCalcModel==Hubbard){
    for(ia=0;ia<X->Check.sdim;ia++){
      i=ia;
      num_up =  tmp_num_up;
      num_down =  tmp_num_down;
      
      for(j=0;j<X->Def.Nsite;j++){
	div_up    = i & X->Def.Tpow[2*j];
	div_up    = div_up/X->Def.Tpow[2*j];
	div_down  = i & X->Def.Tpow[2*j+1];
	div_down  = div_down/X->Def.Tpow[2*j+1];
	num_up += div_up;
	num_down += div_down;
      }
      if(num_up == X->Def.Nup && num_down == X->Def.Ndown){
	list_1[ja+jb]=ia+ib*ihfbit;
      list_2_1[ia]=ja;
      list_2_2[ib]=jb;
      ja+=1;
      } 
    }
  }
  else if(X->Def.iCalcModel==HubbardNConserved){
    for(ia=0;ia<X->Check.sdim;ia++){
      i=ia;
      num_up =  tmp_num_up;
      num_down =  tmp_num_down;
      
      for(j=0;j<X->Def.Nsite;j++){
	div_up    = i & X->Def.Tpow[2*j];
	div_up    = div_up/X->Def.Tpow[2*j];
	div_down  = i & X->Def.Tpow[2*j+1];
	div_down  = div_down/X->Def.Tpow[2*j+1];
	num_up += div_up;
	num_down += div_down;
      }
      if( (num_up+num_down) == X->Def.Ne){
	list_1[ja+jb]=ia+ib*ihfbit;
	list_2_1[ia]=ja;
	list_2_2[ib]=jb;
	ja+=1;
      } 
    }  
  }
  ja=ja-1;    
  return ja; 
}

/** 
 * 
 * 
 * @param ib 
 * @param ihfbit 
 * @param N2 
 * @param X 
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int child_omp_sz_Kondo(long unsigned int ib, long unsigned int ihfbit,int N2,struct BindStruct *X){

  long unsigned int i,j; 
  long unsigned int ia,ja,jb;
  long unsigned int div_down, div_up;
  long unsigned int num_up,num_down;
  long unsigned int tmp_num_up,tmp_num_down;
  int icheck_loc;
    
  jb = list_jb[ib];
  i  = ib*ihfbit;
    
  num_up   = 0;
  num_down = 0;
  icheck_loc=1;
  for(j=X->Def.Nsite/2; j< X->Def.Nsite ;j++){
    div_up    = i & X->Def.Tpow[2*j];
    div_up    = div_up/X->Def.Tpow[2*j];
    div_down  = i & X->Def.Tpow[2*j+1];
    div_down  = div_down/X->Def.Tpow[2*j+1];

    if(X->Def.LocSpn[j] == ITINERANT){
      num_up   += div_up;        
      num_down += div_down;  
    }else{    
      num_up   += div_up;        
      num_down += div_down;
      if(X->Def.Nsite%2==1 && j==(X->Def.Nsite/2)){
	icheck_loc= icheck_loc;
      }
      else{
	icheck_loc   = icheck_loc*(div_up^div_down);// exclude doubllly ocupited site
      }
    }
  }
  
  ja=1;
  tmp_num_up   = num_up;
  tmp_num_down = num_down;
  if(icheck_loc ==1){
    for(ia=0;ia<X->Check.sdim;ia++){
      i=ia;
      num_up =  tmp_num_up;
      num_down =  tmp_num_down;
      icheck_loc=1;
      for(j=0;j<(X->Def.Nsite+1)/2;j++){
	div_up    = i & X->Def.Tpow[2*j];
	div_up    = div_up/X->Def.Tpow[2*j];
	div_down  = i & X->Def.Tpow[2*j+1];
	div_down  = div_down/X->Def.Tpow[2*j+1];
	
	if(X->Def.LocSpn[j] ==  ITINERANT){
	  num_up   += div_up;        
	  num_down += div_down;  
	}else{    
	  num_up   += div_up;        
	  num_down += div_down;  
	  if(X->Def.Nsite%2==1 && j==(X->Def.Nsite/2)){
	    icheck_loc= icheck_loc;
	  }
	  else{
	    icheck_loc   = icheck_loc*(div_up^div_down);// exclude doubllly ocupited site
	  }
	}
      }
      
      if(icheck_loc == 1 && X->Def.LocSpn[X->Def.Nsite/2] != ITINERANT && X->Def.Nsite%2==1){
	div_up    = ia & X->Def.Tpow[X->Def.Nsite-1];
	div_up    = div_up/X->Def.Tpow[X->Def.Nsite-1];
	div_down  = (ib*ihfbit) & X->Def.Tpow[X->Def.Nsite];
	div_down  = div_down/X->Def.Tpow[X->Def.Nsite];
	icheck_loc= icheck_loc*(div_up^div_down);
      }
      
      if(num_up == X->Def.Nup && num_down == X->Def.Ndown && icheck_loc==1){
	list_1[ja+jb]=ia+ib*ihfbit;
	list_2_1[ia]=ja;
	list_2_2[ib]=jb;
	ja+=1;
      }
    }
  }
  ja=ja-1;    
  return ja; 
}

/** 
 * 
 * 
 * @param ib 
 * @param ihfbit 
 * @param N2 
 * @param X 
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int child_omp_sz_KondoGC(long unsigned int ib, long unsigned int ihfbit,int N2,struct BindStruct *X){

  long unsigned int i,j; 
  long unsigned int ia,ja,jb;
  long unsigned int div_down, div_up;
  int icheck_loc;
    
  jb = list_jb[ib];
  i  = ib*ihfbit;
  icheck_loc=1;
  for(j=X->Def.Nsite/2; j< X->Def.Nsite ;j++){
    div_up    = i & X->Def.Tpow[2*j];
    div_up    = div_up/X->Def.Tpow[2*j];
    div_down  = i & X->Def.Tpow[2*j+1];
    div_down  = div_down/X->Def.Tpow[2*j+1];
    if(X->Def.LocSpn[j] !=  ITINERANT){
      if(X->Def.Nsite%2==1 && j==(X->Def.Nsite/2)){
	icheck_loc= icheck_loc;
      }
      else{
	icheck_loc   = icheck_loc*(div_up^div_down);// exclude doubllly ocupited site
      }
    }
  }

  ja=1;
  if(icheck_loc ==1){
    for(ia=0;ia<X->Check.sdim;ia++){
      i=ia;
      icheck_loc =1;
      for(j=0;j<(X->Def.Nsite+1)/2;j++){
	div_up    = i & X->Def.Tpow[2*j];
	div_up    = div_up/X->Def.Tpow[2*j];
	div_down  = i & X->Def.Tpow[2*j+1];
	div_down  = div_down/X->Def.Tpow[2*j+1];	
	if(X->Def.LocSpn[j] !=  ITINERANT){
	  if(X->Def.Nsite%2==1 && j==(X->Def.Nsite/2)){
	    icheck_loc= icheck_loc;
	  }
	  else{
	    icheck_loc   = icheck_loc*(div_up^div_down);// exclude doubllly ocupited site
	  }
	}
      }

      if(icheck_loc == 1 && X->Def.LocSpn[X->Def.Nsite/2] != ITINERANT && X->Def.Nsite%2==1){
	div_up    = ia & X->Def.Tpow[X->Def.Nsite-1];
	div_up    = div_up/X->Def.Tpow[X->Def.Nsite-1];
	div_down  = (ib*ihfbit) & X->Def.Tpow[X->Def.Nsite];
	div_down  = div_down/X->Def.Tpow[X->Def.Nsite];
	icheck_loc= icheck_loc*(div_up^div_down);
      }
      
      if(icheck_loc==1){
	list_1[ja+jb]=ia+ib*ihfbit;
	list_2_1[ia]=ja;
	list_2_2[ib]=jb;
	ja+=1;
      }
    }
  }
  ja=ja-1;
    
  return ja; 
}

/** 
 * 
 * 
 * @param ib 
 * @param ihfbit 
 * @param N 
 * @param X 
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int child_omp_sz_spin(
		      long unsigned int ib, 
		      long unsigned int ihfbit,
		      int N,
		      struct BindStruct *X
		      )
{
  long unsigned int i,j,div; 
  long unsigned int ia,ja,jb;
  long unsigned int num_up;
  int tmp_num_up;
  
  jb = list_jb[ib];
  i  = ib*ihfbit;
  num_up=0;
  for(j=0;j<N;j++){
    div=i & X->Def.Tpow[j];
    div=div/X->Def.Tpow[j];
    num_up+=div;
  }
  ja=1;
  tmp_num_up   = num_up;
  
  for(ia=0;ia<ihfbit;ia++){
    i=ia;
    num_up =  tmp_num_up;
    for(j=0;j<N;j++){
      div=i & X->Def.Tpow[j];
      div=div/X->Def.Tpow[j];
      num_up+=div;
    }

    if(num_up == X->Def.Ne){
      list_1[ja+jb]=ia+ib*ihfbit;
      list_2_1[ia]=ja;
      list_2_2[ib]=jb;
      ja+=1;
    } 
  }
  ja=ja-1;
  return ja; 
}

/** 
 * 
 * 
 * @param ib 
 * @param ihfbit 
 * @param N 
 * @param X 
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int child_omp_sz_GeneralSpin(
		      long unsigned int ib, 
		      long unsigned int ihfbit,
		      int N,
		      struct BindStruct *X
		      )
{
  long unsigned int ia,ja,jb;  
  int list_2_2_Sz_ib=0;
  int tmp_2Sz=0;
  jb = list_jb[ib];
  list_2_2_Sz_ib =list_2_2_Sz[ib];
  ja=1;
  for(ia=0;ia<ihfbit;ia++){
    tmp_2Sz=list_2_1_Sz[ia]+list_2_2_Sz_ib;
    if(tmp_2Sz == X->Def.Total2Sz){
      list_1[ja+jb]=ia+ib*ihfbit;
      list_2_1[ia]=ja;
      list_2_2[ib]=jb;
      ja+=1;
    } 
  }
  ja=ja-1;
  return ja; 
}

/** 
 * 
 * 
 * @param X 
 * @param irght 
 * @param ilft 
 * @param ihfbit 
 * @param i_max 
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int Read_sz
(
 struct BindStruct *X,
 const long unsigned int irght,
 const long unsigned int ilft,
 const long unsigned int ihfbit,
 long unsigned int *i_max
 )
{
  FILE *fp,*fp_err;
  char sdt[D_FileNameMax];
  char buf[D_FileNameMax];
    
  long unsigned int icnt=0; 
  long unsigned int ia,ib;
  long unsigned int ja=0;
  long unsigned int jb=0;
  long unsigned int ibpatn=0;
  long unsigned int dam; 

  TimeKeeper(X,cFileNameSzTimeKeep,cReadSzStart, "a");

  switch(X->Def.iCalcModel){
  case Hubbard:
  case HubbardGC:
  case Spin:
  case SpinGC:
    sprintf(sdt,cFileNameListModel, X->Def.Nsite, X->Def.Nup, X->Def.Ndown);
    break;
  case Kondo:
    sprintf(sdt,"ListForKondo_Ns%d_Ncond%d.dat",X->Def.Nsite,X->Def.Ne);
    break;
  }
  if(childfopenMPI(sdt,"r", &fp)!=0){
    exitMPI(-1);
  }  

  if(fp == NULL){
    if(childfopenMPI(cFileNameErrorSz,"a",&fp_err)!=0){
      exitMPI(-1);
    }
    fprintf(fp_err, cErrSz_NoFile);
    fprintf(stderr, cErrSz_NoFile);
    fprintf(fp_err, cErrSz_NoFile_Show,sdt);
    fprintf(stderr, cErrSz_NoFile_Show, sdt);
    fclose(fp_err);
  }else{
    while(NULL != fgetsMPI(buf,sizeof(buf),fp)){  
      dam=atol(buf);  
      list_1[icnt]=dam;
            
      ia= dam & irght;
      ib= dam & ilft;
      ib=ib/ihfbit; 
            
      if(ib==ibpatn){
	ja=ja+1;
      }else{
	ibpatn=ib;
	ja=1;
	jb=icnt-1;
      }
            
      list_2_1[ia]=ja;
      list_2_2[ib]=jb;
      icnt+=1;
                
    }
    fclose(fp);
    *i_max=icnt-1;
  }

  TimeKeeper(X, cFileNameSzTimeKeep, cReadSzEnd, "a");

  return 0;
}
