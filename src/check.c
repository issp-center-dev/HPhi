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

#include "bitcalc.h"
#include "sz.h"
#include "FileIO.h"
#include "common/setmemory.h"
#include "check.h"
#include "wrapperMPI.h"
#include "CheckMPI.h"

/**
 * @file   check.c
 * @version 0.1, 0.2
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 * @brief  File for giving a function of calculating size of Hilbert space.
 * 
 */


/** 
 * @brief A program to check size of dimension for Hilbert-space.
 * 
 * @param[in,out] X  Common data set used in HPhi.
 * 
 * @retval TRUE normal termination
 * @retval FALSE abnormal termination
 * @retval MPIFALSE CheckMPI abnormal termination
 * @version 0.2
 * @details add function of calculating Hilbert space for canonical ensemble.
 *  
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int check(struct BindStruct *X){
    
  FILE *fp;
  long unsigned int i,tmp_sdim;
  int NLocSpn,NCond,Nup,Ndown;
  long unsigned int u_tmp;
  long unsigned int tmp;
  long unsigned int Ns,comb_1,comb_2,comb_3,comb_sum, comb_up, comb_down;
  int u_loc;
  long int **comb;
  long unsigned int idimmax=0;
  long unsigned int idim=0;
  long unsigned int isite=0;
  int tmp_sz=0;
  int iMinup=0;
  if(X->Def.iCalcModel ==Spin ||X->Def.iCalcModel ==SpinGC )
  {
    X->Def.Ne=X->Def.Nup;
  }

  int iAllup=X->Def.Ne;

  if(X->Def.iFlgScaLAPACK == 0) {
    /*
      Set Site number per MPI process
    */
    if (CheckMPI(X) != TRUE) {
      return MPIFALSE;
    }
  }
  else{
    X->Def.NsiteMPI = X->Def.Nsite;
    X->Def.Total2SzMPI = X->Def.Total2Sz;
  }

  Ns = X->Def.Nsite;

  comb = li_2d_allocate(Ns+1,Ns+1);

  //idim_max
  switch(X->Def.iCalcModel){
  case HubbardGC:
    //comb_sum = 2^(2*Ns)=4^Ns
    comb_sum = 1;
    for(i=0;i<2*X->Def.Nsite;i++){
      comb_sum= 2*comb_sum;     
    }
    break;
  case SpinGC:
    //comb_sum = 2^(Ns)
    comb_sum = 1;
    if(X->Def.iFlgGeneralSpin ==FALSE){
      for(i=0;i<X->Def.Nsite;i++){
        comb_sum= 2*comb_sum;     
      }
    }
    else{
      for(i=0; i<X->Def.Nsite;i++){
        comb_sum=comb_sum*X->Def.SiteToBit[i];
      }
    }
    break;

  case Hubbard:
    comb_up= Binomial(Ns, X->Def.Nup, comb, Ns);
    comb_down= Binomial(Ns, X->Def.Ndown, comb, Ns);
    comb_sum=comb_up*comb_down;
    break;

  case HubbardNConserved:
    comb_sum=0;
    if(X->Def.Ne > X->Def.Nsite){
      iMinup = X->Def.Ne-X->Def.Nsite;
      iAllup = X->Def.Nsite;
    }

    for(i=iMinup; i<= iAllup; i++){
      comb_up= Binomial(Ns, i, comb, Ns);
      comb_down= Binomial(Ns, X->Def.Ne-i, comb, Ns);
      comb_sum +=comb_up*comb_down;
    }
    break;
    
  case Kondo:
    //idim_max
    // calculation of dimension
    // Nup      = u_loc+u_cond
    // Ndown    = d_loc+d_cond
    // NLocSpn  = u_loc+d_loc
    // Ncond    = Nsite-NLocSpn
    // idim_max = \sum_{u_loc=0}^{u_loc=Nup} 
    //              Binomial(NLocSpn,u_loc)
    //             *Binomial(NCond,Nup-u_loc)
    //             *Binomial(NCond,Ndown+u_loc-NLocSpn)
    //comb_1 = Binomial(NLocSpn,u_loc)
    //comb_2 = Binomial(NCond,Nup-u_loc)
    //comb_3 = Binomial(NCond,Ndown+u_loc-NLocSpn)
    Nup     = X->Def.Nup;
    Ndown   = X->Def.Ndown;
    NCond   = X->Def.Nsite-X->Def.NLocSpn;
    NLocSpn = X->Def.NLocSpn;
    comb_sum = 0;
    for(u_loc=0;u_loc<=X->Def.Nup;u_loc++){
      comb_1     = Binomial(NLocSpn,u_loc,comb,Ns);
      comb_2     = Binomial(NCond,Nup-u_loc,comb,Ns);
      comb_3     = Binomial(NCond,Ndown+u_loc-NLocSpn,comb,Ns);
      comb_sum  += comb_1*comb_2*comb_3;
    }
    break;
  case KondoGC:
    comb_sum = 1;
    NCond   = X->Def.Nsite-X->Def.NLocSpn;
    NLocSpn = X->Def.NLocSpn;
    //4^Nc*2^Ns
    for(i=0;i<(2*NCond+NLocSpn);i++){
      comb_sum= 2*comb_sum;     
    }
    break;
  case Spin:

    if(X->Def.iFlgGeneralSpin ==FALSE){
      if(X->Def.Nup+X->Def.Ndown != X->Def.Nsite){
        fprintf(stderr, " 2Sz is incorrect.\n");
        return FALSE;
      }
      //comb_sum= Binomial(Ns, X->Def.Ne, comb, Ns);
      comb_sum= Binomial(Ns, X->Def.Nup, comb, Ns);
    }
    else{
      idimmax = 1;
      X->Def.Tpow[0]=idimmax;
      for(isite=0; isite<X->Def.Nsite;isite++){
        idimmax=idimmax*X->Def.SiteToBit[isite];
        X->Def.Tpow[isite+1]=idimmax;
      }
      comb_sum=0;
#pragma omp parallel for default(none) reduction(+:comb_sum) private(tmp_sz, isite) firstprivate(idimmax, X) 
      for(idim=0; idim<idimmax; idim++){
        tmp_sz=0;
        for(isite=0; isite<X->Def.Nsite;isite++){
          tmp_sz += GetLocal2Sz(isite+1,idim, X->Def.SiteToBit, X->Def.Tpow );          
        }
        if(tmp_sz == X->Def.Total2Sz){
          comb_sum +=1;
        }
      }
      
    }
    
    break;
  default:
    fprintf(stderr, cErrNoModel, X->Def.iCalcModel);
    free_li_2d_allocate(comb);
    return FALSE;
  }  

  //fprintf(stdoutMPI, "Debug: comb_sum= %ld \n",comb_sum);

  X->Check.idim_max = comb_sum;
  switch(X->Def.iCalcType) {
    case Lanczos:
      switch (X->Def.iCalcModel) {
        case Hubbard:
        case HubbardNConserved:
        case Kondo:
        case KondoGC:
        case Spin:
          X->Check.max_mem = 5.5 * X->Check.idim_max * 8.0 / (pow(10, 9));
          break;
        case HubbardGC:
        case SpinGC:
          X->Check.max_mem = 4.5 * X->Check.idim_max * 8.0 / (pow(10, 9));
          break;
      }
      break;
    case CG:
      switch (X->Def.iCalcModel) {
        case Hubbard:
        case HubbardNConserved:
        case Kondo:
        case KondoGC:
        case Spin:
          X->Check.max_mem = (6 * X->Def.k_exct + 2) * X->Check.idim_max * 16.0 / (pow(10, 9));
          break;
        case HubbardGC:
        case SpinGC:
          X->Check.max_mem = (6 * X->Def.k_exct + 1.5) * X->Check.idim_max * 16.0 / (pow(10, 9));
          break;
      }
      break;
    case TPQCalc:
    case cTPQ:
      switch (X->Def.iCalcModel) {
        case Hubbard:
        case HubbardNConserved:
        case Kondo:
        case KondoGC:
        case Spin:
          if (X->Def.iFlgCalcSpec != CALCSPEC_NOT) {
            X->Check.max_mem = (2) * X->Check.idim_max * 16.0 / (pow(10, 9));
          } else {
            X->Check.max_mem = 4.5 * X->Check.idim_max * 16.0 / (pow(10, 9));
          }
          break;
        case HubbardGC:
        case SpinGC:
          if (X->Def.iFlgCalcSpec != CALCSPEC_NOT) {
            X->Check.max_mem = (2) * X->Check.idim_max * 16.0 / (pow(10, 9));
          } else {
            X->Check.max_mem = 3.5 * X->Check.idim_max * 16.0 / (pow(10, 9));
          }
          break;
      }
      break;
    case FullDiag:
      X->Check.max_mem = X->Check.idim_max * 8.0 * X->Check.idim_max * 8.0 / (pow(10, 9));
      break;
    case TimeEvolution:
      X->Check.max_mem = (4 + 2 + 1) * X->Check.idim_max * 16.0 / (pow(10, 9));
      break;
    default:
      return FALSE;
      //break;
  }

  //fprintf(stdoutMPI, "  MAX DIMENSION idim_max=%ld \n",X->Check.idim_max);
  //fprintf(stdoutMPI, "  APPROXIMATE REQUIRED MEMORY  max_mem=%lf GB \n",X->Check.max_mem);
  unsigned long int li_dim_max=MaxMPI_li(X->Check.idim_max);
  fprintf(stdoutMPI, "  MAX DIMENSION idim_max=%ld \n",li_dim_max);
  double dmax_mem=MaxMPI_d(X->Check.max_mem);
  fprintf(stdoutMPI, "  APPROXIMATE REQUIRED MEMORY  max_mem=%lf GB \n",dmax_mem);
  if(childfopenMPI(cFileNameCheckMemory,"w", &fp)!=0){
    free_li_2d_allocate(comb);
    return FALSE;
  }
  fprintf(fp,"  MAX DIMENSION idim_max=%ld \n", li_dim_max);
  fprintf(fp,"  APPROXIMATE REQUIRED MEMORY  max_mem=%lf GB \n", dmax_mem);

  
  /*
  fprintf(fp,"  MAX DIMENSION idim_max=%ld \n",X->Check.idim_max);
  fprintf(fp,"  APPROXIMATE REQUIRED MEMORY  max_mem=%lf GB \n",X->Check.max_mem);
  */
  fclose(fp);

  //sdim 
  tmp=1;
  tmp_sdim=1;

  switch(X->Def.iCalcModel){
  case HubbardGC:
  case KondoGC:
  case HubbardNConserved:
  case Hubbard:
  case Kondo:
    while(tmp <= X->Def.Nsite){
      tmp_sdim=tmp_sdim*2;
      tmp+=1;
    }
    break;
  case Spin:
  case SpinGC:
    if(X->Def.iFlgGeneralSpin==FALSE){ 
      while(tmp <= X->Def.Nsite/2){
        tmp_sdim=tmp_sdim*2;
        tmp+=1;
      }
    }
    else{
      GetSplitBitForGeneralSpin(X->Def.Nsite, &tmp_sdim, X->Def.SiteToBit);
    }
    break;
  default:
    fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
    free_li_2d_allocate(comb);
    return FALSE;
  }  
  X->Check.sdim=tmp_sdim;
  
  if(childfopenMPI(cFileNameCheckSdim,"w", &fp)!=0){
    free_li_2d_allocate(comb);
    return FALSE;
  }

  switch(X->Def.iCalcModel){
  case HubbardGC:
  case KondoGC:
  case HubbardNConserved:
  case Hubbard:
  case Kondo:
    //fprintf(stdoutMPI, "sdim=%ld =2^%d\n",X->Check.sdim,X->Def.Nsite);
    fprintf(fp,"sdim=%ld =2^%d\n",X->Check.sdim,X->Def.Nsite);
    break;
  case Spin:
  case SpinGC:
    if(X->Def.iFlgGeneralSpin==FALSE){
      //fprintf(stdoutMPI, "sdim=%ld =2^%d\n",X->Check.sdim,X->Def.Nsite/2);
      fprintf(fp,"sdim=%ld =2^%d\n",X->Check.sdim,X->Def.Nsite/2);
    }
    break;
  default:
    break;
  }

  free_li_2d_allocate(comb);

  u_tmp=1;
  X->Def.Tpow[0]=u_tmp;
  switch(X->Def.iCalcModel){
  case HubbardGC:
  case KondoGC:
    for(i=1;i<=2*X->Def.Nsite;i++){
      u_tmp=u_tmp*2;
      X->Def.Tpow[i]=u_tmp;
      fprintf(fp,"%ld %ld \n",i,u_tmp);
    }
    break;
  case HubbardNConserved:
  case Hubbard:
  case Kondo:
    for(i=1;i<=2*X->Def.Nsite-1;i++){
      u_tmp=u_tmp*2;
      X->Def.Tpow[i]=u_tmp;
      fprintf(fp,"%ld %ld \n",i,u_tmp);
    }
    break;
 case SpinGC:
   if(X->Def.iFlgGeneralSpin==FALSE){
     for(i=1;i<=X->Def.Nsite;i++){
       u_tmp=u_tmp*2;
       X->Def.Tpow[i]=u_tmp;
       fprintf(fp,"%ld %ld \n",i,u_tmp);
     }
   }
   else{
     X->Def.Tpow[0]=u_tmp;
     fprintf(fp,"%d %ld \n", 0, u_tmp);
     for(i=1;i<X->Def.Nsite;i++){
       u_tmp=u_tmp*X->Def.SiteToBit[i-1];
       X->Def.Tpow[i]=u_tmp;
       fprintf(fp,"%ld %ld \n",i,u_tmp);
     }
   }
   break;
 case Spin:
   if(X->Def.iFlgGeneralSpin==FALSE){
     for(i=1;i<=X->Def.Nsite-1;i++){
       u_tmp=u_tmp*2;
       X->Def.Tpow[i]=u_tmp;
       fprintf(fp,"%ld %ld \n",i,u_tmp);
     }
   }
   else{
     for(i=0;i<X->Def.Nsite;i++){
       fprintf(fp,"%ld %ld \n",i,X->Def.Tpow[i]);
     }
   }     
    break;
  default:
    fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
    free_li_2d_allocate(comb);
    return FALSE;
  }  
  fclose(fp);
  /*
    Print MPI-site information and Modify Tpow 
    in the inter process region.
  */
  CheckMPI_Summary(X);
  
  return TRUE;
}    
