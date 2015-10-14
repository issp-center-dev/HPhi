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

#include "mfmemory.h"
#include "check.h"
#include "wrapperMPI.h"

/** 
 * @brief A program to check size of dimension for hirbert-space.
 * 
 * @param[in,out] X  Common data set used in HPhi.
 * 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @return 
 */
int check(struct BindStruct *X){
    
  FILE *fp;
  long unsigned int i,tmp_sdim;
  int NLocSpn,NCond,Nup,Ndown;
  long unsigned int u_tmp;
  long unsigned int tmp;
  long unsigned int Ns,comb_1,comb_2,comb_3,comb_sum, comb_up, comb_down;
  int u_loc;
  int mfint[7];
  long int **comb;    
  
  Ns = X->Def.Nsite;
  li_malloc2(comb, Ns+1,Ns+1);
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
    for(i=0;i<X->Def.Nsite;i++){
      comb_sum= 2*comb_sum;     
    }
    break;

  case Hubbard:
    comb_up= Binomial(Ns, X->Def.Nup, comb, Ns);
    comb_down= Binomial(Ns, X->Def.Ndown, comb, Ns);
    comb_sum=comb_up*comb_down;
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
    comb_sum= Binomial(Ns, X->Def.Ne, comb, Ns);
    break;
  default:
    fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
    i_free2(comb, Ns+1, Ns+1);
    return -1;
  }  

  fprintf(stdoutMPI, "comb_sum= %ld \n",comb_sum);
  X->Check.idim_max = comb_sum;
  X->Check.max_mem=(3+2+1)*X->Check.idim_max*16.0/(pow(10,9));
  fprintf(stdoutMPI, "MAX DIMENSION idim_max=%ld \n",X->Check.idim_max);
  fprintf(stdoutMPI, "APPROXIMATE REQUIRED MEMORY  max_mem=%lf GB \n",X->Check.max_mem);
  if(childfopenMPI(cFileNameCheckMemory,"w", &fp)!=0){
    i_free2(comb, Ns+1, Ns+1);
    return -1;
  }
  fprintf(fp,"MAX DIMENSION idim_max=%ld \n",X->Check.idim_max);
  fprintf(fp,"APPROXIMATE REQUIRED MEMORY  max_mem=%lf GB \n",X->Check.max_mem);
  fcloseMPI(fp);

  //sdim 
  tmp=1;
  tmp_sdim=1;

  switch(X->Def.iCalcModel){
  case HubbardGC:
  case KondoGC:
  case Hubbard:
  case Kondo:
    while(tmp <= X->Def.Nsite){
      tmp_sdim=tmp_sdim*2;
      tmp+=1;
    }
    break;
  case Spin:
  case SpinGC:
    while(tmp <= X->Def.Nsite/2){
      tmp_sdim=tmp_sdim*2;
      tmp+=1;
    }
    break;
  default:
    fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
    i_free2(comb, Ns+1, Ns+1);
    return -1;
  }  
  X->Check.sdim=tmp_sdim;
  
  if(childfopenMPI(cFileNameCheckSdim,"w", &fp)!=0){
    i_free2(comb, Ns+1, Ns+1);
    return -1;
  }

  switch(X->Def.iCalcModel){
  case HubbardGC:
  case KondoGC:
  case Hubbard:
  case Kondo:
    fprintf(stdoutMPI, "sdim=%ld =2^%d\n",X->Check.sdim,X->Def.Nsite);
    fprintf(fp,"sdim=%ld =2^%d\n",X->Check.sdim,X->Def.Nsite);
    break;
  case Spin:
  case SpinGC:
    fprintf(stdoutMPI, "sdim=%ld =2^%d\n",X->Check.sdim,X->Def.Nsite/2);
    fprintf(fp,"sdim=%ld =2^%d\n",X->Check.sdim,X->Def.Nsite/2);
  default:
    break;
  }
  
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
  case Hubbard:
  case Kondo:
    for(i=1;i<=2*X->Def.Nsite-1;i++){
      u_tmp=u_tmp*2;
      X->Def.Tpow[i]=u_tmp;
      fprintf(fp,"%ld %ld \n",i,u_tmp);
    }
    break;
 case SpinGC:
    for(i=1;i<=X->Def.Nsite;i++){
      u_tmp=u_tmp*2;
      X->Def.Tpow[i]=u_tmp;
      fprintf(fp,"%ld %ld \n",i,u_tmp);
    }
    break;
 case Spin:
    for(i=1;i<=X->Def.Nsite-1;i++){
      u_tmp=u_tmp*2;
      X->Def.Tpow[i]=u_tmp;
      fprintf(fp,"%ld %ld \n",i,u_tmp);
    }
    break;
  default:
    fprintf(stdoutMPI, cErrNoModel, X->Def.iCalcModel);
    i_free2(comb, Ns+1, Ns+1);
    return -1;
  }  
  fcloseMPI(fp);	 
 
  i_free2(comb, Ns+1, Ns+1);
  fprintf(stdoutMPI, "Indices and Parameters of Definition files(*.def) are complete.\n");
  return 0;
}    
