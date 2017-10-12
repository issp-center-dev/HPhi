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
#include "wrapperMPI.h"

static unsigned long int mfint[7];/*for malloc*/

void setmem_HEAD
(
 struct BindStruct *X
 )
{
  X->Def.CDataFileHead = (char*)malloc(D_FileNameMax*sizeof(char));
  X->Def.CParaFileHead = (char*)malloc(D_FileNameMax*sizeof(char));
}

void setmem_def
(
 struct BindStruct *X,
 struct BoostList *xBoost
 )
{
  unsigned long int i=0;
  lui_malloc1(X->Def.Tpow, 2*X->Def.Nsite+2);
  lui_malloc1(X->Def.OrgTpow, 2*X->Def.Nsite+2);
  for(i=0; i<2*X->Def.Nsite+2; i++){
    X->Def.Tpow[i]=0;
    X->Def.OrgTpow[i]=0;
  }
  li_malloc1(X->Def.SiteToBit, X->Def.Nsite+1);
  for(i=0; i<X->Def.Nsite+1; i++){
    X->Def.SiteToBit[i]=0;
  }
  
  i_malloc1(X->Def.LocSpn, X->Def.Nsite);
  d_malloc1(X->Phys.spin_real_cor, X->Def.Nsite*X->Def.Nsite);
  d_malloc1(X->Phys.charge_real_cor, X->Def.Nsite*X->Def.Nsite);
  d_malloc1(X->Phys.loc_spin_z, X->Def.Nsite*X->Def.Nsite);
  
  i_malloc1(X->Def.EDChemi, X->Def.EDNChemi+X->Def.NInterAll+X->Def.NTransfer);
  i_malloc1(X->Def.EDSpinChemi, X->Def.EDNChemi+X->Def.NInterAll+X->Def.NTransfer);
  d_malloc1(X->Def.EDParaChemi, X->Def.EDNChemi+X->Def.NInterAll+X->Def.NTransfer);

  i_malloc2(X->Def.EDGeneralTransfer, X->Def.NTransfer, 4);
  i_malloc2(X->Def.GeneralTransfer, X->Def.NTransfer, 4);
  c_malloc1(X->Def.EDParaGeneralTransfer, X->Def.NTransfer);
  c_malloc1(X->Def.ParaGeneralTransfer, X->Def.NTransfer);  
  
  i_malloc2(X->Def.CoulombIntra, X->Def.NCoulombIntra, 1);
  d_malloc1(X->Def.ParaCoulombIntra, X->Def.NCoulombIntra);
  i_malloc2(X->Def.CoulombInter, X->Def.NCoulombInter+X->Def.NIsingCoupling, 2);
  d_malloc1(X->Def.ParaCoulombInter, X->Def.NCoulombInter+X->Def.NIsingCoupling);
  i_malloc2(X->Def.HundCoupling, X->Def.NHundCoupling+X->Def.NIsingCoupling, 2);
  d_malloc1(X->Def.ParaHundCoupling, X->Def.NHundCoupling+X->Def.NIsingCoupling);
  i_malloc2(X->Def.PairHopping, X->Def.NPairHopping, 2);
  d_malloc1(X->Def.ParaPairHopping, X->Def.NPairHopping); 
  i_malloc2(X->Def.ExchangeCoupling, X->Def.NExchangeCoupling, 2);
  d_malloc1(X->Def.ParaExchangeCoupling, X->Def.NExchangeCoupling);
  i_malloc2(X->Def.PairLiftCoupling, X->Def.NPairLiftCoupling, 2);
  d_malloc1(X->Def.ParaPairLiftCoupling, X->Def.NPairLiftCoupling);

  i_malloc2(X->Def.InterAll, X->Def.NInterAll, 8);
  c_malloc1(X->Def.ParaInterAll, X->Def.NInterAll);
    
  i_malloc2(X->Def.CisAjt, X->Def.NCisAjt, 4);
  i_malloc2(X->Def.CisAjtCkuAlvDC, X->Def.NCisAjtCkuAlvDC, 8);

  i_malloc2(X->Def.SingleExcitationOperator, X->Def.NSingleExcitationOperator, 3);
  c_malloc1(X->Def.ParaSingleExcitationOperator, X->Def.NSingleExcitationOperator);
  i_malloc2(X->Def.PairExcitationOperator, X->Def.NPairExcitationOperator, 5);
  c_malloc1(X->Def.ParaPairExcitationOperator, X->Def.NPairExcitationOperator);

  unsigned int ipivot,iarrayJ,ispin;
  xBoost->list_6spin_star = (int **)malloc(sizeof(int*) * xBoost->R0 * xBoost->num_pivot);
  for (ipivot = 0; ipivot <  xBoost->R0 * xBoost->num_pivot; ipivot++) {
    xBoost->list_6spin_star[ipivot] = (int *)malloc(sizeof(int) * 7);
  }
  
  xBoost->list_6spin_pair = (int ***)malloc(sizeof(int**) * xBoost->R0 * xBoost->num_pivot);
  for (ipivot = 0; ipivot <  xBoost->R0 * xBoost->num_pivot; ipivot++) {
    xBoost->list_6spin_pair[ipivot] = (int **)malloc(sizeof(int*) * 7);
    for (ispin = 0; ispin < 7; ispin++) {
      xBoost->list_6spin_pair[ipivot][ispin] = (int *)malloc(sizeof(int) * 15);
    }
  }

  xBoost->arrayJ = (double complex ***)malloc(sizeof(double complex**) * xBoost->NumarrayJ);
for (iarrayJ = 0; iarrayJ < xBoost->NumarrayJ; iarrayJ++) {
  xBoost->arrayJ[iarrayJ] = (double complex **)malloc(sizeof(double complex*) * 3);
  for (i = 0; i < 3; i++) {
    xBoost->arrayJ[iarrayJ][i] = (double complex *)malloc(sizeof(double complex) * 3);
  }
}

  
}


int setmem_large
(
 struct BindStruct *X
 )
{

  unsigned long int j=0;
  unsigned int idim_maxMPI;
  
  idim_maxMPI = MaxMPI_li(X->Check.idim_max);

   if(GetlistSize(X)==TRUE) {
       lui_malloc1(list_1, X->Check.idim_max + 1);
#ifdef MPI
       lui_malloc1(list_1buf, idim_maxMPI + 1);
     for(j=0; j<X->Check.idim_max+1; j++){
       list_1buf[j]=0;
     }
#endif // MPI
       lui_malloc1(list_2_1, X->Large.SizeOflist_2_1);
       lui_malloc1(list_2_2, X->Large.SizeOflist_2_2);
       if(list_1==NULL
          || list_2_1==NULL
          || list_2_2==NULL
               )
       {
           return -1;
       }
       for(j=0; j<X->Check.idim_max+1; j++){
         list_1[j]=0;
       }
       for(j =0; j<X->Large.SizeOflist_2_1; j++){
           list_2_1[j]=0;
       }
       for(j =0; j<X->Large.SizeOflist_2_2; j++){
           list_2_2[j]=0;
       }
   }

  d_malloc1(list_Diagonal, X->Check.idim_max+1);
  c_malloc1(v0, X->Check.idim_max+1);
  c_malloc1(v1, X->Check.idim_max+1);
  for(j =0; j<X->Check.idim_max+1; j++){
    list_Diagonal[j]=0;
    v0[j]=0;
    v1[j]=0;
  }
#ifdef MPI
  c_malloc1(v1buf, idim_maxMPI + 1);
  for(j =0; j<X->Check.idim_max+1; j++) {
    v1buf[j]=0;
  }
#endif // MPI
  if (X->Def.iCalcType == TPQCalc) {c_malloc1(vg, 1); vg[0]=0;}
  else {
    c_malloc1(vg, X->Check.idim_max+1);
    for(j =0; j<X->Check.idim_max+1; j++) {
      vg[j]=0;
    } 
  }
  d_malloc1(alpha, X->Def.Lanczos_max+1);
  d_malloc1(beta, X->Def.Lanczos_max+1);

  if(
     list_Diagonal==NULL
     || v0==NULL
     || v1==NULL
     || vg==NULL
     ) {
    return -1;
  }


    if(X->Def.iCalcType == TPQCalc || X->Def.iFlgCalcSpec != CALCSPEC_NOT) {
        c_malloc2(vec, X->Def.Lanczos_max + 1, X->Def.Lanczos_max + 1);
    }
    else if(X->Def.iCalcType==Lanczos || X->Def.iCalcType == CG){
        if(X->Def.LanczosTarget>X->Def.nvec){
            c_malloc2(vec,X->Def.LanczosTarget+1, X->Def.Lanczos_max+1);
        }
        else{
            c_malloc2(vec,X->Def.nvec+1, X->Def.Lanczos_max+1);
        }
    }
    
  if(X->Def.iCalcType == FullDiag){
    d_malloc1(X->Phys.all_num_down, X->Check.idim_max+1);
    d_malloc1(X->Phys.all_num_up, X->Check.idim_max+1);
    d_malloc1(X->Phys.all_energy, X->Check.idim_max+1);
    d_malloc1(X->Phys.all_doublon, X->Check.idim_max+1);
    d_malloc1(X->Phys.all_sz, X->Check.idim_max+1);
    d_malloc1(X->Phys.all_s2, X->Check.idim_max+1);
    c_malloc2(Ham, X->Check.idim_max+1,X->Check.idim_max+1);
    c_malloc2(L_vec, X->Check.idim_max+1,X->Check.idim_max+1);

    if(X->Phys.all_num_down == NULL
       ||X->Phys.all_num_up == NULL
       ||X->Phys.all_energy == NULL
       ||X->Phys.all_doublon == NULL
       ||X->Phys.all_s2 ==NULL
       )
      {
	return -1;
      }
    for(j=0; j<X->Check.idim_max+1; j++){
      if(Ham[j]==NULL || L_vec[j]==NULL){
	return -1;
      }
    }
  }
  else if (X->Def.iCalcType == CG) {
    d_malloc1(X->Phys.all_num_down, X->Def.k_exct);
    d_malloc1(X->Phys.all_num_up, X->Def.k_exct);
    d_malloc1(X->Phys.all_energy, X->Def.k_exct);
    d_malloc1(X->Phys.all_doublon, X->Def.k_exct);
    d_malloc1(X->Phys.all_sz, X->Def.k_exct);
    d_malloc1(X->Phys.all_s2, X->Def.k_exct);
  }
  
  fprintf(stdoutMPI, "%s", cProFinishAlloc);
  return 0;
  }

void setmem_IntAll_Diagonal
(
 struct DefineList *X
 )
{
  i_malloc2(X->InterAll_OffDiagonal, X->NInterAll, 8);
  c_malloc1(X->ParaInterAll_OffDiagonal, X->NInterAll);

  i_malloc2(X->InterAll_Diagonal, X->NInterAll, 4);
  d_malloc1(X->ParaInterAll_Diagonal, X->NInterAll);
}


int GetlistSize
        (
                struct BindStruct *X
        )
{
   // unsigned int idim_maxMPI;

//    idim_maxMPI = MaxMPI_li(X->Check.idim_max);

    switch(X->Def.iCalcModel) {
        case Spin:
        case Hubbard:
        case HubbardNConserved:
        case Kondo:
        case KondoGC:
            if (X->Def.iFlgGeneralSpin == FALSE) {
                if (X->Def.iCalcModel == Spin && X->Def.Nsite % 2 == 1) {
                    X->Large.SizeOflist_2_1 = X->Check.sdim * 2 + 2;
                } else {
                    X->Large.SizeOflist_2_1 = X->Check.sdim + 2;
                }
                X->Large.SizeOflist_2_2 = X->Check.sdim + 2;
                X->Large.SizeOflistjb = X->Check.sdim + 2;
            } else {//for spin-canonical general spin
                X->Large.SizeOflist_2_1 = X->Check.sdim + 2;
                X->Large.SizeOflist_2_2 =
                        X->Def.Tpow[X->Def.Nsite - 1] * X->Def.SiteToBit[X->Def.Nsite - 1] / X->Check.sdim + 2;
                X->Large.SizeOflistjb =
                        X->Def.Tpow[X->Def.Nsite - 1] * X->Def.SiteToBit[X->Def.Nsite - 1] / X->Check.sdim + 2;
            }
            break;
        default:
            return FALSE;
    }
    return TRUE;
}

