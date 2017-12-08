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

/**
 * @file   xsetmem.c
 *
 * @brief  Set size of memories to be needed for calculation.
 * @version 2.0
 * @version 1.2
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */


static unsigned long int mfint[7];/*for malloc*/

///
/// \brief Set size of memories headers of output files.
/// \param X [out] BindStruct to get headers of files.\n
/// Output: CDataFileHead, CParaFileHead
/// \version 0.1
void setmem_HEAD
(
 struct BindStruct *X
 )
{
  X->Def.CDataFileHead = (char*)malloc(D_FileNameMax*sizeof(char));
  X->Def.CParaFileHead = (char*)malloc(D_FileNameMax*sizeof(char));
}

///
/// \brief Set size of memories for Def and Phys in BindStruct.
/// \param X [in,out] BindStruct to get information of Def and Phys structs.
/// \param xBoost [in,out] Struct for Boost mode.
/// \version 0.1
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

  i_malloc2(X->Def.GeneralTransfer, X->Def.NTransfer, 4);
  c_malloc1(X->Def.ParaGeneralTransfer, X->Def.NTransfer);

  if(X->Def.iCalcType == TimeEvolution){
    i_malloc2(X->Def.EDGeneralTransfer, X->Def.NTransfer+X->Def.NTETransferMax, 4);
    c_malloc1(X->Def.EDParaGeneralTransfer, X->Def.NTransfer+X->Def.NTETransferMax);

  }
  else {
    i_malloc2(X->Def.EDGeneralTransfer, X->Def.NTransfer, 4);
    c_malloc1(X->Def.EDParaGeneralTransfer, X->Def.NTransfer);
  }
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

  d_malloc1(X->Def.ParaLaser, X->Def.NLaser);

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

  int NInterAllSet;
  NInterAllSet= (X->Def.iCalcType==TimeEvolution) ? X->Def.NInterAll+X->Def.NTEInterAllMax: X->Def.NInterAll;
  i_malloc2(X->Def.InterAll_OffDiagonal, NInterAllSet, 8);
  c_malloc1(X->Def.ParaInterAll_OffDiagonal, NInterAllSet);
  i_malloc2(X->Def.InterAll_Diagonal, NInterAllSet, 4);
  d_malloc1(X->Def.ParaInterAll_Diagonal, NInterAllSet);

  if (X->Def.iCalcType == TimeEvolution){
    d_malloc1(X->Def.TETime, X->Def.NTETimeSteps);
    //Time-dependent Transfer
    ui_malloc1(X->Def.NTETransfer, X->Def.NTETimeSteps);
    ui_malloc1(X->Def.NTETransferDiagonal, X->Def.NTETimeSteps);
    i_malloc3(X->Def.TETransfer, X->Def.NTETimeSteps, X->Def.NTETransferMax, 4);
    i_malloc3(X->Def.TETransferDiagonal, X->Def.NTETimeSteps, X->Def.NTETransferMax, 2);
    c_malloc2(X->Def.ParaTETransfer, X->Def.NTETimeSteps, X->Def.NTETransferMax);
    d_malloc2(X->Def.ParaTETransferDiagonal, X->Def.NTETimeSteps,X->Def.NTETransferMax);
    //Time-dependent InterAll
    ui_malloc1(X->Def.NTEInterAll, X->Def.NTETimeSteps);
    ui_malloc1(X->Def.NTEInterAllDiagonal, X->Def.NTETimeSteps);
    i_malloc3(X->Def.TEInterAll, X->Def.NTETimeSteps, X->Def.NTEInterAllMax, 8);
    i_malloc3(X->Def.TEInterAllDiagonal, X->Def.NTETimeSteps, X->Def.NTEInterAllMax, 2);
    c_malloc2(X->Def.ParaTEInterAll, X->Def.NTETimeSteps, X->Def.NTEInterAllMax);
    d_malloc2(X->Def.ParaTEInterAllDiagonal, X->Def.NTETimeSteps,X->Def.NTEInterAllMax);

    ui_malloc1(X->Def.NTEInterAllOffDiagonal, X->Def.NTETimeSteps);
    i_malloc3(X->Def.TEInterAllOffDiagonal, X->Def.NTETimeSteps, X->Def.NTEInterAllMax, 8);
    c_malloc2(X->Def.ParaTEInterAllOffDiagonal, X->Def.NTETimeSteps,X->Def.NTEInterAllMax);

    //Time-dependent Chemi generated by InterAll diagonal complnents
    ui_malloc1(X->Def.NTEChemi, X->Def.NTETimeSteps);
    i_malloc2(X->Def.TEChemi,  X->Def.NTETimeSteps, X->Def.NTEInterAllMax);
    i_malloc2(X->Def.SpinTEChemi,  X->Def.NTETimeSteps, X->Def.NTEInterAllMax);
    d_malloc2(X->Def.ParaTEChemi, X->Def.NTETimeSteps, X->Def.NTEInterAllMax);
  }
}


///
/// \brief Set size of memories for Hamiltonian (Ham, L_vec), vectors(vg, v0, v1, v2, vec, alpha, beta), lists (list_1, list_2_1, list_2_2, list_Diagonal) and Phys(BindStruct.PhysList) struct in the case of Full Diag mode.
/// \param X [in,out] BindStruct to give information and give size of memories for Hamiltonian, vectors, lists and Phys struct in the case of Full Diag mode.
/// \retval -1 Fail to set memories.
/// \retval 0 Normal to set memories.
/// \version 0.1
int setmem_large
(
 struct BindStruct *X
 ) {

  unsigned long int j = 0;
  unsigned long int idim_maxMPI;

  idim_maxMPI = MaxMPI_li(X->Check.idim_max);

  if (GetlistSize(X) == TRUE) {
    lui_malloc1(list_1, X->Check.idim_max + 1);
#ifdef MPI
    lui_malloc1(list_1buf, idim_maxMPI + 1);
    for (j = 0; j < X->Check.idim_max + 1; j++) {
      list_1buf[j] = 0;
    }
#endif // MPI
    lui_malloc1(list_2_1, X->Large.SizeOflist_2_1);
    lui_malloc1(list_2_2, X->Large.SizeOflist_2_2);
    if (list_1 == NULL
        || list_2_1 == NULL
        || list_2_2 == NULL
            ) {
      return -1;
    }
    for (j = 0; j < X->Check.idim_max + 1; j++) {
      list_1[j] = 0;
    }
    for (j = 0; j < X->Large.SizeOflist_2_1; j++) {
      list_2_1[j] = 0;
    }
    for (j = 0; j < X->Large.SizeOflist_2_2; j++) {
      list_2_2[j] = 0;
    }
  }

  d_malloc1(list_Diagonal, X->Check.idim_max + 1);
  c_malloc1(v0, X->Check.idim_max + 1);
  c_malloc1(v1, X->Check.idim_max + 1);
  for (j = 0; j < X->Check.idim_max + 1; j++) {
    list_Diagonal[j] = 0;
    v0[j] = 0;
    v1[j] = 0;
    if (X->Def.iCalcType == TimeEvolution) {
      c_malloc1(v2, X->Check.idim_max + 1);
    } else {
      c_malloc1(v2, 1);
    }
#ifdef MPI
    c_malloc1(v1buf, idim_maxMPI + 1);
    for (j = 0; j < X->Check.idim_max + 1; j++) {
      v1buf[j] = 0;
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

    if (
            list_Diagonal == NULL
            || v0 == NULL
            || v1 == NULL
            || vg == NULL
            ) {
      return -1;
    }


    if (X->Def.iCalcType == TPQCalc || X->Def.iFlgCalcSpec != CALCSPEC_NOT) {
      c_malloc2(vec, X->Def.Lanczos_max + 1, X->Def.Lanczos_max + 1);
    } else if (X->Def.iCalcType == Lanczos || X->Def.iCalcType == CG) {
      if (X->Def.LanczosTarget > X->Def.nvec) {
        c_malloc2(vec, X->Def.LanczosTarget + 1, X->Def.Lanczos_max + 1);
      } else {
        c_malloc2(vec, X->Def.nvec + 1, X->Def.Lanczos_max + 1);
      }
    }

    if (X->Def.iCalcType == FullDiag) {
      d_malloc1(X->Phys.all_num_down, X->Check.idim_max + 1);
      d_malloc1(X->Phys.all_num_up, X->Check.idim_max + 1);
      d_malloc1(X->Phys.all_energy, X->Check.idim_max + 1);
      d_malloc1(X->Phys.all_doublon, X->Check.idim_max + 1);
      d_malloc1(X->Phys.all_sz, X->Check.idim_max + 1);
      d_malloc1(X->Phys.all_s2, X->Check.idim_max + 1);
      c_malloc2(Ham, X->Check.idim_max + 1, X->Check.idim_max + 1);
      c_malloc2(L_vec, X->Check.idim_max + 1, X->Check.idim_max + 1);

      if (X->Phys.all_num_down == NULL
          || X->Phys.all_num_up == NULL
          || X->Phys.all_energy == NULL
          || X->Phys.all_doublon == NULL
          || X->Phys.all_s2 == NULL
              ) {
        return -1;
      }
      for (j = 0; j < X->Check.idim_max + 1; j++) {
        if (Ham[j] == NULL || L_vec[j] == NULL) {
          return -1;
        }
      }
    } else if (X->Def.iCalcType == CG) {
      d_malloc1(X->Phys.all_num_down, X->Def.k_exct);
      d_malloc1(X->Phys.all_num_up, X->Def.k_exct);
      d_malloc1(X->Phys.all_energy, X->Def.k_exct);
      d_malloc1(X->Phys.all_doublon, X->Def.k_exct);
      d_malloc1(X->Phys.all_sz, X->Def.k_exct);
      d_malloc1(X->Phys.all_s2, X->Def.k_exct);
    }
    fprintf(stdoutMPI, "%s", cProFinishAlloc);
  }
  return 0;
}

///
/// \brief Set the size of memories for InterAllDiagonal and InterAllOffDiagonal arrays.
/// \param InterAllOffDiagonal [in/out] Arrays of cites and spin indexes of off-diagonal parts of InterAll interactions.
/// \param ParaInterAllOffDiagonal [in/out] Arrays of parameters of off-diagonal parts of InterAll interactions.
/// \param InterAllDiagonal [in/out] Arrays of cites and spin indexes of diagonal parts of InterAll interactions.
/// \param ParaInterAllDiagonal [in/out] Arrays of parameters of diagonal parts of InterAll interactions.
/// \param NInterAll [in] Total number of InterAll interactions.
/// \author Kazuyoshi Yoshimi
/// \version 1.2
  void setmem_IntAll_Diagonal
          (
                  int **InterAllOffDiagonal,
                  double complex *ParaInterAllOffDiagonal,
                  int **InterAllDiagonal,
                  double *ParaInterAllDiagonal,
                  const int NInterAll
          ) {
    i_malloc2(InterAllOffDiagonal, NInterAll, 8);
    c_malloc1(ParaInterAllOffDiagonal, NInterAll);
    i_malloc2(InterAllDiagonal, NInterAll, 4);
    d_malloc1(ParaInterAllDiagonal, NInterAll);
  }


///
/// \brief Set size of lists for the canonical ensemble.
/// \param X [in,out] Give the information for getting the list size and get the lists.\n
/// Input: DefineList.iFlgGeneralSpin, DefineList.iCalcModel, DefineList.Nsite, CheckList.sdim, DefineList.Tpow, DefineList.SiteToBit\n
/// Output: LargeList.SizeOflist_2_1, LargeList.SizeOflist_2_2, LargeList.SizeOflistjb
/// \retval TRUE: Normally finished
/// \retval FALSE: Unnormally finished
/// \author Kazuyoshi Yoshimi
/// \version 1.2
  int GetlistSize
          (
                  struct BindStruct *X
          ) {
    // unsigned int idim_maxMPI;

//    idim_maxMPI = MaxMPI_li(X->Check.idim_max);

    switch (X->Def.iCalcModel) {
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
/**
@page page_setmem Malloc vectors

setmem_def()
setmem_large()
mfmemmory.h
*/