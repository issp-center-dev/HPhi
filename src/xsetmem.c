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

/**
 * @file xsetmem.c
 *
 * @brief Memory allocation for HPhi data structures
 *
 * Allocates memory for all major data structures used in HPhi calculations.
 * Memory allocation is organized by structure type:
 *
 * setmem_HEAD: File name headers
 * setmem_def: Definition struct (Hamiltonian parameters)
 *   - Tpow[]: Bit position powers (2^i for site i)
 *   - Transfer, CoulombIntra/Inter, Exchange, etc.
 *
 * setmem_large: Large vectors and temporary buffers
 *   - v0, v1: Lanczos vectors (size: idim_max)
 *   - v1buf: MPI receive buffer
 *   - list_1, list_2_1, list_2_2: Hilbert space indexing
 *   - list_Diagonal: Pre-computed diagonal elements
 *
 * Key arrays:
 * - Tpow[i]: bit mask for the i-th degree of freedom. = 2^i for Hubbard
 *   (2 bits/site; i = 2*site + spin) and for Spin / SpinlessFermion
 *   (1 bit/site; i = site). For general spin, Tpow is the running
 *   product of preceding sites' local Hilbert dimensions.
 * - list_1[j]: Maps restricted index j to full bit representation
 * - list_Diagonal[j]: Diagonal Hamiltonian element for state j
 *
 * @version 2.0, 1.2, 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
#include "Common.h"
#include "common/setmemory.h"
#include "xsetmem.h"
#include "wrapperMPI.h"

/**
 * @brief Allocate memory for output file name headers
 *
 * Allocates CDataFileHead and CParaFileHead strings used as
 * prefixes for output file names (e.g., "zvo" -> "zvo_energy.dat").
 *
 * @param X BindStruct to receive allocated headers [out]
 */
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
 ) {
  X->Def.Tpow = lui_1d_allocate(2 * X->Def.Nsite + 2);
  X->Def.OrgTpow = lui_1d_allocate(2 * X->Def.Nsite + 2);
  X->Def.SiteToBit = li_1d_allocate(X->Def.Nsite + 1);
  X->Def.LocSpn = i_1d_allocate(X->Def.Nsite);
  X->Phys.spin_real_cor = d_1d_allocate(X->Def.Nsite * X->Def.Nsite);
  X->Phys.charge_real_cor = d_1d_allocate(X->Def.Nsite * X->Def.Nsite);
  X->Phys.loc_spin_z = d_1d_allocate(X->Def.Nsite * X->Def.Nsite);

  X->Def.EDChemi = i_1d_allocate(X->Def.EDNChemi + X->Def.NInterAll + X->Def.NTransfer);
  X->Def.EDSpinChemi = i_1d_allocate(X->Def.EDNChemi + X->Def.NInterAll + X->Def.NTransfer);
  X->Def.EDParaChemi = d_1d_allocate(X->Def.EDNChemi + X->Def.NInterAll + X->Def.NTransfer);
  X->Def.GeneralTransfer = i_2d_allocate(X->Def.NTransfer, 4);
  X->Def.ParaGeneralTransfer = cd_1d_allocate(X->Def.NTransfer);

  if (X->Def.iCalcType == TimeEvolution) {
    X->Def.EDGeneralTransfer = i_2d_allocate(X->Def.NTransfer + X->Def.NTETransferMax, 4);
    X->Def.EDParaGeneralTransfer = cd_1d_allocate(X->Def.NTransfer + X->Def.NTETransferMax);
  } else {
    X->Def.EDGeneralTransfer = i_2d_allocate(X->Def.NTransfer, 4);
    X->Def.EDParaGeneralTransfer = cd_1d_allocate(X->Def.NTransfer);
  }

  X->Def.CoulombIntra = i_2d_allocate(X->Def.NCoulombIntra, 1);
  X->Def.ParaCoulombIntra = d_1d_allocate(X->Def.NCoulombIntra);
  X->Def.CoulombInter = i_2d_allocate(X->Def.NCoulombInter + X->Def.NIsingCoupling, 2);
  X->Def.ParaCoulombInter = d_1d_allocate(X->Def.NCoulombInter + X->Def.NIsingCoupling);
  X->Def.HundCoupling = i_2d_allocate(X->Def.NHundCoupling + X->Def.NIsingCoupling, 2);
  X->Def.ParaHundCoupling = d_1d_allocate(X->Def.NHundCoupling + X->Def.NIsingCoupling);
  X->Def.PairHopping = i_2d_allocate(X->Def.NPairHopping, 2);
  X->Def.ParaPairHopping = d_1d_allocate(X->Def.NPairHopping);
  X->Def.ExchangeCoupling = i_2d_allocate(X->Def.NExchangeCoupling, 2);
  X->Def.ParaExchangeCoupling = d_1d_allocate(X->Def.NExchangeCoupling);
  X->Def.PairLiftCoupling = i_2d_allocate(X->Def.NPairLiftCoupling, 2);
  X->Def.ParaPairLiftCoupling = d_1d_allocate(X->Def.NPairLiftCoupling);

  X->Def.InterAll = i_2d_allocate(X->Def.NInterAll, 8);
  X->Def.ParaInterAll = cd_1d_allocate(X->Def.NInterAll);

  {
    const unsigned int nbody_terms = X->Def.NNBodyInterAll > 0 ? X->Def.NNBodyInterAll : 1;
    const unsigned int nbody_factors = X->Def.NBodyInterAll_TotalFactors > 0 ? X->Def.NBodyInterAll_TotalFactors : 1;
    X->Def.NBodyInterAll_N = ui_1d_allocate(nbody_terms);
    X->Def.NBodyInterAll_Offset = ui_1d_allocate(nbody_terms);
    X->Def.NBodyInterAll_Factors = i_2d_allocate(nbody_factors, 4);
    X->Def.ParaNBodyInterAll = cd_1d_allocate(nbody_terms);
    X->Def.NBodyInterAll_CanonicalN = ui_1d_allocate(nbody_terms);
    X->Def.NBodyInterAll_CanonicalOffset = ui_1d_allocate(nbody_terms);
    X->Def.NBodyInterAll_CanonicalFactors = i_2d_allocate(nbody_factors, 4);
    X->Def.NBodyInterAll_DiagonalIndex = ui_1d_allocate(nbody_terms);
    X->Def.NBodyInterAll_OffDiagonalIndex = ui_1d_allocate(nbody_terms);
  }

  {
    const unsigned int nbodyg_terms = X->Def.NNBodyG > 0 ? X->Def.NNBodyG : 1;
    const unsigned int nbodyg_factors = X->Def.NBodyG_TotalFactors > 0 ? X->Def.NBodyG_TotalFactors : 1;
    X->Def.NBodyG_N = ui_1d_allocate(nbodyg_terms);
    X->Def.NBodyG_Offset = ui_1d_allocate(nbodyg_terms);
    X->Def.NBodyG_Factors = i_2d_allocate(nbodyg_factors, 4);
    X->Def.NBodyG_IsZero = i_1d_allocate(nbodyg_terms);
    X->Def.NBodyG_CanonicalN = ui_1d_allocate(nbodyg_terms);
    X->Def.NBodyG_CanonicalOffset = ui_1d_allocate(nbodyg_terms);
    X->Def.NBodyG_CanonicalFactors = i_2d_allocate(nbodyg_factors, 4);
  }

  {
    const unsigned int anomalous_terms = X->Def.NAnomalousTerm > 0 ? X->Def.NAnomalousTerm : 1;
    const unsigned int anomalousg_terms = X->Def.NAnomalousG > 0 ? X->Def.NAnomalousG : 1;
    X->Def.AnomalousTerm = i_2d_allocate(anomalous_terms, 5);
    X->Def.ParaAnomalousTerm = cd_1d_allocate(anomalous_terms);
    X->Def.AnomalousG = i_2d_allocate(anomalousg_terms, 5);
  }

  X->Def.CisAjt = i_2d_allocate(X->Def.NCisAjt, 4);
  X->Def.CisAjtCkuAlvDC = i_2d_allocate(X->Def.NCisAjtCkuAlvDC, 8);

  X->Def.TBody = i_2d_allocate(X->Def.NTBody, 12);
  X->Def.FBody = i_2d_allocate(X->Def.NFBody, 16);
  X->Def.SBody = i_2d_allocate(X->Def.NSBody, 24);

  X->Def.SingleExcitationOperator = i_2d_allocate(X->Def.NSingleExcitationOperator, 3);
  X->Def.ParaSingleExcitationOperator = cd_1d_allocate(X->Def.NSingleExcitationOperator);
  X->Def.PairExcitationOperator = i_2d_allocate(X->Def.NPairExcitationOperator, 5);
  X->Def.ParaPairExcitationOperator = cd_1d_allocate(X->Def.NPairExcitationOperator);

  X->Def.SingleExcitationOperatorBra = i_2d_allocate(X->Def.NSingleExcitationOperatorBra, 3);
  X->Def.ParaSingleExcitationOperatorBra = cd_1d_allocate(X->Def.NSingleExcitationOperatorBra);
  X->Def.PairExcitationOperatorBra = i_2d_allocate(X->Def.NPairExcitationOperatorBra, 5);
  X->Def.ParaPairExcitationOperatorBra = cd_1d_allocate(X->Def.NPairExcitationOperatorBra);

  X->Def.ParaLaser = d_1d_allocate(X->Def.NLaser);

  xBoost->list_6spin_star = i_2d_allocate(xBoost->R0 * xBoost->num_pivot, 7);
  xBoost->list_6spin_pair = i_3d_allocate(xBoost->R0 * xBoost->num_pivot, 7, 15);
  xBoost->arrayJ = cd_3d_allocate(xBoost->NumarrayJ, 3, 3);

  int NInterAllSet;
  NInterAllSet = (X->Def.iCalcType == TimeEvolution) ? X->Def.NInterAll + X->Def.NTEInterAllMax : X->Def.NInterAll;
  X->Def.InterAll_OffDiagonal = i_2d_allocate(NInterAllSet, 8);
  X->Def.ParaInterAll_OffDiagonal = cd_1d_allocate(NInterAllSet);
  X->Def.InterAll_Diagonal = i_2d_allocate(NInterAllSet, 4);
  X->Def.ParaInterAll_Diagonal = d_1d_allocate(NInterAllSet);

  if (X->Def.iCalcType == TimeEvolution) {
    X->Def.TETime = d_1d_allocate(X->Def.NTETimeSteps);
    //Time-dependent Transfer
    X->Def.NTETransfer = ui_1d_allocate(X->Def.NTETimeSteps);
    X->Def.NTETransferDiagonal = ui_1d_allocate(X->Def.NTETimeSteps);
    X->Def.TETransfer = i_3d_allocate(X->Def.NTETimeSteps, X->Def.NTETransferMax, 4);
    X->Def.TETransferDiagonal = i_3d_allocate(X->Def.NTETimeSteps, X->Def.NTETransferMax, 2);
    X->Def.ParaTETransfer = cd_2d_allocate(X->Def.NTETimeSteps, X->Def.NTETransferMax);
    X->Def.ParaTETransferDiagonal = d_2d_allocate(X->Def.NTETimeSteps, X->Def.NTETransferMax);
    //Time-dependent InterAll
    X->Def.NTEInterAll = ui_1d_allocate(X->Def.NTETimeSteps);
    X->Def.NTEInterAllDiagonal = ui_1d_allocate(X->Def.NTETimeSteps);
    X->Def.TEInterAll = i_3d_allocate(X->Def.NTETimeSteps, X->Def.NTEInterAllMax, 8);
    X->Def.TEInterAllDiagonal = i_3d_allocate(X->Def.NTETimeSteps, X->Def.NTEInterAllMax, 4);
    X->Def.ParaTEInterAll = cd_2d_allocate(X->Def.NTETimeSteps, X->Def.NTEInterAllMax);
    X->Def.ParaTEInterAllDiagonal = d_2d_allocate(X->Def.NTETimeSteps, X->Def.NTEInterAllMax);
    X->Def.NTEInterAllOffDiagonal = ui_1d_allocate(X->Def.NTETimeSteps);
    X->Def.TEInterAllOffDiagonal = i_3d_allocate(X->Def.NTETimeSteps, X->Def.NTEInterAllMax, 8);
    X->Def.ParaTEInterAllOffDiagonal = cd_2d_allocate(X->Def.NTETimeSteps, X->Def.NTEInterAllMax);
    //Time-dependent Chemi generated by InterAll diagonal components
    X->Def.NTEChemi = ui_1d_allocate(X->Def.NTETimeSteps);
    X->Def.TEChemi = i_2d_allocate(X->Def.NTETimeSteps, X->Def.NTEInterAllMax);
    X->Def.SpinTEChemi = i_2d_allocate(X->Def.NTETimeSteps, X->Def.NTEInterAllMax);
    X->Def.ParaTEChemi = d_2d_allocate(X->Def.NTETimeSteps, X->Def.NTEInterAllMax);
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
      list_1 = lui_1d_allocate(X->Check.idim_max + 1);
#ifdef MPI
      list_1buf = lui_1d_allocate(idim_maxMPI + 1);
#endif // MPI
      list_2_1 = lui_1d_allocate(X->Large.SizeOflist_2_1);
      list_2_2 = lui_1d_allocate(X->Large.SizeOflist_2_2);
      if (list_1 == NULL
          || list_2_1 == NULL
          || list_2_2 == NULL
              ) {
          return -1;
      }
  }

    list_Diagonal = d_1d_allocate(X->Check.idim_max + 1);
    v0 = cd_1d_allocate(X->Check.idim_max + 1);
    v1 = cd_1d_allocate(X->Check.idim_max + 1);
  if (X->Def.iCalcType == TimeEvolution || X->Def.iCalcType == cTPQ) {
      v2 = cd_1d_allocate(X->Check.idim_max + 1);
  } else {
      v2 = cd_1d_allocate(1);
  }
#ifdef MPI
    v1buf = cd_1d_allocate(idim_maxMPI + 1);
#endif // MPI
  if (X->Def.iCalcType == TPQCalc || X->Def.iCalcType == cTPQ) {
      vg = cd_1d_allocate(1);
  } else {
      vg = cd_1d_allocate(X->Check.idim_max + 1);
  }
    alpha = d_1d_allocate(X->Def.Lanczos_max + 1);
    beta = d_1d_allocate(X->Def.Lanczos_max + 1);

  if (
          list_Diagonal == NULL
          || v0 == NULL
          || v1 == NULL
          || vg == NULL
          ) {
    return -1;
  }

  if (X->Def.iCalcType == TPQCalc || X->Def.iCalcType == cTPQ || X->Def.iFlgCalcSpec != CALCSPEC_NOT) {
      vec = cd_2d_allocate(X->Def.Lanczos_max + 1, X->Def.Lanczos_max + 1);
  } else if (X->Def.iCalcType == Lanczos || X->Def.iCalcType == CG) {
    if (X->Def.LanczosTarget > X->Def.nvec) {
        vec = cd_2d_allocate(X->Def.LanczosTarget + 1, X->Def.Lanczos_max + 1);
    } else {
        vec = cd_2d_allocate(X->Def.nvec + 1, X->Def.Lanczos_max + 1);
    }
  }

  if (X->Def.iCalcType == FullDiag) {
      X->Phys.all_num_down = d_1d_allocate(X->Check.idim_max + 1);
      X->Phys.all_num_up = d_1d_allocate( X->Check.idim_max + 1);
      X->Phys.all_energy = d_1d_allocate(X->Check.idim_max + 1);
      X->Phys.all_doublon = d_1d_allocate(X->Check.idim_max + 1);
      X->Phys.all_sz = d_1d_allocate(X->Check.idim_max + 1);
      X->Phys.all_s2 = d_1d_allocate(X->Check.idim_max + 1);
#ifdef _ELPA
      if (X->Def.iSolver == SOLVER_ELPA && nproc > 1) {
        /* Distributed-panel mode (design doc sec. 3 phase 2):
           each rank stores only its owned 1D column block.
           NC = ceil(N/P); rank p owns 1-based columns
           [p*NC+1, min((p+1)*NC, N)]. Ham/L_vec stay unallocated. */
        long int NN = X->Check.idim_max;
        long int NC = (NN + nproc - 1) / nproc;
        long int jb = (long int)myrank * NC + 1;
        long int je = ((long int)myrank + 1) * NC;
        if (je > NN) je = NN;
        if (jb > NN) { jb = 1; je = 0; } /* rank owns no column */
        HamColBegin = jb;
        HamColEnd = je;
        HamPanelLd = NN;
        iHamPanelActive = 1;
        {
          long int ncols = (je >= jb) ? (je - jb + 1) : 0;
          long int nelem = NN * ncols;
          if (nelem < 1) nelem = 1;
          Ham_local = (double complex *)malloc(nelem * sizeof(double complex));
          if (Ham_local == NULL) return -1;
          for (long int k = 0; k < nelem; k++) Ham_local[k] = 0.0;
        }
      } else
#endif
      {
        Ham = cd_2d_allocate(X->Check.idim_max + 1, X->Check.idim_max + 1);
        L_vec = cd_2d_allocate(X->Check.idim_max + 1, X->Check.idim_max + 1);
      }

    if (X->Phys.all_num_down == NULL
        || X->Phys.all_num_up == NULL
        || X->Phys.all_energy == NULL
        || X->Phys.all_doublon == NULL
        || X->Phys.all_s2 == NULL
            ) {
      return -1;
    }
    if (!iHamPanelActive) {
      for (j = 0; j < X->Check.idim_max + 1; j++) {
        if (Ham[j] == NULL || L_vec[j] == NULL) {
          return -1;
        }
      }
    }
  } else if (X->Def.iCalcType == CG) {
      X->Phys.all_num_down = d_1d_allocate(X->Def.k_exct);
      X->Phys.all_num_up = d_1d_allocate(X->Def.k_exct);
      X->Phys.all_energy = d_1d_allocate(X->Def.k_exct);
      X->Phys.all_doublon = d_1d_allocate(X->Def.k_exct);
      X->Phys.all_sz = d_1d_allocate(X->Def.k_exct);
      X->Phys.all_s2 = d_1d_allocate( X->Def.k_exct);
  }
  fprintf(stdoutMPI, "%s", cProFinishAlloc);
  return 0;
}

///
/// \brief Set the size of memories for InterAllDiagonal and InterAllOffDiagonal arrays.
/// \param InterAllOffDiagonal [in,out] Arrays of cites and spin indexes of off-diagonal parts of InterAll interactions.
/// \param ParaInterAllOffDiagonal [in,out] Arrays of parameters of off-diagonal parts of InterAll interactions.
/// \param InterAllDiagonal [in,out] Arrays of cites and spin indexes of diagonal parts of InterAll interactions.
/// \param ParaInterAllDiagonal [in,out] Arrays of parameters of diagonal parts of InterAll interactions.
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
    InterAllOffDiagonal = i_2d_allocate(NInterAll, 8);
    ParaInterAllOffDiagonal = cd_1d_allocate(NInterAll);
    InterAllDiagonal = i_2d_allocate(NInterAll, 4);
    ParaInterAllDiagonal = d_1d_allocate(NInterAll);
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
      case SpinlessFermion:
      case Hubbard:
      case HubbardNConserved:
      case Kondo:
      case KondoNConserved:
      case KondoGC:
      case tJ:
      case tJNConserved:
      case tJGC:
        if (X->Def.iFlgGeneralSpin == FALSE) {
          if ((X->Def.iCalcModel == Spin || X->Def.iCalcModel == SpinlessFermion) && X->Def.Nsite % 2 == 1) {
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

 To set memory, we use ```?malloc%``` function defined in @c mfmemmory.h,
 where ```?``` indicates the type of the array and ```%``` means the dimension.

 For example, ```char_malloc2(X, N1, N2)``` function sets the size of memories N1@f$ \times @f$ N2 characters to two dimensional array X.

 To set memories to global arrays, we prepare two functions, setmem_def() and setmem_large() functions.

 - setmem_def()

    In this function, the memories of the arrays which do not have large memory are stored.

    Arrays for defining interactions and correlation functions are mainly defined.

 - setmem_large()

    In this function, the memories of the arrays which have large memory are stored.

    Arrays for defining Hamiltonian and vectors are mainly defined.

@sa setmem_def(), setmem_large()
*/
