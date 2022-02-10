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

#include <bitcalc.h>
#include "mltplyCommon.h"
#include "mltplyHubbardCore.h"
#include "mltplySpinCore.h"
#include "makeHam_medium.h"
#include "wrapperMPI.h"

/**
 * @file   makeHam.c
 * 
 * @brief  Making Hamiltonian for the full diagonalization method.
 * 
 * @version 0.2
 * @details add function to treat the case of generalspin
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)

 */


/** 
 * @brief Making Hamiltonian for the full diagonalization method.\n
 * The Hamiltonian is stored in the two dimensional array @f$ \verb|Ham| @f$.
 * 
 * @param X [in] Struct for getting the information of the operators
 * 
 * @retval 0  normally finished
 * @retval -1 unnormally finished
 * 
 * @version 0.2
 * @details add function to treat the case of generalspin
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int makeHam_medium(struct BindStruct *X) {

  long unsigned int i, j;
  long unsigned int is1_spin;
  long unsigned int irght, ilft, ihfbit;
  double complex dmv;
  double num1;
  long unsigned int off;
  long unsigned int isite1, isite2, isite3, isite4;
  int sigma1, sigma2, sigma3, sigma4;
  long unsigned int isA_up, isB_up;
  double complex tmp_trans, tmp_V;
  long unsigned int Asum, Bsum, Adiff, Bdiff;
  long unsigned int tmp_off, tmp_off_2;
  int tmp_sgn;
  off = 0;
  tmp_off = 0;
  tmp_off_2 = 0;
  long unsigned int i_max;
  i_max = X->Check.idim_max;
  int ihermite = 0;
  int idx = 0;

  if (GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit) != 0) {
    return -1;
  }
  X->Large.i_max = i_max;
  X->Large.irght = irght;
  X->Large.ilft = ilft;
  X->Large.ihfbit = ihfbit;
  X->Large.prdct = 0.0;
  X->Large.mode = M_Ham;
  unsigned long int NHam_offdiagonal = 0;
  NHam_offdiagonal = X->Def.EDNTransfer+X->Def.NPairHopping+X->Def.NExchangeCoupling+X->Def.NPairLiftCoupling+X->Def.NInterAll;

  for (i = 0; i <= NHam_offdiagonal; i++) {
#pragma omp parallel for default(none) firstprivate(i, i_max) shared(Ham)
      for (j = 0; j <= i_max; j++) {
      Ham[i][j] = 0.0;
    }
  }

#pragma omp parallel for default(none) firstprivate(i_max) private(j) shared(Ham, list_Diagonal, v0, v1)
  for (j = 1; j <= i_max; j++) {
    Ham[0][j] += list_Diagonal[j];
    v0[j] = 1.0;
    v1[j] = 1.0;
  }
  long unsigned int count_idx;
  count_idx = 1; // 0 is used for diagonal component

  switch (X->Def.iCalcModel) {
    case HubbardGC:
      //Transfer
      for (i = 0; i < X->Def.EDNTransfer / 2; i++) {
        for (ihermite = 0; ihermite < 2; ihermite++) {
          idx = 2 * i + ihermite;
          isite1 = X->Def.EDGeneralTransfer[idx][0] + 1;
          isite2 = X->Def.EDGeneralTransfer[idx][2] + 1;
          sigma1 = X->Def.EDGeneralTransfer[idx][1];
          sigma2 = X->Def.EDGeneralTransfer[idx][3];

          if (general_hopp_GetInfo(X, isite1, isite2, sigma1, sigma2) != 0) {
            return -1;
          }
          tmp_trans = -X->Def.EDParaGeneralTransfer[idx];

#pragma omp parallel for default(none) firstprivate(X, tmp_trans, count_idx) private(j, tmp_off, dmv) shared(Ham, lui_counter_vec, v0, v1)
          for (j = 1; j <= X->Large.i_max; j++) {
            dmv = tmp_trans *
                  GC_CisAjt(j, v0, v1, X, X->Large.is1_spin, X->Large.is2_spin, X->Large.isA_spin, X->Large.A_spin,
                            tmp_trans, &tmp_off);
            Ham[count_idx][j] += dmv;
            lui_counter_vec[count_idx][j] = tmp_off+1;
          }
          count_idx++;
        }
      }


      for (i = 0; i < X->Def.NInterAll_OffDiagonal / 2; i++) {
          for (ihermite = 0; ihermite < 2; ihermite++) {
              idx = 2 * i + ihermite;
              isite1 = X->Def.InterAll_OffDiagonal[idx][0] + 1;
              isite2 = X->Def.InterAll_OffDiagonal[idx][2] + 1;
              isite3 = X->Def.InterAll_OffDiagonal[idx][4] + 1;
              isite4 = X->Def.InterAll_OffDiagonal[idx][6] + 1;
              sigma1 = X->Def.InterAll_OffDiagonal[idx][1];
              sigma2 = X->Def.InterAll_OffDiagonal[idx][3];
              sigma3 = X->Def.InterAll_OffDiagonal[idx][5];
              sigma4 = X->Def.InterAll_OffDiagonal[idx][7];
              tmp_V = X->Def.ParaInterAll_OffDiagonal[idx];
              general_int_GetInfo(
                      i,
                      X,
                      isite1,
                      isite2,
                      isite3,
                      isite4,
                      sigma1,
                      sigma2,
                      sigma3,
                      sigma4,
                      tmp_V
              );

              i_max = X->Large.i_max;
              isite1 = X->Large.is1_spin;
              isite2 = X->Large.is2_spin;
              Asum = X->Large.isA_spin;
              Adiff = X->Large.A_spin;

              isite3 = X->Large.is3_spin;
              isite4 = X->Large.is4_spin;
              Bsum = X->Large.isB_spin;
              Bdiff = X->Large.B_spin;

              tmp_V = X->Large.tmp_V;

              if (isite1 == isite2 && isite3 == isite4) { //diagonal
                  for (j = 1; j <= i_max; j++) {
                      dmv = GC_CisAisCisAis_element(j, isite1, isite3, tmp_V, v0, v1, X, &tmp_off);
                      Ham[0][j] += dmv;
                  }
              } else {
                  if (isite1 == isite2 && isite3 != isite4) {
                      for (j = 1; j <= i_max; j++) {
                          dmv = GC_CisAisCjtAku_element(j, isite1, isite3, isite4, Bsum, Bdiff, tmp_V, v0, v1, X,
                                                        &tmp_off);
                          Ham[count_idx][j] += dmv;
                          lui_counter_vec[count_idx][j] = tmp_off + 1;
                      }
                  } else if (isite1 != isite2 && isite3 == isite4) {
                      for (j = 1; j <= i_max; j++) {
                          dmv = GC_CisAjtCkuAku_element(j, isite1, isite2, isite3, Asum, Adiff, tmp_V, v0, v1, X,
                                                        &tmp_off);
                          Ham[count_idx][j] += dmv;
                          lui_counter_vec[count_idx][j] = tmp_off + 1;
                      }
                  } else if (isite1 != isite2 && isite3 != isite4) {
                      for (j = 1; j <= i_max; j++) {
                          dmv = GC_CisAjtCkuAlv_element(j, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff,
                                                        tmp_V,
                                                        v0, v1, X, &tmp_off_2);
                          Ham[count_idx][j] += dmv;
                          lui_counter_vec[count_idx][j] = tmp_off_2 + 1;
                      }
                  }
                  count_idx++;
              }
          }
      }
      //Pair hopping
      for (i = 0; i < X->Def.NPairHopping / 2; i++) {
        for (ihermite = 0; ihermite < 2; ihermite++) {
          idx = 2 * i + ihermite;
          pairhopp_GetInfo(idx, X);
          for (j = 1; j <= X->Large.i_max; j++) {
            dmv = GC_pairhopp_element(j, v0, v1, X, &tmp_off);
            Ham[count_idx][j] += dmv;
            lui_counter_vec[count_idx][j] = tmp_off+1;
          }
          count_idx++;
        }
      }
      //Exchange
      for (i = 0; i < X->Def.NExchangeCoupling; i++) {
        exchange_GetInfo(i, X);
        for (j = 1; j <= X->Large.i_max; j++) {
          dmv = GC_exchange_element(j, v0, v1, X, &tmp_off);
          Ham[count_idx][j] += dmv;
          lui_counter_vec[count_idx][j] = tmp_off+1;
        }
        count_idx++;
      }
      break;
    case KondoGC:
    case Hubbard:
    case Kondo:
      //Transfer
      for (i = 0; i < X->Def.EDNTransfer / 2; i++) {
        for (ihermite = 0; ihermite < 2; ihermite++) {
          idx = 2 * i + ihermite;

          isite1 = X->Def.EDGeneralTransfer[idx][0] + 1;
          isite2 = X->Def.EDGeneralTransfer[idx][2] + 1;
          sigma1 = X->Def.EDGeneralTransfer[idx][1];
          sigma2 = X->Def.EDGeneralTransfer[idx][3];

          if (general_hopp_GetInfo(X, isite1, isite2, sigma1, sigma2) != 0) {
            return -1;
          }
          tmp_trans = -X->Def.EDParaGeneralTransfer[idx];

#pragma omp parallel for default(none) firstprivate(X, tmp_trans, count_idx) private(j, tmp_off, dmv) shared(Ham, lui_counter_vec, list_1)
          for (j = 1; j <= X->Large.i_max; j++) {
            dmv = tmp_trans *
                  child_CisAjt(list_1[j], X, X->Large.is1_spin, X->Large.is2_spin, X->Large.isA_spin, X->Large.A_spin,
                           &tmp_off);
            Ham[count_idx][j] += dmv;
            lui_counter_vec[count_idx][j] = tmp_off;
          }
          count_idx++;
        }
      }

      //InterAll
      for (i = 0; i < X->Def.NInterAll_OffDiagonal / 2; i++) {
        for (ihermite = 0; ihermite < 2; ihermite++) {
          idx = 2 * i + ihermite;
          isite1 = X->Def.InterAll_OffDiagonal[idx][0] + 1;
          isite2 = X->Def.InterAll_OffDiagonal[idx][2] + 1;
          isite3 = X->Def.InterAll_OffDiagonal[idx][4] + 1;
          isite4 = X->Def.InterAll_OffDiagonal[idx][6] + 1;
          sigma1 = X->Def.InterAll_OffDiagonal[idx][1];
          sigma2 = X->Def.InterAll_OffDiagonal[idx][3];
          sigma3 = X->Def.InterAll_OffDiagonal[idx][5];
          sigma4 = X->Def.InterAll_OffDiagonal[idx][7];
          tmp_V = X->Def.ParaInterAll_OffDiagonal[idx];
          if (isite1 == 1 && sigma1 == 0 && isite2 == 4 && sigma2 == 0 && isite3 == 17 && sigma3 == 0 && isite4 == 19 &&
              sigma4 == 0) {
            tmp_V = tmp_V * 1.0;
          }
//  fprintf(stdoutMPI, "Debug: %d, %d, %d, %d, %d, %d, %d, %d\n ", isite1, sigma1,isite2, sigma2,isite3, sigma3,isite4, sigma4);
          general_int_GetInfo(
                  i,
                  X,
                  isite1,
                  isite2,
                  isite3,
                  isite4,
                  sigma1,
                  sigma2,
                  sigma3,
                  sigma4,
                  tmp_V
          );

          i_max = X->Large.i_max;
          isite1 = X->Large.is1_spin;
          isite2 = X->Large.is2_spin;
          Asum = X->Large.isA_spin;
          Adiff = X->Large.A_spin;

          isite3 = X->Large.is3_spin;
          isite4 = X->Large.is4_spin;
          Bsum = X->Large.isB_spin;
          Bdiff = X->Large.B_spin;

          tmp_V = X->Large.tmp_V;

          if (isite1 == isite2 && isite3 == isite4) { //diagonal
            for (j = 1; j <= i_max; j++) {
              dmv = CisAisCisAis_element(j, isite1, isite3, tmp_V, v0, v1, X, &tmp_off);
              Ham[0][j] += dmv;
            }
          }
          else {
              if (isite1 == isite2 && isite3 != isite4) {
#pragma omp parallel for default(none) firstprivate(X, tmp_V, isite1, isite3, isite4, Bsum, Bdiff, count_idx, i_max) private(j, tmp_off, dmv) shared(Ham, lui_counter_vec, list_1, v0, v1)
                  for (j = 1; j <= i_max; j++) {
                      dmv = CisAisCjtAku_element(j, isite1, isite3, isite4, Bsum, Bdiff, tmp_V, v0, v1, X, &tmp_off);
                      Ham[count_idx][j] += dmv;
                      lui_counter_vec[count_idx][j] = tmp_off;
                  }
              } else if (isite1 != isite2 && isite3 == isite4) {
#pragma omp parallel for default(none) firstprivate(X, isite1, isite2, isite3, isite4, Asum, Adiff,tmp_V, count_idx, i_max) private(j,  tmp_off, dmv) shared(Ham, lui_counter_vec, list_1, v0, v1)
                  for (j = 1; j <= i_max; j++) {
                      dmv = CisAjtCkuAku_element(j, isite1, isite2, isite3, Asum, Adiff, tmp_V, v0, v1, X, &tmp_off);
                      Ham[count_idx][j] += dmv;
                      lui_counter_vec[count_idx][j] = tmp_off;
                  }
              } else if (isite1 != isite2 && isite3 != isite4) {
#pragma omp parallel for default(none) firstprivate(X, tmp_V, count_idx, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff, i_max) private(j, tmp_off, tmp_off_2, dmv) shared(Ham, lui_counter_vec, list_1, v0, v1)
                  for (j = 1; j <= i_max; j++) {
                      dmv = CisAjtCkuAlv_element(j, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff, tmp_V, v0,
                                                 v1, X, &tmp_off_2);
                      Ham[count_idx][j] += dmv;
                      lui_counter_vec[count_idx][j] = tmp_off_2;
                  }
              }
              count_idx++;
          }
        }
      }

      //Pair hopping
      for (i = 0; i < X->Def.NPairHopping / 2; i++) {
        for (ihermite = 0; ihermite < 2; ihermite++) {
          idx = 2 * i + ihermite;
          pairhopp_GetInfo(idx, X);
#pragma omp parallel for default(none) firstprivate(X, count_idx) private(j, tmp_off, dmv) shared(Ham, lui_counter_vec, list_1, v0, v1)
            for (j = 1; j <= X->Large.i_max; j++) {
            dmv = pairhopp_element(j, v0, v1, X, &tmp_off);
            Ham[count_idx][j] += dmv;
            lui_counter_vec[count_idx][j] = tmp_off;
          }
          count_idx++;
        }
      }
      //Exchange
      for (i = 0; i < X->Def.NExchangeCoupling; i++) {
          exchange_GetInfo(i, X);
#pragma omp parallel for default(none) firstprivate(X, count_idx) private(j, tmp_off, dmv) shared(Ham, lui_counter_vec, v0, v1)
          for (j = 1; j <= X->Large.i_max; j++) {
              dmv = exchange_element(j, v0, v1, X, &tmp_off);
              Ham[count_idx][j] += dmv;
              lui_counter_vec[count_idx][j] = tmp_off;
          }
          count_idx++;
      }
      break;

    case SpinGC:
      if (X->Def.iFlgGeneralSpin == FALSE) {
        //Transfer
        for (i = 0; i < X->Def.EDNTransfer / 2; i++) {
          for (ihermite = 0; ihermite < 2; ihermite++) {
            idx = 2 * i + ihermite;
            isite1 = X->Def.EDGeneralTransfer[idx][0] + 1;
            isite2 = X->Def.EDGeneralTransfer[idx][2] + 1;
            sigma1 = X->Def.EDGeneralTransfer[idx][1];
            sigma2 = X->Def.EDGeneralTransfer[idx][3];
            tmp_trans = -X->Def.EDParaGeneralTransfer[idx];

            if (general_hopp_GetInfo(X, isite1, isite2, sigma1, sigma2) != 0) {
              return -1;
            }

            if (isite1 == isite2) {
              is1_spin = X->Def.Tpow[isite1 - 1];
              if (sigma1 == sigma2) {
                // longitudinal magnetic field
#pragma omp parallel for default(none) firstprivate(X, tmp_trans, is1_spin, sigma1, i_max) private(j) shared(Ham)
                  for (j = 1; j <= i_max; j++) {
                  Ham[0][j] += tmp_trans * child_Spin_CisAis(j, X, is1_spin, sigma1);
                }
              } else {
                // transverse magnetic field
                is1_spin = X->Def.Tpow[isite1 - 1];
#pragma omp parallel for default(none) firstprivate(X, tmp_trans, is1_spin, sigma2, i_max, count_idx) private(j, off, dmv) shared(Ham, lui_counter_vec)
                  for (j = 1; j <= i_max; j++) {
                    dmv = tmp_trans * child_SpinGC_CisAit(j, X, is1_spin, sigma2, &off);
                    Ham[count_idx][j] += dmv;
                    lui_counter_vec[count_idx][j] = off+1;
                }
                count_idx++;
              }
            } else {
              // hopping is not allowed in localized spin system
              return -1;
            }
          }
        }

        //InterAll
        for (i = 0; i < X->Def.NInterAll_OffDiagonal / 2; i++) {
            for (ihermite = 0; ihermite < 2; ihermite++) {
                idx = 2 * i + ihermite;
                isite1 = X->Def.InterAll_OffDiagonal[idx][0] + 1;
                isite2 = X->Def.InterAll_OffDiagonal[idx][4] + 1;
                sigma1 = X->Def.InterAll_OffDiagonal[idx][1];
                sigma2 = X->Def.InterAll_OffDiagonal[idx][3];
                sigma3 = X->Def.InterAll_OffDiagonal[idx][5];
                sigma4 = X->Def.InterAll_OffDiagonal[idx][7];
                tmp_V = X->Def.ParaInterAll_OffDiagonal[idx];

                general_int_spin_GetInfo(X, isite1, isite2, sigma1, sigma2, sigma3, sigma4, tmp_V);
                isA_up = X->Def.Tpow[isite1 - 1];
                isB_up = X->Def.Tpow[isite2 - 1];

                if (sigma1 == sigma2 && sigma3 == sigma4) { //diagonal
#pragma omp parallel for default(none) firstprivate(X, isA_up, isB_up, sigma2, sigma4, tmp_V, i_max) private(j, tmp_off, dmv) shared(Ham, v0, v1)
                    for (j = 1; j <= i_max; j++) {
                        dmv = GC_CisAisCisAis_spin_element(j, isA_up, isB_up, sigma2, sigma4, tmp_V, v0, v1, X);
                        Ham[0][j] += dmv;
                    }
                } else { //off-diagonal
                    if (sigma1 == sigma2 && sigma3 != sigma4) {
#pragma omp parallel for default(none) firstprivate(X, sigma2, sigma4, isA_up, isB_up, tmp_V, i_max, count_idx) private(j, tmp_off, dmv) shared(Ham, lui_counter_vec, v0, v1)
                        for (j = 1; j <= i_max; j++) {
                            dmv = GC_CisAisCitAiu_spin_element(j, sigma2, sigma4, isA_up, isB_up, tmp_V, v0, v1, X,
                                                               &tmp_off);
                            Ham[count_idx][j] += dmv;
                            lui_counter_vec[count_idx][j] = tmp_off + 1;
                        }
                    } else if (sigma1 != sigma2 && sigma3 == sigma4) {
#pragma omp parallel for default(none) firstprivate(X, sigma2, sigma4, isA_up, isB_up, tmp_V, i_max, count_idx) private(j, tmp_off, dmv) shared(Ham, lui_counter_vec, v0, v1)
                        for (j = 1; j <= i_max; j++) {
                            dmv = GC_CisAitCiuAiu_spin_element(j, sigma2, sigma4, isA_up, isB_up, tmp_V, v0, v1, X,
                                                               &tmp_off);
                            Ham[count_idx][j] += dmv;
                            lui_counter_vec[count_idx][j] = tmp_off + 1;
                        }
                    } else if (sigma1 != sigma2 && sigma3 != sigma4) {
#pragma omp parallel for default(none) firstprivate(X, sigma2, sigma4, isA_up, isB_up, tmp_V, i_max, count_idx) private(j, tmp_off_2, dmv) shared(Ham, lui_counter_vec, v0, v1)
                        for (j = 1; j <= i_max; j++) {
                            dmv = GC_CisAitCiuAiv_spin_element(j, sigma2, sigma4, isA_up, isB_up, tmp_V, v0, v1, X,
                                                               &tmp_off_2);
                            Ham[count_idx][j] += dmv;
                            lui_counter_vec[count_idx][j] = tmp_off_2 + 1;
                        }
                    }
                    count_idx++;
                }
            }
        }
        //Exchange
        for (i = 0; i < X->Def.NExchangeCoupling; i++) {
            exchange_spin_GetInfo(i, X);
#pragma omp parallel for default(none) firstprivate(X,count_idx) private(j, tmp_off, dmv) shared(Ham, lui_counter_vec, v0, v1)
            for (j = 1; j <= X->Large.i_max; j++) {
            dmv = GC_exchange_spin_element(j, v0, v1, X, &tmp_off);
            Ham[count_idx][j] += dmv;
            lui_counter_vec[count_idx][j] = tmp_off+1;
          }
          count_idx++;
        }

        //PairLift
        for (i = 0; i < X->Def.NPairLiftCoupling / 2; i++) {
          for (ihermite = 0; ihermite < 2; ihermite++) {
            idx = 2 * i + ihermite;
            pairlift_spin_GetInfo(idx, X);
#pragma omp parallel for default(none) firstprivate(X,count_idx) private(j, tmp_off, dmv) shared(Ham, lui_counter_vec, v0, v1)
            for (j = 1; j <= X->Large.i_max; j++) {
              dmv = GC_pairlift_spin_element(j, v0, v1, X, &tmp_off);
              Ham[count_idx][j] += dmv;
              lui_counter_vec[count_idx][j] = tmp_off+1;
            }
            count_idx++;
          }
        }
      } else { //For General spin
        for (i = 0; i < X->Def.EDNTransfer / 2; i++) {
          for (ihermite = 0; ihermite < 2; ihermite++) {
            idx = 2 * i + ihermite;
            isite1 = X->Def.EDGeneralTransfer[idx][0] + 1;
            isite2 = X->Def.EDGeneralTransfer[idx][2] + 1;
            sigma1 = X->Def.EDGeneralTransfer[idx][1];
            sigma2 = X->Def.EDGeneralTransfer[idx][3];
            tmp_trans = -X->Def.EDParaGeneralTransfer[idx];

            if (isite1 == isite2) {
              // longitudinal magnetic field is absorbed in diagonal calculation.
              // transverse magnetic field
#pragma omp parallel for default(none) firstprivate(X, isite1, sigma2, sigma1, tmp_trans, count_idx, i_max) private(j, off, num1) shared(Ham, lui_counter_vec, v0, v1)
              for (j = 1; j <= i_max; j++) {
                num1 = GetOffCompGeneralSpin(j - 1, isite1, sigma2, sigma1, &off, X->Def.SiteToBit, X->Def.Tpow);
                Ham[count_idx][j] += tmp_trans * num1;
                lui_counter_vec[count_idx][j] = off+1;
              }
              count_idx++;
            } else {
              // hopping is not allowed in localized spin system
              return -1;
            }
          }
        }

        //InterAll
        for (i = 0; i < X->Def.NInterAll_OffDiagonal / 2; i++) {
          for (ihermite = 0; ihermite < 2; ihermite++) {
            idx = 2 * i + ihermite;
            isite1 = X->Def.InterAll_OffDiagonal[idx][0] + 1;
            isite2 = X->Def.InterAll_OffDiagonal[idx][4] + 1;
            sigma1 = X->Def.InterAll_OffDiagonal[idx][1];
            sigma2 = X->Def.InterAll_OffDiagonal[idx][3];
            sigma3 = X->Def.InterAll_OffDiagonal[idx][5];
            sigma4 = X->Def.InterAll_OffDiagonal[idx][7];
            tmp_V = X->Def.ParaInterAll_OffDiagonal[idx];
#pragma omp parallel for default(none) firstprivate(X, isite1, isite2, sigma4, sigma3, sigma2, sigma1, tmp_V, i_max, count_idx) private(j, tmp_off, off, num1) shared(Ham, lui_counter_vec)
              for (j = 1; j <= i_max; j++) {
              num1 = GetOffCompGeneralSpin(j - 1, isite1, sigma2, sigma1, &tmp_off, X->Def.SiteToBit, X->Def.Tpow);
              if (num1 != 0) {
                num1 = GetOffCompGeneralSpin(tmp_off, isite2, sigma4, sigma3, &off, X->Def.SiteToBit, X->Def.Tpow);
                if (num1 != 0) {
                    Ham[count_idx][j] += tmp_V * num1;
                    lui_counter_vec[count_idx][j] = off+1;
                }
              }
            }
            count_idx++;
          }
        }
      }
      break;

    case Spin:
      if (X->Def.iFlgGeneralSpin == FALSE) {
        //Transfer is abosrbed in diagonal term.
        //InterAll
        for (i = 0; i < X->Def.NInterAll_OffDiagonal / 2; i++) {
          for (ihermite = 0; ihermite < 2; ihermite++) {
            idx = 2 * i + ihermite;
            isite1 = X->Def.InterAll_OffDiagonal[idx][0] + 1;
            isite2 = X->Def.InterAll_OffDiagonal[idx][4] + 1;
            sigma1 = X->Def.InterAll_OffDiagonal[idx][1];
            sigma2 = X->Def.InterAll_OffDiagonal[idx][3];
            sigma3 = X->Def.InterAll_OffDiagonal[idx][5];
            sigma4 = X->Def.InterAll_OffDiagonal[idx][7];
            tmp_V = X->Def.ParaInterAll_OffDiagonal[idx];

            general_int_spin_GetInfo(X, isite1, isite2, sigma1, sigma2, sigma3, sigma4, tmp_V);
            isA_up = X->Large.is1_up;
            isB_up = X->Large.is2_up;

#pragma omp parallel for default(none) firstprivate(X, isA_up, isB_up, sigma2, sigma4, tmp_V, i_max,count_idx) private(j, tmp_sgn, tmp_off, dmv) shared(Ham, lui_counter_vec)
            for (j = 1; j <= i_max; j++) {
              tmp_sgn = child_exchange_spin_element(j, X, isA_up, isB_up, sigma2, sigma4, &tmp_off);
              dmv = tmp_sgn * tmp_V;
              Ham[count_idx][j] += dmv;
              lui_counter_vec[count_idx][j] = tmp_off;
            }
            count_idx++;
          }
        }

        //Exchange
        for (i = 0; i < X->Def.NExchangeCoupling; i++) {
          exchange_spin_GetInfo(i, X);
#pragma omp parallel for default(none) firstprivate(X, count_idx, i_max) private(j, tmp_off, dmv) shared(Ham, lui_counter_vec, v0, v1)
            for (j = 1; j <= i_max; j++) {
            dmv = exchange_spin_element(j, v0, v1, X, &tmp_off);
            Ham[count_idx][j] += dmv;
            lui_counter_vec[count_idx][j] = tmp_off;
          }
          count_idx++;
        }

      } else { //For General spin
        //Transfer absorbed in Diagonal term.

        //InterAll
        for (i = 0; i < X->Def.NInterAll_OffDiagonal / 2; i++) {
          for (ihermite = 0; ihermite < 2; ihermite++) {
            idx = 2 * i + ihermite;
            isite1 = X->Def.InterAll_OffDiagonal[idx][0] + 1;
            isite2 = X->Def.InterAll_OffDiagonal[idx][4] + 1;
            sigma1 = X->Def.InterAll_OffDiagonal[idx][1];
            sigma2 = X->Def.InterAll_OffDiagonal[idx][3];
            sigma3 = X->Def.InterAll_OffDiagonal[idx][5];
            sigma4 = X->Def.InterAll_OffDiagonal[idx][7];
            tmp_V = X->Def.ParaInterAll_OffDiagonal[idx];

#pragma omp parallel for default(none) firstprivate(X, i_max, isite1, sigma2, sigma1, isite2, sigma4, sigma3, tmp_V, count_idx) private(j, tmp_off, off, num1) shared(Ham, lui_counter_vec, list_1)
              for (j = 1; j <= i_max; j++) {
              num1 = GetOffCompGeneralSpin(list_1[j], isite1, sigma2, sigma1, &tmp_off, X->Def.SiteToBit, X->Def.Tpow);
              if (num1 != 0) {
                num1 = GetOffCompGeneralSpin(tmp_off, isite2, sigma4, sigma3, &off, X->Def.SiteToBit, X->Def.Tpow);
                if (num1 != 0) {
                  ConvertToList1GeneralSpin(off, X->Check.sdim, &tmp_off);
                  Ham[count_idx][j] += tmp_V*num1;
                  lui_counter_vec[count_idx][j] = tmp_off;
                }
              }
            }
            count_idx++;
          }
        }
      }
      break;
  }

  return 0;
}
