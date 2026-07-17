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
#include "makeHam.h"
#include "wrapperMPI.h"
#include "nbody_interall.h"
#include "anomalous_pair.h"
#include "hamstore.h"

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
int makeHam(struct BindStruct *X) {

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
  /* Owned-column bounds for pure column-scatter loops below (design doc
     sec. 3 phase 2): in replicated mode (iHamPanelActive == 0) these
     always evaluate to [1, i_max], so the loop range is unchanged and
     the replicated path stays bit-identical. */
  long int hs_jb, hs_je;

  if (GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit) != 0) {
    return -1;
  }
  X->Large.i_max = i_max;
  X->Large.irght = irght;
  X->Large.ilft = ilft;
  X->Large.ihfbit = ihfbit;
  X->Large.prdct = 0.0;
  X->Large.mode = M_Ham;

  if (!iHamPanelActive) {
    for (i = 0; i <= i_max; i++) {
      for (j = 0; j <= i_max; j++) {
        Ham[i][j] = 0;
      }
    }
  } else {
    /* Ham_local is zero-initialized at allocation (xsetmem), but that
       only covers the first call. Zero it here too so a hypothetical
       second makeHam() call (e.g. future re-entrant use) cannot
       accumulate onto stale values from a previous call. */
    if (HamColEnd >= HamColBegin) {
      long int hs_ncols = HamColEnd - HamColBegin + 1;
      long int hs_nelem = HamPanelLd * hs_ncols;
      long int hs_k;
      for (hs_k = 0; hs_k < hs_nelem; hs_k++) Ham_local[hs_k] = 0;
    }
  }
#pragma omp parallel for default(none) firstprivate(i_max) private(j) shared(Ham, list_Diagonal, v0, v1, iHamPanelActive, HamColBegin, HamColEnd, HamPanelLd, Ham_local)
  for (j = 1; j <= i_max; j++) {
    if (HAM_OWNED_COL(j)) AddHamElem(j, j, list_Diagonal[j]);
    v0[j] = 1.0;
    v1[j] = 1.0;
    //printf("%ld, %f\n", j, list_Diagonal[j]);
  }
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

          hs_jb = iHamPanelActive ? HamColBegin : 1;
          hs_je = iHamPanelActive ? HamColEnd : (long int)X->Large.i_max;
          for (j = hs_jb; j <= hs_je; j++) {
            dmv = tmp_trans *
                  GC_CisAjt(j, v0, v1, X, X->Large.is1_spin, X->Large.is2_spin, X->Large.isA_spin, X->Large.A_spin,
                            tmp_trans, &tmp_off);
            AddHamElem(tmp_off + 1, j, dmv);
          }
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

          hs_jb = iHamPanelActive ? HamColBegin : 1;
          hs_je = iHamPanelActive ? HamColEnd : (long int)i_max;
          if (isite1 == isite2 && isite3 == isite4) {

            for (j = hs_jb; j <= hs_je; j++) {
              dmv = GC_CisAisCisAis_element(j, isite1, isite3, tmp_V, v0, v1, X, &tmp_off);
              AddHamElem(j, j, dmv);
            }
          } else if (isite1 == isite2 && isite3 != isite4) {

            for (j = hs_jb; j <= hs_je; j++) {
              dmv = GC_CisAisCjtAku_element(j, isite1, isite3, isite4, Bsum, Bdiff, tmp_V, v0, v1, X, &tmp_off);
              AddHamElem(tmp_off + 1, j, dmv);
            }
          } else if (isite1 != isite2 && isite3 == isite4) {

            for (j = hs_jb; j <= hs_je; j++) {
              dmv = GC_CisAjtCkuAku_element(j, isite1, isite2, isite3, Asum, Adiff, tmp_V, v0, v1, X, &tmp_off);
              AddHamElem(tmp_off + 1, j, dmv);
            }
          } else if (isite1 != isite2 && isite3 != isite4) {

            for (j = hs_jb; j <= hs_je; j++) {
              dmv = GC_CisAjtCkuAlv_element(j, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff, tmp_V,
                                                  v0, v1, X, &tmp_off_2);
              AddHamElem(tmp_off_2 + 1, j, dmv);
            }
          }
        }
      }
      if (X->Def.NNBodyInterAll_OffDiagonal > 0) {
        if (AddNBodyInterAllToHamHubbardGC(X) != 0) {
          return -1;
        }
      }
      if (X->Def.NAnomalousTerm > 0) {
        if (AddAnomalousTermToHamHubbardGC(X) != 0) {
          return -1;
        }
      }
      //Pair hopping
      for (i = 0; i < X->Def.NPairHopping / 2; i++) {
        for (ihermite = 0; ihermite < 2; ihermite++) {
          idx = 2 * i + ihermite;
          pairhopp_GetInfo(idx, X);
          hs_jb = iHamPanelActive ? HamColBegin : 1;
          hs_je = iHamPanelActive ? HamColEnd : (long int)X->Large.i_max;
          for (j = hs_jb; j <= hs_je; j++) {
            dmv = GC_pairhopp_element(j, v0, v1, X, &tmp_off);
            AddHamElem(tmp_off + 1, j, dmv);
          }
        }
      }
      //Exchange
      for (i = 0; i < X->Def.NExchangeCoupling; i++) {
        exchange_GetInfo(i, X);
        hs_jb = iHamPanelActive ? HamColBegin : 1;
        hs_je = iHamPanelActive ? HamColEnd : (long int)X->Large.i_max;
        for (j = hs_jb; j <= hs_je; j++) {
          dmv = GC_exchange_element(j, v0, v1, X, &tmp_off);
          AddHamElem(tmp_off + 1, j, dmv);
        }
      }
      break;
    case Hubbard:
    case tJ:
    case tJGC:
    case Kondo:
    case KondoGC:
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

          hs_jb = iHamPanelActive ? HamColBegin : 1;
          hs_je = iHamPanelActive ? HamColEnd : (long int)X->Large.i_max;
          for (j = hs_jb; j <= hs_je; j++) {
            dmv = tmp_trans *
                  child_CisAjt(list_1[j], X, X->Large.is1_spin, X->Large.is2_spin, X->Large.isA_spin, X->Large.A_spin,
                           &tmp_off);
            /* child_CisAjt sets *tmp_off=0 (and returns sgn 0 -> dmv 0) on an
               annihilated hop, and a valid 1-based row (>=1) with dmv!=0 only
               on a surviving hop. Skipping tmp_off==0 avoids the row-0 panel
               underrun and is numerically identical to the replicated path,
               which added dmv==0 into the unused Ham[0][j]. */
            if (tmp_off > 0) AddHamElem(tmp_off, j, dmv);
          }
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

          hs_jb = iHamPanelActive ? HamColBegin : 1;
          hs_je = iHamPanelActive ? HamColEnd : (long int)i_max;
          if (isite1 == isite2 && isite3 == isite4) {

            for (j = hs_jb; j <= hs_je; j++) {
              dmv = CisAisCisAis_element(j, isite1, isite3, tmp_V, v0, v1, X, &tmp_off);
              AddHamElem(j, j, dmv);
            }
          } else if (isite1 == isite2 && isite3 != isite4) {

            for (j = hs_jb; j <= hs_je; j++) {
              dmv = CisAisCjtAku_element(j, isite1, isite3, isite4, Bsum, Bdiff, tmp_V, v0, v1, X, &tmp_off);
              /* Element helper returns dam_pr==0 on every annihilated path;
                 dmv!=0 implies its internal child_CisAjt (which couples
                 *tmp_off=0 with return 0) survived, so tmp_off is a valid
                 1-based row. Guarding on dmv is exact (adding 0 is a no-op)
                 and robust to a stale tmp_off left by a dead gate. */
              if (dmv != 0.0) AddHamElem(tmp_off, j, dmv);
            }
          } else if (isite1 != isite2 && isite3 == isite4) {

            for (j = hs_jb; j <= hs_je; j++) {
              dmv = CisAjtCkuAku_element(j, isite1, isite2, isite3, Asum, Adiff, tmp_V, v0, v1, X, &tmp_off);
              /* child_CisAis gate runs BEFORE child_CisAjt here, so a failed
                 gate leaves tmp_off STALE (untouched) -> tmp_off could be 0
                 (row-0 underrun) or a wrong prior row. dmv==0 on every dead
                 path, so guard on dmv (exact: skipped write was +=0). */
              if (dmv != 0.0) AddHamElem(tmp_off, j, dmv);
            }
          } else if (isite1 != isite2 && isite3 != isite4) {

            for (j = hs_jb; j <= hs_je; j++) {
              dmv = CisAjtCkuAlv_element(j, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff, tmp_V, v0,
                                               v1, X, &tmp_off_2);
              /* Same stale-out-param hazard: a failed child_GC_CisAjt
                 intermediate leaves tmp_off_2 untouched. dmv==0 on all dead
                 paths; guard on dmv (exact, robust to stale tmp_off_2). */
              if (dmv != 0.0) AddHamElem(tmp_off_2, j, dmv);
            }
          }
        }
      }

      if ((X->Def.iCalcModel == Hubbard ||
           X->Def.iCalcModel == tJ ||
           X->Def.iCalcModel == tJGC ||
           X->Def.iCalcModel == Kondo ||
           X->Def.iCalcModel == KondoGC) &&
          X->Def.NNBodyInterAll_OffDiagonal > 0) {
        if (AddNBodyInterAllToHamHubbard(X) != 0) {
          return -1;
        }
      }

      //Pair hopping
      for (i = 0; i < X->Def.NPairHopping / 2; i++) {
        for (ihermite = 0; ihermite < 2; ihermite++) {
          idx = 2 * i + ihermite;
          pairhopp_GetInfo(idx, X);
          hs_jb = iHamPanelActive ? HamColBegin : 1;
          hs_je = iHamPanelActive ? HamColEnd : (long int)X->Large.i_max;
          for (j = hs_jb; j <= hs_je; j++) {
            dmv = pairhopp_element(j, v0, v1, X, &tmp_off);
            /* pairhopp_element leaves tmp_off STALE on its dead branch (and
               returns 0 without setting it on GetOffComp failure). dmv==0 on
               every no-op; guard on dmv (exact, robust to stale tmp_off). */
            if (dmv != 0.0) AddHamElem(tmp_off, j, dmv);
          }
        }
      }
      //Exchange
      for (i = 0; i < X->Def.NExchangeCoupling; i++) {
        exchange_GetInfo(i, X);
        hs_jb = iHamPanelActive ? HamColBegin : 1;
        hs_je = iHamPanelActive ? HamColEnd : (long int)X->Large.i_max;
        for (j = hs_jb; j <= hs_je; j++) {
          dmv = exchange_element(j, v0, v1, X, &tmp_off);
          /* exchange_element leaves tmp_off STALE on its dead branch. dmv==0
             on every no-op; guard on dmv (exact, robust to stale tmp_off). */
          if (dmv != 0.0) AddHamElem(tmp_off, j, dmv);
        }
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
              hs_jb = iHamPanelActive ? HamColBegin : 1;
              hs_je = iHamPanelActive ? HamColEnd : (long int)i_max;
              if (sigma1 == sigma2) {
                // longitudinal magnetic field
                for (j = hs_jb; j <= hs_je; j++) {
                  AddHamElem(j, j, tmp_trans * child_Spin_CisAis(j, X, is1_spin, sigma1));
                }
              } else {
                // transverse magnetic field
                is1_spin = X->Def.Tpow[isite1 - 1];

                for (j = hs_jb; j <= hs_je; j++) {
                  /* child_SpinGC_CisAit writes *off as an out-param; keep the
                     function call in its own statement (as the rest of this
                     file does via tmp_off) so off is fully updated before
                     AddHamElem reads it for the row index — evaluating both
                     in one expression would make the row index depend on
                     unspecified evaluation order between off's read (for the
                     index) and its write (inside the call). */
                  dmv = tmp_trans * child_SpinGC_CisAit(j, X, is1_spin, sigma2, &off);
                  AddHamElem(off + 1, j, dmv);
                }
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

            hs_jb = iHamPanelActive ? HamColBegin : 1;
            hs_je = iHamPanelActive ? HamColEnd : (long int)i_max;
            if (sigma1 == sigma2 && sigma3 == sigma4) { //diagonal
              for (j = hs_jb; j <= hs_je; j++) {
                dmv = GC_CisAisCisAis_spin_element(j, isA_up, isB_up, sigma2, sigma4, tmp_V, v0, v1, X);
                AddHamElem(j, j, dmv);
              }
            } else if (sigma1 == sigma2 && sigma3 != sigma4) {
              for (j = hs_jb; j <= hs_je; j++) {
                dmv = GC_CisAisCitAiu_spin_element(j, sigma2, sigma4, isA_up, isB_up, tmp_V, v0, v1, X, &tmp_off);
                AddHamElem(tmp_off + 1, j, dmv);
              }
            } else if (sigma1 != sigma2 && sigma3 == sigma4) {
              for (j = hs_jb; j <= hs_je; j++) {
                dmv = GC_CisAitCiuAiu_spin_element(j, sigma2, sigma4, isA_up, isB_up, tmp_V, v0, v1, X, &tmp_off);
                AddHamElem(tmp_off + 1, j, dmv);
              }
            } else if (sigma1 != sigma2 && sigma3 != sigma4) {
              for (j = hs_jb; j <= hs_je; j++) {
                dmv = GC_CisAitCiuAiv_spin_element(j, sigma2, sigma4, isA_up, isB_up, tmp_V, v0, v1, X,
                                                         &tmp_off_2);
                AddHamElem(tmp_off_2 + 1, j, dmv);
              }
            }
          }
        }
        if (X->Def.NNBodyInterAll_OffDiagonal > 0) {
          if (AddNBodyInterAllToHamSpinGC(X) != 0) {
            return -1;
          }
        }
        //Exchange
        for (i = 0; i < X->Def.NExchangeCoupling; i++) {
          exchange_spin_GetInfo(i, X);
          hs_jb = iHamPanelActive ? HamColBegin : 1;
          hs_je = iHamPanelActive ? HamColEnd : (long int)X->Large.i_max;
          for (j = hs_jb; j <= hs_je; j++) {
            dmv = GC_exchange_spin_element(j, v0, v1, X, &tmp_off);
            AddHamElem(tmp_off + 1, j, dmv);
          }
        }

        //PairLift
        for (i = 0; i < X->Def.NPairLiftCoupling / 2; i++) {
          for (ihermite = 0; ihermite < 2; ihermite++) {
            idx = 2 * i + ihermite;
            pairlift_spin_GetInfo(idx, X);

            hs_jb = iHamPanelActive ? HamColBegin : 1;
            hs_je = iHamPanelActive ? HamColEnd : (long int)X->Large.i_max;
            for (j = hs_jb; j <= hs_je; j++) {
              dmv = GC_pairlift_spin_element(j, v0, v1, X, &tmp_off);
              AddHamElem(tmp_off + 1, j, dmv);
            }
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
              hs_jb = iHamPanelActive ? HamColBegin : 1;
              hs_je = iHamPanelActive ? HamColEnd : (long int)i_max;
              for (j = hs_jb; j <= hs_je; j++) {
                num1 = GetOffCompGeneralSpin(j - 1, isite1, sigma2, sigma1, &off, X->Def.SiteToBit, X->Def.Tpow);
                AddHamElem(off + 1, j, tmp_trans * num1);
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
            hs_jb = iHamPanelActive ? HamColBegin : 1;
            hs_je = iHamPanelActive ? HamColEnd : (long int)i_max;
            for (j = hs_jb; j <= hs_je; j++) {
              num1 = GetOffCompGeneralSpin(j - 1, isite1, sigma2, sigma1, &tmp_off, X->Def.SiteToBit, X->Def.Tpow);
              if (num1 != 0) {
                num1 = GetOffCompGeneralSpin(tmp_off, isite2, sigma4, sigma3, &off, X->Def.SiteToBit, X->Def.Tpow);
                if (num1 != 0) {
                  AddHamElem(off + 1, j, tmp_V * num1);
                }
              }
            }
          }
        }

        if (X->Def.NNBodyInterAll_OffDiagonal > 0) {
          if (AddNBodyInterAllToHamSpinGC(X) != 0) {
            return -1;
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

            hs_jb = iHamPanelActive ? HamColBegin : 1;
            hs_je = iHamPanelActive ? HamColEnd : (long int)i_max;
            for (j = hs_jb; j <= hs_je; j++) {
              tmp_sgn = child_exchange_spin_element(j, X, isA_up, isB_up, sigma2, sigma4, &tmp_off);
              dmv = tmp_sgn * tmp_V;
              /* child_exchange_spin_element sets *tmp_off=0 both on its dead
                 branch AND (via an unchecked GetOffComp) if the exchanged
                 state is out of the sector -- but in the latter case it still
                 returns tmp_sgn=1, so dmv!=0. A dmv-guard would therefore
                 write into row 0. Guard on tmp_off (row-0 skip matches the
                 replicated Ham[0][j], which the eigensolver never reads). */
              if (tmp_off > 0) AddHamElem(tmp_off, j, dmv);
            }
          }
        }

        if (X->Def.NNBodyInterAll_OffDiagonal > 0) {
          if (AddNBodyInterAllToHamSpinGC(X) != 0) {
            return -1;
          }
        }

        //Exchange
        for (i = 0; i < X->Def.NExchangeCoupling; i++) {
          exchange_spin_GetInfo(i, X);
          hs_jb = iHamPanelActive ? HamColBegin : 1;
          hs_je = iHamPanelActive ? HamColEnd : (long int)X->Large.i_max;
          for (j = hs_jb; j <= hs_je; j++) {
            dmv = exchange_spin_element(j, v0, v1, X, &tmp_off);
            /* exchange_spin_element leaves tmp_off STALE on its dead branch,
               and on its live branch sets *tmp_off via an unchecked GetOffComp
               (=0 if out of sector) while dmv!=0. So neither a stale value nor
               a dmv-guard is safe; guard on tmp_off (row-0 skip matches the
               replicated unused Ham[0][j]). */
            if (tmp_off > 0) AddHamElem(tmp_off, j, dmv);
          }
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

            hs_jb = iHamPanelActive ? HamColBegin : 1;
            hs_je = iHamPanelActive ? HamColEnd : (long int)i_max;
            for (j = hs_jb; j <= hs_je; j++) {
              num1 = GetOffCompGeneralSpin(list_1[j], isite1, sigma2, sigma1, &tmp_off, X->Def.SiteToBit, X->Def.Tpow);
              if (num1 != 0) {
                num1 = GetOffCompGeneralSpin(tmp_off, isite2, sigma4, sigma3, &off, X->Def.SiteToBit, X->Def.Tpow);
                if (num1 != 0) {
                  ConvertToList1GeneralSpin(off, X->Check.sdim, &tmp_off);
                  /* ConvertToList1GeneralSpin sets *tmp_off=0 when the state
                     is not in list_1. The amplitude here is the constant
                     tmp_V (no dmv signal), so guard on tmp_off; row-0 skip
                     matches the replicated unused Ham[0][j]. */
                  if (tmp_off > 0) AddHamElem(tmp_off, j, tmp_V);
                }
              }
            }
          }
        }

        if (X->Def.NNBodyInterAll_OffDiagonal > 0) {
          if (AddNBodyInterAllToHamSpinGC(X) != 0) {
            return -1;
          }
        }
      }

      break;
  }
  return 0;
}
