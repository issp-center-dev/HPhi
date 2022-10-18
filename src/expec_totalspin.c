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
#include "mltplySpinCore.h"
#include "wrapperMPI.h"
#include "mltplyMPISpin.h"
#include "mltplyMPISpinCore.h"
#include "expec_totalspin.h"
/**
 * @file   expec_totalspin.c
 *
 * @brief  File for calculating total spin
 *
 * @version 0.2
 * @details modify to treat the case of general spin
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 *
 */
/**
 * @brief function of calculating totalspin for Hubbard model
 *
 * @param[in,out] X data list of calculation parameters
 * @param vec eigenvector
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void totalspin_Hubbard(
  struct BindStruct *X,
  int nstate,
  double complex **vec
) {
  long unsigned int j;
  long unsigned int irght, ilft, ihfbit;
  long unsigned int isite1, isite2;
  long unsigned int is1_up, is2_up, is1_down, is2_down;
  long unsigned int iexchg, off;
  int num1_up, num2_up, istate;
  int num1_down, num2_down;
  long unsigned int ibit1_up, ibit2_up, ibit1_down, ibit2_down;
  double complex tmp_spn_z;
  long unsigned i_max;
  i_max = X->Check.idim_max;

  GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit);
  for (istate = 0; istate < nstate; istate++) {
    X->Phys.s2[istate] = 0.0;
    X->Phys.Sz[istate] = 0.0;
  }
  for (isite1 = 1; isite1 <= X->Def.NsiteMPI; isite1++) {
    is1_up = X->Def.Tpow[2 * isite1 - 2];
    is1_down = X->Def.Tpow[2 * isite1 - 1];

    for (isite2 = 1; isite2 <= X->Def.NsiteMPI; isite2++) {
      is2_up = X->Def.Tpow[2 * isite2 - 2];
      is2_down = X->Def.Tpow[2 * isite2 - 1];

      for (j = 1; j <= i_max; j++) {

        ibit1_up = list_1[j] & is1_up;
        num1_up = ibit1_up / is1_up;
        ibit2_up = list_1[j] & is2_up;
        num2_up = ibit2_up / is2_up;

        ibit2_down = list_1[j] & is2_down;
        num2_down = ibit2_down / is2_down;
        ibit1_down = list_1[j] & is1_down;
        num1_down = ibit1_down / is1_down;

        tmp_spn_z = (num1_up - num1_down) * (num2_up - num2_down);
        for (istate = 0; istate < nstate; istate++)
          X->Phys.s2[istate] += conj(vec[j][istate]) * vec[j][istate] * tmp_spn_z / 4.0;
        if (isite1 == isite2) {
          for (istate = 0; istate < nstate; istate++) {
            X->Phys.s2[istate] += conj(vec[j][istate]) * vec[j][istate] * (num1_up + num1_down - 2 * num1_up * num1_down) / 2.0;
            X->Phys.Sz[istate] += conj(vec[j][istate]) * vec[j][istate] * (num1_up - num1_down) / 2.0;
          }
        }
        else {
          if (ibit1_up != 0 && ibit1_down == 0 && ibit2_up == 0 && ibit2_down != 0) {
            iexchg = list_1[j] - (is1_up + is2_down);
            iexchg += (is2_up + is1_down);
            GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off);
            for (istate = 0; istate < nstate; istate++)
              X->Phys.s2[istate] += conj(vec[j][istate]) * vec[off][istate] / 2.0;
          }
          else if (ibit1_up == 0 && ibit1_down != 0 && ibit2_up != 0 && ibit2_down == 0) {
            iexchg = list_1[j] - (is1_down + is2_up);
            iexchg += (is2_down + is1_up);
            GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off);
            for (istate = 0; istate < nstate; istate++)
              X->Phys.s2[istate] += conj(vec[j][istate]) * vec[off][istate] / 2.0;
          }
        }
      }
    }
  }
  SumMPI_dv(nstate, X->Phys.s2);
  SumMPI_dv(nstate, X->Phys.Sz);
}
/**
 * @brief function of calculating totalspin for Hubbard model in grand canonical ensemble
 *
 * @param[in,out] X data list of calculation parameters
 * @param vec eigenvector
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void totalspin_HubbardGC(
  struct BindStruct *X,
  int nstate,
  double complex **vec
) {
  long unsigned int j;
  long unsigned int isite1, isite2;
  long unsigned int is1_up, is2_up, is1_down, is2_down;
  long unsigned int iexchg, off;
  int num1_up, num2_up, istate;
  int num1_down, num2_down;
  long unsigned int ibit1_up, ibit2_up, ibit1_down, ibit2_down, list_1_j;
  double complex tmp_spn_z;
  long unsigned int i_max;

  i_max = X->Check.idim_max;

  for (istate = 0; istate < nstate; istate++) {
    X->Phys.s2[istate] = 0.0;
    X->Phys.Sz[istate] = 0.0;
  }
  for (isite1 = 1; isite1 <= X->Def.NsiteMPI; isite1++) {
    for (isite2 = 1; isite2 <= X->Def.NsiteMPI; isite2++) {
      is1_up = X->Def.Tpow[2 * isite1 - 2];
      is1_down = X->Def.Tpow[2 * isite1 - 1];
      is2_up = X->Def.Tpow[2 * isite2 - 2];
      is2_down = X->Def.Tpow[2 * isite2 - 1];

      for (j = 1; j <= i_max; j++) {
        list_1_j = j - 1;
        ibit1_up = list_1_j & is1_up;
        num1_up = ibit1_up / is1_up;
        ibit2_up = list_1_j & is2_up;
        num2_up = ibit2_up / is2_up;

        ibit1_down = list_1_j & is1_down;
        num1_down = ibit1_down / is1_down;
        ibit2_down = list_1_j & is2_down;
        num2_down = ibit2_down / is2_down;

        tmp_spn_z = (num1_up - num1_down) * (num2_up - num2_down);
        for (istate = 0; istate < nstate; istate++)
          X->Phys.s2[istate] += conj(vec[j][istate]) * vec[j][istate] * tmp_spn_z / 4.0;
        if (isite1 == isite2) {
          X->Phys.s2[istate] += conj(vec[j][istate]) * vec[j][istate] * (num1_up + num1_down - 2 * num1_up * num1_down) / 2.0;
          X->Phys.Sz[istate] += conj(vec[j][istate]) * vec[j][istate] * (num1_up - num1_down) / 2.0;
        }
        else {
          if (ibit1_up != 0 && ibit1_down == 0 && ibit2_up == 0 && ibit2_down != 0) {
            iexchg = list_1_j - (is1_up + is2_down);
            iexchg += (is2_up + is1_down);
            off = iexchg + 1;
            for (istate = 0; istate < nstate; istate++)
              X->Phys.s2[istate] += conj(vec[j][istate]) * vec[off][istate] / 2.0;
          }
          else if (ibit1_up == 0 && ibit1_down != 0 && ibit2_up != 0 && ibit2_down == 0) {
            iexchg = list_1_j - (is1_down + is2_up);
            iexchg += (is2_down + is1_up);
            off = iexchg + 1;
            for (istate = 0; istate < nstate; istate++)
              X->Phys.s2[istate] += conj(vec[j][istate]) * vec[off][istate] / 2.0;
          }
        }
      }
    }
  }
  SumMPI_dv(nstate, X->Phys.s2);
  SumMPI_dv(nstate, X->Phys.Sz);
}
/**
 * @brief function of calculating totalspin for spin model
 *
 * @param[in,out] X data list of calculation parameters
 * @param vec eigenvector
 * @version 0.2
 * @details modify for hybrid parallel
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void totalspin_Spin(
  struct BindStruct *X, 
  int nstate,
  double complex **vec
) {
  long unsigned int j;
  long unsigned int irght, ilft, ihfbit;
  long unsigned int isite1, isite2;
  long unsigned int tmp_isite1, tmp_isite2;

  long unsigned int is1_up, is2_up;
  long unsigned int iexchg, off, off_2;

  int num1_up, num2_up;
  int num1_down, num2_down;
  int sigma_1, sigma_2, istate;
  long unsigned int ibit1_up, ibit2_up, ibit_tmp, is_up;
  double complex spn_z = 0.0;
  double complex spn_z1 = 0.0;
  double complex spn_z2 = 0.0;
  double complex spn_zd = 0.0;
  double complex spn = 0.0;
  long unsigned int i_max;

  i_max = X->Check.idim_max;
  for (istate = 0; istate < nstate; istate++) {
    X->Phys.s2[istate] = 0.0;
    X->Phys.Sz[istate] = 0.0;
  }
  if (X->Def.iFlgGeneralSpin == FALSE) {
    GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit);
    for (isite1 = 1; isite1 <= X->Def.NsiteMPI; isite1++) {
      for (isite2 = 1; isite2 <= X->Def.NsiteMPI; isite2++) {

        if (isite1 > X->Def.Nsite && isite2 > X->Def.Nsite) {
#ifdef MPI
          is1_up = X->Def.Tpow[isite1 - 1];
          is2_up = X->Def.Tpow[isite2 - 1];
          is_up = is1_up + is2_up;
          num1_up = child_SpinGC_CisAis((unsigned long int) myrank + 1, X, is1_up, 1);
          num1_down = 1 - num1_up;
          num2_up = child_SpinGC_CisAis((unsigned long int) myrank + 1, X, is2_up, 1);
          num2_down = 1 - num2_up;
          spn_z = (num1_up - num1_down) * (num2_up - num2_down);

          for (j = 1; j <= i_max; j++) {
            for (istate = 0; istate < nstate; istate++)
              X->Phys.s2[istate] += conj(vec[j][istate]) * vec[j][istate] * spn_z / 4.0;
          }
          if (isite1 == isite2) {
            for (j = 1; j <= i_max; j++) {
              for (istate = 0; istate < nstate; istate++)
                X->Phys.s2[istate] += conj(vec[j][istate]) * vec[j][istate] / 2.0;
            }
          } else {//off diagonal
            spn += child_general_int_spin_TotalS_MPIdouble(isite1 - 1, isite2 - 1, X, nstate, vec, vec);
          }
#endif
        }
        else if (isite1 > X->Def.Nsite || isite2 > X->Def.Nsite) {
#ifdef MPI
          if (isite1 < isite2) {
            tmp_isite1 = isite1;
            tmp_isite2 = isite2;
          }
          else {
            tmp_isite1 = isite2;
            tmp_isite2 = isite1;
          }

          is1_up = X->Def.Tpow[tmp_isite1 - 1];
          is2_up = X->Def.Tpow[tmp_isite2 - 1];
          num2_up = child_SpinGC_CisAis((unsigned long int) myrank + 1, X, is2_up, 1);
          num2_down = 1 - num2_up;

          //diagonal
          for (j = 1; j <= i_max; j++) {
            ibit1_up = list_1[j] & is1_up;
            num1_up = ibit1_up / is1_up;
            num1_down = 1 - num1_up;
            spn_z = (num1_up - num1_down) * (num2_up - num2_down);
            for (istate = 0; istate < nstate; istate++)
              X->Phys.s2[istate] += conj(vec[j][istate]) * vec[j][istate] * spn_z / 4.0;
          }
          if (isite1 < isite2) {
            spn += child_general_int_spin_MPIsingle(isite1 - 1, 0, 1, isite2 - 1, 1, 0, 1.0, X, nstate, vec, vec);
          } else {
            spn += conj(child_general_int_spin_MPIsingle(isite2 - 1, 1, 0, isite1 - 1, 0, 1, 1.0, X, nstate, vec, vec));
          }
#endif
        }//isite1 > Nsite || isite2 > Nsite
        else {
          is1_up = X->Def.Tpow[isite1 - 1];
          is2_up = X->Def.Tpow[isite2 - 1];
          is_up = is1_up + is2_up;

          for (j = 1; j <= i_max; j++) {
            ibit1_up = list_1[j] & is1_up;
            num1_up = ibit1_up / is1_up;
            num1_down = 1 - num1_up;
            ibit2_up = list_1[j] & is2_up;
            num2_up = ibit2_up / is2_up;
            num2_down = 1 - num2_up;

            spn_z = (num1_up - num1_down) * (num2_up - num2_down);
            for (istate = 0; istate < nstate; istate++)
              X->Phys.s2[istate] += conj(vec[j][istate]) * vec[j][istate] * spn_z / 4.0;

            if (isite1 == isite2) {
              for (istate = 0; istate < nstate; istate++)
                X->Phys.s2[istate] += conj(vec[j][istate]) * vec[j][istate] / 2.0;
            }
            else {
              ibit_tmp = (num1_up) ^ (num2_up);
              if (ibit_tmp != 0) {
                iexchg = list_1[j] ^ (is_up);
                GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off);
                for (istate = 0; istate < nstate; istate++)
                  X->Phys.s2[istate] += conj(vec[j][istate]) * vec[off][istate] / 2.0;
              }
            }
          }// j
        }
      }//isite2
    }//isite1
  }//generalspin=FALSE
  else {
    double S1 = 0;
    double S2 = 0;
    for (isite1 = 1; isite1 <= X->Def.NsiteMPI; isite1++) {
      for (isite2 = 1; isite2 <= X->Def.NsiteMPI; isite2++) {
        S1 = 0.5 * (X->Def.SiteToBit[isite1 - 1] - 1);
        S2 = 0.5 * (X->Def.SiteToBit[isite2 - 1] - 1);
        if (isite1 == isite2) {
          for (j = 1; j <= i_max; j++) {
            spn_z1 = 0.5 * GetLocal2Sz(isite1, list_1[j], X->Def.SiteToBit, X->Def.Tpow);
            for (istate = 0; istate < nstate; istate++) {
              X->Phys.s2[istate] += conj(vec[j][istate]) * vec[j][istate] * S1 * (S1 + 1.0);
              X->Phys.Sz[istate] += conj(vec[j][istate]) * vec[j][istate] * spn_z1;
            }
          }
        }
        else {
          for (j = 1; j <= i_max; j++) {
            spn_z1 = 0.5 * GetLocal2Sz(isite1, list_1[j], X->Def.SiteToBit, X->Def.Tpow);
            spn_z2 = 0.5 * GetLocal2Sz(isite2, list_1[j], X->Def.SiteToBit, X->Def.Tpow);
            for (istate = 0; istate < nstate; istate++)
              X->Phys.s2[istate] += conj(vec[j][istate]) * vec[j][istate] * spn_z1 * spn_z2;

            sigma_1 = GetBitGeneral(isite1, list_1[j], X->Def.SiteToBit, X->Def.Tpow);
            sigma_2 = GetBitGeneral(isite2, list_1[j], X->Def.SiteToBit, X->Def.Tpow);

            ibit_tmp = GetOffCompGeneralSpin(list_1[j], isite2, sigma_2, sigma_2 + 1, &off, X->Def.SiteToBit,
              X->Def.Tpow);
            if (ibit_tmp == TRUE) {
              ibit_tmp = GetOffCompGeneralSpin(off, isite1, sigma_1, sigma_1 - 1, &off_2, X->Def.SiteToBit,
                X->Def.Tpow);
              if (ibit_tmp == TRUE) {
                ConvertToList1GeneralSpin(off_2, X->Check.sdim, &off);
                for (istate = 0; istate < nstate; istate++)
                  X->Phys.s2[istate] += conj(vec[j][istate]) * vec[off][istate]
                  * sqrt(S2 * (S2 + 1) - spn_z2 * (spn_z2 + 1)) 
                  * sqrt(S1 * (S1 + 1) - spn_z1 * (spn_z1 - 1)) / 2.0;
              }
            }

            ibit_tmp = GetOffCompGeneralSpin(list_1[j], isite2, sigma_2, sigma_2 - 1, &off, X->Def.SiteToBit,
              X->Def.Tpow);
            if (ibit_tmp == TRUE) {
              ibit_tmp = GetOffCompGeneralSpin(off, isite1, sigma_1, sigma_1 + 1, &off_2, X->Def.SiteToBit,
                X->Def.Tpow);
              if (ibit_tmp == TRUE) {
                ConvertToList1GeneralSpin(off_2, X->Check.sdim, &off);
                for (istate = 0; istate < nstate; istate++)
                  X->Phys.s2[istate] += conj(vec[j][istate]) * vec[off][istate]
                  * sqrt(S2 * (S2 + 1) - spn_z2 * (spn_z2 - 1.0)) 
                  * sqrt(S1 * (S1 + 1) - spn_z1 * (spn_z1 + 1)) / 2.0;
              }
            }
          }
        }
      }
    }
  }
  SumMPI_dv(nstate, X->Phys.s2);
  SumMPI_dv(nstate, X->Phys.Sz);
}
/**
 * @brief function of calculating totalspin for spin model in grand canonical ensemble
 *
 * @param[in,out] X data list of calculation parameters
 * @param vec eigenvector
 * @version 0.2
 * @details add function to treat a calculation of total spin for general spin
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void totalspin_SpinGC(
  struct BindStruct *X,
  int nstate,
  double complex **vec
) {
  long unsigned int j;
  long unsigned int isite1, isite2, tmp_isite1, tmp_isite2;
  long unsigned int is1_up, is2_up;
  long unsigned int iexchg, off, off_2;
  int num1_up, num2_up, istate;
  int num1_down, num2_down;
  int sigma_1, sigma_2;
  long unsigned int ibit1_up, ibit2_up, ibit_tmp, is_up;
  double complex spn_z1, spn_z2;
  long unsigned int list_1_j;
  long unsigned int i_max;

  i_max = X->Check.idim_max;
  X->Large.mode = M_TOTALS;
  for (istate = 0; istate < nstate; istate++) {
    X->Phys.s2[istate] = 0.0;
    X->Phys.Sz[istate] = 0.0;
  }
  if (X->Def.iFlgGeneralSpin == FALSE) {
    for (isite1 = 1; isite1 <= X->Def.NsiteMPI; isite1++) {
      if (isite1 > X->Def.Nsite) {
        is1_up = X->Def.Tpow[isite1 - 1];
        ibit1_up = myrank & is1_up;
        num1_up = ibit1_up / is1_up;
        num1_down = 1 - num1_up;
        for (j = 1; j <= i_max; j++) {
          for (istate = 0; istate < nstate; istate++)
            X->Phys.Sz[istate] += conj(vec[j][istate])*vec[j][istate] * (num1_up - num1_down) / 2.0;
        }
      }
      else {
        is1_up = X->Def.Tpow[isite1 - 1];
        for (j = 1; j <= i_max; j++) {
          list_1_j = j - 1;
          ibit1_up = list_1_j & is1_up;
          num1_up = ibit1_up / is1_up;
          num1_down = 1 - num1_up;
          for (istate = 0; istate < nstate; istate++)
            X->Phys.Sz[istate] += conj(vec[j][istate])*vec[j][istate] * (num1_up - num1_down) / 2.0;
        }
      }
      for (isite2 = 1; isite2 <= X->Def.NsiteMPI; isite2++) {

        if (isite1 > X->Def.Nsite && isite2 > X->Def.Nsite) {
          is1_up = X->Def.Tpow[isite1 - 1];
          is2_up = X->Def.Tpow[isite2 - 1];
          num1_up = child_SpinGC_CisAis((unsigned long int)myrank + 1, X, is1_up, 1);
          num1_down = 1 - num1_up;
          num2_up = child_SpinGC_CisAis((unsigned long int)myrank + 1, X, is2_up, 1);
          num2_down = 1 - num2_up;
          spn_z2 = (num1_up - num1_down)*(num2_up - num2_down) / 4.0;
          for (j = 1; j <= i_max; j++) {
            for (istate = 0; istate < nstate; istate++)
              X->Phys.s2[istate] += conj(vec[j][istate])*vec[j][istate] * spn_z2;
          }
          if (isite1 == isite2) {
            for (j = 1; j <= i_max; j++) {
              for (istate = 0; istate < nstate; istate++)
                X->Phys.s2[istate] += conj(vec[j][istate])*vec[j][istate] / 2.0;
            }
          }//isite1 = isite2
          else {//off diagonal
            //debug spn += X_GC_child_CisAitCiuAiv_spin_MPIdouble(
            //debug   isite1 - 1, 0, 1, isite2 - 1, 1, 0, 1.0, X, nstate, vec, vec) / 2.0;
          }
        }
        else if (isite1 > X->Def.Nsite || isite2 > X->Def.Nsite) {
          if (isite1 < isite2) {
            tmp_isite1 = isite1;
            tmp_isite2 = isite2;
          }
          else {
            tmp_isite1 = isite2;
            tmp_isite2 = isite1;
          }
          is1_up = X->Def.Tpow[tmp_isite1 - 1];
          is2_up = X->Def.Tpow[tmp_isite2 - 1];
          num2_up = child_SpinGC_CisAis((unsigned long int)myrank + 1, X, is2_up, 1);
          num2_down = 1 - num2_up;
          //diagonal
          for (j = 1; j <= i_max; j++) {
            list_1_j = j - 1;
            ibit1_up = list_1_j & is1_up;
            num1_up = ibit1_up / is1_up;
            num1_down = 1 - num1_up;
            spn_z2 = (num1_up - num1_down)*(num2_up - num2_down);
            for (istate = 0; istate < nstate; istate++)
              X->Phys.s2[istate] += conj(vec[j][istate])*vec[j][istate] * spn_z2 / 4.0;
          }
          if (isite1 < isite2) {
            //debug spn += X_GC_child_CisAitCiuAiv_spin_MPIsingle(isite1 - 1, 0, 1, isite2 - 1, 1, 0, 1.0, X, nstate, vec, vec) / 2.0;
          }
          else {
            //debug spn += conj(X_GC_child_CisAitCiuAiv_spin_MPIsingle(isite2 - 1, 1, 0, isite1 - 1, 0, 1, 1.0, X, nstate, vec, vec)) / 2.0;
          }
        }
        else {
          is2_up = X->Def.Tpow[isite2 - 1];
          is_up = is1_up + is2_up;
          for (j = 1; j <= i_max; j++) {
            list_1_j = j - 1;
            ibit1_up = list_1_j & is1_up;
            num1_up = ibit1_up / is1_up;
            num1_down = 1 - num1_up;
            ibit2_up = list_1_j & is2_up;
            num2_up = ibit2_up / is2_up;
            num2_down = 1 - num2_up;

            spn_z2 = (num1_up - num1_down)*(num2_up - num2_down);
            for (istate = 0; istate < nstate; istate++)
              X->Phys.s2[istate] += conj(vec[j][istate])*vec[j][istate] * spn_z2 / 4.0;

            if (isite1 == isite2) {
              for (istate = 0; istate < nstate; istate++)
                X->Phys.s2[istate] += conj(vec[j][istate])*vec[j][istate] / 2.0;
            }
            else {
              ibit_tmp = (num1_up) ^ (num2_up);
              if (ibit_tmp != 0) {
                iexchg = list_1_j ^ (is_up);
                off = iexchg + 1;
                for (istate = 0; istate < nstate; istate++)
                  X->Phys.s2[istate] += conj(vec[j][istate])*vec[off][istate] / 2.0;
              }
            }
          }//j  
        }//else
      }
    }
  }
  else {//general spin
    double S1 = 0;
    double S2 = 0;
    for (isite1 = 1; isite1 <= X->Def.NsiteMPI; isite1++) {
      S1 = 0.5*(X->Def.SiteToBit[isite1 - 1] - 1);
      if (isite1 > X->Def.Nsite) {
        spn_z1 = 0.5*GetLocal2Sz(isite1, (unsigned long int) myrank, X->Def.SiteToBit, X->Def.Tpow);
        for (j = 1; j <= i_max; j++) {
          for (istate = 0; istate < nstate; istate++) {
            X->Phys.s2[istate] += conj(vec[j][istate])*vec[j][istate] * S1*(S1 + 1.0);
            X->Phys.Sz[istate] += conj(vec[j][istate])*vec[j][istate] * spn_z1;
          }
        }
      }
      else {
        for (j = 1; j <= i_max; j++) {
          spn_z1 = 0.5*GetLocal2Sz(isite1, j - 1, X->Def.SiteToBit, X->Def.Tpow);
          for (istate = 0; istate < nstate; istate++) {
            X->Phys.s2[istate] += conj(vec[j][istate])*vec[j][istate] * S1*(S1 + 1.0);
            X->Phys.Sz[istate] += conj(vec[j][istate])*vec[j][istate] * spn_z1;
          }
        }
      }
      for (isite2 = 1; isite2 <= X->Def.NsiteMPI; isite2++) {
        if (isite1 == isite2) continue;
        S2 = 0.5*(X->Def.SiteToBit[isite2 - 1] - 1);
        if (isite1 > X->Def.Nsite && isite2 > X->Def.Nsite) {
        }
        else if (isite1 > X->Def.Nsite || isite2 > X->Def.Nsite) {
        }
        else { //inner-process
          for (j = 1; j <= i_max; j++) {
            spn_z1 = 0.5*GetLocal2Sz(isite1, j - 1, X->Def.SiteToBit, X->Def.Tpow);
            spn_z2 = 0.5*GetLocal2Sz(isite2, j - 1, X->Def.SiteToBit, X->Def.Tpow);
            for (istate = 0; istate < nstate; istate++)
              X->Phys.s2[istate] += conj(vec[j][istate])*vec[j][istate] * spn_z1*spn_z2;

            sigma_1 = GetBitGeneral(isite1, j - 1, X->Def.SiteToBit, X->Def.Tpow);
            sigma_2 = GetBitGeneral(isite2, j - 1, X->Def.SiteToBit, X->Def.Tpow);

            ibit_tmp = GetOffCompGeneralSpin(j - 1, isite2, sigma_2, sigma_2 + 1, &off, X->Def.SiteToBit, X->Def.Tpow);
            if (ibit_tmp != 0) {
              ibit_tmp = GetOffCompGeneralSpin(off, isite1, sigma_1, sigma_1 - 1, &off_2, X->Def.SiteToBit, X->Def.Tpow);
              if (ibit_tmp != 0) {
                for (istate = 0; istate < nstate; istate++)
                  X->Phys.s2[istate] += conj(vec[j][istate])*vec[off_2 + 1][istate]
                  * sqrt(S2*(S2 + 1) - spn_z2 * (spn_z2 + 1))
                  * sqrt(S1*(S1 + 1) - spn_z1 * (spn_z1 - 1)) / 2.0;
              }
            }
            ibit_tmp = GetOffCompGeneralSpin(j - 1, isite2, sigma_2, sigma_2 - 1, &off, X->Def.SiteToBit, X->Def.Tpow);
            if (ibit_tmp != 0) {
              ibit_tmp = GetOffCompGeneralSpin(off, isite1, sigma_1, sigma_1 + 1, &off_2, X->Def.SiteToBit, X->Def.Tpow);
              if (ibit_tmp != 0) {
                for (istate = 0; istate < nstate; istate++)
                  X->Phys.s2[istate] += conj(vec[j][istate])*vec[off_2 + 1][istate]
                  * sqrt(S2*(S2 + 1) - spn_z2 * (spn_z2 - 1.0))
                  * sqrt(S1*(S1 + 1) - spn_z1 * (spn_z1 + 1)) / 2.0;
              }
            }
          }//j  
        }//inner-process
      }//isite2
    }//isite1
  }
  SumMPI_dv(nstate, X->Phys.s2);
  SumMPI_dv(nstate, X->Phys.Sz);
}
/**
 * @brief Parent function of calculation of total spin
 *
 * @param[in,out] X data list of calculation parameters
 * @param[in] vec eigenvector
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @retval 0 calculation is normally finished
 */
int expec_totalspin
(
  struct BindStruct *X,
  int nstate,
  double complex **vec
)
{
  int istate;

  X->Large.mode = M_TOTALS;
  switch (X->Def.iCalcModel) {
  case Spin:
    totalspin_Spin(X, nstate, vec);
    for (istate = 0; istate < nstate; istate++)
      X->Phys.Sz[istate] = X->Def.Total2SzMPI / 2.;
    break;
  case SpinGC:
    totalspin_SpinGC(X, nstate, vec);
    break;
  case Hubbard:
  case Kondo:
    totalspin_Hubbard(X, nstate, vec);
    break;
  case HubbardGC:
  case KondoGC:
    totalspin_HubbardGC(X, nstate, vec);
    break;
  default:
    for (istate = 0; istate < nstate; istate++) {
      X->Phys.s2[istate] = 0.0;
      X->Phys.Sz[istate] = 0.0;
    }
  }
  return 0;
}

int expec_totalSz(
  struct BindStruct* X,
  double complex** vec
) {
  X->Large.mode = M_TOTALS;
  switch (X->Def.iCalcModel) {
  case Spin:
    X->Phys.Sz[0] = X->Def.Total2SzMPI / 2.;
    break;
  case SpinGC:
    totalSz_SpinGC(X, vec);
    break;
  case Hubbard:
  case Kondo:
    X->Phys.Sz[0] = X->Def.Total2SzMPI / 2.;

    break;
  case HubbardGC:
  case KondoGC:
    totalSz_HubbardGC(X, vec);
    break;
  default:
    X->Phys.Sz[0] = 0.0;
  }

  return 0;
}


/**
 * @brief function of calculating totalSz for Hubbard model in grand canonical ensemble
 *
 * @param[in,out] X data list of calculation parameters
 * @param vec eigenvector
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void totalSz_HubbardGC
(
  struct BindStruct* X,
  double complex** vec
) {

  long unsigned int j;
  long unsigned int isite1;
  long unsigned int is1_up, is1_down;
  int num1_up, num1_down, num1_sz;
  long unsigned int ibit1_up, ibit1_down, list_1_j;
  double complex spn_z;
  long unsigned int i_max;

  i_max = X->Check.idim_max;
  spn_z = 0.0;
  for (isite1 = 1; isite1 <= X->Def.NsiteMPI; isite1++) {
    if (isite1 > X->Def.Nsite) {
#ifdef MPI
      is1_up = X->Def.Tpow[2 * isite1 - 2];
      is1_down = X->Def.Tpow[2 * isite1 - 1];
      ibit1_up = (unsigned long int) myrank & is1_up;
      num1_up = ibit1_up / is1_up;
      ibit1_down = (unsigned long int) myrank & is1_down;
      num1_down = ibit1_down / is1_down;
      num1_sz = num1_up - num1_down;
#pragma omp parallel for reduction(+:spn_z) default(none) firstprivate(i_max) private(j) shared(num1_sz,vec)
      for (j = 1; j <= i_max; j++) {
        spn_z += (num1_sz) / 2.0 * conj(vec[j][0]) * vec[j][0];
      }
#endif
    }
    else {//isite1 > X->Def.Nsite
      is1_up = X->Def.Tpow[2 * isite1 - 2];
      is1_down = X->Def.Tpow[2 * isite1 - 1];
#pragma omp parallel for reduction(+:spn_z) default(none) firstprivate(i_max, is1_up, is1_down, isite1) \
  private(list_1_j, ibit1_up, num1_up, ibit1_down, num1_down) shared(vec)
      for (j = 1; j <= i_max; j++) {
        list_1_j = j - 1;
        ibit1_up = list_1_j & is1_up;
        num1_up = ibit1_up / is1_up;
        ibit1_down = list_1_j & is1_down;
        num1_down = ibit1_down / is1_down;
        spn_z += conj(vec[j][0]) * vec[j][0] * (num1_up - num1_down) / 2.0;
      }
    }
  }
  spn_z = SumMPI_dc(spn_z);
  X->Phys.Sz[0] = creal(spn_z);
}

/**
 * @brief function of calculating totalSz for Spin model in grand canonical ensemble
 *
 * @param[in,out] X data list of calculation parameters
 * @param vec eigenvector
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void totalSz_SpinGC
(
  struct BindStruct* X,
  double complex** vec
) {
  long unsigned int j, list_1_j;
  long unsigned int isite1;
  long unsigned int is1_up;
  int num1_up;
  int num1_down;
  long unsigned int ibit1_up;
  double complex spn_z, spn_z1;
  long unsigned int i_max;
  i_max = X->Check.idim_max;
  X->Large.mode = M_TOTALS;
  spn_z = 0.0;
  if (X->Def.iFlgGeneralSpin == FALSE) {
    for (isite1 = 1; isite1 <= X->Def.NsiteMPI; isite1++) {
      if (isite1 > X->Def.Nsite) {
#ifdef MPI
        is1_up = X->Def.Tpow[isite1 - 1];
        ibit1_up = myrank & is1_up;
        num1_up = ibit1_up / is1_up;
        num1_down = 1 - num1_up;
#pragma omp parallel for reduction(+: spn_z) default(none) firstprivate(i_max, is1_up, num1_up, num1_down) shared(vec)
        for (j = 1; j <= i_max; j++) {
          spn_z += conj(vec[j][0]) * vec[j][0] * (num1_up - num1_down) / 2.0;
        }
#endif
      }
      else {
        is1_up = X->Def.Tpow[isite1 - 1];
#pragma omp parallel for reduction(+: spn_z) default(none) firstprivate(i_max, is1_up) private(list_1_j, ibit1_up, num1_up, num1_down) shared(vec)
        for (j = 1; j <= i_max; j++) {
          list_1_j = j - 1;
          ibit1_up = list_1_j & is1_up;
          num1_up = ibit1_up / is1_up;
          num1_down = 1 - num1_up;
          spn_z += conj(vec[j][0]) * vec[j][0] * (num1_up - num1_down) / 2.0;
        }
      }//else
    }
  }
  else {//general spin
    spn_z = 0.0;
    for (isite1 = 1; isite1 <= X->Def.NsiteMPI; isite1++) {
      if (isite1 > X->Def.Nsite) {
        spn_z1 = 0.5 * GetLocal2Sz(isite1, (unsigned long int) myrank, X->Def.SiteToBit, X->Def.Tpow);
#pragma omp parallel for reduction(+: spn_z) default(none) firstprivate(spn_z1, i_max) shared(vec)
        for (j = 1; j <= i_max; j++) {
          spn_z += conj(vec[j][0]) * vec[j][0] * spn_z1;
        }
      }
      else {
#pragma omp parallel for reduction(+:  spn_z) default(none) firstprivate(i_max, isite1, X) private(spn_z1) shared(vec)
        for (j = 1; j <= i_max; j++) {
          spn_z1 = 0.5 * GetLocal2Sz(isite1, j - 1, X->Def.SiteToBit, X->Def.Tpow);
          spn_z += conj(vec[j][0]) * vec[j][0] * spn_z1;
        }
      }
    }//isite1
  }
  spn_z = SumMPI_dc(spn_z);
  X->Phys.Sz[0] = creal(spn_z);
}
