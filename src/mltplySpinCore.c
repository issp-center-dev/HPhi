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
 * @file mltplySpinCore.c
 *
 * @brief Core functions for spin model Hamiltonian operations
 *
 * This module provides building blocks for spin Hamiltonians including
 * exchange, pair-lift, and general spin interactions.
 *
 * Bit representation (for S=1/2 spin):
 * - Each site has 1 bit: bit[i] = spin state (0=down, 1=up)
 * - Tpow[i] = 2^i is the bit mask for site i
 * - State |...s_2 s_1 s_0> stored as integer where s_i in {0,1}
 *
 * For general spin S:
 * - Each site has ceil(log2(2S+1)) bits
 * - SiteToBit[i] gives bits needed for site i
 * - State stored as product of local spin states
 *
 * Key operators:
 * - Exchange: S+_i S-_j + S-_i S+_j (flip two antiparallel spins)
 * - Pair-lift: S+_i S+_j (grand canonical only, changes total Sz)
 * - Ising: Sz_i Sz_j (diagonal, no state change)
 *
 * Naming convention:
 * - is1_up, is2_up: Bit masks for sites
 * - isA_spin: Combined mask for two-site operations
 * - tmp_J: Coupling constant for current term
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */

#include <bitcalc.h>
#include "xsetmem.h"
#include "wrapperMPI.h"
#include "mltplyCommon.h"
#include "mltplySpinCore.h"
#include "expec_trace_internal.h"

/******************************************************************************/
//[s] GetInfo functions
/******************************************************************************/

/**
 * @brief Set up masks for spin exchange term J(S+_i S-_j + S-_i S+_j)
 *
 * The exchange interaction flips two antiparallel spins:
 *   |...up_i...down_j...> <-> |...down_i...up_j...>
 *
 * Prepares X->Large with:
 * - is1_up, is2_up: Bit masks for checking spin states at sites i,j
 * - isA_spin: Combined mask for XOR operation (flips both spins)
 * - tmp_J: Exchange coupling constant
 *
 * The exchange operation is: state' = state XOR isA_spin
 * (valid only when exactly one of the two sites has spin up)
 *
 * @param iExchange Index into ExchangeCoupling array [in]
 * @param X Struct to store computed masks [inout]
 *
 * @return 0 (always succeeds)
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int exchange_spin_GetInfo(
  int iExchange,//!<[in] Counter of exchange interaction
  struct BindStruct *X//!<[inout]
) {
  int isite1 = X->Def.ExchangeCoupling[iExchange][0] + 1;
  int isite2 = X->Def.ExchangeCoupling[iExchange][1] + 1;
  /**
   Set the exchange coupling constant (LargeList::tmp_J)
  */
  X->Large.tmp_J = X->Def.ParaExchangeCoupling[iExchange];
  /**
  Set the bit mask for computing spin state of both site
  (LargeList::is1_up, LargeList::is2_up)
  */
  X->Large.is1_up = X->Def.Tpow[isite1 - 1];
  X->Large.is2_up = X->Def.Tpow[isite2 - 1];
  /**
  Set the bit mask for exchange 2 spins (LargeList::isA_spin)
  */
  X->Large.isA_spin = X->Large.is1_up + X->Large.is2_up;
  return 0;
}/*int exchange_spin_GetInfo*/
/**
@brief Set parameters for the bit operation of spin-pairlift term
@return Always return 0
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int pairlift_spin_GetInfo(
  int iPairLift,
  struct BindStruct *X
) {
  int isite1 = X->Def.PairLiftCoupling[iPairLift][0] + 1;
  int isite2 = X->Def.PairLiftCoupling[iPairLift][1] + 1;
  /**
  Set the pairlift coupling constant (LargeList::tmp_J)
  */
  X->Large.tmp_J = X->Def.ParaPairLiftCoupling[iPairLift];
  /**
  Set the bit mask for computing spin state of both site
  (LargeList::is1_up, LargeList::is2_up)
  */
  X->Large.is1_up = X->Def.Tpow[isite1 - 1];
  X->Large.is2_up = X->Def.Tpow[isite2 - 1];
  /**
  Set the bit mask for exchange 2 spins (LargeList::isA_spin)
  */
  X->Large.isA_spin = X->Large.is1_up + X->Large.is2_up;
  return 0;
}/*int pairlift_spin_GetInfo*/
/**
@brief Set parameters for the bit operation of spin-general interaction term
@return Always return 0
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int general_int_spin_GetInfo(
  struct BindStruct *X,//!<[inout]
  long unsigned int isite1,//!<[in] Site 1
  long unsigned int isite2,//!<[in] Site 2
  long unsigned int sigma1,//!<[in] Sigma 1, final state of site 1
  long unsigned int sigma2,//!<[in] Sigma 3, initial state of site 1
  long unsigned int sigma3,//!<[in] Sigma 3, final state of site 2
  long unsigned int sigma4,//!<[in] Sigma 4, initial state of site 2
  double complex tmp_V//!<[in] General interaction constatnt
) {
  /**
  Set the pairlift coupling constant (LargeList::tmp_J)
  */
  X->Large.tmp_V = tmp_V;
  X->Large.isite1 = isite1;
  X->Large.isite2 = isite2;
  /**
  Set the bit mask for computing spin state of both site
  (LargeList::is1_up, LargeList::is2_up)
  */
  X->Large.is1_up = X->Def.Tpow[isite1 - 1];
  X->Large.is2_up = X->Def.Tpow[isite2 - 1];
  /**
  Set the bit mask for general interaction 
  (LargeList::is1_spin, LargeList::is2_spin, LargeList::is3_spin, LargeList::is4_spin)
  */
  X->Large.is1_spin = sigma1;
  X->Large.is2_spin = sigma2;
  X->Large.is3_spin = sigma3;
  X->Large.is4_spin = sigma4;
  return 0;
}/*int general_int_spin_GetInfo*/

/******************************************************************************/
//[e] GetInfo functions
/******************************************************************************/

/******************************************************************************/
//[s] core routines
/******************************************************************************/

/**
@brief Compute index of final wavefunction by @f$c_{is}^\dagger c_{it}@f$ term.
@return 1 if bit-mask is1_spin is sigma2, 0 for the other
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int child_Spin_CisAit(
  long unsigned int j,//!<[in] Index of initial wavefunction
  struct BindStruct *X,//!<[inout]
  long unsigned int is1_spin,//!<[in] Bit mask for computing spin state
  long unsigned int sigma2,//!<[in] Spin state at site 2
  long unsigned int *list_1_Org_,//!<[in] Similar to ::list_1
  long unsigned int *list_2_1_,//!<[in] Similar to ::list_2_1
  long unsigned int *list_2_2_,//!<[in] Similar to ::list_2_2
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  long unsigned int list_1_j;
  long unsigned int off;
  list_1_j = list_1_Org_[j];
  if (child_SpinGC_CisAit(list_1_j + 1, X, is1_spin, sigma2, &off) != 0) {
    GetOffComp(list_2_1_, list_2_2_, off, X->Large.irght, X->Large.ilft, X->Large.ihfbit, tmp_off);
    return 1;
  }
  else {
    *tmp_off = 0;
    return 0;
  }
}/*int child_Spin_CisAit*/
/**
@brief Compute the spin state with bit mask is1_spin
@return 1 if the spin state with bit mask is1_spin is sigma1, 0 for the other.
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int child_Spin_CisAis(
  long unsigned int j,//!<[in] Index of wavefunction
  struct BindStruct *X,//!<[inout]
  long unsigned int is1_spin,//!<[in] Bit mask
  long unsigned int sigma1//!<[in] Target spin state
) {
  int A_ibit_tmp;
  // off = j
  A_ibit_tmp = ((list_1[j] & is1_spin) / is1_spin) ^ (1 - sigma1);
  return A_ibit_tmp;
}/*int child_Spin_CisAis*/
/**
@brief Compute the grandcanonical spin state with bit mask is1_spin
@return 1 if the spin state with bit mask is1_spin is sigma1, 0 for the other.
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int child_SpinGC_CisAis(
  long unsigned int j,//!<[in] Index of wavefunction
  struct BindStruct *X,//!<[inout]
  long unsigned int is1_spin,//!<[in] Bit mask
  long unsigned int sigma1//!<[in] Target spin state
) {
  int A_ibit_tmp;
  long unsigned int list_1_j;
  // off = j
  list_1_j = j - 1;
  A_ibit_tmp = ((list_1_j & is1_spin) / is1_spin) ^ (1 - sigma1);
  return A_ibit_tmp;
}/*int child_SpinGC_CisAis*/
/**
@brief Compute index of final wavefunction by @f$c_{is}^\dagger c_{it}@f$ term
(grandcanonical).
@return 1 if bit-mask is1_spin is sigma2, 0 for the other
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int child_SpinGC_CisAit(
  long unsigned int j,//!<[in] Index of wavefunction
  struct BindStruct *X,//!<[inout]
  long unsigned int is1_spin,//!<[in] Bit mask for computing spin state
  long unsigned int sigma2,//!<[in] Spin state at site 2
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  long unsigned int list_1_j, ibit_tmp_1;

  list_1_j = j - 1;

  ibit_tmp_1 = list_1_j & is1_spin;
  if (ibit_tmp_1 == 0 && sigma2 == 0) {    // down -> up
    *tmp_off = list_1_j + is1_spin;
    return 1;
  }
  else if (ibit_tmp_1 != 0 && sigma2 == 1) { // up -> down
    *tmp_off = list_1_j - is1_spin;
    return 1;
  }
  else {
    *tmp_off = 0;
    return 0;
  }
}/*int child_SpinGC_CisAit*/

/**
@brief Trace-probe adapter for ::child_Spin_CisAis (canonical half-spin diagonal
one-body). The underlying function is already a pure map, so the probe simply
reports it: diagonal, always returns 1 with *kprime_out = j-1 and *amp_out the
0/1 spin match.
*/
int child_Spin_CisAis_TraceProbe(
  long unsigned int j, struct BindStruct *X,
  long unsigned int is1_spin, long unsigned int sigma1,
  long int *kprime_out, double complex *amp_out
) {
  *kprime_out = (long int)j - 1;
  *amp_out = (double complex)child_Spin_CisAis(j, X, is1_spin, sigma1);
  return 1;
}/*int child_Spin_CisAis_TraceProbe*/
/**
@brief Trace-probe adapter for ::child_SpinGC_CisAis (grand-canonical half-spin
diagonal one-body).
*/
int child_SpinGC_CisAis_TraceProbe(
  long unsigned int j, struct BindStruct *X,
  long unsigned int is1_spin, long unsigned int sigma1,
  long int *kprime_out, double complex *amp_out
) {
  *kprime_out = (long int)j - 1;
  *amp_out = (double complex)child_SpinGC_CisAis(j, X, is1_spin, sigma1);
  return 1;
}/*int child_SpinGC_CisAis_TraceProbe*/
/**
@brief Trace-probe adapter for ::child_SpinGC_CisAit (grand-canonical half-spin
transverse off-diagonal one-body). Bare-bit flip: on a surviving flip *kprime_out
= tmp_off (0-based) and *amp_out = the sign (+1); otherwise annihilated.
*/
int child_SpinGC_CisAit_TraceProbe(
  long unsigned int j, struct BindStruct *X,
  long unsigned int is1_spin, long unsigned int sigma2,
  long int *kprime_out, double complex *amp_out
) {
  long unsigned int tmp_off = 0;
  int sgn = child_SpinGC_CisAit(j, X, is1_spin, sigma2, &tmp_off);
  if (sgn != 0) {
    *kprime_out = (long int)tmp_off;
    *amp_out = (double complex)sgn;
    return 1;
  }
  *kprime_out = -1;
  *amp_out = 0.0;
  return 0;
}/*int child_SpinGC_CisAit_TraceProbe*/

/******************************************************************************/
//[e] core routines
/******************************************************************************/

/**
@brief Compute index of final wavefunction associated to spin-exchange term
@return 1 if spin of site 1 is sigmaA and spin of site 2 is sigmaB. 0 if not. 
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int child_exchange_spin_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  struct BindStruct *X,//!<[inout]
  long unsigned int isA_up,//!<[in] Bit mask for spin 1
  long unsigned int isB_up,//!<[in] Bit mask for spin 2
  long unsigned int sigmaA,//!<[in] Target of spin 1
  long unsigned int sigmaB,//!<[in] Target of spin 2
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  long unsigned int iexchg, off;
  long unsigned int irght = X->Large.irght;
  long unsigned int ilft = X->Large.ilft;
  long unsigned int ihfbit = X->Large.ihfbit;
  long unsigned int ibit_tmp_A, ibit_tmp_B;

  ibit_tmp_A = ((list_1[j] & isA_up) / isA_up);
  ibit_tmp_B = ((list_1[j] & isB_up) / isB_up);
  if (ibit_tmp_A == sigmaA && ibit_tmp_B == sigmaB) {
    iexchg = list_1[j] ^ (isA_up + isB_up);
    GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off);
    *tmp_off = off;
    return 1;
  }
  else {
    *tmp_off = 0; // just tentative
    return 0;
  }
}/*int child_exchange_spin_element*/
/**
@brief Trace-probe adapter for ::child_exchange_spin_element (canonical half-spin
two-body EXCHANGE branch). The underlying function is already a pure map
(GetOffComp yields a 1-based canonical index), so the probe reports
*kprime_out = tmp_off-1 and *amp_out = the 0/1 exchange sign. NOTE: the Mode-1
caller (expec_cisajscktalt_SpinHalf's exchange branch) multiplies this
contribution by the raw sign ONLY, deliberately NOT by tmp_V -- so this probe
reports the bare sign and the extraction driver does not reintroduce tmp_V.
See docs/superpowers/specs/2026-07-11-expec-call-inventory.md 2c.
*/
int child_exchange_spin_element_TraceProbe(
  long unsigned int j, struct BindStruct *X,
  long unsigned int isA_up, long unsigned int isB_up,
  long unsigned int sigmaA, long unsigned int sigmaB,
  long int *kprime_out, double complex *amp_out
) {
  long unsigned int tmp_off = 0;
  int flag = child_exchange_spin_element(j, X, isA_up, isB_up, sigmaA, sigmaB, &tmp_off);
  if (flag != 0) {
    *kprime_out = (long int)tmp_off - 1;
    *amp_out = (double complex)flag;
    return 1;
  }
  *kprime_out = -1;
  *amp_out = 0.0;
  return 0;
}/*int child_exchange_spin_element_TraceProbe*/
/**
@brief Multiply Hamiltonian of exchange term of canonical spin system
@return @f$\langle v_1 | H_{\rm this}| v_1 \rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex exchange_spin_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  double complex *tmp_v0,//!<[out] Resulting wavefunction
  double complex *tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  long unsigned int off;
  double complex dmv;
  long unsigned int iexchg;
  long unsigned int is_up = X->Large.isA_spin;
  long unsigned int irght = X->Large.irght;
  long unsigned int ilft = X->Large.ilft;
  long unsigned int ihfbit = X->Large.ihfbit;
  double complex tmp_J = X->Large.tmp_J;
  int mode = X->Large.mode;
  double complex dam_pr = 0;
  long unsigned int ibit_tmp;

  ibit_tmp = (list_1[j] & is_up);
  if (ibit_tmp == 0 || ibit_tmp == is_up) {
    return dam_pr;
  }
  else {
    iexchg = list_1[j] ^ is_up;
    GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off);
    *tmp_off = off;
    dmv = tmp_J * tmp_v1[j];
    if (mode == M_MLTPLY) {
      tmp_v0[off] += dmv;
    }
    dam_pr += dmv * conj(tmp_v1[off]);
    return dam_pr;
  }
}/*double complex exchange_spin_element*/
/**
@brief Multiply Hamiltonian of exchange term of grandcanonical spin system
@return @f$\langle v_1 | H_{\rm this}| v_1 \rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex GC_exchange_spin_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  double complex *tmp_v0,//!<[out] Resulting wavefunction
  double complex *tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  double complex dmv;
  long unsigned int is_up = X->Large.isA_spin;
  double complex tmp_J = X->Large.tmp_J;
  int mode = X->Large.mode;
  long unsigned int list_1_j, list_1_off;

  double complex dam_pr = 0;
  list_1_j = j - 1;

  long unsigned int ibit_tmp;
  ibit_tmp = (list_1_j & is_up);
  if (ibit_tmp == 0 || ibit_tmp == is_up) {
    return dam_pr;
  }
  else {
    list_1_off = list_1_j ^ is_up;
    *tmp_off = list_1_off;
    dmv = tmp_J * tmp_v1[j];
    if (mode == M_MLTPLY) {
      tmp_v0[list_1_off + 1] += dmv;
    }
    dam_pr += dmv * conj(tmp_v1[list_1_off + 1]);
    return dam_pr;
  }
}/*double complex GC_exchange_spin_element*/
/**
@brief Multiply Hamiltonian of pairlift term of grandcanonical spin system
@return @f$\langle v_1 | H_{\rm this}| v_1 \rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex GC_pairlift_spin_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  double complex *tmp_v0,//!<[out] Resulting wavefunction
  double complex *tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  double complex dmv;
  long unsigned int is_up = X->Large.isA_spin;
  double complex tmp_J = X->Large.tmp_J;
  int mode = X->Large.mode;
  double complex dam_pr = 0;
  long unsigned int list_1_off;
  long unsigned int list_1_j = j - 1;
  long unsigned int ibit_tmp;
  //ibit_tmp = ((list_1_j & is1_up) / is1_up) ^ ((list_1_j & is2_up) / is2_up);
  ibit_tmp = (list_1_j & is_up);
  if (ibit_tmp == 0 || ibit_tmp == is_up) {
    list_1_off = list_1_j ^ is_up; //Change: ++ -> -- or -- -> ++
    *tmp_off = list_1_off;
    dmv = tmp_J * tmp_v1[j];//* ibit_tmp;
    if (mode == M_MLTPLY) {
      tmp_v0[list_1_off + 1] += dmv;
    }
    dam_pr += dmv * conj(tmp_v1[list_1_off + 1]);
    return dam_pr;
  }
  else {
    return dam_pr;
  }
}/*double complex GC_pairlift_spin_element*/

//[s]Spin
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{is}^\dagger c_{is}@f$ term of 
canonical spsin system
@return Fragment of @f$\langle v_1 | H_{\rm this}|v_1 \rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
/**
@brief Mapping core of ::CisAisCisAis_spin_element (canonical half-spin diagonal).
@return signed spin-match product (0/1); kprime_out = j-1 (always valid).
*/
static int CisAisCisAis_spin_element_map(
  long unsigned int j, long unsigned int isA_up, long unsigned int isB_up,
  long unsigned int org_sigma2, long unsigned int org_sigma4,
  struct BindStruct *X, long int *kprime_out
) {
  int tmp_sgn;
  tmp_sgn = child_Spin_CisAis(j, X, isB_up, org_sigma4);
  tmp_sgn *= child_Spin_CisAis(j, X, isA_up, org_sigma2);
  *kprime_out = (long int)j - 1;
  return tmp_sgn;
}
double complex CisAisCisAis_spin_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int isA_up,//!<[in] Bit mask for spin 1
  long unsigned int isB_up,//!<[in] Bit mask for spin 2
  long unsigned int org_sigma2,//!<[in] Target for spin 1
  long unsigned int org_sigma4,//!<[in] Target for spin 2
  double complex tmp_V,//!<[in] Coupling constatnt
  double complex *tmp_v0,//!<[in] Resulting wavefunction
  double complex *tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X//!<[inout]
) {
  int tmp_sgn;
  long int kprime;
  double complex dmv;
  double complex dam_pr = 0;

  tmp_sgn = CisAisCisAis_spin_element_map(j, isA_up, isB_up, org_sigma2, org_sigma4, X, &kprime);
  dmv = tmp_v1[j] * tmp_sgn * tmp_V;
  if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
    tmp_v0[kprime + 1] += dmv;
  }
  dam_pr = conj(tmp_v1[kprime + 1]) * dmv;
  return dam_pr;
}/*double complex CisAisCisAis_spin_element*/
/**
@brief Trace-probe adapter for ::CisAisCisAis_spin_element (canonical half-spin
diagonal two-body). Diagonal: always returns 1, *amp_out = tmp_V * signed match.
*/
int CisAisCisAis_spin_element_TraceProbe(
  long unsigned int j, long unsigned int isA_up, long unsigned int isB_up,
  long unsigned int org_sigma2, long unsigned int org_sigma4,
  double complex tmp_V, struct BindStruct *X,
  long int *kprime_out, double complex *amp_out
) {
  int tmp_sgn = CisAisCisAis_spin_element_map(j, isA_up, isB_up, org_sigma2, org_sigma4, X, kprime_out);
  *amp_out = tmp_V * (double complex)tmp_sgn;
  return 1;
}/*int CisAisCisAis_spin_element_TraceProbe*/

//[e]Spin

//[s]GC Spin
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{is}^\dagger c_{is}@f$ term of 
grandcanonical spsin system
@return Fragment of @f$\langle v_1 | H_{\rm this}|v_1 \rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
/**
@brief Mapping core of ::GC_CisAisCisAis_spin_element (grand-canonical diagonal).
@return signed spin-match product (0/1); kprime_out = j-1 (always valid).
*/
static int GC_CisAisCisAis_spin_element_map(
  long unsigned int j, long unsigned int isA_up, long unsigned int isB_up,
  long unsigned int org_sigma2, long unsigned int org_sigma4,
  struct BindStruct *X, long int *kprime_out
) {
  int tmp_sgn;
  tmp_sgn = child_SpinGC_CisAis(j, X, isB_up, org_sigma4);
  tmp_sgn *= child_SpinGC_CisAis(j, X, isA_up, org_sigma2);
  *kprime_out = (long int)j - 1;
  return tmp_sgn;
}
double complex GC_CisAisCisAis_spin_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int isA_up,//!<[in] Bit mask for spin 1
  long unsigned int isB_up,//!<[in] Bit mask for spin 2
  long unsigned int org_sigma2,//!<[in] Target for spin 1
  long unsigned int org_sigma4,//!<[in] Target for spin 2
  double complex tmp_V,//!<[in] Coupling constatnt
  double complex *tmp_v0,//!<[in] Resulting wavefunction
  double complex *tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X//!<[inout]
) {
  int tmp_sgn;
  long int kprime;
  double complex dmv = 0;
  double complex dam_pr = 0;

  tmp_sgn = GC_CisAisCisAis_spin_element_map(j, isA_up, isB_up, org_sigma2, org_sigma4, X, &kprime);
  if (tmp_sgn != 0) {
    dmv = tmp_v1[j] * tmp_sgn * tmp_V;
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
      tmp_v0[kprime + 1] += dmv;
      dam_pr = conj(tmp_v1[kprime + 1]) * dmv;
    }else if(X->Large.mode == H_CORR){
      dam_pr = conj(tmp_v0[kprime + 1]) * dmv;
    }else{
      dam_pr = conj(tmp_v1[kprime + 1]) * dmv;
    }
  }
  return dam_pr;
}/*double complex GC_CisAisCisAis_spin_element*/
/**
@brief Trace-probe adapter for ::GC_CisAisCisAis_spin_element (SpinGC-half
diagonal two-body). Diagonal: always returns 1, *amp_out = tmp_V * signed match.
The H_CORR branch is never used by extraction (M_CORR only).
*/
int GC_CisAisCisAis_spin_element_TraceProbe(
  long unsigned int j, long unsigned int isA_up, long unsigned int isB_up,
  long unsigned int org_sigma2, long unsigned int org_sigma4,
  double complex tmp_V, struct BindStruct *X,
  long int *kprime_out, double complex *amp_out
) {
  int tmp_sgn = GC_CisAisCisAis_spin_element_map(j, isA_up, isB_up, org_sigma2, org_sigma4, X, kprime_out);
  *amp_out = tmp_V * (double complex)tmp_sgn;
  return 1;
}/*int GC_CisAisCisAis_spin_element_TraceProbe*/
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{it}^\dagger c_{iu}@f$ term of 
grandcanonical spsin system
@return Fragment of @f$\langle v_1 | H_{\rm this}|v_1 \rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
/**
@brief Mapping core of ::GC_CisAisCitAiu_spin_element. child_SpinGC_CisAit
(bare-bit *tmp_off) gated by child_SpinGC_CisAis. @return signed result
(0 = dead); kprime_out = *tmp_off on survival, else -1.
*/
static int GC_CisAisCitAiu_spin_element_map(
  long unsigned int j, long unsigned int org_sigma2, long unsigned int org_sigma4,
  long unsigned int isA_up, long unsigned int isB_up,
  struct BindStruct *X, long unsigned int *tmp_off, long int *kprime_out
) {
  int tmp_sgn;
  tmp_sgn = child_SpinGC_CisAit(j, X, isB_up, org_sigma4, tmp_off);
  if (tmp_sgn != 0) {
    tmp_sgn *= child_SpinGC_CisAis((*tmp_off + 1), X, isA_up, org_sigma2);
    if (tmp_sgn != 0) {
      *kprime_out = (long int)(*tmp_off);
      return tmp_sgn;
    }
  }
  *kprime_out = -1;
  return 0;
}
double complex GC_CisAisCitAiu_spin_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int org_sigma2,//!<[in] Target for spin 1
  long unsigned int org_sigma4,//!<[in] Target for spin 2
  long unsigned int isA_up,//!<[in] Bit mask for spin 1
  long unsigned int isB_up,//!<[in] Bit mask for spin 2
  double complex tmp_V,//!<[in] Coupling constatnt
  double complex *tmp_v0,//!<[in] Resulting wavefunction
  double complex *tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  long int kprime;
  double complex dmv;
  double complex dam_pr = 0 + 0 * I;
  tmp_sgn = GC_CisAisCitAiu_spin_element_map(j, org_sigma2, org_sigma4, isA_up, isB_up, X, tmp_off, &kprime);
  if (tmp_sgn != 0) {
    dmv = tmp_v1[j] * tmp_sgn * tmp_V;
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
      tmp_v0[kprime + 1] += dmv;
      dam_pr = conj(tmp_v1[kprime + 1]) * dmv;
    }else if(X->Large.mode == H_CORR){
      dam_pr = conj(tmp_v0[kprime + 1]) * dmv;
    }else{
      dam_pr = conj(tmp_v1[kprime + 1]) * dmv;
    }
  }/*if (tmp_sgn != 0)*/
  return dam_pr;
}/*double complex GC_CisAisCitAiu_spin_element*/
/**
@brief Trace-probe adapter for ::GC_CisAisCitAiu_spin_element (SpinGC-half
two-body).
*/
int GC_CisAisCitAiu_spin_element_TraceProbe(
  long unsigned int j, long unsigned int org_sigma2, long unsigned int org_sigma4,
  long unsigned int isA_up, long unsigned int isB_up,
  double complex tmp_V, struct BindStruct *X,
  long int *kprime_out, double complex *amp_out
) {
  long unsigned int tmp_off = 0;
  int tmp_sgn = GC_CisAisCitAiu_spin_element_map(j, org_sigma2, org_sigma4, isA_up, isB_up, X, &tmp_off, kprime_out);
  if (tmp_sgn != 0) {
    *amp_out = tmp_V * (double complex)tmp_sgn;
    return 1;
  }
  *amp_out = 0.0;
  return 0;
}/*int GC_CisAisCitAiu_spin_element_TraceProbe*/
/**
@brief Compute @f$c_{is}^\dagger c_{it} c_{iu}^\dagger c_{iu}@f$ term of 
grandcanonical spsin system
@return Fragment of @f$\langle v_1 | H_{\rm this}|v_1 \rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
/**
@brief Mapping core of ::GC_CisAitCiuAiu_spin_element. child_SpinGC_CisAis gate
then child_SpinGC_CisAit (bare-bit *tmp_off). @return signed result (0 = dead);
kprime_out = *tmp_off on survival, else -1.
*/
static int GC_CisAitCiuAiu_spin_element_map(
  long unsigned int j, long unsigned int org_sigma2, long unsigned int org_sigma4,
  long unsigned int isA_up, long unsigned int isB_up,
  struct BindStruct *X, long unsigned int *tmp_off, long int *kprime_out
) {
  int tmp_sgn;
  tmp_sgn = child_SpinGC_CisAis(j, X, isB_up, org_sigma4);
  if (tmp_sgn != 0) {
    tmp_sgn *= child_SpinGC_CisAit(j, X, isA_up, org_sigma2, tmp_off);
    if (tmp_sgn != 0) {
      *kprime_out = (long int)(*tmp_off);
      return tmp_sgn;
    }
  }
  *kprime_out = -1;
  return 0;
}
double complex GC_CisAitCiuAiu_spin_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int org_sigma2,//!<[in] Target for spin 1
  long unsigned int org_sigma4,//!<[in] Target for spin 2
  long unsigned int isA_up,//!<[in] Bit mask for spin 1
  long unsigned int isB_up,//!<[in] Bit mask for spin 2
  double complex tmp_V,//!<[in] Coupling constatnt
  double complex *tmp_v0,//!<[in] Resulting wavefunction
  double complex *tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  long int kprime;
  double complex dmv;
  double complex dam_pr = 0 + 0 * I;
  tmp_sgn = GC_CisAitCiuAiu_spin_element_map(j, org_sigma2, org_sigma4, isA_up, isB_up, X, tmp_off, &kprime);
  if (tmp_sgn != 0) {
    dmv = tmp_v1[j] * tmp_sgn * tmp_V;
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
      tmp_v0[kprime + 1] += dmv;
      dam_pr = conj(tmp_v1[kprime + 1]) * dmv;
    }else if(X->Large.mode == H_CORR){
      dam_pr = conj(tmp_v0[kprime + 1]) * dmv;
    }else{
      dam_pr = conj(tmp_v1[kprime + 1]) * dmv;
    }
  }/*if (tmp_sgn != 0)*/
  return dam_pr;
}/*double complex GC_CisAitCiuAiu_spin_element*/
/**
@brief Trace-probe adapter for ::GC_CisAitCiuAiu_spin_element (SpinGC-half
two-body).
*/
int GC_CisAitCiuAiu_spin_element_TraceProbe(
  long unsigned int j, long unsigned int org_sigma2, long unsigned int org_sigma4,
  long unsigned int isA_up, long unsigned int isB_up,
  double complex tmp_V, struct BindStruct *X,
  long int *kprime_out, double complex *amp_out
) {
  long unsigned int tmp_off = 0;
  int tmp_sgn = GC_CisAitCiuAiu_spin_element_map(j, org_sigma2, org_sigma4, isA_up, isB_up, X, &tmp_off, kprime_out);
  if (tmp_sgn != 0) {
    *amp_out = tmp_V * (double complex)tmp_sgn;
    return 1;
  }
  *amp_out = 0.0;
  return 0;
}/*int GC_CisAitCiuAiu_spin_element_TraceProbe*/
/**
@brief Compute @f$c_{is}^\dagger c_{it} c_{iu}^\dagger c_{iv}@f$ term of
grandcanonical spsin system
@return Fragment of @f$\langle v_1 | H_{\rm this}|v_1 \rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
/**
@brief Mapping core of ::GC_CisAitCiuAiv_spin_element. Two chained
child_SpinGC_CisAit (intermediate bare-bit tmp_off_1, final bare-bit
*tmp_off_2). @return signed result (0 = dead); kprime_out = *tmp_off_2 on
survival, else -1.
*/
static int GC_CisAitCiuAiv_spin_element_map(
  long unsigned int j, long unsigned int org_sigma2, long unsigned int org_sigma4,
  long unsigned int isA_up, long unsigned int isB_up,
  struct BindStruct *X, long unsigned int *tmp_off_2, long int *kprime_out
) {
  int tmp_sgn;
  long unsigned int tmp_off_1;
  tmp_sgn = child_SpinGC_CisAit(j, X, isB_up, org_sigma4, &tmp_off_1);
  if (tmp_sgn != 0) {
    tmp_sgn *= child_SpinGC_CisAit((tmp_off_1 + 1), X, isA_up, org_sigma2, tmp_off_2);
    if (tmp_sgn != 0) {
      *kprime_out = (long int)(*tmp_off_2);
      return tmp_sgn;
    }
  }
  *kprime_out = -1;
  return 0;
}
double complex GC_CisAitCiuAiv_spin_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int org_sigma2,//!<[in] Target for spin 1
  long unsigned int org_sigma4,//!<[in] Target for spin 2
  long unsigned int isA_up,//!<[in] Bit mask for spin 1
  long unsigned int isB_up,//!<[in] Bit mask for spin 2
  double complex tmp_V,//!<[in] Coupling constatnt
  double complex *tmp_v0,//!<[in] Resulting wavefunction
  double complex *tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off_2//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  long int kprime;
  double complex dmv;
  double complex dam_pr = 0 + 0 * I;
  tmp_sgn = GC_CisAitCiuAiv_spin_element_map(j, org_sigma2, org_sigma4, isA_up, isB_up, X, tmp_off_2, &kprime);
  if (tmp_sgn != 0) {
    dmv = tmp_v1[j] * tmp_sgn * tmp_V;
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
      tmp_v0[kprime + 1] += dmv;
      dam_pr = conj(tmp_v1[kprime + 1]) * dmv;
    }else if(X->Large.mode == H_CORR){
      dam_pr = conj(tmp_v0[kprime + 1]) * dmv;
    }else{
      dam_pr = conj(tmp_v1[kprime + 1]) * dmv;
    }
  }/*if (tmp_sgn != 0)*/
  return dam_pr;
}/*double complex GC_CisAitCiuAiv_spin_element*/
/**
@brief Trace-probe adapter for ::GC_CisAitCiuAiv_spin_element (SpinGC-half
two-body).
*/
int GC_CisAitCiuAiv_spin_element_TraceProbe(
  long unsigned int j, long unsigned int org_sigma2, long unsigned int org_sigma4,
  long unsigned int isA_up, long unsigned int isB_up,
  double complex tmp_V, struct BindStruct *X,
  long int *kprime_out, double complex *amp_out
) {
  long unsigned int tmp_off_2 = 0;
  int tmp_sgn = GC_CisAitCiuAiv_spin_element_map(j, org_sigma2, org_sigma4, isA_up, isB_up, X, &tmp_off_2, kprime_out);
  if (tmp_sgn != 0) {
    *amp_out = tmp_V * (double complex)tmp_sgn;
    return 1;
  }
  *amp_out = 0.0;
  return 0;
}/*int GC_CisAitCiuAiv_spin_element_TraceProbe*/
//[e]GC Spin
