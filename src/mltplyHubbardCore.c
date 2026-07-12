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
 * @file mltplyHubbardCore.c
 *
 * @brief Core functions for Hubbard model Hamiltonian operations
 *
 * This module provides the fundamental building blocks for Hubbard-type
 * Hamiltonians, including hopping, on-site interaction, exchange, and
 * pair-hopping terms.
 *
 * Key concepts:
 *
 * Bit representation (for Hubbard):
 * - Each site has 2 bits: bit[2*i] = up-spin, bit[2*i+1] = down-spin
 *   (see CheckMPI.c: SpinNum==1 (binary 01) -> Nup; SpinNum==2 (binary 10) -> Ndown)
 * - Tpow[2*i+sigma] = 2^(2*i+sigma) is the bit mask for (site i, spin sigma)
 * - State |...n_{1,down} n_{1,up} n_{0,down} n_{0,up}> stored as integer
 *   (left = MSB; so the rightmost pair (n_{0,down} n_{0,up}) is bits 1,0)
 *
 * Mask operations (set by GetInfo functions):
 * - is1_spin, is2_spin: Bit masks for checking occupations
 * - A_spin: Mask for computing fermion sign (bits between sites)
 * - isA_spin = is1_spin + is2_spin: Combined mask for hopping
 *
 * Fermion sign:
 * - When c†_i c_j acts, sign = (-1)^(number of particles between i and j)
 * - Computed via SgnBit(state & A_spin)
 *
 * Naming convention:
 * - Cis = c†_{i,sigma} (creation)
 * - Ajt = c_{j,tau} (annihilation)
 * - CisAis = n_{i,sigma} = c†_{i,sigma} c_{i,sigma} (number operator)
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
#include <bitcalc.h>
#include "xsetmem.h"
#include "wrapperMPI.h"
#include "mltplyCommon.h"
#include "mltplyHubbardCore.h"
#include "expec_trace_internal.h"

/******************************************************************************/
//[s] GetInfo functions
/******************************************************************************/

/**
 * @brief Set up bit masks for hopping term c†_{i1,s1} c_{i2,s2}
 *
 * Prepares X->Large with masks needed for hopping operations:
 * - is1_spin: Mask for site i1, spin s1 (= Tpow[2*i1-2+s1])
 * - is2_spin: Mask for site i2, spin s2 (= Tpow[2*i2-2+s2])
 * - A_spin: Mask for bits between (i1,s1) and (i2,s2), used for fermion sign
 * - isA_spin: Combined mask = is1_spin + is2_spin
 *
 * The fermion sign for hopping from (i2,s2) to (i1,s1) is:
 *   sign = SgnBit(state & A_spin)
 *
 * @param X Struct to store computed masks [inout]
 * @param isite1 Site index (1-based) [in]
 * @param isite2 Site index (1-based) [in]
 * @param sigma1 Spin index (0=down, 1=up) [in]
 * @param sigma2 Spin index (0=down, 1=up) [in]
 *
 * @return 0 (always succeeds)
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int general_hopp_GetInfo(
  struct BindStruct *X,//!<[inout]
  unsigned long int isite1,//!<[in] Site index
  unsigned long int isite2,//!<[in] Site index
  unsigned long int sigma1,//!<[in] Spin index
  unsigned long int sigma2//!<[in] Spin index
) {
  /**
  Compute mask for checking occupations of @f$(i_1,\sigma_1)@f$ (LargeList::is1_spin)
  and @f$(i_2,\sigma_2)@f$ (LargeList::is2_spin)
  */
  X->Large.is1_spin = X->Def.Tpow[2 * isite1 - 2 + sigma1];
  X->Large.is2_spin = X->Def.Tpow[2 * isite2 - 2 + sigma2];
  /**
  Compute mask for Fermion sign (LargeList::A_spin)
  */
  if (isite1 > isite2) {
    X->Large.A_spin = (X->Def.Tpow[2 * isite1 - 2 + sigma1] - X->Def.Tpow[2 * isite2 - 1 + sigma2]);
  }
  else if (isite1 < isite2) {
    X->Large.A_spin = (X->Def.Tpow[2 * isite2 - 2 + sigma2] - X->Def.Tpow[2 * isite1 - 1 + sigma1]);
  }
  else {
    if (sigma1 > sigma2) {
      X->Large.A_spin = (X->Def.Tpow[2 * isite1 - 2 + sigma1] - X->Def.Tpow[2 * isite2 - 1 + sigma2]);
    }
    else {
      X->Large.A_spin = (X->Def.Tpow[2 * isite2 - 2 + sigma2] - X->Def.Tpow[2 * isite1 - 1 + sigma1]);
    }
  }
  /**
  Compute mask for hopping (LargeList::isA_spin)
  */
  X->Large.isA_spin = X->Large.is1_spin + X->Large.is2_spin;
  return 0;
}/*int general_hopp_GetInfo*/
/**
@brief Compute mask for bit operation of general interaction term.
@return Error-code, always return 0
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int general_int_GetInfo(
  int iInterAll,//!<[in] It is not used
  struct BindStruct *X,//!<[inout]
  long unsigned int isite1,//!<[in] Site index
  long unsigned int isite2,//!<[in] Site index
  long unsigned int isite3,//!<[in] Site index
  long unsigned int isite4,//!<[in] Site index
  long unsigned int sigma1,//!<[in] Spin index
  long unsigned int sigma2,//!<[in] Spin index
  long unsigned int sigma3,//!<[in] Spin index
  long unsigned int sigma4,//!<[in] Spin index
  double complex tmp_V//!<[in] Coupling constant
) {
  long unsigned int is1_spin, is2_spin, is3_spin, is4_spin;
  long unsigned int A_spin, B_spin;
  long unsigned int isA_spin, isB_spin;
  /**
  Compute mask for checking occupations of @f$(i_1,\sigma_1)@f$ (LargeList::is1_spin)
  and @f$(i_2,\sigma_2)@f$ (LargeList::is2_spin)
  */
  is1_spin = X->Def.Tpow[2 * isite1 - 2 + sigma1];
  is2_spin = X->Def.Tpow[2 * isite2 - 2 + sigma2];
  /**
  Compute mask for Fermion sign (LargeList::A_spin)
  */
  if (isite1 > isite2) {
    A_spin = (X->Def.Tpow[2 * isite1 - 2 + sigma1] - X->Def.Tpow[2 * isite2 - 1 + sigma2]);
  }
  else if (isite2 > isite1) {
    A_spin = (X->Def.Tpow[2 * isite2 - 2 + sigma2] - X->Def.Tpow[2 * isite1 - 1 + sigma1]);
  }
  else {//isite1=isite2
    if (sigma1 > sigma2) {
      A_spin = (X->Def.Tpow[2 * isite1 - 2 + sigma1] - X->Def.Tpow[2 * isite2 - 1 + sigma2]);
    }
    else {
      A_spin = (X->Def.Tpow[2 * isite2 - 2 + sigma2] - X->Def.Tpow[2 * isite1 - 1 + sigma1]);
    }
  }
  /**
  Compute mask for checking occupations of @f$(i_3,\sigma_3)@f$ (LargeList::is3_spin)
  and @f$(i_4,\sigma_4)@f$ (LargeList::is4_spin)
  */
  is3_spin = X->Def.Tpow[2 * isite3 - 2 + sigma3];
  is4_spin = X->Def.Tpow[2 * isite4 - 2 + sigma4];
  /**
  Compute mask for Fermion sign (LargeList::B_spin)
  */
  if (isite3 > isite4) {
    B_spin = (X->Def.Tpow[2 * isite3 - 2 + sigma3] - X->Def.Tpow[2 * isite4 - 1 + sigma4]);
  }
  else if (isite3 < isite4) {
    B_spin = (X->Def.Tpow[2 * isite4 - 2 + sigma4] - X->Def.Tpow[2 * isite3 - 1 + sigma3]);
  }
  else {//isite3=isite4
    if (sigma3 > sigma4) {
      B_spin = (X->Def.Tpow[2 * isite3 - 2 + sigma3] - X->Def.Tpow[2 * isite4 - 1 + sigma4]);
    }
    else {
      B_spin = (X->Def.Tpow[2 * isite4 - 2 + sigma4] - X->Def.Tpow[2 * isite3 - 1 + sigma3]);
    }
  }
  /**
  Compute mask for hopping (LargeList::isA_spin, LargeList::isB_spin)
  */
  isA_spin = is1_spin + is2_spin;
  isB_spin = is3_spin + is4_spin;

  X->Large.is1_spin = is1_spin;
  X->Large.is2_spin = is2_spin;
  X->Large.is3_spin = is3_spin;
  X->Large.is4_spin = is4_spin;
  X->Large.isA_spin = isA_spin;
  X->Large.isB_spin = isB_spin;
  X->Large.A_spin = A_spin;
  X->Large.B_spin = B_spin;
  /**
  Copy coupling constant (LargeList::tmp_V)
  */
  X->Large.tmp_V = tmp_V;
  X->Large.isite1 = isite1;
  X->Large.isite2 = isite2;
  X->Large.isite3 = isite3;
  X->Large.isite4 = isite4;

  return 0;
}/*int general_int_GetInfo*/
/**
@brief Compute mask for bit operation of pairhop term.
@return Error-code, always return 0
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int pairhopp_GetInfo(
  int iPairHopp,//!<[in] Index of pairhopp interaction
  struct BindStruct *X//!<[inout]
) {
  int isite1 = X->Def.PairHopping[iPairHopp][0] + 1;
  int isite2 = X->Def.PairHopping[iPairHopp][1] + 1;
  /**
  Copy coupling constant (LargeList::tmp_J)
  */
  X->Large.tmp_J = X->Def.ParaPairHopping[iPairHopp];
  /**
  Compute mask for checking occupations of 
  @f$(i_1,\uparrow)@f$ (LargeList::is1_up), @f$(i_1,\downarrow)@f$ (LargeList::is1_down)
  @f$(i_2,\uparrow)@f$ (LargeList::is2_up), @f$(i_2,\downarrow)@f$ (LargeList::is2_down)
  */
  X->Large.is1_up = X->Def.Tpow[2 * isite1 - 2];
  X->Large.is1_down = X->Def.Tpow[2 * isite1 - 1];
  X->Large.is2_up = X->Def.Tpow[2 * isite2 - 2];
  X->Large.is2_down = X->Def.Tpow[2 * isite2 - 1];

  return 0;
}/*int pairhopp_GetInfo*/
/**
@brief Compute mask for bit operation of exchange term.
@return Error-code, always return 0
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int exchange_GetInfo(
  int iExchange,//!<[in] Index of exchange interaction
  struct BindStruct *X//!<[inout]
) {
  int isite1 = X->Def.ExchangeCoupling[iExchange][0] + 1;
  int isite2 = X->Def.ExchangeCoupling[iExchange][1] + 1;
  /**
  Copy coupling constant (LargeList::tmp_J)
  */
  X->Large.tmp_J = -X->Def.ParaExchangeCoupling[iExchange];
  /**
  Compute mask for checking occupations of
  @f$(i_1,\uparrow)@f$ (LargeList::is1_up), @f$(i_1,\downarrow)@f$ (LargeList::is1_down)
  @f$(i_2,\uparrow)@f$ (LargeList::is2_up), @f$(i_2,\downarrow)@f$ (LargeList::is2_down)
  */
  X->Large.is1_up = X->Def.Tpow[2 * isite1 - 2];
  X->Large.is1_down = X->Def.Tpow[2 * isite1 - 1];
  X->Large.is2_up = X->Def.Tpow[2 * isite2 - 2];
  X->Large.is2_down = X->Def.Tpow[2 * isite2 - 1];

  return 0;
}/*int exchange_GetInfo*/

/******************************************************************************/
//[e] GetInfo functions
/******************************************************************************/

/******************************************************************************/
//[s] core routines
/******************************************************************************/
/**
@brief Operation of @f$t c_{i\sigma}^\dagger c_{i\sigma}@f$ (Grandcanonical)
@return Fragment of @f$\langle v_1|{\hat H}|v_1\rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
 */
/**
@brief Mapping core of ::GC_CisAis shared 3-way (original / probe). Diagonal:
returns the occupation (0/1) of @f$(i\sigma)@f$ for bare-bit state j-1.
*/
static long unsigned int GC_CisAis_map(
  long unsigned int j,
  long unsigned int is1_spin
) {
  return ((j - 1) & is1_spin) / is1_spin;
}
double complex GC_CisAis(
  long unsigned int j,//!<[in] Index of element of wavefunction
  double complex *tmp_v0,//!<[inout] Result vector
  double complex *tmp_v1,//!<[in] Input producted vector
  struct BindStruct *X,//!<[inout]
  long unsigned int is1_spin,//!<[in] Mask for occupation of @f$(i \sigma)@f$
  double complex tmp_trans//!<[in] Transfer integral
) {
  long unsigned int A_ibit_tmp;
  double complex dmv;
  double complex dam_pr;

  A_ibit_tmp = GC_CisAis_map(j, is1_spin);
  dmv = tmp_v1[j] * A_ibit_tmp;
  if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
    tmp_v0[j] += dmv * tmp_trans;
  }/*if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC)*/
  dam_pr = dmv * conj(tmp_v1[j]);
  return dam_pr;
}/*double complex GC_CisAis*/
/**
@brief Trace-probe adapter for ::GC_CisAis (grand-canonical diagonal one-body).
Diagonal: always returns 1 with *kprime_out = j-1 and *amp_out = occupation.
*/
int GC_CisAis_TraceProbe(
  long unsigned int j, struct BindStruct *X, long unsigned int is1_spin,
  long int *kprime_out, double complex *amp_out
) {
  (void)X;
  *kprime_out = (long int)j - 1;
  *amp_out = (double complex)GC_CisAis_map(j, is1_spin);
  return 1;
}/*int GC_CisAis_TraceProbe*/
/**
@brief Operation of @f$t c_{i\sigma} c_{i\sigma}^\dagger@f$ (Grandcanonical)
@return Fragment of @f$\langle v_1|{\hat H}|v_1\rangle@f$
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex GC_AisCis(
  long unsigned int j,//!<[in] Index of element of wavefunction
  double complex *tmp_v0,//!<[inout] Result vector
  double complex *tmp_v1,//!<[in] Input producted vector
  struct BindStruct *X,//!<[inout]
  long unsigned int is1_spin,//!<[in] Mask for occupation of @f$(i \sigma)@f$
  double complex tmp_trans//!<[in] Transfer integral
) {
  long unsigned int A_ibit_tmp;
  long unsigned int list_1_j;
  double complex dmv;
  double complex dam_pr;

  list_1_j = j - 1;
  A_ibit_tmp = (list_1_j & is1_spin) / is1_spin;
  dmv = tmp_v1[j] * (1 - A_ibit_tmp);
  if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
    tmp_v0[j] += dmv * tmp_trans;
  }/*if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC)*/
  dam_pr = dmv * conj(tmp_v1[j]);
  return dam_pr;
}/*double complex GC_AisCis*/
/**
@brief @f$c_{is}\\dagger c_{is}@f$ term in Hubbard (canonical) 
@return
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int child_CisAis(
  long unsigned int list_1_j,
  struct BindStruct *X,
  long unsigned int is1_spin
) {
  int A_ibit_tmp;

  // off = j
  A_ibit_tmp = (list_1_j & is1_spin) / is1_spin;
  return A_ibit_tmp;
}
/**
@brief @f$c_{is}^\dagger c_{jt}@f$ term for canonical Hubbard
@return @f$\langle v_1|{\hat H}_{\rm this}|v_1\rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
/**
@brief Mapping core of ::CisAjt shared 3-way (original / probe). Computes the
canonical hopping destination and Fermion sign WITHOUT touching any vector.
@return 1 if the hop survives (kprime_out = 0-based destination in the basis,
i.e. GetOffComp's 1-based off minus 1; sgn_out = +-1); 0 if annihilated or the
exchanged bit pattern is not in the basis (kprime_out = -1).
*/
static int CisAjt_map(
  long unsigned int j, struct BindStruct *X,
  long unsigned int is1_spin, long unsigned int is2_spin,
  long unsigned int sum_spin, long unsigned int diff_spin,
  long int *kprime_out, int *sgn_out
) {
  long unsigned int ibit_tmp_1, ibit_tmp_2;
  long unsigned int bit, iexchg, off;
  int sgn;

  ibit_tmp_1 = (list_1[j] & is1_spin);
  ibit_tmp_2 = (list_1[j] & is2_spin);
  if (ibit_tmp_1 == 0 && ibit_tmp_2 != 0) {
    bit = list_1[j] & diff_spin;
    SgnBit(bit, &sgn); // Fermion sign
    iexchg = list_1[j] ^ sum_spin;
    if(GetOffComp(list_2_1, list_2_2, iexchg, X->Large.irght, X->Large.ilft, X->Large.ihfbit, &off)==FALSE){
      *kprime_out = -1;
      return 0;
    }
    *kprime_out = (long int)off - 1;
    *sgn_out = sgn;
    return 1;
  }
  *kprime_out = -1;
  return 0;
}
double complex CisAjt(
  long unsigned int j,//!<[in] Index of wavefunction
  double complex *tmp_v0,//!<[inout] @f$v_0 = H v_1@f$
  double complex *tmp_v1,//!<[in] Vector to be producted
  struct BindStruct *X,//!<[inout]
  long unsigned int is1_spin,//!<[in] Mask for occupation of (is)
  long unsigned int is2_spin,//!<[in] Mask for occupation of (jt)
  long unsigned int sum_spin,//!<[in] Mask for hopping
  long unsigned int diff_spin,//!<[in] Mask for Fermion sign
  double complex tmp_V//!<[in] Hopping integral
) {
  long int kprime;
  int sgn;
  double complex dmv, dam_pr;

  if (CisAjt_map(j, X, is1_spin, is2_spin, sum_spin, diff_spin, &kprime, &sgn)) {
    dmv = sgn * tmp_v1[j];
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
      tmp_v0[kprime + 1] += tmp_V * dmv;
    }
    dam_pr = dmv * conj(tmp_v1[kprime + 1]);
    return dam_pr;
  }
  else {
    return 0;
  }
}
/**
@brief Trace-probe adapter for ::CisAjt (canonical off-diagonal one-body).
One-body operator: implicit coupling 1.0, so *amp_out is the bare Fermion sign.
*/
int CisAjt_TraceProbe(
  long unsigned int j, struct BindStruct *X,
  long unsigned int is1_spin, long unsigned int is2_spin,
  long unsigned int sum_spin, long unsigned int diff_spin,
  long int *kprime_out, double complex *amp_out
) {
  int sgn;
  if (CisAjt_map(j, X, is1_spin, is2_spin, sum_spin, diff_spin, kprime_out, &sgn)) {
    *amp_out = (double complex)sgn;
    return 1;
  }
  *amp_out = 0.0;
  return 0;
}/*int CisAjt_TraceProbe*/
/**
@brief @f$c_{is}^\dagger c_{jt}@f$ term for grandcanonical Hubbard
@return @f$\langle v_1|{\hat H}_{\rm this}|v_1\rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
/**
@brief Mapping core of ::GC_CisAjt shared 3-way (original / probe). Bare-bit
grand-canonical hop: no vector access. @return 1 if the hop survives
(kprime_out = 0-based destination = bare-bit pattern list_1_off; sgn_out = +-1);
0 if annihilated (kprime_out = -1).
*/
static int GC_CisAjt_map(
  long unsigned int j,
  long unsigned int is1_spin, long unsigned int is2_spin,
  long unsigned int sum_spin, long unsigned int diff_spin,
  long int *kprime_out, int *sgn_out
) {
  long unsigned int list_1_j, list_1_off;
  long unsigned int ibit_tmp_1, ibit_tmp_2;
  long unsigned int bit;
  int sgn;

  list_1_j = j - 1;
  ibit_tmp_1 = (list_1_j & is1_spin);
  ibit_tmp_2 = (list_1_j & is2_spin);
  if (ibit_tmp_1 == 0 && ibit_tmp_2 != 0) {
    bit = list_1_j & diff_spin;
    SgnBit(bit, &sgn); // Fermion sign
    list_1_off = list_1_j ^ sum_spin;
    *kprime_out = (long int)list_1_off;
    *sgn_out = sgn;
    return 1;
  }
  *kprime_out = -1;
  return 0;
}
double complex GC_CisAjt(
  long unsigned int j,//!<[in] Index of wavefunction
  double complex *tmp_v0,//!<[in] @f$v_0 = H v_1@f$
  double complex *tmp_v1,//!<[in]Vector to be producted
  struct BindStruct *X,//!<[inout]
  long unsigned int is1_spin,//!<[in] Mask for occupation of (is)
  long unsigned int is2_spin,//!<[in] Mask for occupation of (jt)
  long unsigned int sum_spin,//!<[in] Mask for hopping
  long unsigned int diff_spin,//!<[in] Mask for Fermion sign
  double complex tmp_V,//!<[in] Hopping
  long unsigned int *tmp_off//!<[in] Index of wavefunction of final state
) {
  long int kprime;
  int sgn;
  double complex dmv, dam_pr;

  *tmp_off = 0;
  if (GC_CisAjt_map(j, is1_spin, is2_spin, sum_spin, diff_spin, &kprime, &sgn)) {
    *tmp_off = (long unsigned int)kprime;
    dmv = sgn * tmp_v1[j];
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
      tmp_v0[kprime + 1] += dmv * tmp_V;
    }
    dam_pr = dmv * conj(tmp_v1[kprime + 1]);
    return dam_pr;
  }
  else {
    return 0;
  }
}/*double complex GC_CisAjt*/
/**
@brief Trace-probe adapter for ::GC_CisAjt (grand-canonical off-diagonal
one-body). One-body operator: implicit coupling 1.0, *amp_out = Fermion sign.
*/
int GC_CisAjt_TraceProbe(
  long unsigned int j, struct BindStruct *X,
  long unsigned int is1_spin, long unsigned int is2_spin,
  long unsigned int sum_spin, long unsigned int diff_spin,
  long int *kprime_out, double complex *amp_out
) {
  int sgn;
  (void)X;
  if (GC_CisAjt_map(j, is1_spin, is2_spin, sum_spin, diff_spin, kprime_out, &sgn)) {
    *amp_out = (double complex)sgn;
    return 1;
  }
  *amp_out = 0.0;
  return 0;
}/*int GC_CisAjt_TraceProbe*/
/**
@brief Compute index of wavefunction of final state
@return Fermion sign
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int child_CisAjt(
  long unsigned int list_1_j,//!<[in] Similer to ::list_1 ?
  struct BindStruct *X,//!<[in]
  long unsigned int is1_spin,//!<[in] Mask for occupation of (is)
  long unsigned int is2_spin,//!<[in] Mask for occupation of (jt)
  long unsigned int sum_spin,//!<[in] Mask for hopping
  long unsigned int diff_spin,//!<[in] Mask for Fermion sign
  long unsigned int *tmp_off//!<[in] Index of wavefunction of final state
) {
  long unsigned int off;
  int sgn = 1;

  sgn = child_GC_CisAjt(list_1_j, X, is1_spin, is2_spin, sum_spin, diff_spin, tmp_off);
  if (sgn != 0) {
    if(GetOffComp(list_2_1, list_2_2, *tmp_off, X->Large.irght, X->Large.ilft, X->Large.ihfbit, &off)!=TRUE){
      *tmp_off = 0;
      return 0;
    }
    *tmp_off = off;
    return sgn;
  }
  else {
    *tmp_off = 0;
    return 0;
  }
}/*int child_CisAjt*/
/**
@brief Compute index of wavefunction of final state
@return Fermion sign
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
int child_GC_CisAjt(
  long unsigned int list_1_j,//!<[in] ::list_1 ?
  struct BindStruct *X,//!<[in]
  long unsigned int is1_spin,//!<[in] Mask for occupation of (is)
  long unsigned int is2_spin,//!<[in] Mask for occupation of (jt)
  long unsigned int sum_spin,//!<[in] Mask for hopping
  long unsigned int diff_spin,//!<[in] Mask for Fermion sign
  long unsigned int *tmp_off//!<[out] Index of wavefunction of final state
) {
  long unsigned int ibit_tmp_1, ibit_tmp_2;
  long unsigned int bit, off;
  int sgn = 1;

  ibit_tmp_1 = (list_1_j & is1_spin);
  ibit_tmp_2 = (list_1_j & is2_spin);

  if (ibit_tmp_1 == 0 && ibit_tmp_2 != 0) {
    bit = list_1_j & diff_spin;
    SgnBit(bit, &sgn); // Fermion sign
    off = list_1_j ^ sum_spin;
    *tmp_off = off;
    return sgn; // pm 1
  }
  else {
    *tmp_off = 0;
    return 0;
  }
}
/******************************************************************************/
//[e] core routines
/******************************************************************************/

/******************************************************************************/
//[s] child element functions
/******************************************************************************/
/**
@brief Compute exchange term of canonical-Hubbard
@return @return @f$\langle v_1|{\hat H}_{\rm this}|v_1\rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex exchange_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  double complex *tmp_v0,//!<[inout] @f$v_0 = H v_1@f$
  double complex *tmp_v1,//!<[in] Vector to be producted
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[off] Index of wavefunction of final state
) {
  long unsigned int off;
  long unsigned int ibit1_up, ibit2_up, ibit1_down, ibit2_down;
  double complex dmv;
  long unsigned int iexchg;
  long unsigned int is1_up = X->Large.is1_up;
  long unsigned int is2_up = X->Large.is2_up;
  long unsigned int is1_down = X->Large.is1_down;
  long unsigned int is2_down = X->Large.is2_down;
  long unsigned int irght = X->Large.irght;
  long unsigned int ilft = X->Large.ilft;
  long unsigned int ihfbit = X->Large.ihfbit;
  double complex tmp_J = X->Large.tmp_J;
  int mode = X->Large.mode;
  double complex dam_pr = 0;

  ibit1_up = list_1[j] & is1_up;
  ibit2_up = list_1[j] & is2_up;
  ibit1_down = list_1[j] & is1_down;
  ibit2_down = list_1[j] & is2_down;

  if (ibit1_up == 0 && ibit1_down != 0 && ibit2_up != 0 && ibit2_down == 0) {
    iexchg = list_1[j] - (is1_down + is2_up);
    iexchg += (is1_up + is2_down);
    if(GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off)!=TRUE){
      return 0;
    }
    *tmp_off = off;
    dmv = tmp_J * tmp_v1[j];
    if (mode == M_MLTPLY) {
      tmp_v0[off] += dmv;
    }
    dam_pr += dmv * conj(tmp_v1[off]);
  }
  else if (ibit1_up != 0 && ibit1_down == 0 && ibit2_up == 0 && ibit2_down != 0) {
    iexchg = list_1[j] - (is1_up + is2_down);
    iexchg += (is1_down + is2_up);
    if(GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off)!=TRUE){
      return 0;
    }
    *tmp_off = off;
    dmv = tmp_J * tmp_v1[j];
    if (mode == M_MLTPLY) {
      tmp_v0[off] += dmv;
    }
    dam_pr += dmv * conj(tmp_v1[off]);
  }
  return dam_pr;
}/*double complex exchange_element*/
/**
@brief Compute pairhopp term of canonical Hubbard system
@return Fragment of @f$\langle v_1 | H_{\rm this} | v_1 \rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex pairhopp_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  long unsigned int off;
  long unsigned int ibit1_up, ibit2_up, ibit1_down, ibit2_down;
  double complex dmv;
  long unsigned int iexchg;
  long unsigned int is1_up = X->Large.is1_up;
  long unsigned int is2_up = X->Large.is2_up;
  long unsigned int is1_down = X->Large.is1_down;
  long unsigned int is2_down = X->Large.is2_down;
  long unsigned int irght = X->Large.irght;
  long unsigned int ilft = X->Large.ilft;
  long unsigned int ihfbit = X->Large.ihfbit;
  double complex tmp_J = X->Large.tmp_J;
  int mode = X->Large.mode;
  double complex dam_pr = 0;

  ibit1_up = list_1[j] & is1_up;
  ibit2_up = list_1[j] & is2_up;
  ibit1_down = list_1[j] & is1_down;
  ibit2_down = list_1[j] & is2_down;

  if (ibit1_up == 0 && ibit1_down == 0 && ibit2_up != 0 && ibit2_down != 0) {
    iexchg = list_1[j] - (is2_up + is2_down);
    iexchg += (is1_up + is1_down);

    if(GetOffComp(list_2_1, list_2_2, iexchg, irght, ilft, ihfbit, &off)!=TRUE){
      return 0;
    }
    *tmp_off = off;
    dmv = tmp_J * tmp_v1[j];
    if (mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
      tmp_v0[off] += dmv;
    }
    dam_pr += dmv * conj(tmp_v1[off]);
  }
  return dam_pr;
}/*double complex pairhopp_element*/
/**
@brief Compute exchange term of grandcanonical Hubbard system
@return Fragment of @f$\langle v_1 | H_{\rm this} | v_1 \rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex GC_exchange_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  long unsigned int ibit1_up, ibit2_up, ibit1_down, ibit2_down;
  double complex dmv;
  long unsigned int iexchg;
  long unsigned int is1_up = X->Large.is1_up;
  long unsigned int is2_up = X->Large.is2_up;
  long unsigned int is1_down = X->Large.is1_down;
  long unsigned int is2_down = X->Large.is2_down;
  long unsigned int list_1_j, list_1_off;
  double complex tmp_J = X->Large.tmp_J;
  int mode = X->Large.mode;
  double complex dam_pr = 0;

  list_1_j = j - 1;
  ibit1_up = list_1_j & is1_up;
  ibit2_up = list_1_j & is2_up;
  ibit1_down = list_1_j & is1_down;
  ibit2_down = list_1_j & is2_down;

  if (ibit1_up == 0 && ibit1_down != 0 && ibit2_up != 0 && ibit2_down == 0) {

    iexchg = list_1_j - (is1_down + is2_up);
    iexchg += (is1_up + is2_down);
    list_1_off = iexchg;
    *tmp_off = list_1_off;

    dmv = tmp_J * tmp_v1[j];
    if (mode == M_MLTPLY) {
      tmp_v0[list_1_off + 1] += dmv;
    }
    dam_pr += dmv * conj(tmp_v1[list_1_off + 1]);
  }
  else if (ibit1_up != 0 && ibit1_down == 0 && ibit2_up == 0 && ibit2_down != 0) {
    iexchg = list_1_j - (is1_up + is2_down);
    iexchg += (is1_down + is2_up);
    list_1_off = iexchg;
    *tmp_off = list_1_off;

    dmv = tmp_J * tmp_v1[j];
    if (mode == M_MLTPLY) {
      tmp_v0[list_1_off + 1] += dmv;
    }
    dam_pr += dmv * conj(tmp_v1[list_1_off + 1]);
  }
  return dam_pr;
}/*double complex GC_exchange_element*/
/**
@brief Compute pairhopp term of grandcanonical Hubbard system
@return Fragment of @f$\langle v_1 | H_{\rm this} | v_1 \rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
double complex GC_pairhopp_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  long unsigned int ibit1_up, ibit2_up, ibit1_down, ibit2_down;
  double complex dmv;
  long unsigned int iexchg;
  long unsigned int is1_up = X->Large.is1_up;
  long unsigned int is2_up = X->Large.is2_up;
  long unsigned int is1_down = X->Large.is1_down;
  long unsigned int is2_down = X->Large.is2_down;
  long unsigned int list_1_j, list_1_off;
  double complex tmp_J = X->Large.tmp_J;
  int mode = X->Large.mode;

  double complex dam_pr = 0 + 0 * I;
  list_1_j = j - 1;

  ibit1_up = list_1_j & is1_up;

  ibit2_up = list_1_j & is2_up;

  ibit1_down = list_1_j & is1_down;

  ibit2_down = list_1_j & is2_down;

  if (ibit1_up == 0 && ibit1_down == 0 && ibit2_up != 0 && ibit2_down != 0) {
    iexchg = list_1_j - (is2_up + is2_down);
    iexchg += (is1_up + is1_down);
    list_1_off = iexchg;
    *tmp_off = list_1_off;
    dmv = tmp_J * tmp_v1[j];
    if (mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) {
      tmp_v0[list_1_off + 1] += dmv;
    }
    dam_pr += dmv * conj(tmp_v1[list_1_off + 1]);
  }
  return dam_pr;
}
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{is}^\dagger c_{is}@f$
term of canonical Hubbard system
@return Fragment of @f$\langle v_1 | H_{\rm this} | v_1 \rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
/**
@brief Mapping core of ::CisAisCisAis_element (diagonal). @return the signed
occupation product (0/1); kprime_out = j-1 (diagonal, always valid).
*/
static int CisAisCisAis_element_map(
  long unsigned int j, long unsigned int isite1, long unsigned int isite3,
  struct BindStruct *X, long int *kprime_out
) {
  int tmp_sgn;
  tmp_sgn = child_CisAis(list_1[j], X, isite3);
  tmp_sgn *= child_CisAis(list_1[j], X, isite1);
  *kprime_out = (long int)j - 1;
  return tmp_sgn;
}
double complex CisAisCisAis_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int isite1,//!<[in] Site 1
  long unsigned int isite3,//!<[in] Site 3
  double complex tmp_V,//!<[in] Coupling constant
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  long int kprime;
  double complex dmv;
  double complex dam_pr = 0 + 0 * I;
  (void)tmp_off;
  tmp_sgn = CisAisCisAis_element_map(j, isite1, isite3, X, &kprime);
  dmv = tmp_V * tmp_v1[j] * tmp_sgn;
  if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
    tmp_v0[kprime + 1] += dmv;
  }
  dam_pr = conj(tmp_v1[kprime + 1]) * dmv;
  return dam_pr;
}/*double complex CisAisCisAis_element*/
/**
@brief Trace-probe adapter for ::CisAisCisAis_element (canonical diagonal
two-body). Diagonal: always returns 1, *amp_out = tmp_V * signed occupation.
*/
int CisAisCisAis_element_TraceProbe(
  long unsigned int j, long unsigned int isite1, long unsigned int isite3,
  double complex tmp_V, struct BindStruct *X,
  long int *kprime_out, double complex *amp_out
) {
  int tmp_sgn = CisAisCisAis_element_map(j, isite1, isite3, X, kprime_out);
  *amp_out = tmp_V * (double complex)tmp_sgn;
  return 1;
}/*int CisAisCisAis_element_TraceProbe*/
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{jt}^\dagger c_{ku}@f$
term of canonical Hubbard system
@return Fragment of @f$\langle v_1 | H_{\rm this} | v_1 \rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
/**
@brief Mapping core of ::CisAisCjtAku_element. Chains child_CisAjt (canonical
1-based destination via *tmp_off) then child_CisAis. @return signed result
(0 = annihilated); kprime_out = *tmp_off-1 on survival, else -1.
*/
static int CisAisCjtAku_element_map(
  long unsigned int j, long unsigned int isite1,
  long unsigned int isite3, long unsigned int isite4,
  long unsigned int Bsum, long unsigned int Bdiff,
  struct BindStruct *X, long unsigned int *tmp_off, long int *kprime_out
) {
  int tmp_sgn;
  tmp_sgn = child_CisAjt(list_1[j], X, isite3, isite4, Bsum, Bdiff, tmp_off);
  if (tmp_sgn != 0) {
    tmp_sgn *= child_CisAis(list_1[*tmp_off], X, isite1);
    if (tmp_sgn != 0) {
      *kprime_out = (long int)(*tmp_off) - 1;
      return tmp_sgn;
    }
  }
  *kprime_out = -1;
  return 0;
}
double complex CisAisCjtAku_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int isite1,//!<[in] Site 1
  long unsigned int isite3,//!<[in] Site 3
  long unsigned int isite4,//!<[in] Site 4
  long unsigned int Bsum,//!<[in] Bit mask for hopping
  long unsigned int Bdiff,//!<[in] Bit mask for Fermion sign
  double complex tmp_V,//!<[in] Coupling constant
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  long int kprime;
  double complex dmv;
  double complex dam_pr = 0 + 0 * I;
  tmp_sgn = CisAisCjtAku_element_map(j, isite1, isite3, isite4, Bsum, Bdiff, X, tmp_off, &kprime);
  if (tmp_sgn != 0) {
    dmv = tmp_V * tmp_v1[j] * tmp_sgn;
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
      tmp_v0[kprime + 1] += dmv;
    }
    dam_pr = conj(tmp_v1[kprime + 1]) * dmv;
  }
  return dam_pr;
}/*double complex CisAisCjtAku_element*/
/**
@brief Trace-probe adapter for ::CisAisCjtAku_element (canonical two-body).
*/
int CisAisCjtAku_element_TraceProbe(
  long unsigned int j, long unsigned int isite1,
  long unsigned int isite3, long unsigned int isite4,
  long unsigned int Bsum, long unsigned int Bdiff,
  double complex tmp_V, struct BindStruct *X,
  long int *kprime_out, double complex *amp_out
) {
  long unsigned int tmp_off;
  int tmp_sgn = CisAisCjtAku_element_map(j, isite1, isite3, isite4, Bsum, Bdiff, X, &tmp_off, kprime_out);
  if (tmp_sgn != 0) {
    *amp_out = tmp_V * (double complex)tmp_sgn;
    return 1;
  }
  *amp_out = 0.0;
  return 0;
}/*int CisAisCjtAku_element_TraceProbe*/
/**
@brief Compute @f$c_{is}^\dagger c_{jt} c_{ku}^\dagger c_{ku}@f$
term of canonical Hubbard system
@return Fragment of @f$\langle v_1 | H_{\rm this} | v_1 \rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
/**
@brief Mapping core of ::CisAjtCkuAku_element. child_CisAis gate then
child_CisAjt (canonical 1-based *tmp_off). @return signed result (0 = dead);
kprime_out = *tmp_off-1 on survival, else -1.
*/
static int CisAjtCkuAku_element_map(
  long unsigned int j, long unsigned int isite1, long unsigned int isite2,
  long unsigned int isite3, long unsigned int Asum, long unsigned int Adiff,
  struct BindStruct *X, long unsigned int *tmp_off, long int *kprime_out
) {
  int tmp_sgn;
  tmp_sgn = child_CisAis(list_1[j], X, isite3);
  if (tmp_sgn != 0) {
    tmp_sgn *= child_CisAjt(list_1[j], X, isite1, isite2, Asum, Adiff, tmp_off);
    if (tmp_sgn != 0) {
      *kprime_out = (long int)(*tmp_off) - 1;
      return tmp_sgn;
    }
  }
  *kprime_out = -1;
  return 0;
}
double complex CisAjtCkuAku_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int isite1,//!<[in] Site 1
  long unsigned int isite2,//!<[in] Site 2
  long unsigned int isite3,//!<[in] Site 3
  long unsigned int Asum,//!<[in] Bit mask for hopping
  long unsigned int Adiff,//!<[in] Bit mask for Fermion sign
  double complex tmp_V,//!<[in] Coupling constant
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  long int kprime;
  double complex dmv;
  double complex dam_pr;
  dam_pr = 0;
  tmp_sgn = CisAjtCkuAku_element_map(j, isite1, isite2, isite3, Asum, Adiff, X, tmp_off, &kprime);
  if (tmp_sgn != 0) {
    dmv = tmp_V * tmp_v1[j] * tmp_sgn;
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
      tmp_v0[kprime + 1] += dmv;
    }
    dam_pr = conj(tmp_v1[kprime + 1]) * dmv;
  }
  return dam_pr;
}/*double complex CisAjtCkuAku_element*/
/**
@brief Trace-probe adapter for ::CisAjtCkuAku_element (canonical two-body).
*/
int CisAjtCkuAku_element_TraceProbe(
  long unsigned int j, long unsigned int isite1, long unsigned int isite2,
  long unsigned int isite3, long unsigned int Asum, long unsigned int Adiff,
  double complex tmp_V, struct BindStruct *X,
  long int *kprime_out, double complex *amp_out
) {
  long unsigned int tmp_off;
  int tmp_sgn = CisAjtCkuAku_element_map(j, isite1, isite2, isite3, Asum, Adiff, X, &tmp_off, kprime_out);
  if (tmp_sgn != 0) {
    *amp_out = tmp_V * (double complex)tmp_sgn;
    return 1;
  }
  *amp_out = 0.0;
  return 0;
}/*int CisAjtCkuAku_element_TraceProbe*/
/**
@brief Compute @f$c_{is}^\dagger c_{jt} c_{ku}^\dagger c_{lv}@f$
term of canonical Hubbard system
@return Fragment of @f$\langle v_1 | H_{\rm this} | v_1 \rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
/**
@brief Mapping core of ::CisAjtCkuAlv_element. child_GC_CisAjt (bare-bit
intermediate tmp_off_1) then child_CisAjt (canonical 1-based *tmp_off_2).
@return signed result (0 = dead); kprime_out = *tmp_off_2-1 on survival, else -1.
*/
static int CisAjtCkuAlv_element_map(
  long unsigned int j, long unsigned int isite1, long unsigned int isite2,
  long unsigned int isite3, long unsigned int isite4,
  long unsigned int Asum, long unsigned int Adiff,
  long unsigned int Bsum, long unsigned int Bdiff,
  struct BindStruct *X, long unsigned int *tmp_off_2, long int *kprime_out
) {
  int tmp_sgn;
  long unsigned int tmp_off_1;
  tmp_sgn = child_GC_CisAjt(list_1[j], X, isite3, isite4, Bsum, Bdiff, &tmp_off_1);
  if (tmp_sgn != 0) {
    tmp_sgn *= child_CisAjt(tmp_off_1, X, isite1, isite2, Asum, Adiff, tmp_off_2);
    if (tmp_sgn != 0) {
      *kprime_out = (long int)(*tmp_off_2) - 1;
      return tmp_sgn;
    }
  }
  *kprime_out = -1;
  return 0;
}
double complex CisAjtCkuAlv_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int isite1,//!<[in] Site 1
  long unsigned int isite2,//!<[in] Site 2
  long unsigned int isite3,//!<[in] Site 3
  long unsigned int isite4,//!<[in] Site 4
  long unsigned int Asum,//!<[in] Bit mask for hopping
  long unsigned int Adiff,//!<[in] Bit mask for Fermion sign
  long unsigned int Bsum,//!<[in] Bit mask for hopping
  long unsigned int Bdiff,//!<[in] Bit mask for Fermion sign
  double complex tmp_V,//!<[in] Coupling constant
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off_2//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  long int kprime;
  double complex dmv;
  double complex dam_pr = 0;
  tmp_sgn = CisAjtCkuAlv_element_map(j, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff, X, tmp_off_2, &kprime);
  if (tmp_sgn != 0) {
    dmv = tmp_V * tmp_v1[j] * tmp_sgn;
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
      tmp_v0[kprime + 1] += dmv;
    }
    dam_pr = conj(tmp_v1[kprime + 1]) * dmv;
  }
  return dam_pr;
}/*double complex CisAjtCkuAlv_element*/
/**
@brief Trace-probe adapter for ::CisAjtCkuAlv_element (canonical two-body).
*/
int CisAjtCkuAlv_element_TraceProbe(
  long unsigned int j, long unsigned int isite1, long unsigned int isite2,
  long unsigned int isite3, long unsigned int isite4,
  long unsigned int Asum, long unsigned int Adiff,
  long unsigned int Bsum, long unsigned int Bdiff,
  double complex tmp_V, struct BindStruct *X,
  long int *kprime_out, double complex *amp_out
) {
  long unsigned int tmp_off_2;
  int tmp_sgn = CisAjtCkuAlv_element_map(j, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff, X, &tmp_off_2, kprime_out);
  if (tmp_sgn != 0) {
    *amp_out = tmp_V * (double complex)tmp_sgn;
    return 1;
  }
  *amp_out = 0.0;
  return 0;
}/*int CisAjtCkuAlv_element_TraceProbe*/
//[s] Grand Canonical
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{is}^\dagger c_{is}@f$
term of grandcanonical Hubbard system
@return Fragment of @f$\langle v_1 | H_{\rm this} | v_1 \rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
/**
@brief Mapping core of ::GC_CisAisCisAis_element (grand-canonical diagonal).
@return signed occupation product (0/1); kprime_out = j-1 (always valid).
*/
static int GC_CisAisCisAis_element_map(
  long unsigned int j, long unsigned int isite1, long unsigned int isite3,
  struct BindStruct *X, long int *kprime_out
) {
  int tmp_sgn;
  tmp_sgn = child_CisAis(j - 1, X, isite3);
  tmp_sgn *= child_CisAis(j - 1, X, isite1);
  *kprime_out = (long int)j - 1;
  return tmp_sgn;
}
double complex GC_CisAisCisAis_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int isite1,//!<[in] Site 1
  long unsigned int isite3,//!<[in] Site 3
  double complex tmp_V,//!<[in] Coupling constant
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  long int kprime;
  double complex dmv = 0.0;
  double complex dam_pr = 0;
  tmp_sgn = GC_CisAisCisAis_element_map(j, isite1, isite3, X, &kprime);
  if (tmp_sgn != 0) {
    dmv = tmp_V * tmp_v1[j] * tmp_sgn;
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
      tmp_v0[kprime + 1] += dmv;
    }
    dam_pr = conj(tmp_v1[kprime + 1]) * dmv;
  }
  return dam_pr;
}/*double complex GC_CisAisCisAis_element*/
/**
@brief Trace-probe adapter for ::GC_CisAisCisAis_element (grand-canonical
diagonal two-body). Diagonal: always returns 1, *amp_out = tmp_V * signed occ.
*/
int GC_CisAisCisAis_element_TraceProbe(
  long unsigned int j, long unsigned int isite1, long unsigned int isite3,
  double complex tmp_V, struct BindStruct *X,
  long int *kprime_out, double complex *amp_out
) {
  int tmp_sgn = GC_CisAisCisAis_element_map(j, isite1, isite3, X, kprime_out);
  *amp_out = tmp_V * (double complex)tmp_sgn;
  return 1;
}/*int GC_CisAisCisAis_element_TraceProbe*/
/**
@brief Compute @f$c_{is}^\dagger c_{is} c_{jt}^\dagger c_{ku}@f$
term of grandcanonical Hubbard system
@return Fragment of @f$\langle v_1 | H_{\rm this} | v_1 \rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
/**
@brief Mapping core of ::GC_CisAisCjtAku_element. child_GC_CisAjt (bare-bit
*tmp_off) gated by child_CisAis. @return signed result (0 = dead);
kprime_out = *tmp_off on survival, else -1.
*/
static int GC_CisAisCjtAku_element_map(
  long unsigned int j, long unsigned int isite1,
  long unsigned int isite3, long unsigned int isite4,
  long unsigned int Bsum, long unsigned int Bdiff,
  struct BindStruct *X, long unsigned int *tmp_off, long int *kprime_out
) {
  int tmp_sgn;
  tmp_sgn = child_GC_CisAjt((j - 1), X, isite3, isite4, Bsum, Bdiff, tmp_off);
  if (tmp_sgn != 0) {
    tmp_sgn *= child_CisAis(*tmp_off, X, isite1);
    if (tmp_sgn != 0) {
      *kprime_out = (long int)(*tmp_off);
      return tmp_sgn;
    }
  }
  *kprime_out = -1;
  return 0;
}
double complex GC_CisAisCjtAku_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int isite1,//!<[in] Site 1
  long unsigned int isite3,//!<[in] Site 3
  long unsigned int isite4,//!<[in] Site 4
  long unsigned int Bsum,//!<[in] Bit mask for hopping
  long unsigned int Bdiff,//!<[in] Bit mask for Fermion sign
  double complex tmp_V,//!<[in] Coupling constant
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  long int kprime;
  double complex dmv;
  double complex dam_pr = 0 + 0 * I;
  tmp_sgn = GC_CisAisCjtAku_element_map(j, isite1, isite3, isite4, Bsum, Bdiff, X, tmp_off, &kprime);
  if (tmp_sgn != 0) {
    dmv = tmp_V * tmp_v1[j] * tmp_sgn;
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
      tmp_v0[kprime + 1] += dmv;
    }
    dam_pr = conj(tmp_v1[kprime + 1]) * dmv;
  }
  return dam_pr;
}/*double complex GC_CisAisCjtAku_element*/
/**
@brief Trace-probe adapter for ::GC_CisAisCjtAku_element (grand-canonical
two-body).
*/
int GC_CisAisCjtAku_element_TraceProbe(
  long unsigned int j, long unsigned int isite1,
  long unsigned int isite3, long unsigned int isite4,
  long unsigned int Bsum, long unsigned int Bdiff,
  double complex tmp_V, struct BindStruct *X,
  long int *kprime_out, double complex *amp_out
) {
  long unsigned int tmp_off;
  int tmp_sgn = GC_CisAisCjtAku_element_map(j, isite1, isite3, isite4, Bsum, Bdiff, X, &tmp_off, kprime_out);
  if (tmp_sgn != 0) {
    *amp_out = tmp_V * (double complex)tmp_sgn;
    return 1;
  }
  *amp_out = 0.0;
  return 0;
}/*int GC_CisAisCjtAku_element_TraceProbe*/
/**
@brief Compute @f$c_{is}^\dagger c_{jt} c_{ku}^\dagger c_{ku}@f$
term of grandcanonical Hubbard system
@return Fragment of @f$\langle v_1 | H_{\rm this} | v_1 \rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
/**
@brief Mapping core of ::GC_CisAjtCkuAku_element. child_CisAis gate then
child_GC_CisAjt (bare-bit *tmp_off). @return signed result (0 = dead);
kprime_out = *tmp_off on survival, else -1.
*/
static int GC_CisAjtCkuAku_element_map(
  long unsigned int j, long unsigned int isite1, long unsigned int isite2,
  long unsigned int isite3, long unsigned int Asum, long unsigned int Adiff,
  struct BindStruct *X, long unsigned int *tmp_off, long int *kprime_out
) {
  int tmp_sgn;
  tmp_sgn = child_CisAis((j - 1), X, isite3);
  if (tmp_sgn != 0) {
    tmp_sgn *= child_GC_CisAjt((j - 1), X, isite1, isite2, Asum, Adiff, tmp_off);
    if (tmp_sgn != 0) {
      *kprime_out = (long int)(*tmp_off);
      return tmp_sgn;
    }
  }
  *kprime_out = -1;
  return 0;
}
double complex GC_CisAjtCkuAku_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int isite1,//!<[in] Site 1
  long unsigned int isite2,//!<[in] Site 2
  long unsigned int isite3,//!<[in] Site 3
  long unsigned int Asum,//!<[in] Bit mask for hopping
  long unsigned int Adiff,//!<[in] Bit mask for Fermion sign
  double complex tmp_V,//!<[in] Coupling constant
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  long int kprime;
  double complex dmv;
  double complex dam_pr = 0 + 0 * I;
  tmp_sgn = GC_CisAjtCkuAku_element_map(j, isite1, isite2, isite3, Asum, Adiff, X, tmp_off, &kprime);
  if (tmp_sgn != 0) {
    dmv = tmp_V * tmp_v1[j] * tmp_sgn;
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
      tmp_v0[kprime + 1] += dmv;
    }
    dam_pr = conj(tmp_v1[kprime + 1]) * dmv;
  }
  return dam_pr;
}/*double complex GC_CisAjtCkuAku_element*/
/**
@brief Trace-probe adapter for ::GC_CisAjtCkuAku_element (grand-canonical
two-body).
*/
int GC_CisAjtCkuAku_element_TraceProbe(
  long unsigned int j, long unsigned int isite1, long unsigned int isite2,
  long unsigned int isite3, long unsigned int Asum, long unsigned int Adiff,
  double complex tmp_V, struct BindStruct *X,
  long int *kprime_out, double complex *amp_out
) {
  long unsigned int tmp_off;
  int tmp_sgn = GC_CisAjtCkuAku_element_map(j, isite1, isite2, isite3, Asum, Adiff, X, &tmp_off, kprime_out);
  if (tmp_sgn != 0) {
    *amp_out = tmp_V * (double complex)tmp_sgn;
    return 1;
  }
  *amp_out = 0.0;
  return 0;
}/*int GC_CisAjtCkuAku_element_TraceProbe*/
/**
@brief Compute @f$c_{is}^\dagger c_{jt} c_{ku}^\dagger c_{lv}@f$
term of grandcanonical Hubbard system
@return Fragment of @f$\langle v_1 | H_{\rm this} | v_1 \rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
*/
/**
@brief Mapping core of ::GC_CisAjtCkuAlv_element. Two chained child_GC_CisAjt
(intermediate bare-bit tmp_off_1, final bare-bit *tmp_off_2). @return signed
result (0 = dead); kprime_out = *tmp_off_2 on survival, else -1.
*/
static int GC_CisAjtCkuAlv_element_map(
  long unsigned int j, long unsigned int isite1, long unsigned int isite2,
  long unsigned int isite3, long unsigned int isite4,
  long unsigned int Asum, long unsigned int Adiff,
  long unsigned int Bsum, long unsigned int Bdiff,
  struct BindStruct *X, long unsigned int *tmp_off_2, long int *kprime_out
) {
  int tmp_sgn;
  long unsigned int tmp_off_1;
  tmp_sgn = child_GC_CisAjt((j - 1), X, isite3, isite4, Bsum, Bdiff, &tmp_off_1);
  if (tmp_sgn != 0) {
    tmp_sgn *= child_GC_CisAjt(tmp_off_1, X, isite1, isite2, Asum, Adiff, tmp_off_2);
    if (tmp_sgn != 0) {
      *kprime_out = (long int)(*tmp_off_2);
      return tmp_sgn;
    }
  }
  *kprime_out = -1;
  return 0;
}
double complex GC_CisAjtCkuAlv_element(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int isite1,//!<[in] Site 1
  long unsigned int isite2,//!<[in] Site 2
  long unsigned int isite3,//!<[in] Site 3
  long unsigned int isite4,//!<[in] Site 4
  long unsigned int Asum,//!<[in] Bit mask for hopping
  long unsigned int Adiff,//!<[in] Bit mask for Fermion sign
  long unsigned int Bsum,//!<[in] Bit mask for hopping
  long unsigned int Bdiff,//!<[in] Bit mask for Fermion sign
  double complex tmp_V,//!<[in] Coupling constant
  double complex *tmp_v0,//!<[inout] Resulting wavefunction
  double complex *tmp_v1,//!<[in] Wavefunction to be multiplied
  struct BindStruct *X,//!<[inout]
  long unsigned int *tmp_off_2//!<[out] Index of final wavefunction
) {
  int tmp_sgn;
  long int kprime;
  double complex dmv;
  double complex dam_pr = 0 + 0 * I;

  tmp_sgn = GC_CisAjtCkuAlv_element_map(j, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff, X, tmp_off_2, &kprime);
  if (tmp_sgn != 0) {
    dmv = tmp_V * tmp_v1[j] * tmp_sgn;
    if (X->Large.mode == M_MLTPLY || X->Large.mode == M_CALCSPEC) { // for multply
      tmp_v0[kprime + 1] += dmv;
    }
    dam_pr = conj(tmp_v1[kprime + 1]) * dmv;
  }
  return dam_pr;
}/*double complex GC_CisAjtCkuAlv_element*/
/**
@brief Trace-probe adapter for ::GC_CisAjtCkuAlv_element (grand-canonical
two-body).
*/
int GC_CisAjtCkuAlv_element_TraceProbe(
  long unsigned int j, long unsigned int isite1, long unsigned int isite2,
  long unsigned int isite3, long unsigned int isite4,
  long unsigned int Asum, long unsigned int Adiff,
  long unsigned int Bsum, long unsigned int Bdiff,
  double complex tmp_V, struct BindStruct *X,
  long int *kprime_out, double complex *amp_out
) {
  long unsigned int tmp_off_2;
  int tmp_sgn = GC_CisAjtCkuAlv_element_map(j, isite1, isite2, isite3, isite4, Asum, Adiff, Bsum, Bdiff, X, &tmp_off_2, kprime_out);
  if (tmp_sgn != 0) {
    *amp_out = tmp_V * (double complex)tmp_sgn;
    return 1;
  }
  *amp_out = 0.0;
  return 0;
}/*int GC_CisAjtCkuAlv_element_TraceProbe*/
//[e] Grand Canonical
/**
@brief Compute @f$c_{is}^\dagger@f$
term of grandcanonical Hubbard system
@return Fragment of @f$\langle v_1 | H_{\rm this} | v_1 \rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
@author Youhei Yamaji (The University of Tokyo)
*/
double complex GC_Cis(
  long unsigned int j,//!<[in] Index of initial wavefunction
  double complex *tmp_v0,//!<[in] Resulting wavefunction
  double complex *tmp_v1,//!<[in] Wavefunction to be multiplied
  long unsigned int is1_spin,//!<[in] Bit mask 
  double complex tmp_V,//!<[in] Coupling constant
  long unsigned int *tmp_off//!<[in] Index of final wavefunction
) {
  long unsigned int list_1_j, list_1_off;
  long unsigned int ibit_tmp_1;
  long unsigned int bit;
  int sgn, ipsgn;
  double complex dmv, dam_pr;

  list_1_j = j - 1;

  ibit_tmp_1 = (list_1_j & is1_spin);
  // is1_spin >= 1
  // is1_spin = Tpow[2*isite + ispin]

  *tmp_off = 0;

  if (ibit_tmp_1 == 0) {
    // able to create an electron at the is1_spin state
    bit = list_1_j - (list_1_j & (2 * is1_spin - 1));
    SgnBit(bit, &sgn); // Fermion sign
    ipsgn = 1;
#ifdef MPI
    SgnBit(myrank, &ipsgn); // Fermion sign
#endif
    list_1_off = list_1_j | is1_spin; // OR
    *tmp_off = list_1_off;
    dmv = ipsgn * sgn * tmp_v1[j];
    //if (X->Large.mode == M_MLTPLY) { // for multply
    tmp_v0[list_1_off + 1] += dmv * tmp_V;
    //}
    dam_pr = dmv * conj(tmp_v1[list_1_off + 1]);
    return dam_pr;
  }
  else {
    return 0;
  }
}/*double complex GC_Cis*/
/**
@brief Compute @f$c_{jt}@f$
term of grandcanonical Hubbard system
@return Fragment of @f$\langle v_1 | H_{\rm this} | v_1 \rangle@f$
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
@author Youhei Yamaji (The University of Tokyo)
*/
double complex GC_Ajt(
  long unsigned int j,//!<[in] Index of initial wavefunction
  double complex *tmp_v0,//!<[in] Resulting wavefunction
  double complex *tmp_v1,//!<[in] Wavefunction to be multiplied
  long unsigned int is1_spin,//!<[in] Bit mask
  double complex tmp_V,//!<[in] Coupling constant
  long unsigned int *tmp_off//!<[in] Index of final wavefunction
) {
  long unsigned int list_1_j, list_1_off;
  long unsigned int ibit_tmp_1;
  long unsigned int bit;
  int sgn, ipsgn;
  double complex dmv, dam_pr;

  list_1_j = j - 1;

  ibit_tmp_1 = (list_1_j & is1_spin);
  // is1_spin >= 1

  *tmp_off = 0;

  if (ibit_tmp_1 == is1_spin) {
    // able to create an electron at the is1_spin state
    bit = list_1_j - (list_1_j & (2 * is1_spin - 1));
    SgnBit(bit, &sgn); // Fermion sign
    ipsgn = 1;
#ifdef MPI
    SgnBit(myrank, &ipsgn); // Fermion sign
#endif
    list_1_off = list_1_j ^ is1_spin;
    *tmp_off = list_1_off;
    dmv = ipsgn * sgn * tmp_v1[j];
    //if (X->Large.mode == M_MLTPLY) { // for multply
    tmp_v0[list_1_off + 1] += dmv * tmp_V;
    //}
    dam_pr = dmv * conj(tmp_v1[list_1_off + 1]);
    return dam_pr;
  }
  else {
    return 0;
  }
}/*double complex GC_Ajt*/
/**
@brief Compute index of final wavefunction associatesd to @f$c_{is}^\dagger@f$
term of canonical Hubbard system
@return 1 if electron (is) absent, 0 if not.
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
@author Youhei Yamaji (The University of Tokyo)
*/
int child_Cis(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int is1_spin,//!<[in] Bit mask
  long unsigned int *tmp_off,//!<[out] Index of final wavefunction
  long unsigned int *list_1_org,//!<[in] Similar to ::list_1
  long unsigned int *list_2_1_target,//!<[in] Similar to ::list_2_1
  long unsigned int *list_2_2_target,//!<[in] Similar to ::list_2_2
  long unsigned int _irght,//!<[in] Similar to LargeList::irght
  long unsigned int _ilft,//!<[in] Similar to LargeList::ilft
  long unsigned int _ihfbit//!<[in] Similar to LargeList::ihfbit
) {
  long unsigned int list_1_j, list_1_off;
  long unsigned int ibit_tmp_1;
  long unsigned int bit;
  int sgn, ipsgn;

  list_1_j = list_1_org[j];

  ibit_tmp_1 = (list_1_j & is1_spin);
  // is1_spin >= 1
  // is1_spin = Tpow[2*isite + ispin]

  *tmp_off = 0;

  if (ibit_tmp_1 == 0) {
    // able to create an electron at the is1_spin state
    bit = list_1_j - (list_1_j & (2 * is1_spin - 1));
    SgnBit(bit, &sgn); // Fermion sign
    ipsgn = 1;
#ifdef MPI
    SgnBit(myrank, &ipsgn); // Fermion sign
#endif
    list_1_off = list_1_j | is1_spin; // OR

    if(GetOffComp(list_2_1_target, list_2_2_target, list_1_off, _irght, _ilft, _ihfbit, tmp_off)!=TRUE){
      *tmp_off=0;
      return 0;
    }
    sgn *= ipsgn;
    return (sgn);
  }
  else {
    *tmp_off = 0;
    return 0;
  }
}/*int child_Cis*/
/**
@brief Compute index of final wavefunction associatesd to @f$c_{jt}@f$
term of canonical Hubbard system
@return 1 if electron (jt) exist, 0 if not.
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
@author Youhei Yamaji (The University of Tokyo)
*/
double complex child_Ajt(
  long unsigned int j,//!<[in] Index of initial wavefunction
  long unsigned int is1_spin,//!<[in] Bit mask
  long unsigned int *tmp_off,//!<[out] Index of final wavefunction
  long unsigned int *list_1_org,//!<[in] Similar to ::list_1
  long unsigned int *list_2_1_target,//!<[in] Similar to ::list_2_1
  long unsigned int *list_2_2_target,//!<[in] Similar to ::list_2_2
  long unsigned int _irght,//!<[in] Similar to LargeList::irght
  long unsigned int _ilft,//!<[in] Similar to LargeList::ilft
  long unsigned int _ihfbit//!<[in] Similar to LargeList::ihfbit
) {
  long unsigned int list_1_j, list_1_off;
  long unsigned int ibit_tmp_1;
  long unsigned int bit;
  int sgn, ipsgn;

  list_1_j = list_1_org[j];

  ibit_tmp_1 = (list_1_j & is1_spin);
// is1_spin >= 1
// is1_spin = Tpow[2*isite + ispin]

  *tmp_off = 0;
  if (ibit_tmp_1 != 0) {
    // able to delete an electron at the is1_spin state
    bit = list_1_j - (list_1_j & (2 * is1_spin - 1));
    SgnBit(bit, &sgn); // Fermion sign
    ipsgn = 1;
#ifdef MPI
    SgnBit(myrank, &ipsgn); // Fermion sign
#endif
    list_1_off = list_1_j ^ is1_spin;
    if(GetOffComp(list_2_1_target, list_2_2_target, list_1_off, _irght, _ilft, _ihfbit, tmp_off)!=TRUE){
      *tmp_off=0;
      return 0;
    }
    sgn *= ipsgn;
    return(sgn);
  }
  else {
    *tmp_off = 0;
    return 0;
  }
}/*double complex child_Ajt*/

/******************************************************************************/
//[e] child element functions
/******************************************************************************/
