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
#include "wrapperMPI.h"

/**
 * @file   bitcalc.c
 * @version 0.1, 0.2
 * @author Takahiro Misawa (The University of Tokyo) 
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 * 
 * @brief  File for giving functions of treating bits on the target of Hilbert space.
 * 
 * 
 */


/** 
 * 
 * @brief function of getting right, left and half bits corresponding to a original Hilbert space.
 }
 * @param Nsite a total number of sites
 * @param irght a bit to split original Hilbert space into right space @f$2^{(Ns+2)/2}-1@f$
 * @param ilft a bit to split original Hilbert space into left space 
 * @param ihfbit a half bit to split original Hilbert space @f$2^{(Ns+2)/2}@f$
 * 
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo) 
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 * @return 
 */
int GetSplitBit(
                const int Nsite, //!<[in]
                long unsigned int *irght, //!<[out]
                long unsigned int *ilft, //!<[out]
                long unsigned int *ihfbit//!<[out]
){
  if(Nsite<1){
    fprintf(stderr, "%s", cErrSiteNumber);
    return -1;
  }
  *ihfbit=1;
  *ihfbit=(*ihfbit<<(unsigned long int)((Nsite+1)/2));
  *irght = *ihfbit-1;
  *ilft=1;
  *ilft = (*ilft<<(unsigned long int)Nsite)-1;
  *ilft= *ilft ^ *irght;
  return 0;
}

/** 
 * @brief function of splitting original bit into right and left  spaces.
 * 
 * @param Nsite a total number of sites
 * @param iCalcModel Calc model defined in CalcMode file
 * @param irght a bit to split original  space into right space
 * @param ilft a bit to split original  space into left space
 * @param ihfbit a half bit to split original  space
 * 
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo) 
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 * @return 
 */
int GetSplitBitByModel(
                       const int Nsite,  //!<[in]
                       const int iCalcModel,  //!<[in]
                       long unsigned int *irght,  //!<[out]
                       long unsigned int *ilft,  //!<[out]
                       long unsigned int *ihfbit  //!<[out]
                       )
{
  int tmpNsite=Nsite;
  switch(iCalcModel){    
  case HubbardGC:
  case KondoGC:
  case HubbardNConserved:
  case Hubbard:
  case Kondo:
    tmpNsite *= 2;
    break;
  case Spin:
  case SpinGC:   
    break;
  default:
    fprintf(stderr, cErrNoModel, iCalcModel);
    return -1;
  }

  if(GetSplitBit(tmpNsite, irght, ilft, ihfbit)!=0){
    return -1;
  }
  
  return 0;
}


/** 
 * 
 * @brief function of getting right, left and half bits corresponding to a original  space.
 * @param Nsite a total number of sites
 * @param ihfbit a bit to split original  space
 * 
 * @retval 0 normally finished
 * @retval -1 unnormally finished
 *
 * @version 0.2
 * @author Takahiro Misawa (The University of Tokyo) 
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 */
int GetSplitBitForGeneralSpin(
                const int Nsite,  //!<[in]
                long unsigned int *ihfbit, //!<[out]
                const long int *SiteToBit  //!<[in]
){
  int isite=0;
  long int isqrtMaxDim=1;
  long int tmpbit=1;
  
  if(Nsite<1){
    fprintf(stderr, "%s", cErrSiteNumber);
    return -1;
  }
  
  for(isite=1; isite<=Nsite ; isite++){
    isqrtMaxDim *= SiteToBit[isite-1];
  }
  isqrtMaxDim =(long int)sqrt(isqrtMaxDim);

  for(isite=1; isite<=Nsite ; isite++){
    tmpbit *= SiteToBit[isite-1];
    if(tmpbit >= isqrtMaxDim) break;
  }
  *ihfbit=tmpbit;
  return 0;
}


/** 
 * 
 * @brief function of splitting a original bit to right and left spaces
 * 
 * @param ibit a original bit
 * @param irght a bit to split original  space into right space
 * @param ilft a bit to split original  space into left space
 * @param ihfbit a half bit to split original  space
 * @param isplited_Bit_right a splitted bit reflected on right space
 * @param isplited_Bit_left a splitted bit reflected on left space
 * @version 0.1 
 * @author Takahiro Misawa (The University of Tokyo) 
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 */
void SplitBit(
              const long unsigned int ibit,  //!<[in]
              const long unsigned int irght,  //!<[in]
              const long unsigned int ilft,  //!<[in]
              const long unsigned int ihfbit,  //!<[in]
              long unsigned int *isplited_Bit_right,  //!<[out]
              long unsigned int *isplited_Bit_left  //!<[out]
)
{
  *isplited_Bit_right=ibit & irght;
  *isplited_Bit_left=ibit & ilft;
  *isplited_Bit_left=*isplited_Bit_left/ihfbit;
}

/** 
 * 
 * @brief function of getting off-diagonal component
 * 
 * @param _list_2_1 list to right space
 * @param _list_2_2 list to left space
 * @param _ibit a original bit 
 * @param _irght a bit to split original  space into right space
 * @param _ilft a bit to split original  space into left space
 * @param _ihfbit a half bit to split original  space
 * @param _ioffComp an off diagonal component
 * @version 0.1 
 * @author Takahiro Misawa (The University of Tokyo) 
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 */
int GetOffComp(
               long unsigned int *_list_2_1,  //!<[in]
               long unsigned int *_list_2_2, //!<[in]
               long unsigned int _ibit,  //!<[in]
               const long unsigned int _irght, //!<[in]
               const long unsigned int _ilft,  //!<[in]
               const long unsigned int _ihfbit,  //!<[in]
               long unsigned int *_ioffComp  //!<[out]
)
{
  long unsigned int ia, ib;
  SplitBit(_ibit, _irght, _ilft, _ihfbit, &ia, &ib);
/*
  *_ioffComp =_list_2_1[ia];
  *_ioffComp+=_list_2_2[ib];
*/

  //if(myrank==1)
  //printf( "DEGBUG:_ibit=%ld, _list_2_1=%ld, _list_2_2=%ld\n", _ibit, _list_2_1[ia], _list_2_2[ib]);

  if(_list_2_1[ia]*_list_2_2[ib]==0){
    *_ioffComp=0;
    return FALSE;
  }
  *_ioffComp =_list_2_1[ia]-1;
  *_ioffComp+=_list_2_2[ib]-1;

  return TRUE;
}


/** 
 * 
 * @brief function of getting off-diagonal component for general spin
 * 
 * @param org_ibit a original bit
 * @param org_isite a target site 
 * @param org_ispin a target spin to delete.
 * @param off_ispin a target spin to create.
 * @param _ioffComp a generated bit 
 * @param SiteToBit List for getting bit at a site
 * @param Tpow List for getting total bit at a site before
 * @retval FALSE off-diagonal component does not exist
 * @retval TRUE off-diagonal component exists
 * 
 * @version 0.2
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 */
int GetOffCompGeneralSpin(
                          const long unsigned int org_ibit,  //!<[in]
                          const int org_isite,  //!<[in]
                          const int org_ispin, //!<[in]
                          const int off_ispin, //!<[in]
                          long  unsigned int *_ioffComp,  //!<[out]
                          const long int *SiteToBit,  //!<[in]
                          const long unsigned int *Tpow  //!<[in]
)
{
  if(off_ispin>SiteToBit[org_isite-1]-1 ||
     off_ispin<0                      ||
     org_ispin>SiteToBit[org_isite-1]-1 ||
     org_ispin <0){
    *_ioffComp=0;
    return FALSE;
  }
  if(BitCheckGeneral(org_ibit, org_isite, org_ispin, SiteToBit, Tpow) == FALSE){
    *_ioffComp=0;
    return FALSE;
  }
  
  //delete org_ispin and create off_ispin
  long int tmp_off=0;
  tmp_off=(long int)(off_ispin-org_ispin);
  tmp_off *=Tpow[org_isite-1];
  tmp_off +=org_ibit;
  *_ioffComp =tmp_off;
  return TRUE;
}

/** 
 * 
 * @brief function of converting component to list_1
 * 
 * @param org_ibit a original bit
 * @param ihlfbit a split bit for general spin
 * @param _ilist1Comp a component converted to list_1
 * 
 * @version 0.2
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 */
int ConvertToList1GeneralSpin(
    const long unsigned int org_ibit,  //!<[in]
    const long unsigned int ihlfbit,    //!<[in]
    long unsigned int *_ilist1Comp     //!<[out]
)
{
  long unsigned int ia, ib;
  ia=org_ibit%ihlfbit;
  ib=org_ibit/ihlfbit;
  if(list_2_1[ia]*list_2_2[ib]==0){
    *_ilist1Comp=0;
    return FALSE;
  }
  *_ilist1Comp = list_2_1[ia] + list_2_2[ib] - 2;
  return TRUE;
}

/** 
 * 
 * @brief function of getting fermion signs (for 32bit)
 * 
 * @param org_bit an original bit
 * @param sgn fermion sign 
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo) 
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 */
void SgnBit_old( 
      const long unsigned int org_bit,  //!<[in]
                  int *sgn  //!<[out]
)
{
   long unsigned int bit;

   bit     = org_bit^(org_bit>>1);
   bit     = (bit^(bit>>2) ) & 0x11111111;
   bit     = bit*0x11111111;
   *sgn    = 1-2*((bit>>28) & 1); // sgn = pm 1
}


// for 64 bit
/** 
 * 
 * @brief function of getting fermion sign (64 bit)
 * 
 * @param org_bit an original bit
 * @param sgn fermion sign 
 * @version 0.1
 *
 * @author Takahiro Misawa (The University of Tokyo) 
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 */
void SgnBit( 
      const long unsigned int org_bit,  //!<[in]
                  int *sgn  //!<[out]
)
{
   long unsigned int bit;

   bit =  org_bit^(org_bit>>1);
   bit =  bit^(bit>>2);
   bit =  bit^(bit>>4);
   bit =  bit^(bit>>8);
   bit =  bit^(bit>>16);
   bit =  bit^(bit>>32);
   *sgn    = 1-2*(bit & 1); // sgn = pm 1
}

/** 
 * 
 * @brief bit check function
 *
 * @param org_bit original bit to check
 * @param target_bit target bit to check
 * @retval 1
 * @retval 0 
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo) 
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 */
int BitCheck( 
       const long unsigned int org_bit,  //!<[in]
       const long unsigned int target_bit  //!<[in]
)
{
   return  (org_bit >> target_bit) &1;
   // (org_bit & (2^target_bit))/2^target_bit
}



/** 
 * 
 * @brief bit check function for general spin 
 *
 * @param org_bit original bit to check
 * @param org_isite site index (org_isite >= 1)
 * @param target_ispin target spin to check 
 * @param SiteToBit List for getting bit at a site
 * @param Tpow List for getting total bit at a site before
 * @retval 0 bit does not exists
 * @retval 1 bit exists
 * 
 * @version 0.2
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 */
int BitCheckGeneral(
       const long unsigned int org_bit,  //!<[in]
       const unsigned int org_isite,  //!<[in]
       const unsigned int target_ispin,  //!<[in]
       const long int *SiteToBit,  //!<[in]
       const long unsigned int *Tpow  //!<[in]
)
{

  if(GetBitGeneral(org_isite, org_bit, SiteToBit, Tpow) !=target_ispin){
    return FALSE;
  }
  return TRUE;
}


/** 
 * 
 * @brief get bit at a site for general spin  
 *
 * @param isite site index (isite >= 1)
 * @param org_bit original bit to check 
 * @param SiteToBit List for getting bit at a site
 * @param Tpow List for getting total bit at a site before
 * @return bit at a site
 * 
 * @version 0.2
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 */
int GetBitGeneral( 
       const unsigned int isite,  //!<[in]
       const long unsigned int org_bit,  //!<[in]
       const long int *SiteToBit,  //!<[in]
       const long unsigned int *Tpow  //!<[in]
)
{
  long unsigned int tmp_bit=(org_bit/Tpow[isite-1])%SiteToBit[isite-1] ;
  return (tmp_bit);
}


/** 
 * 
 * @brief get 2sz at a site for general spin  
 *
 * @param isite site index (isite >= 1)
 * @param org_bit original bit to check 
 * @param SiteToBit List for getting bit at a site
 * @param Tpow List for getting total bit at a site before
 * @return 2sz at isite
 * 
 * @version 0.2
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 */

int GetLocal2Sz
(
 const unsigned int isite,  //!<[in]
 const long unsigned int org_bit,  //!<[in]
 const long int *SiteToBit,  //!<[in]
 const long unsigned int *Tpow  //!<[in]
 )
{
  int TwiceSz=0;
  int bitAtSite=0;
  //get bit
  bitAtSite=GetBitGeneral(isite, org_bit, SiteToBit, Tpow);
  TwiceSz=-(SiteToBit[isite-1]-1)+2*bitAtSite; //-2S^{total}_i+2Sz_i
  return TwiceSz;
}
/** 
 * 
 * @brief "finding the next higher number after a given number that has the same number of 1-bits" 
 * This method is introduced in S.H. Warren, Hacker’s Delight, second ed., Addison-Wesley, ISBN: 0321842685, 2012.
 *
 * @param x 
 * 
 * @version 2.0
 * @author Takahiro Misawa (The University of Tokyo) 
 */

unsigned long int snoob(unsigned long int x){
  unsigned long int smallest, ripple, ones;
  smallest = x &(-x);
  ripple   = x+ smallest;
  ones     = x ^ ripple;
  ones     = (ones>>2)/smallest;
  return   ripple|ones;
}
/** 
 * 
 * @brief calculating number of 1-bits in x (32 bit)
 * This method is introduced in S.H. Warren, Hacker’s Delight, second ed., Addison-Wesley, ISBN: 0321842685, 2012.
 *
 * @param x 
 * 
 * @version 2.0
 * @author Takahiro Misawa (The University of Tokyo) 
 */
int pop(unsigned int x){
  x = x - ((x>>1) & 0x55555555);
  x = (x & 0x33333333)+ ((x>>2)& 0x33333333);
  x = (x+(x>>4)) & 0x0F0F0F0F;
  x = x+ (x>>8);
  x = x+ (x>>16);
  return  x & 0x0000003F;
}
