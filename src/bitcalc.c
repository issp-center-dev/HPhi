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
 * @brief function of getting right, left and half bits corresponding to a original hilbert space.
 }
 * @param Nsite a total number of sites
 * @param irght a bit to split original hilbert space into right space
 * @param ilft a bit to split original hilbert space into left space
 * @param ihfbit a half bit to split original hilbert space
 * 
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo) 
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 * @return 
 */
int GetSplitBit(
		const int Nsite,
		long unsigned int *irght,
		long unsigned int *ilft,
		long unsigned int *ihfbit
){
  if(Nsite<1){
    fprintf(stderr, "%s", cErrSiteNumber);
    return -1;
  }
  *irght=pow(2,((Nsite+1)/2))-1;
  *ilft=pow(2,Nsite)-1;
  *ilft=*ilft ^ *irght; 
  *ihfbit=pow(2,((Nsite+1)/2));
  return 0;
}

/** 
 * @brief function of splitting original bit into right and left hilbert spaces.
 * 
 * @param Nsite a total number of sites
 * @param iCalcModel Calc model defined in CalcMode file
 * @param irght a bit to split original hilbert space into right space
 * @param ilft a bit to split original hilbert space into left space
 * @param ihfbit a half bit to split original hilbert space
 * 
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo) 
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 * @return 
 */
int GetSplitBitByModel(
		       const int Nsite,
		       const int iCalcModel,
		       long unsigned int *irght,
		       long unsigned int *ilft,
		       long unsigned int *ihfbit
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

  if(!GetSplitBit(tmpNsite, irght, ilft, ihfbit)==0){
    return -1;
  }
  
  return 0;
}


/** 
 * 
 * @brief function of getting right, left and half bits corresponding to a original hilbert space.
 * @param Nsite a total number of sites
 * @param ihfbit a bit to split original hilbert space
 * 
 * @retval 0 normally finished
 * @retval -1 unnormally finished
 *
 * @version 0.2
 * @author Takahiro Misawa (The University of Tokyo) 
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 */
int GetSplitBitForGeneralSpin(
		const int Nsite,
		long unsigned int *ihfbit,
		const long int *SiteToBit
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
 * @param irght a bit to split original hilbert space into right space
 * @param ilft a bit to split original hilbert space into left space
 * @param ihfbit a half bit to split original hilbert space
 * @param isplited_Bit_right a splitted bit reflected on right space
 * @param isplited_Bit_left a splitted bit reflected on left space
 * @version 0.1 
 * @author Takahiro Misawa (The University of Tokyo) 
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 */
void SplitBit(
		  const long unsigned int ibit,
		  const long unsigned int irght,
		  const long unsigned int ilft,
		  const long unsigned int ihfbit,
		  long unsigned int *isplited_Bit_right,
		  long unsigned int *isplited_Bit_left
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
 * @param irght a bit to split original hilbert space into right space
 * @param ilft a bit to split original hilbert space into left space
 * @param ihfbit a half bit to split original hilbert space
 * @param _ioffComp an off diagonal component
 * @version 0.1 
 * @author Takahiro Misawa (The University of Tokyo) 
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 */
int GetOffComp(
	    long unsigned int *_list_2_1,
	    long unsigned int *_list_2_2,
		long unsigned int _ibit,
		const long unsigned int _irght,
		const long unsigned int _ilft,
		const long unsigned int _ihfbit,
		long unsigned int *_ioffComp
)
{
  long unsigned int ia, ib;
  SplitBit(_ibit, _irght, _ilft, _ihfbit, &ia, &ib);
  *_ioffComp =_list_2_1[ia];
  *_ioffComp+=_list_2_2[ib];
  if(*_ioffComp !=0) return TRUE;
  else return FALSE;
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
 * @param _SiteToBit List for getting bit at a site
 * @param _Tpow List for getting total bit at a site before
 * @retval FALSE off-diagonal component does not exist
 * @retval TRUE off-diagonal component exists
 * 
 * @version 0.2
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 */
int GetOffCompGeneralSpin(
		const long unsigned int org_ibit,
		const int org_isite,
		const int org_ispin,
		const int off_ispin,
		long  unsigned int *_ioffComp,
		const long int *SiteToBit,
		const long int *Tpow
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
void ConvertToList1GeneralSpin(
		const long unsigned int org_ibit,
		const long unsigned int ihlfbit,
		long unsigned int *_ilist1Comp
)
{
  long unsigned int ia, ib;
  ia=org_ibit%ihlfbit;
  ib=org_ibit/ihlfbit;
  *_ilist1Comp=list_2_1[ia]+list_2_2[ib];
}

/** 
 * 
 * @brief function of getting fermion sign (for 32bit)
 * 
 * @param org_bit an original bit
 * @param _sgn fermion sign 
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo) 
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 */
void SgnBit_old( 
		  const long unsigned int org_bit,
                  int *sgn
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
 * @brief function of getting fermion sign
 * 
 * @param org_bit an original bit
 * @param _sgn fermion sign 
 * @version 0.1
 *
 * @author Takahiro Misawa (The University of Tokyo) 
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 */
void SgnBit( 
		  const long unsigned int org_bit,
                  int *sgn
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
	     const long unsigned int org_bit,
	     const long unsigned int target_bit
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
 * @param _SiteToBit List for getting bit at a site
 * @param _TPow List for getting total bit at a site before
 * @retval 0 bit does not exists
 * @retval 1 bit exists
 * 
 * @version 0.2
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 */
int BitCheckGeneral(
	     const long unsigned int org_bit,
	     const unsigned int org_isite,
	     const unsigned int target_ispin,
	     const long int *SiteToBit,
	     const long int *Tpow
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
 * @param _SiteToBit List for getting bit at a site
 * @param _Tpow List for getting total bit at a site before
 * @return bit at a site
 * 
 * @version 0.2
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 */
int GetBitGeneral( 
	     const unsigned int isite,
	     const long unsigned int org_bit,
	     const long int *SiteToBit,
	     const long int *Tpow
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
 * @param _SiteToBit List for getting bit at a site
 * @param _Tpow List for getting total bit at a site before
 * @return 2sz at isite
 * 
 * @version 0.2
 * @author Kazuyoshi Yoshimi (The University of Tokyo) 
 */

int GetLocal2Sz
(
 const unsigned int isite,
 const long unsigned int org_bit,
 const long int *SiteToBit,
 const long unsigned int *Tpow
 )
{
  int TwiceSz=0;
  int bitAtSite=0;
  //get bit
  bitAtSite=GetBitGeneral(isite, org_bit, SiteToBit, Tpow);
  TwiceSz=-(SiteToBit[isite-1]-1)+2*bitAtSite; //-2S^{total}_i+2Sz_i
  return TwiceSz;
}

unsigned long int snoob(unsigned long int x){
  unsigned long int smallest, ripple, ones;
  smallest = x &(-x);
  ripple   = x+ smallest;
  ones     = x ^ ripple;
  ones     = (ones>>2)/smallest;
  return   ripple|ones;
}
