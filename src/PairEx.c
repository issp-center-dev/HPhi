#include "bitcalc.h"
#include "wrapperMPI.h"
#include "mltplyMPI.h"
#include "mltply.h"
#include "PairEx.h"
#ifdef MPI
#include "mpi.h"
#endif

int GetPairExcitedState
(
 struct BindStruct *X,
 double complex *tmp_v0, /**< [out] Result v0 = H v1*/
 double complex *tmp_v1 /**< [in] v0 = H v1*/
 )
{

  FILE *fp;
  char sdt[D_FileNameMax];

  long unsigned int i,j;
  long unsigned int irght,ilft,ihfbit;
  long unsigned int isite1;
  long unsigned int org_isite1,org_isite2,org_sigma1,org_sigma2;
  long unsigned int tmp_off=0;
  double complex tmp_trans=0;
  long int i_max;
  int tmp_sgn, num1;
  int idx=0;
  int ihermite=0;
  long int ibit1, ibit;
  long unsigned int is1_up, is;
  //For TPQ
  int step=0;
  int rand_i=0;

  i_max = X->Check.idim_max;
  if(GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit)!=0){
    return -1;
  }
  X->Large.i_max    = i_max;
  X->Large.irght    = irght;
  X->Large.ilft     = ilft;
  X->Large.ihfbit   = ihfbit;
  X->Large.mode     = M_MLTPLY;

  switch(X->Def.iCalcModel){
  case HubbardGC:
    for(i=0;i<X->Def.NPairExcitationOperator;i++){
      org_isite1 = X->Def.PairExcitationOperator[i][0]+1;
      org_isite2 = X->Def.PairExcitationOperator[i][2]+1;
      org_sigma1 = X->Def.PairExcitationOperator[i][1];
      org_sigma2 = X->Def.PairExcitationOperator[i][3];
      tmp_trans = X->Def.ParaPairExcitationOperator[i];
      if (org_isite1  > X->Def.Nsite &&
          org_isite2  > X->Def.Nsite) {
        if(org_isite1==org_isite2 && org_sigma1==org_sigma2){
          if(org_sigma1==0){
            is   = X->Def.Tpow[2 * org_isite1 - 2];
          }
          else{
            is = X->Def.Tpow[2 * org_isite1 - 1];
          }
          ibit = (unsigned long int)myrank & is;
          if (ibit == is) {
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1)	\
  firstprivate(i_max, tmp_trans) private(j)
            for (j = 1; j <= i_max; j++) tmp_v0[j] += tmp_trans*tmp_v1[j];
          }
        }
        else{
          X_GC_child_general_hopp_MPIdouble(org_isite1-1, org_sigma1, org_isite2-1, org_sigma2, -tmp_trans, X, tmp_v0, tmp_v1);
        }
      }
      else if (org_isite2  > X->Def.Nsite || org_isite1  > X->Def.Nsite){
        if(org_isite1<org_isite2){
          X_GC_child_general_hopp_MPIsingle(org_isite1-1, org_sigma1, org_isite2-1, org_sigma2, -tmp_trans, X, tmp_v0, tmp_v1);
        }
        else{
          X_GC_child_general_hopp_MPIsingle(org_isite2-1, org_sigma2, org_isite1-1, org_sigma1, -conj(tmp_trans), X, tmp_v0, tmp_v1);
        }
      }
      else{
        if(child_general_hopp_GetInfo( X,org_isite1,org_isite2,org_sigma1,org_sigma2)!=0){
          return -1;
        }
        GC_child_general_hopp(tmp_v0, tmp_v1, X, tmp_trans);
      }
    }
    break;

  case KondoGC:
  case Hubbard:
  case Kondo:
    for(i=0;i<X->Def.NPairExcitationOperator;i++){
      org_isite1 = X->Def.PairExcitationOperator[i][0]+1;
      org_isite2 = X->Def.PairExcitationOperator[i][2]+1;
      org_sigma1 = X->Def.PairExcitationOperator[i][1];
      org_sigma2 = X->Def.PairExcitationOperator[i][3];
      tmp_trans = X->Def.ParaPairExcitationOperator[i];

      if(X->Def.iFlgSzConserved ==TRUE){
        if(org_sigma1 != org_sigma2){
          continue;
        }
      }
      if (org_isite1  > X->Def.Nsite &&
          org_isite2  > X->Def.Nsite) {
        if(org_isite1==org_isite2 && org_sigma1==org_sigma2){//diagonal

          is   = X->Def.Tpow[2 * org_isite1 - 2+org_sigma1];
          ibit = (unsigned long int)myrank & is;
          if (ibit == is) {
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1)	\
  firstprivate(i_max, tmp_trans) private(j)
            for (j = 1; j <= i_max; j++) tmp_v0[j] += tmp_trans*tmp_v1[j];
          }

        }
        else{
          X_child_general_hopp_MPIdouble(org_isite1-1, org_sigma1, org_isite2-1, org_sigma2, -tmp_trans, X, tmp_v0, tmp_v1);
        }
      }
      else if (org_isite2  > X->Def.Nsite || org_isite1  > X->Def.Nsite){
        if(org_isite1 < org_isite2){
          X_child_general_hopp_MPIsingle(org_isite1-1, org_sigma1,org_isite2-1, org_sigma2, -tmp_trans, X, tmp_v0, tmp_v1);
        }
        else{
          X_child_general_hopp_MPIsingle(org_isite2-1, org_sigma2, org_isite1-1, org_sigma1, -conj(tmp_trans), X, tmp_v0, tmp_v1);
        }
      }
      else{
        if(child_general_hopp_GetInfo( X,org_isite1,org_isite2,org_sigma1,org_sigma2)!=0){
          return -1;
        }
        if(org_isite1==org_isite2 && org_sigma1==org_sigma2){

          is   = X->Def.Tpow[2 * org_isite1 - 2 + org_sigma1];

#pragma omp parallel for default(none) shared(list_1, tmp_v0, tmp_v1) firstprivate(i_max, is, tmp_trans) private(num1, ibit)
          for(j = 1;j <= i_max;j++){
            ibit = list_1[j]&is;
            num1  = ibit/is;
            tmp_v0[j] += tmp_trans*num1*tmp_v1[j];
          }
        }
        else{
          child_general_hopp(tmp_v0, tmp_v1,X,tmp_trans);
        }
      }
    }
    break;

  case Spin: // for the Sz-conserved spin system

    if(X->Def.iFlgGeneralSpin==FALSE){
      for(i=0;i<X->Def.NPairExcitationOperator;i++) {
        org_isite1 = X->Def.PairExcitationOperator[i][0] + 1;
        org_isite2 = X->Def.PairExcitationOperator[i][2] + 1;
        org_sigma1 = X->Def.PairExcitationOperator[i][1];
        org_sigma2 = X->Def.PairExcitationOperator[i][3];
        tmp_trans = X->Def.ParaPairExcitationOperator[i];
        if (org_sigma1 == org_sigma2) {
          if (org_isite1 == org_isite2) {
            if (org_isite1 > X->Def.Nsite) {
              is1_up = X->Def.Tpow[org_isite1 - 1];
              ibit1 = X_SpinGC_CisAis((unsigned long int) myrank + 1, X, is1_up, org_sigma1);
              if (ibit1 != 0) {
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1)	\
  firstprivate(i_max, tmp_trans) private(j)
                for (j = 1; j <= i_max; j++) tmp_v0[j] += tmp_trans * tmp_v1[j];
              }
            }// org_isite1 > X->Def.Nsite
            else {
              isite1 = X->Def.Tpow[org_isite1 - 1];
#pragma omp parallel for default(none) private(j) firstprivate(i_max, isite1, org_sigma1, X, tmp_trans) shared(tmp_v0, tmp_v1)
              for (j = 1; j <= i_max; j++) {
                tmp_v0[j] += X_Spin_CisAis(j, X, isite1, org_sigma1) * tmp_v1[j] * tmp_trans;
              }
            }
          } else {
            // for the canonical case
          }
        }
      }
    }//FlgGeneralSpin=FALSE
    else{
      for(i=0;i<X->Def.NPairExcitationOperator;i++){
        org_isite1 = X->Def.PairExcitationOperator[i][0]+1;
        org_isite2 = X->Def.PairExcitationOperator[i][2]+1;
        org_sigma1 = X->Def.PairExcitationOperator[i][1];
        org_sigma2 = X->Def.PairExcitationOperator[i][3];
        tmp_trans = X->Def.ParaPairExcitationOperator[i];
        if(org_isite1 == org_isite2){
          if(org_isite1 >X->Def.Nsite){
            if(org_sigma1==org_sigma2){
              // longitudinal magnetic field
              num1 = BitCheckGeneral((unsigned long int)myrank,
                                     org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
              if (num1 != 0) {
#pragma omp parallel for default(none) private(j) firstprivate(i_max, tmp_trans) shared(tmp_v0,tmp_v1)
                for(j=1;j<=i_max;j++){
                  tmp_v0[j]+= tmp_trans*tmp_v1[j];
                }
              }
            }
          }
          else {//org_isite1 <= X->Def.Nsite
            if(org_sigma1==org_sigma2){
              // longitudinal magnetic field
#pragma omp parallel for default(none) private(j, num1) firstprivate(i_max, org_isite1, org_sigma1, X, tmp_trans) shared(tmp_v0,tmp_v1, list_1)
              for(j=1;j<=i_max;j++){
                num1 = BitCheckGeneral(list_1[j], org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
                tmp_v0[j]+= tmp_trans*tmp_v1[j]*num1;
              }
            }
          }
        }else{
          // hopping is not allowed in localized spin system
        }//org_isite1 != org_isite2
      }
    }//general spin

    break;

  case SpinGC:

    if(X->Def.iFlgGeneralSpin==FALSE){
      for(i=0;i<X->Def.NPairExcitationOperator;i++){
        org_isite1 = X->Def.PairExcitationOperator[i][0]+1;
        org_isite2 = X->Def.PairExcitationOperator[i][2]+1;
        org_sigma1 = X->Def.PairExcitationOperator[i][1];
        org_sigma2 = X->Def.PairExcitationOperator[i][3];
        tmp_trans = X->Def.ParaPairExcitationOperator[i];

        if(org_isite1 == org_isite2){
          if(org_isite1 > X->Def.Nsite){
            if(org_sigma1==org_sigma2){  // longitudinal magnetic field
              X_GC_child_CisAis_spin_MPIdouble(org_isite1-1, org_sigma1, tmp_trans, X, tmp_v0, tmp_v1);
            }
            else{  // transverse magnetic field
              X_GC_child_CisAit_spin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, tmp_trans, X, tmp_v0, tmp_v1);
            }
          }else{
            isite1 = X->Def.Tpow[org_isite1-1];
            if(org_sigma1==org_sigma2){
              // longitudinal magnetic field
#pragma omp parallel for default(none) private(j, tmp_sgn) firstprivate(i_max, isite1, org_sigma1, X,tmp_trans) shared(tmp_v0, tmp_v1)
              for(j=1;j<=i_max;j++){
                tmp_v0[j] += X_SpinGC_CisAis(j, X, isite1, org_sigma1)*tmp_v1[j]*tmp_trans;
              }
            }else{
              // transverse magnetic field
#pragma omp parallel for default(none) private(j, tmp_sgn, tmp_off) firstprivate(i_max, isite1, org_sigma2, X, tmp_trans) shared(tmp_v0, tmp_v1)
              for(j=1;j<=i_max;j++){
                tmp_sgn  =  X_SpinGC_CisAit(j,X, isite1,org_sigma2,&tmp_off);
                if(tmp_sgn !=0){
                  tmp_v0[j]+= tmp_sgn*tmp_v1[j]*tmp_trans;
                }
              }
            }
          }
        }else{
          // hopping is not allowed in localized spin system
        }
      }

    }//FlgGeneralSpin=FALSE
    else{
      for(i=0;i<X->Def.NPairExcitationOperator;i++){
        org_isite1 = X->Def.PairExcitationOperator[i][0]+1;
        org_isite2 = X->Def.PairExcitationOperator[i][2]+1;
        org_sigma1 = X->Def.PairExcitationOperator[i][1];
        org_sigma2 = X->Def.PairExcitationOperator[i][3];
        tmp_trans = X->Def.ParaPairExcitationOperator[i];
        if(org_isite1 == org_isite2){
          if(org_isite1 > X->Def.Nsite){
            if(org_sigma1==org_sigma2){
              // longitudinal magnetic field
              X_GC_child_CisAis_GeneralSpin_MPIdouble(org_isite1-1, org_sigma1, tmp_trans, X, tmp_v0, tmp_v1);
            }else{
              // transverse magnetic field
              X_GC_child_CisAit_GeneralSpin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, tmp_trans, X, tmp_v0, tmp_v1);
            }
          }
          else{//org_isite1 <= X->Def.Nsite
            if(org_sigma1==org_sigma2){
              // longitudinal magnetic field
#pragma omp parallel for default(none) private(j, num1) firstprivate(i_max, org_isite1, org_sigma1, X, tmp_trans) shared(tmp_v0, tmp_v1)
              for(j=1;j<=i_max;j++){
                num1 = BitCheckGeneral(j-1, org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
                tmp_v0[j] += tmp_trans* tmp_v1[j]*num1;
              }
            }else{
              // transverse magnetic field
#pragma omp parallel for default(none) private(j, num1) firstprivate(i_max, org_isite1, org_sigma1, org_sigma2, X,tmp_off, tmp_trans) shared(tmp_v0, tmp_v1)
              for(j=1;j<=i_max;j++){
                num1 = GetOffCompGeneralSpin(j-1, org_isite1, org_sigma2, org_sigma1, &tmp_off, X->Def.SiteToBit, X->Def.Tpow);
                if(num1 !=0){
                  tmp_v0[j]  +=  tmp_trans*tmp_v1[j]*num1;
                }
              }
            }
          }
        }else{
          // hopping is not allowed in localized spin system
        }
      }
    }
    break;

  default:
    return FALSE;
  }

  return TRUE;
}
