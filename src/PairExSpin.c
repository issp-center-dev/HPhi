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
/*-------------------------------------------------------------*/
#include "PairExSpin.h"
#include "bitcalc.h"
#include "wrapperMPI.h"
#include "mltplyMPISpinCore.h"
#include "mltplySpinCore.h"
#ifdef MPI
#include "common/setmemory.h"
#endif

//
/// \brief Calculation of pair excited state for Spin Grand canonical system
/// \param X [in,out] define list to get and put information of calculation
/// \param tmp_v0 [out] Result v0 = H v1
/// \param tmp_v1 [in] v0 = H v1
/// \returns TRUE: Normally finished
/// \returns FALSE: Abnormally finished
/// \author Kazuyoshi Yoshimi
/// \version 1.2
int GetPairExcitedStateSpinGC(
        struct BindStruct *X,/**< [in,out] define list to get and put information of calculation*/
        double complex *tmp_v0, /**< [out] Result v0 = H v1*/
        double complex *tmp_v1 /**< [in] v0 = H v1*/

){

    int iret=0;
    if (X->Def.iFlgGeneralSpin == FALSE) {
        iret=GetPairExcitedStateHalfSpinGC(X, tmp_v0, tmp_v1);
    }
    else{
        iret=GetPairExcitedStateGeneralSpinGC(X, tmp_v0, tmp_v1);
    }
    return iret;
}


//
/// Calculation of pair excited state for Half Spin Grand canonical system
/// \param X [in,out] define list to get and put information of calculation
/// \param tmp_v0 [out] Result v0 = H v1
/// \param tmp_v1 [in] v0 = H v1
/// \returns TRUE: Normally finished
/// \returns FALSE: Abnormally finished
/// \author Kazuyoshi Yoshimi
/// \version 1.2
int GetPairExcitedStateHalfSpinGC(
        struct BindStruct *X,/**< [in,out] define list to get and put information of calculation*/
        double complex *tmp_v0, /**< [out] Result v0 = H v1*/
        double complex *tmp_v1 /**< [in] v0 = H v1*/

){
    long unsigned int i,j;
    long unsigned int isite1;
    long unsigned int org_isite1,org_isite2,org_sigma1,org_sigma2;
    long unsigned int tmp_off=0;

    double complex tmp_trans=0;
    long int i_max;
    int tmp_sgn;
    i_max = X->Check.idim_maxOrg;

    for(i=0;i<X->Def.NPairExcitationOperator;i++){
        org_isite1 = X->Def.PairExcitationOperator[i][0]+1;
        org_isite2 = X->Def.PairExcitationOperator[i][2]+1;
        org_sigma1 = X->Def.PairExcitationOperator[i][1];
        org_sigma2 = X->Def.PairExcitationOperator[i][3];
        tmp_trans = X->Def.ParaPairExcitationOperator[i];
        if(org_isite1 == org_isite2){
            if(org_isite1 > X->Def.Nsite){
                if(org_sigma1==org_sigma2){  // longitudinal magnetic field
                    if(X->Def.PairExcitationOperator[i][4]==0) {
                        child_GC_AisCis_spin_MPIdouble(org_isite1 - 1, org_sigma1, -tmp_trans, X, tmp_v0, tmp_v1);
                    }
                    else{
                        child_GC_CisAis_spin_MPIdouble(org_isite1 - 1, org_sigma1, tmp_trans, X, tmp_v0, tmp_v1);
                    }
                }
                else{  // transverse magnetic field
                    //fprintf(stdoutMPI, "Debug: test, org_isite1=%d, org_sigma1=%d, orgsima_2=%d\n", org_isite1, org_sigma1, org_sigma2);
                    child_GC_CisAit_spin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, tmp_trans, X, tmp_v0, tmp_v1);
                }
            }else{
                isite1 = X->Def.Tpow[org_isite1-1];
                if(org_sigma1==org_sigma2) {
                    if (X->Def.PairExcitationOperator[i][4] == 0) {
                        // longitudinal magnetic field
#pragma omp parallel for default(none) private(j, tmp_sgn) firstprivate(i_max, isite1, org_sigma1, X,tmp_trans) shared(tmp_v0, tmp_v1)
                        for (j = 1; j <= i_max; j++) {
                            tmp_v0[j] += (1.0-child_SpinGC_CisAis(j, X, isite1, org_sigma1)) * tmp_v1[j] * (-tmp_trans);
                        }
                    }
                    else {
                        // longitudinal magnetic field
#pragma omp parallel for default(none) private(j, tmp_sgn) firstprivate(i_max, isite1, org_sigma1, X,tmp_trans) shared(tmp_v0, tmp_v1)
                        for (j = 1; j <= i_max; j++) {
                            tmp_v0[j] += child_SpinGC_CisAis(j, X, isite1, org_sigma1) * tmp_v1[j] * tmp_trans;
                        }
                    }
                }else{
                    // transverse magnetic field
                    // fprintf(stdoutMPI, "Debug: isite1=%d, org_sigma2=%d\n", isite1, org_sigma2);
#pragma omp parallel for default(none) private(j, tmp_sgn, tmp_off) firstprivate(i_max, isite1, org_sigma2, X, tmp_trans) shared(tmp_v0, tmp_v1)
                    for(j=1;j<=i_max;j++){
                        tmp_sgn  =  child_SpinGC_CisAit(j,X, isite1,org_sigma2,&tmp_off);
                        if(tmp_sgn !=0){
                            tmp_v0[tmp_off+1]+= tmp_sgn*tmp_v1[j]*tmp_trans;
                        }
                    }
                }
            }
        }else{
            fprintf(stdoutMPI, "ERROR: hopping is not allowed in localized spin system\n");
            return FALSE;
        }
    }
    return TRUE;
}

//
/// Calculation of pair excited state for general Spin Grand canonical system
/// \param X [in,out] define list to get and put information of calculation
/// \param tmp_v0 [out] Result v0 = H v1
/// \param tmp_v1 [in] v0 = H v1
/// \returns TRUE: Normally finished
/// \returns FALSE: Abnormally finished
/// \author Kazuyoshi Yoshimi
/// \version 1.2
int GetPairExcitedStateGeneralSpinGC(
        struct BindStruct *X,/**< [in,out] define list to get and put information of calculation*/
        double complex *tmp_v0, /**< [out] Result v0 = H v1*/
        double complex *tmp_v1 /**< [in] v0 = H v1*/

) {
    long unsigned int i, j;
    int num1;
    long unsigned int org_isite1, org_isite2, org_sigma1, org_sigma2;
    long unsigned int tmp_off = 0;

    double complex tmp_trans = 0;
    long int i_max;
    i_max = X->Check.idim_maxOrg;

    for(i=0;i<X->Def.NPairExcitationOperator;i++){
        org_isite1 = X->Def.PairExcitationOperator[i][0]+1;
        org_isite2 = X->Def.PairExcitationOperator[i][2]+1;
        org_sigma1 = X->Def.PairExcitationOperator[i][1];
        org_sigma2 = X->Def.PairExcitationOperator[i][3];
        tmp_trans = X->Def.ParaPairExcitationOperator[i];
        if(org_isite1 == org_isite2){
            if(org_isite1 > X->Def.Nsite){
                if(org_sigma1==org_sigma2){
                    if(X->Def.PairExcitationOperator[i][4]==0) {
                        // longitudinal magnetic field
                        child_GC_AisCis_GeneralSpin_MPIdouble(org_isite1 - 1, org_sigma1, -tmp_trans, X, tmp_v0, tmp_v1);
                    }
                    else{
                        child_GC_CisAis_GeneralSpin_MPIdouble(org_isite1 - 1, org_sigma1, tmp_trans, X, tmp_v0, tmp_v1);
                    }
                }else{
                    // transverse magnetic field
                    child_GC_CisAit_GeneralSpin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, tmp_trans, X, tmp_v0, tmp_v1);
                }
            }
            else{//org_isite1 <= X->Def.Nsite
                if(org_sigma1==org_sigma2){
                    if(X->Def.PairExcitationOperator[i][4]==0) {
                        // longitudinal magnetic field
#pragma omp parallel for default(none) private(j, num1) firstprivate(i_max, org_isite1, org_sigma1, X, tmp_trans) shared(tmp_v0, tmp_v1)
                        for (j = 1; j <= i_max; j++) {
                            num1 = BitCheckGeneral(j - 1, org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
                            tmp_v0[j] += -tmp_trans * tmp_v1[j] * (1.0-num1);
                        }
                    }
                    else{
                        // longitudinal magnetic field
#pragma omp parallel for default(none) private(j, num1) firstprivate(i_max, org_isite1, org_sigma1, X, tmp_trans) shared(tmp_v0, tmp_v1)
                        for (j = 1; j <= i_max; j++) {
                            num1 = BitCheckGeneral(j - 1, org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
                            tmp_v0[j] += tmp_trans * tmp_v1[j] * num1;
                        }
                    }
                }else{
                    // transverse magnetic field
#pragma omp parallel for default(none) private(j, num1) firstprivate(i_max, org_isite1, org_sigma1, org_sigma2, X,tmp_off, tmp_trans) shared(tmp_v0, tmp_v1)
                    for(j=1;j<=i_max;j++){
                        num1 = GetOffCompGeneralSpin(j-1, org_isite1, org_sigma2, org_sigma1, &tmp_off, X->Def.SiteToBit, X->Def.Tpow);
                        if(num1 !=0){
                            tmp_v0[tmp_off+1]  +=  tmp_trans*tmp_v1[j]*num1;
                        }
                    }
                }
            }
        }else{
            fprintf(stdoutMPI, "ERROR: hopping is not allowed in localized spin system\n");
            return FALSE;
        }
    }
    return TRUE;
}

//
/// Calculation of pair excited state for Spin canonical system
/// \param X [in,out] define list to get and put information of calculation
/// \param tmp_v0 [out] Result v0 = H v1
/// \param tmp_v1 [in] v0 = H v1
/// \returns TRUE: Normally finished
/// \returns FALSE: Abnormally finished
/// \author Kazuyoshi Yoshimi
/// \version 1.2
int GetPairExcitedStateSpin(
        struct BindStruct *X,/**< [in,out] define list to get and put information of calculation*/
        double complex *tmp_v0, /**< [out] Result v0 = H v1*/
        double complex *tmp_v1 /**< [in] v0 = H v1*/

){
    int iret=0;
    if (X->Def.iFlgGeneralSpin == FALSE) {
        iret=GetPairExcitedStateHalfSpin(X, tmp_v0, tmp_v1);
    }
    else{
        iret=GetPairExcitedStateGeneralSpin(X, tmp_v0, tmp_v1);
    }
    return iret;
}

//
/// Calculation of pair excited state for Half Spin canonical system
/// \param X [in,out] define list to get and put information of calculation
/// \param tmp_v0 [out] Result v0 = H v1
/// \param tmp_v1 [in] v0 = H v1
/// \returns TRUE: Normally finished
/// \returns FALSE: Abnormally finished
/// \author Kazuyoshi Yoshimi
/// \version 1.2
int GetPairExcitedStateHalfSpin(
        struct BindStruct *X,/**< [in,out] define list to get and put information of calculation*/
        double complex *tmp_v0, /**< [out] Result v0 = H v1*/
        double complex *tmp_v1 /**< [in] v0 = H v1*/

)
{
    long unsigned int i,j, idim_maxMPI;
    long unsigned int isite1;
    long unsigned int org_isite1,org_isite2,org_sigma1,org_sigma2;
    long unsigned int tmp_off=0;

    double complex tmp_trans=0;
    long int i_max;
    int num1;
    long int ibit1;
    long unsigned int is1_up;

    i_max = X->Check.idim_maxOrg;

    double complex *tmp_v1bufOrg;
    //set size
#ifdef MPI
    idim_maxMPI = MaxMPI_li(X->Check.idim_maxOrg);
    tmp_v1bufOrg=cd_1d_allocate(idim_maxMPI + 1);
#endif // MPI

    for (i = 0; i < X->Def.NPairExcitationOperator; i++) {
        org_isite1 = X->Def.PairExcitationOperator[i][0] + 1;
        org_isite2 = X->Def.PairExcitationOperator[i][2] + 1;
        org_sigma1 = X->Def.PairExcitationOperator[i][1];
        org_sigma2 = X->Def.PairExcitationOperator[i][3];
        tmp_trans = X->Def.ParaPairExcitationOperator[i];
        if (org_sigma1 == org_sigma2) {
            if (org_isite1 == org_isite2) {
                if (org_isite1 > X->Def.Nsite) {
                    is1_up = X->Def.Tpow[org_isite1 - 1];
                    ibit1 = child_SpinGC_CisAis((unsigned long int) myrank + 1, X, is1_up, org_sigma1);
                    if (X->Def.PairExcitationOperator[i][4] == 0) {
                        if (ibit1 == 0) {
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1)	\
  firstprivate(i_max, tmp_trans) private(j)
                            for (j = 1; j <= i_max; j++) tmp_v0[j] += -tmp_trans * tmp_v1[j];
                        }
                    } else {
                        if (ibit1 != 0) {
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1)	\
  firstprivate(i_max, tmp_trans) private(j)
                            for (j = 1; j <= i_max; j++) tmp_v0[j] += tmp_trans * tmp_v1[j];
                        }
                    }
                }// org_isite1 > X->Def.Nsite
                else {
                    isite1 = X->Def.Tpow[org_isite1 - 1];
                    if (org_isite1 == org_isite2 && org_sigma1 == org_sigma2 &&
                        X->Def.PairExcitationOperator[i][4] == 0) {
#pragma omp parallel for default(none) private(j) firstprivate(i_max, isite1, org_sigma1, X, tmp_trans) shared(tmp_v0, tmp_v1)
                        for (j = 1; j <= i_max; j++) {
                            tmp_v0[j] += (1.0 - child_Spin_CisAis(j, X, isite1, org_sigma1)) * tmp_v1[j] * (-tmp_trans);
                        }
                    } else {
#pragma omp parallel for default(none) private(j) firstprivate(i_max, isite1, org_sigma1, X, tmp_trans) shared(tmp_v0, tmp_v1)
                        for (j = 1; j <= i_max; j++) {
                            tmp_v0[j] += child_Spin_CisAis(j, X, isite1, org_sigma1) * tmp_v1[j] * tmp_trans;
                        }
                    }
                }
            } else {
                fprintf(stdoutMPI, "Error: isite1 must be equal to isite2 for Spin system. \n");
                return FALSE;
            }
        } else { //org_sigma1 != org_sigma2             // for the canonical case
            if (org_isite1 > X->Def.Nsite) {//For MPI
                child_CisAit_spin_MPIdouble(org_isite1-1, org_sigma2, tmp_trans, X, tmp_v0, tmp_v1, tmp_v1bufOrg, i_max, X->Def.Tpow,list_1_org, list_1buf_org, list_2_1, list_2_2, X->Large.irght, X->Large.ilft,X->Large.ihfbit);

            } else {
                isite1 = X->Def.Tpow[org_isite1 - 1];
#pragma omp parallel for default(none) private(j, tmp_off, num1)        \
  firstprivate(i_max, isite1, org_sigma2, X, tmp_trans, list_1_org, list_1, list_2_1, list_2_2) shared(tmp_v0, tmp_v1)
                for (j = 1; j <= i_max; j++) {
                    num1=child_Spin_CisAit(j, X, isite1, org_sigma2, list_1_org, list_2_1, list_2_2, &tmp_off);
                    if (num1 != 0) tmp_v0[tmp_off] +=  tmp_v1[j] * tmp_trans*(double)num1;
                }
            }
        }
    }
#ifdef MPI
    free_cd_1d_allocate(tmp_v1bufOrg);
#endif
    return TRUE;
}


//
/// Calculation of pair excited state for general Spin canonical system
/// \param X [in,out] define list to get and put information of calculation
/// \param tmp_v0 [out] Result v0 = H v1
/// \param tmp_v1 [in] v0 = H v1
/// \returns TRUE: Normally finished
/// \returns FALSE: Abnormally finished
/// \author Kazuyoshi Yoshimi
/// \version 1.2
int GetPairExcitedStateGeneralSpin(
        struct BindStruct *X,/**< [in,out] define list to get and put information of calculation*/
        double complex *tmp_v0, /**< [out] Result v0 = H v1*/
        double complex *tmp_v1 /**< [in] v0 = H v1*/

)
{
    long unsigned int i,j, idim_maxMPI;
    long unsigned int org_isite1,org_isite2,org_sigma1,org_sigma2;
    long unsigned int tmp_off=0;
    long unsigned int off=0;

    double complex tmp_trans=0;
    long int i_max;
    int tmp_sgn, num1;
    i_max = X->Check.idim_maxOrg;

    double complex *tmp_v1bufOrg;
    //set size
#ifdef MPI
    idim_maxMPI = MaxMPI_li(X->Check.idim_maxOrg);
    tmp_v1bufOrg = cd_1d_allocate(idim_maxMPI + 1);
#endif // MPI

    for(i=0;i<X->Def.NPairExcitationOperator;i++) {
        org_isite1 = X->Def.PairExcitationOperator[i][0] + 1;
        org_isite2 = X->Def.PairExcitationOperator[i][2] + 1;
        org_sigma1 = X->Def.PairExcitationOperator[i][1];
        org_sigma2 = X->Def.PairExcitationOperator[i][3];
        tmp_trans = X->Def.ParaPairExcitationOperator[i];
        if (org_isite1 == org_isite2) {
            if (org_isite1 > X->Def.Nsite) {
                if (org_sigma1 == org_sigma2) {
                    // longitudinal magnetic field
                    num1 = BitCheckGeneral((unsigned long int) myrank,
                                           org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
                    if (X->Def.PairExcitationOperator[i][4] == 0) {
                        if (num1 == 0) {
#pragma omp parallel for default(none) private(j) firstprivate(i_max, tmp_trans) shared(tmp_v0, tmp_v1)
                            for (j = 1; j <= i_max; j++) {
                                tmp_v0[j] += -tmp_trans * tmp_v1[j];
                            }
                        }
                    } else {
                        if (num1 != 0) {
#pragma omp parallel for default(none) private(j) firstprivate(i_max, tmp_trans) shared(tmp_v0, tmp_v1)
                            for (j = 1; j <= i_max; j++) {
                                tmp_v0[j] += tmp_trans * tmp_v1[j];
                            }
                        }
                    }
                }//org_sigma1=org_sigma2
                else {//org_sigma1 != org_sigma2
                    child_CisAit_GeneralSpin_MPIdouble(org_isite1 - 1, org_sigma1, org_sigma2, tmp_trans, X, tmp_v0,
                                                         tmp_v1, tmp_v1bufOrg, i_max, list_1_org, list_1buf_org,
                                                         X->Large.ihfbit);
                }
            } else {//org_isite1 <= X->Def.Nsite
                if (org_sigma1 == org_sigma2) {
                    // longitudinal magnetic field
                    if (X->Def.PairExcitationOperator[i][4] == 0) {
#pragma omp parallel for default(none) private(j, num1) firstprivate(i_max, org_isite1, org_sigma1, X, tmp_trans) shared(tmp_v0, tmp_v1, list_1)
                        for (j = 1; j <= i_max; j++) {
                            num1 = BitCheckGeneral(list_1[j], org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
                            tmp_v0[j] += -tmp_trans * tmp_v1[j] * (1.0 - num1);
                        }
                    } else {
#pragma omp parallel for default(none) private(j, num1) firstprivate(i_max, org_isite1, org_sigma1, X, tmp_trans) shared(tmp_v0, tmp_v1, list_1)
                        for (j = 1; j <= i_max; j++) {
                            num1 = BitCheckGeneral(list_1[j], org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
                            tmp_v0[j] += tmp_trans * tmp_v1[j] * num1;
                        }
                    }
                }//org_sigma1=org_sigma2
                else {//org_sigma1 != org_sigma2
#pragma omp parallel for default(none) private(j, tmp_sgn, tmp_off)      \
  firstprivate(i_max, org_isite1, org_sigma1, org_sigma2, X, off, tmp_trans, myrank) \
  shared(tmp_v0, tmp_v1, list_1_org, list_1)
                    for (j = 1; j <= i_max; j++) {
                        tmp_sgn = GetOffCompGeneralSpin(list_1_org[j], org_isite1, org_sigma2, org_sigma1, &off,
                                                        X->Def.SiteToBit, X->Def.Tpow);
                        if (tmp_sgn != FALSE) {
                            ConvertToList1GeneralSpin(off, X->Large.ihfbit, &tmp_off);
#ifdef _DEBUG
                            printf("rank=%d, org=%ld, tmp_off=%ld, list_1=%ld, ihfbit=%ld\n",myrank, list_1_org[j], off, list_1[tmp_off], X->Large.ihfbit);
#endif
                            tmp_v0[tmp_off] += tmp_v1[j] * tmp_trans;
                        }
                    }

                }
            }
        } else {
            fprintf(stdoutMPI, "ERROR: hopping is not allowed in localized spin system\n");
            return FALSE;
        }//org_isite1 != org_isite2
    }
#ifdef MPI
    free_cd_1d_allocate(tmp_v1bufOrg);
#endif // MPI

    return TRUE;
}
