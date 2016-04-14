/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 Takahiro Misawa, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Youhei Yamaji, Synge Todo, Naoki Kawashima */

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
#include "mltply.h"
#include "bitcalc.h"
#include "CalcSpectrum.h"
#include "CalcSpectrumByLanczos.h"
#include "wrapperMPI.h"
#include "mltplyMPI.h"

/**
 * @file   CalcSpectrum.c
 * @version 1.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 * @brief  File for givinvg functions of calculating spectrum 
 * 
 * 
 */

/** 
 * @brief A main function to calculate spectrum 
 * 
 * @param[in,out] X CalcStruct list for getting and pushing calculation information 
 * @retval 0 normally finished
 * @retval -1 unnormally finished
 *
 * @version 1.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 */
int CalcSpectrum(		 
		 struct EDMainCalStruct *X
				 )
{
  char sdt[D_FileNameMax];
  double diff_ene,var;
  double complex cdnorm;
  unsigned long int i;
  unsigned long int i_max=0;
  FILE *fp;
  double dnorm;

  //input eigen vector
  fprintf(stdoutMPI, "An Eigenvector is inputted.\n");
  sprintf(sdt, cFileNameInputEigen, X->Bind.Def.CDataFileHead, X->Bind.Def.k_exct-1, myrank);
  fp = fopen(sdt, "rb");
  if(fp==NULL){
    fprintf(stderr, "Error: A file of Inputvector does not exist.\n");
    fclose(fp);
    exitMPI(-1);
  }
  fread(&i_max, sizeof(long int), 1, fp);
  if(i_max != X->Bind.Check.idim_max){
    fprintf(stderr, "Error: A file of Inputvector is incorrect.\n");
    fclose(fp);
    exitMPI(-1);
  }
  fread(v1, sizeof(complex double),X->Bind.Check.idim_max+1, fp);
  fclose(fp);
  
  //mltply Operator by using mltply.c (not yet modified) 
  mltply(&(X->Bind), v0, v1);
  
  //calculate norm
  dnorm = NormMPI_dc(i_max, v0);
  
  //normalize vector
#pragma omp parallel for default(none) private(i) shared(v1, v0) firstprivate(i_max, dnorm)
  for(i=1;i<=i_max;i++){
    v1[i] = v0[i]/dnorm;
  }

  int CalcSpecByLanczos=0;
  int iCalcSpecType=CalcSpecByLanczos;
  int iret=TRUE;
  switch (iCalcSpecType){
  case 0:
    iret = CalcSpectrumByLanczos(X, v1, dnorm);
    if(iret != TRUE){
      //Error Message will be added.
      return FALSE;
    }    
    break;    
    // case CalcSpecByShiftedKlyrov will be added
  default:
    break;
  }

  return TRUE;
}

int GetExcitedState
(
 struct BindStruct *X,
 double complex *tmp_v0
)
{
    double complex *tmp_v1 = v1;

	if(!GetSingleExcitedState(X,tmp_v0)==TRUE){
		return FALSE;
	}

	if(!GetPairExcitedState(X,tmp_v1, tmp_v0)==TRUE){
		return FALSE;
	}

    tmp_v0=tmp_v1;
  return TRUE;
}

int GetSingleExcitedState
		(
				struct BindStruct *X,
				double complex *vec
		){
	return TRUE;
}


int GetPairExcitedState
		(
				struct BindStruct *X,
                double complex *tmp_v0, /**< [out] Result v0 = H v1*/
                double complex *tmp_v1 /**< [in] v0 = H v1*/
		)
{
/*
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
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(tmp_v0, tmp_v1) \
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
#pragma omp parallel for default(none) reduction(+:dam_pr) shared(tmp_v0, tmp_v1)	\
  firstprivate(i_max) private(j)
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

#pragma omp parallel for default(none) shared(list_1, tmp_v0, tmp_v1) reduction(+:dam_pr) firstprivate(i_max, is, tmp_trans) private(num1, ibit)
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
				for(i=0;i<X->Def.NPairExcitationOperator;i++){
					org_isite1 = X->Def.PairExcitationOperator[i][0]+1;
					org_isite2 = X->Def.PairExcitationOperator[i][2]+1;
					org_sigma1 = X->Def.PairExcitationOperator[i][1];
					org_sigma2 = X->Def.PairExcitationOperator[i][3];
                    tmp_trans = X->Def.ParaPairExcitationOperator[i];
					if(org_sigma1 == org_sigma2){
						if(org_isite1==org_isite2){
							if(org_isite1 > X->Def.Nsite){
								is1_up = X->Def.Tpow[org_isite1 - 1];
								ibit1 = X_SpinGC_CisAis((unsigned long int)myrank + 1, X, is1_up, org_sigma1);
								if(ibit1 !=0){
#pragma omp parallel for reduction(+:dam_pr)default(none) shared(tmp_v0, tmp_v1)	\
  firstprivate(i_max, tmp_trans) private(j)
									for (j = 1; j <= i_max; j++) tmp_v0[j] += tmp_trans*tmp_v1[j];
								}
							}// org_isite1 > X->Def.Nsite
							else{
								isite1     = X->Def.Tpow[org_isite1-1];
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max, isite1, org_sigma1, X, tmp_trans) shared(tmp_v0, tmp_v1)
								for(j=1;j<=i_max;j++){
									tmp_v0[j] += X_Spin_CisAis(j,X, isite1,org_sigma1)*tmp_v1[j]*tmp_trans;
							}
						}
					}else{
						// for the canonical case
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
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j) firstprivate(i_max, tmp_trans) shared(tmp_v0,tmp_v1)
									for(j=1;j<=i_max;j++){
										tmp_v0[j]+= tmp_trans*tmp_v1[j];
									}
								}
							}
						}
						else {//org_isite1 <= X->Def.Nsite
							if(org_sigma1==org_sigma2){
								// longitudinal magnetic field
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, num1) firstprivate(i_max, org_isite1, org_sigma1, X, tmp_trans) shared(tmp_v0,tmp_v1, list_1)
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

					if(org_isite1 == org_isite2){
						if(org_isite1 > X->Def.Nsite){
							if(org_sigma1==org_sigma2){  // longitudinal magnetic field
								dam_pr += X_GC_child_CisAis_spin_MPIdouble(org_isite1-1, org_sigma1, 1.0, X, vec, vec);
							}
							else{  // transverse magnetic field
								dam_pr += X_GC_child_CisAit_spin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, 1.0, X, vec, vec);
							}
						}else{
							isite1 = X->Def.Tpow[org_isite1-1];

							if(org_sigma1==org_sigma2){
								// longitudinal magnetic field
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn) firstprivate(i_max, isite1, org_sigma1, X) shared(vec)
								for(j=1;j<=i_max;j++){
									dam_pr += X_SpinGC_CisAis(j, X, isite1, org_sigma1)*conj(vec[j])*vec[j];
								}
							}else{
								// transverse magnetic field
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, tmp_sgn, tmp_off) firstprivate(i_max, isite1, org_sigma2, X) shared(vec)
								for(j=1;j<=i_max;j++){
									tmp_sgn  =  X_SpinGC_CisAit(j,X, isite1,org_sigma2,&tmp_off);
									if(tmp_sgn !=0){
										dam_pr  +=  tmp_sgn*conj(vec[tmp_off+1])*vec[j];
									}
								}
							}
						}
					}else{
						// hopping is not allowed in localized spin system
						dam_pr=0.0;
					}
					dam_pr= X->Def.ParaPairExcitationOperator[i]*SumMPI_dc(dam_pr);
				}

			}//FlgGeneralSpin=FALSE
			else{
				for(i=0;i<X->Def.NPairExcitationOperator;i++){
					org_isite1 = X->Def.PairExcitationOperator[i][0]+1;
					org_isite2 = X->Def.PairExcitationOperator[i][2]+1;
					org_sigma1 = X->Def.PairExcitationOperator[i][1];
					org_sigma2 = X->Def.PairExcitationOperator[i][3];
					if(org_isite1 == org_isite2){
						if(org_isite1 > X->Def.Nsite){
							if(org_sigma1==org_sigma2){
								// longitudinal magnetic field
								dam_pr=X_GC_child_CisAis_GeneralSpin_MPIdouble(org_isite1-1, org_sigma1, 1.0, X, vec, vec);
							}else{
								// transverse magnetic field
								dam_pr=X_GC_child_CisAit_GeneralSpin_MPIdouble(org_isite1-1, org_sigma1, org_sigma2, 1.0, X, vec, vec);
							}
						}
						else{//org_isite1 <= X->Def.Nsite
							if(org_sigma1==org_sigma2){
								// longitudinal magnetic field
								dam_pr=0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, num1) firstprivate(i_max, org_isite1, org_sigma1, X) shared(vec)
								for(j=1;j<=i_max;j++){
									num1 = BitCheckGeneral(j-1, org_isite1, org_sigma1, X->Def.SiteToBit, X->Def.Tpow);
									dam_pr+=conj(vec[j])*vec[j]*num1;
								}
							}else{
								// transverse magnetic field
								dam_pr=0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, num1) firstprivate(i_max, org_isite1, org_sigma1, org_sigma2, X,tmp_off) shared(vec)
								for(j=1;j<=i_max;j++){
									num1 = GetOffCompGeneralSpin(j-1, org_isite1, org_sigma2, org_sigma1, &tmp_off, X->Def.SiteToBit, X->Def.Tpow);
									if(num1 !=0){
										dam_pr  +=  conj(vec[tmp_off+1])*vec[j]*num1;
									}
								}
							}
						}
					}else{
						// hopping is not allowed in localized spin system
						dam_pr=0.0;
					}
					dam_pr= X->Def.ParaPairExcitationOperator[i]*SumMPI_dc(dam_pr);
				}
			}
			break;

		default:
			return FALSE;
	}
*/
	return TRUE;
}