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
#include "CalcSpectrumByLanczos.h"
#include "Lanczos_EigenValue.h"
#include "FileIO.h"
#include "wrapperMPI.h"

/**
 * @file   CalcSpectrumByLanczos.c
 * @version 1.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 * @brief  File for givinvg functions of calculating spectrum by Lanczos
 * 
 * 
 */

int ReadTMComponents(
        struct EDMainCalStruct *X,
        double *_alpha,
        double *_beta,
        double *_dnorm,
        unsigned long int *i_max
);


int OutputTMComponents(
        struct EDMainCalStruct *X,
        double *_alpha,
        double *_beta,
        double _dnorm,
        int liLanczosStp
);

/** 
 * @brief A main function to calculate spectrum by Lanczos
 * 
 * @param[in,out] X CalcStruct list for getting and pushing calculation information 
 * @retval 0 normally finished
 * @retval -1 unnormally finished
 *
 * @version 1.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 */
int CalcSpectrumByLanczos(		 
			  struct EDMainCalStruct *X,
			  double complex *tmp_v1,
			  double dnorm,
  int Nomega,
  double complex *dcSpectrum,
  double complex *dcomega
				 ) {
    char sdt[D_FileNameMax];
    unsigned long int i, i_max;
    FILE *fp;
    int iret = TRUE;
    unsigned long int liLanczosStp = X->Bind.Def.Lanczos_max;

    //Read diagonal components
    if(X->Bind.Def.iFlgRecalcSpec == RECALC_FROM_TMComponents || X->Bind.Def.iFlgRecalcSpec ==RECALC_FROM_TMComponents_VEC){
        iret=ReadTMComponents(X, alpha, beta, &dnorm, &i_max);
        if(!iret ==TRUE){
            fprintf(stdoutMPI, "  Error: Fail to read TMcomponents\n");
            return FALSE;
        }

        if(X->Bind.Def.iFlgRecalcSpec == RECALC_FROM_TMComponents){
            liLanczosStp=i_max;
        }
    }

    // calculate ai, bi
    if (X->Bind.Def.iFlgRecalcSpec == RECALC_NOT ||
            X->Bind.Def.iFlgRecalcSpec == RECALC_FROM_TMComponents_VEC) {
        fprintf(stdoutMPI, "    Start: Calculate tridiagonal matrix components.\n");
        TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_GetTridiagonalStart, "a");
        // Functions in Lanczos_EigenValue
        iret = Lanczos_GetTridiagonalMatrixComponents(&(X->Bind), alpha, beta, tmp_v1, &(liLanczosStp));
        if (iret != TRUE) {
            //Error Message will be added.
            return FALSE;
        }
        fprintf(stdoutMPI, "    End:   Calculate tridiagonal matrix components.\n\n");
        TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_GetTridiagonalEnd, "a");
        OutputTMComponents(X, alpha,beta, dnorm, liLanczosStp);
    }//X->Bind.Def.iFlgRecalcSpec == RECALC_NOT || RECALC_FROM_TMComponents_VEC

    fprintf(stdoutMPI, "    Start: Caclulate spectrum from tridiagonal matrix components.\n");
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_CalcSpectrumFromTridiagonalStart, "a");
    for( i = 0 ; i < Nomega; i++) {
        iret = GetSpectrumByTridiagonalMatrixComponents(alpha, beta, dnorm, dcomega[i], &dcSpectrum[i], liLanczosStp);
        if (iret != TRUE) {
            //ToDo: Error Message will be added.
            //ReAlloc alpha, beta and Set alpha_start and beta_start in Lanczos_EigenValue
            return FALSE;
        }
        dcSpectrum[i] = dnorm * dcSpectrum[i];
    }
    fprintf(stdoutMPI, "    End:   Caclulate spectrum from tridiagonal matrix components.\n\n");
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_CalcSpectrumFromTridiagonalEnd, "a");

  return TRUE;
}


int GetSpectrumByTridiagonalMatrixComponents(
		double *tmp_alpha,
		double *tmp_beta,
        double dnorm,
		double complex dcomega,
		double complex *dcSpectrum,
        unsigned long int ilLanczosStp
		)
{
    unsigned long int istp=2;
    double complex dcDn;
    double complex dcb0;
    double complex dcbn, dcan;
    double complex dcDeltahn;
    double complex dch;

    if(ilLanczosStp < 1){
        //TODO: Add error message
        return FALSE;
    }

    dcb0 = dcomega-tmp_alpha[1];
    if(ilLanczosStp ==1) {
        if(cabs(dcb0)<eps_Energy){
            dcb0=eps_Energy;
        }
        *dcSpectrum = dnorm / (dcb0);
        return TRUE;
    }

    dcbn = dcomega-tmp_alpha[2];
    dcan = -pow(tmp_beta[1],2);
    dcDn=1.0/dcbn;
    dcDeltahn = dcan*dcDn;
    dch=dcb0+dcDeltahn;

    for(istp=2; istp<=ilLanczosStp; istp++){
        dcbn = dcomega-tmp_alpha[istp+1];
        dcan =-pow(tmp_beta[istp],2);
        dcDn = (dcbn+dcan*dcDn);
        if(cabs(dcDn)<eps_Energy){
            dcDn=eps_Energy;
        }
        dcDn =1.0/dcDn;
        dcDeltahn = (dcbn*dcDn-1.0)*dcDeltahn;
        dch += dcDeltahn;
        if(cabs(dcDeltahn/dch)<cabs(dcb0)*eps) break;
    }
    *dcSpectrum = dnorm/(dch);
  return TRUE;
}

int ReadTMComponents(
        struct EDMainCalStruct *X,
        double *_alpha,
        double *_beta,
        double *_dnorm,
        unsigned long int *_i_max
){
    char sdt[D_FileNameMax];
    char ctmp[256], ctmp2[256];

    unsigned long int i, idx;
    unsigned long int i_max;
    double dnorm;
    double alpha, beta;
    FILE *fp;
    idx=1;
    sprintf(sdt, cFileNameTridiagonalMatrixComponents, X->Bind.Def.CDataFileHead);
    childfopenMPI(sdt,"r", &fp);

    fgetsMPI(ctmp, sizeof(ctmp)/sizeof(char), fp);
    sscanf(ctmp,"%ld \n", &i_max);
    fgetsMPI(ctmp, sizeof(ctmp)/sizeof(char), fp);
    sscanf(ctmp,"%lf \n", &dnorm);
    while(fgetsMPI(ctmp, sizeof(ctmp)/sizeof(char), fp) != NULL){
        sscanf(ctmp,"%lf %lf \n", &_alpha[idx], &_beta[idx]);
        idx++;
    }
    fclose(fp);
    *_dnorm=dnorm;
    *_i_max=i_max;
}

int OutputTMComponents(
        struct EDMainCalStruct *X,
        double *_alpha,
        double *_beta,
        double _dnorm,
        int liLanczosStp
)
{
    char sdt[D_FileNameMax];
    unsigned long int i;
    FILE *fp;

    sprintf(sdt, cFileNameTridiagonalMatrixComponents, X->Bind.Def.CDataFileHead);
    childfopenMPI(sdt,"w", &fp);
    fprintf(fp, "%ld \n",liLanczosStp);
    fprintf(fp, "%.10lf \n",_dnorm);
    for( i = 1 ; i <= liLanczosStp; i++){
        fprintf(fp,"%.10lf %.10lf \n", alpha[i], beta[i]);
    }
    fclose(fp);

}
