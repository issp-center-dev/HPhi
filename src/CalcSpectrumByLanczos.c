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
#include "mfmemory.h"

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
    int iret;
    unsigned long int liLanczosStp = X->Bind.Def.Lanczos_max;
    unsigned long int liLanczosStp_vec=0;

    if(X->Bind.Def.iFlgRecalcSpec == RECALC_FROM_TMComponents_VEC ||
       X->Bind.Def.iFlgRecalcSpec == RECALC_INOUT_TMComponents_VEC) {
        fprintf(stdoutMPI, "  Start: Input vectors for recalculation.\n");
        TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_InputSpectrumRecalcvecStart, "a");

        sprintf(sdt, cFileNameOutputRestartVec, X->Bind.Def.CDataFileHead, myrank);
        if (childfopenALL(sdt, "rb", &fp) != 0) {
            exitMPI(-1);
        }
        fread(&liLanczosStp_vec, sizeof(liLanczosStp_vec),1,fp);
        fread(&i_max, sizeof(long int), 1, fp);
        if(i_max != X->Bind.Check.idim_max){
            fprintf(stderr, "Error: A size of Inputvector is incorrect.\n");
            exitMPI(-1);
        }
        fread(v0, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
        fread(v1, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
        fclose(fp);
        fprintf(stdoutMPI, "  End:   Input vectors for recalculation.\n");
        TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_InputSpectrumRecalcvecEnd, "a");
    }

    //Read diagonal components
    if(X->Bind.Def.iFlgRecalcSpec == RECALC_FROM_TMComponents ||
       X->Bind.Def.iFlgRecalcSpec ==RECALC_FROM_TMComponents_VEC||
       X->Bind.Def.iFlgRecalcSpec == RECALC_INOUT_TMComponents_VEC)
    {
        iret=ReadTMComponents(X, &dnorm, &liLanczosStp);
        if(!iret ==TRUE){
            fprintf(stdoutMPI, "  Error: Fail to read TMcomponents\n");
            return FALSE;
        }

        if(X->Bind.Def.iFlgRecalcSpec == RECALC_FROM_TMComponents){
            X->Bind.Def.Lanczos_restart=0;
        }
        else if(X->Bind.Def.iFlgRecalcSpec == RECALC_INOUT_TMComponents_VEC||
                X->Bind.Def.iFlgRecalcSpec == RECALC_FROM_TMComponents_VEC){
            if(liLanczosStp_vec !=liLanczosStp){
                fprintf(stdoutMPI, "  Error: Input files for vector and TMcomponents are incoorect.\n");
                fprintf(stdoutMPI, "  Error: Input vector %ld th stps, TMcomponents  %ld th stps.\n", liLanczosStp_vec, liLanczosStp);
                return FALSE;
            }
            X->Bind.Def.Lanczos_restart=liLanczosStp;
            liLanczosStp = liLanczosStp+X->Bind.Def.Lanczos_max;
        }
    }

    // calculate ai, bi
    if (X->Bind.Def.iFlgRecalcSpec == RECALC_NOT ||
        X->Bind.Def.iFlgRecalcSpec == RECALC_OUTPUT_TMComponents_VEC ||
        X->Bind.Def.iFlgRecalcSpec == RECALC_FROM_TMComponents_VEC ||
        X->Bind.Def.iFlgRecalcSpec == RECALC_INOUT_TMComponents_VEC
        )
    {
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


    //output vectors for recalculation
    if(X->Bind.Def.iFlgRecalcSpec==RECALC_OUTPUT_TMComponents_VEC ||
       X->Bind.Def.iFlgRecalcSpec==RECALC_INOUT_TMComponents_VEC){
        fprintf(stdoutMPI, "    Start: Output vectors for recalculation.\n");
        TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_OutputSpectrumRecalcvecStart, "a");

        sprintf(sdt, cFileNameOutputRestartVec, X->Bind.Def.CDataFileHead, myrank);
        if(childfopenALL(sdt, "wb", &fp)!=0){
            exitMPI(-1);
        }
        fwrite(&liLanczosStp, sizeof(liLanczosStp),1,fp);
        fwrite(&X->Bind.Check.idim_max, sizeof(X->Bind.Check.idim_max),1,fp);
        fwrite(v0, sizeof(complex double),X->Bind.Check.idim_max+1, fp);
        fwrite(v1, sizeof(complex double),X->Bind.Check.idim_max+1, fp);
        fclose(fp);

        fprintf(stdoutMPI, "    End:   Output vectors for recalculation.\n");
        TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_OutputSpectrumRecalcvecEnd, "a");
    }

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
        double *_dnorm,
        unsigned long int *_i_max
){
    char sdt[D_FileNameMax];
    char ctmp[256], ctmp2[256];

    unsigned long int i, idx;
    unsigned long int i_max;
    double dnorm;
    FILE *fp;
    idx=1;
    sprintf(sdt, cFileNameTridiagonalMatrixComponents, X->Bind.Def.CDataFileHead);
    childfopenMPI(sdt,"r", &fp);

    fgetsMPI(ctmp, sizeof(ctmp)/sizeof(char), fp);
    sscanf(ctmp,"%ld \n", &i_max);
    if(X->Bind.Def.iFlgRecalcSpec == RECALC_INOUT_TMComponents_VEC||
       X->Bind.Def.iFlgRecalcSpec == RECALC_FROM_TMComponents_VEC) {
        alpha=(double*)realloc(alpha, sizeof(double)*(i_max + X->Bind.Def.Lanczos_max + 1));
        beta=(double*)realloc(beta, sizeof(double)*(i_max + X->Bind.Def.Lanczos_max + 1));
    }
    else if(X->Bind.Def.iFlgRecalcSpec==RECALC_FROM_TMComponents){
        alpha=(double*)realloc(alpha, sizeof(double)*(i_max + 1));
        beta=(double*)realloc(beta, sizeof(double)*(i_max + 1));
    }
    fgetsMPI(ctmp, sizeof(ctmp)/sizeof(char), fp);
    sscanf(ctmp,"%lf \n", &dnorm);
    while(fgetsMPI(ctmp, sizeof(ctmp)/sizeof(char), fp) != NULL){
        sscanf(ctmp,"%lf %lf \n", &alpha[idx], &beta[idx]);
        idx++;
    }
    fclose(fp);
    *_dnorm=dnorm;
    *_i_max=i_max;
    return TRUE;
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
