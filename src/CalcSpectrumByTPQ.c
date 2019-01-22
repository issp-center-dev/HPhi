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
#include "CalcSpectrumByTPQ.h"
#include "Lanczos_EigenValue.h"
#include "FileIO.h"
#include "wrapperMPI.h"
#include "vec12.h"
#include "common/setmemory.h"
/**
 * @file   CalcSpectrumByTPQ.c
 * @version 1.2
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 * @brief  Calculate spectrum function for the TPQ state. \n
 * Note: This method is trial and cannot be used in the release mode.
 * 
 * 
 */

/// \brief Read TPQ data at "X->Bind.Large.itr" step in SS_rand file.
/// \param [in] X CalcStruct list for getting and pushing calculation information
/// \param [out] ene energy
/// \param [out] temp temperature
/// \param [out] specificHeat specific heat
/// \retval TRUE succeed to read data
/// \retval FALSE fail to read data
int ReadTPQData(
        struct EDMainCalStruct *X,
        double* ene,
        double* temp,
        double* specificHeat
){
    FILE *fp;
    char sdt[D_FileNameMax];
    char ctmp2[256];
    double dinv_temp;
    double dene, dHvar, dn, ddoublon;
    int istp;
    sprintf(sdt, cFileNameSSRand, X->Bind.Def.irand);
    childfopenMPI(sdt, "r", &fp);
    if(fp==NULL){
        fprintf(stderr, "  Error:  SS_rand%d.dat does not exist.\n", X->Bind.Def.irand);
        fclose(fp);
        return FALSE;
    }
    fgetsMPI(ctmp2, 256, fp);
    while(fgetsMPI(ctmp2, 256, fp) != NULL) {
        sscanf(ctmp2, "%lf %lf %lf %lf %lf %d\n",
               &dinv_temp,
               &dene,
               &dHvar,
               &dn,
               &ddoublon,
               &istp
        );
        if(istp==X->Bind.Large.itr) break;
    }
    fclose(fp);

    *ene = dene;
    *temp = 1.0/dinv_temp;
    *specificHeat = (dHvar-dene*dene)*(dinv_temp*dinv_temp);

    return TRUE;
}

///
/// \brief Calculate spectrum function from the TPQ state.
/// \param dcomega [in] Target frequencies.
/// \param dtemp [in] Temperature corresponding to the target TPQ state.
/// \param dspecificheat [in] Specific heat.
/// \param ene [in] Energy for the target TPQ state.
/// \param tmp_E [in] Energies included in the excited TPQ state obtained by the continued fraction expansions.
/// \param Nsite [in] Total number of sites.
/// \param idim_max [in] Dimension of the Hilbert space.
/// \param dc_tmpSpec [out] Calculated spectrum.
/// \retval FALSE fail to calculate spectrum.
/// \retval TRUE sucsceed to calculate spectrum.
int GetCalcSpectrumTPQ(double complex dcomega, double dtemp, double dspecificheat,
                       double ene, double *tmp_E, int Nsite, int idim_max, double complex * dc_tmpSpec)
{
    int l;
    double tmp_dcSpec;
    double factor, pre_factor;
    pre_factor=2.0*dtemp*dtemp*dspecificheat;
    factor=M_PI*pre_factor;
    factor=1.0/sqrt(factor);
    tmp_dcSpec=0;

    if(cimag(dcomega)>0) {
        for (l = 1; l <= idim_max; l++) {
            //TODO: Check omega is real ?
            //fprintf(stdoutMPI, "Debug: %lf, %lf\n", creal(dcomega) - tmp_E[l] + ene, pre_factor);
            tmp_dcSpec += (double)(vec[l][1] * conj(vec[l][1])) * exp(-pow((creal(dcomega) - tmp_E[l] + ene),2)/(pre_factor));
        }
    }
    else{
        fprintf(stderr, " an imarginary part of omega must be positive.\n");
        return FALSE;
    }
    tmp_dcSpec  *=factor;
    *dc_tmpSpec = tmp_dcSpec;
    return TRUE;
}

/// \brief A main function to calculate spectrum by TPQ (Note: This method is trial)
/// \param X [in,out] CalcStruct list for getting and pushing calculation information
/// \param tmp_v1 [in] Normalized excited state.
/// \param dnorm [in] Norm of the excited state before normalization.
/// \param Nomega [in] Total number of frequencies.
/// \param dcSpectrum [out] Calculated spectrum.
/// \param dcomega [in] Target frequencies.
/// \retval 0 normally finished
/// \retval -1 unnormally finished
/// \version 1.2
/// \author Kazuyoshi Yoshimi (The University of Tokyo)
int CalcSpectrumByTPQ(
        struct EDMainCalStruct *X,
        double complex *tmp_v1,
        double dnorm,
        int Nomega,
        double complex *dcSpectrum,
        double complex *dcomega
)
{
    char sdt[D_FileNameMax];
    unsigned long int i, i_max;
    FILE *fp;
    int iret;
    unsigned long int liLanczosStp = X->Bind.Def.Lanczos_max;
    unsigned long int liLanczosStp_vec=0;
    double dene, dtemp, dspecificHeat;
    double *tmp_E;
    double complex dctmp_Spectrum;
    int stp;
    size_t byte_size;

    //Read Ene, temp, C
    if(ReadTPQData(X, &dene, &dtemp, &dspecificHeat)!=TRUE){
        return FALSE;
    }

    if(X->Bind.Def.iFlgCalcSpec == RECALC_FROM_TMComponents_VEC ||
       X->Bind.Def.iFlgCalcSpec == RECALC_INOUT_TMComponents_VEC) {
        fprintf(stdoutMPI, "  Start: Input vectors for recalculation.\n");
        TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_InputSpectrumRecalcvecStart, "a");

        sprintf(sdt, cFileNameOutputRestartVec, X->Bind.Def.CDataFileHead, myrank);
        if (childfopenALL(sdt, "rb", &fp) != 0) {
            exitMPI(-1);
        }
        byte_size = fread(&liLanczosStp_vec, sizeof(liLanczosStp_vec),1,fp);
        byte_size = fread(&i_max, sizeof(long int), 1, fp);
        if(i_max != X->Bind.Check.idim_max){
            fprintf(stderr, "Error: A size of Inputvector is incorrect.\n");
            exitMPI(-1);
        }
        byte_size = fread(v0, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
        byte_size = fread(v1, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
        fclose(fp);
        fprintf(stdoutMPI, "  End:   Input vectors for recalculation.\n");
        TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_InputSpectrumRecalcvecEnd, "a");
        if (byte_size == 0) printf("byte_size: %d \n", (int)byte_size);
    }

    //Read diagonal components
    if(X->Bind.Def.iFlgCalcSpec == RECALC_FROM_TMComponents ||
       X->Bind.Def.iFlgCalcSpec ==RECALC_FROM_TMComponents_VEC||
       X->Bind.Def.iFlgCalcSpec == RECALC_INOUT_TMComponents_VEC)
    {
        int iFlgTMComp=0;
        if(X->Bind.Def.iFlgCalcSpec == RECALC_INOUT_TMComponents_VEC ||
           X->Bind.Def.iFlgCalcSpec ==  RECALC_FROM_TMComponents_VEC)
        {
            iFlgTMComp=0;
        }
        else{
            iFlgTMComp=1;
        }
        iret=ReadTMComponents(&(X->Bind), &dnorm, &liLanczosStp, iFlgTMComp);
        if(iret !=TRUE){
            fprintf(stdoutMPI, "  Error: Fail to read TMcomponents\n");
            return FALSE;
        }

        if(X->Bind.Def.iFlgCalcSpec == RECALC_FROM_TMComponents){
            X->Bind.Def.Lanczos_restart=0;
        }
        else if(X->Bind.Def.iFlgCalcSpec == RECALC_INOUT_TMComponents_VEC||
                X->Bind.Def.iFlgCalcSpec == RECALC_FROM_TMComponents_VEC){
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
    if (X->Bind.Def.iFlgCalcSpec == RECALC_NOT ||
        X->Bind.Def.iFlgCalcSpec == RECALC_OUTPUT_TMComponents_VEC ||
        X->Bind.Def.iFlgCalcSpec == RECALC_FROM_TMComponents_VEC ||
        X->Bind.Def.iFlgCalcSpec == RECALC_INOUT_TMComponents_VEC
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
        OutputTMComponents(&(X->Bind), alpha,beta, dnorm, liLanczosStp);
    }//X->Bind.Def.iFlgCalcSpec == RECALC_NOT || RECALC_FROM_TMComponents_VEC

    stp=liLanczosStp;
    tmp_E = d_1d_allocate(stp + 1);
    X->Bind.Def.nvec= stp;
    vec12(alpha,beta,stp,tmp_E, &(X->Bind));
    fprintf(stdoutMPI, "    Start: Caclulate spectrum from tridiagonal matrix components.\n");
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_CalcSpectrumFromTridiagonalStart, "a");
    for( i = 0 ; i < Nomega; i++) {
        dctmp_Spectrum=0;
        iret = GetCalcSpectrumTPQ(dcomega[i], dtemp, dspecificHeat, dene, tmp_E, X->Bind.Def.Nsite, stp, &dctmp_Spectrum);
        if (iret != TRUE) {
            //ReAlloc alpha, beta and Set alpha_start and beta_start in Lanczos_EigenValue
            return FALSE;
        }
        dcSpectrum[i] = dnorm * dctmp_Spectrum;
    }
    fprintf(stdoutMPI, "    End:   Caclulate spectrum from tridiagonal matrix components.\n\n");
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_CalcSpectrumFromTridiagonalEnd, "a");

    free_d_1d_allocate(tmp_E);
    //output vectors for recalculation
    if(X->Bind.Def.iFlgCalcSpec==RECALC_OUTPUT_TMComponents_VEC ||
       X->Bind.Def.iFlgCalcSpec==RECALC_INOUT_TMComponents_VEC){
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
