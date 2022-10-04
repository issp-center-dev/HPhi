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
#include "common/setmemory.h"
#include "CalcTime.h"
/**
 * @file   CalcSpectrumByLanczos.c
 * @brief Get the spectrum function by continued fraction expansions.\n
 * Ref. E.Dagotto, Rev. Mod. Phys. 66 (1994), 763.
 * @version 1.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 */

///
/// \brief A main function to calculate spectrum by continued fraction expansions.
/// \param X [in,out] Struct for getting and giving calculation information
/// \param tmp_v1 [in] Normalized excited state.
/// \param dnorm [in] Norm of the excited state before normalization.
/// \param Nomega [in] Total number of frequencies.
/// \param dcSpectrum [out] Calculated spectrum.
/// \param dcomega [in] Target frequencies.
/// \retval 0 normally finished
/// \retval -1 unnormally finished
/// \version 1.1
/// \author Kazuyoshi Yoshimi (The University of Tokyo)
int CalcSpectrumByLanczos(
        struct EDMainCalStruct *X,
        double complex *tmp_v1,
        double dnorm,
        int Nomega,
        double complex *dcSpectrum,
        double complex *dcomega
)
{
    unsigned long int i;
    int iret;
    unsigned long int liLanczosStp = X->Bind.Def.Lanczos_max;
    unsigned long int liLanczosStp_vec=0;

    if(X->Bind.Def.iFlgCalcSpec == RECALC_FROM_TMComponents_VEC ||
       X->Bind.Def.iFlgCalcSpec == RECALC_INOUT_TMComponents_VEC) {

        StartTimer(6201);
        if(ReadInitialVector( &(X->Bind), v0, v1, &liLanczosStp_vec)!=0){
            StopTimer(6201);
            exitMPI(-1);
        }
        StopTimer(6201);
    }

    //Read diagonal components
    if(X->Bind.Def.iFlgCalcSpec == RECALC_FROM_TMComponents ||
       X->Bind.Def.iFlgCalcSpec ==RECALC_FROM_TMComponents_VEC||
       X->Bind.Def.iFlgCalcSpec == RECALC_INOUT_TMComponents_VEC)
    {
        StartTimer(6202);
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
        StopTimer(6202);
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
        StartTimer(6203);
        iret = Lanczos_GetTridiagonalMatrixComponents(&(X->Bind), alpha, beta, tmp_v1, &(liLanczosStp));
        StopTimer(6203);
        if (iret != TRUE) {
            //Error Message will be added.
            return FALSE;
        }
        fprintf(stdoutMPI, "    End:   Calculate tridiagonal matrix components.\n\n");
        TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_GetTridiagonalEnd, "a");
        StartTimer(6204);
        OutputTMComponents(&(X->Bind), alpha,beta, dnorm, liLanczosStp);
        StopTimer(6204);
    }//X->Bind.Def.iFlgCalcSpec == RECALC_NOT || RECALC_FROM_TMComponents_VEC

    fprintf(stdoutMPI, "    Start: Caclulate spectrum from tridiagonal matrix components.\n");
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_CalcSpectrumFromTridiagonalStart, "a");
    StartTimer(6205);
    for( i = 0 ; i < Nomega; i++) {
        iret = GetSpectrumByTridiagonalMatrixComponents(alpha, beta, dnorm, dcomega[i], &dcSpectrum[i], liLanczosStp);
        if (iret != TRUE) {
            //ToDo: Error Message will be added.
            //ReAlloc alpha, beta and Set alpha_start and beta_start in Lanczos_EigenValue
            return FALSE;
        }
        dcSpectrum[i] = dnorm * dcSpectrum[i];
    }
    StopTimer(6205);
    fprintf(stdoutMPI, "    End:   Caclulate spectrum from tridiagonal matrix components.\n\n");
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_CalcSpectrumFromTridiagonalEnd, "a");

    //output vectors for recalculation
    if(X->Bind.Def.iFlgCalcSpec==RECALC_OUTPUT_TMComponents_VEC ||
       X->Bind.Def.iFlgCalcSpec==RECALC_INOUT_TMComponents_VEC){
      StartTimer(6206);
      if(OutputLanczosVector(&(X->Bind), v0, v1, liLanczosStp)!=0){
        StopTimer(6206);
        exitMPI(-1);
      }
        StopTimer(6206);
    }

    return TRUE;
}

///
/// \brief Calculate the spectrum by using tridiagonal matrix components obtained by the Lanczos_GetTridiagonalMatrixComponents function.
/// \param tmp_alpha [in] Tridiagonal matrix components.
/// \param tmp_beta [in] Tridiagonal matrix components.
/// \param dnorm [in] Norm for the excited state.
/// \param dcomega [in] Target frequency.
/// \param dcSpectrum [out] Spectrum at dcomega.
/// \param ilLanczosStp [in] Lanczos step required to get tridiagonal matrix components.
/// \retval FALSE Fail to get the spectrum
/// \retval TRUE  Success to get the spectrum
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
        fprintf(stdoutMPI, "Error: LanczosStep must be greater than 1.\n");
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
