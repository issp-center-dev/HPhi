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
#include "CalcSpectrumByFullDiag.h"
#include "SingleEx.h"
#include "PairEx.h"
#include "wrapperMPI.h"
#include "FileIO.h"
#include "mfmemory.h"

/**
 * @file   CalcSpectrum.c
 * @version 1.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 *
 * @brief  File for givinvg functions of calculating spectrum
 *
 *
 */

int OutputSpectrum(
  struct EDMainCalStruct *X,
  int Nomega,
  double complex *dcSpectrum,
  double complex *dcomega) 
{
  FILE *fp;
  char sdt[D_FileNameMax];
  int i;

  //output spectrum
  sprintf(sdt, cFileNameCalcDynamicalGreen, X->Bind.Def.CDataFileHead);
  childfopenMPI(sdt, "w", &fp);

  for (i = 0; i < Nomega; i++) {
    fprintf(fp, "%.10lf %.10lf %.10lf %.10lf \n",
      creal(dcomega[i]), cimag(dcomega[i]),
      creal(dcSpectrum[i]), cimag(dcSpectrum[i]));
  }/*for (i = 0; i < Nomega; i++)*/

  fclose(fp);
}/*int OutputSpectrum*/

/**
 * @brief A main function to calculate spectrum
 *
 * @param[in,out] X CalcStruct list for getting and pushing calculation information
 * @retval 0 normally finished
 * @retval -1 unnormally finished
 *
 * @version 1.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Youhei Yamaji (The University of Tokyo)
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

  //ToDo: Nomega should be given as a parameter
  int Nomega;
  double complex OmegaMax, OmegaMin;
  double complex *dcSpectrum;
  double complex *dcomega;

  //set omega
  if(SetOmega(&(X->Bind.Def)) != TRUE){
    fprintf(stderr, "Error: Fail to set Omega.\n");
    exitMPI(-1);
  }
  else{
    if(X->Bind.Def.iFlgSpecOmegaIm == FALSE){
      X->Bind.Def.dOmegaIm = (X->Bind.Def.dOmegaMax - X->Bind.Def.dOmegaMin)/(double) X->Bind.Def.iNOmega;
    }
  }
  /*
   Set & malloc omega grid
  */
  Nomega = X->Bind.Def.iNOmega;
  c_malloc1(dcSpectrum, Nomega);
  c_malloc1(dcomega, Nomega);
  OmegaMax = X->Bind.Def.dOmegaMax + X->Bind.Def.dOmegaIm*I;
  OmegaMin = X->Bind.Def.dOmegaMin + X->Bind.Def.dOmegaIm*I;
  for (i = 0; i< Nomega; i++) {
    dcomega[i] = (OmegaMax - OmegaMin) / Nomega*i + OmegaMin;
  }

  if(X->Bind.Def.NSingleExcitationOperator == 0 && X->Bind.Def.NPairExcitationOperator == 0){
      fprintf(stderr, "Error: Any excitation operators are not defined.\n");
      exitMPI(-1);
  }

  //input eigen vector
  fprintf(stdoutMPI, "  Start: An Eigenvector is inputted in CalcSpectrum.\n");
  TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_InputEigenVectorStart, "a");
  sprintf(sdt, cFileNameInputEigen, X->Bind.Def.CDataFileHead, X->Bind.Def.k_exct-1, myrank);
  childfopenALL(sdt, "rb", &fp);
  if(fp==NULL){
    fprintf(stderr, "Error: A file of Inputvector does not exist.\n");
    exitMPI(-1);
  }
  //  fscanf(fp, "%ld\n", &i_max);
  fread(&i_max, sizeof(i_max), 1 ,fp);
  if(i_max != X->Bind.Check.idim_max){
    fprintf(stderr, "Error: myrank=%d, i_max=%ld\n", myrank, i_max);
    fprintf(stderr, "Error: A file of Inputvector is incorrect.\n");
    //exitMPI(-1);
    return -1;
  }
  fread(v1, sizeof(complex double), X->Bind.Check.idim_max+1, fp);
  fclose(fp);
  for (i = 1; i <= i_max; i++) {
    v0[i]=0;
  }
  fprintf(stdoutMPI, "  End:   An Eigenvector is inputted in CalcSpectrum.\n\n");
  TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_InputEigenVectorEnd, "a");

  TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_CalcExcitedStateStart, "a");
  fprintf(stdoutMPI, "  Start: Calculating an excited Eigenvector.\n");
  //mltply Operator
  GetExcitedState( &(X->Bind), v0, v1);

  //calculate norm
  dnorm = NormMPI_dc(i_max, v0);
  if(fabs(dnorm)<pow(10.0, -15)){
    fprintf(stderr, "Error: Excitation vector is illegal; norm becomes 0.\n");
    return -1;
  }

  //normalize vector
#pragma omp parallel for default(none) private(i) shared(v1, v0) firstprivate(i_max, dnorm)
  for(i=1;i<=i_max;i++){
    v1[i] = v0[i]/dnorm;
  }
  fprintf(stdoutMPI, "  End:   Calculating an excited Eigenvector.\n\n");

  TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_CalcExcitedStateEnd, "a");

  int iret=TRUE;

  fprintf(stdoutMPI, "  Start: Caclulating a spectrum.\n\n");
  TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_CalcSpectrumStart, "a");
  switch (X->Bind.Def.iCalcType) {

  case Spectrum:
    iret = CalcSpectrumByLanczos(X, v1, dnorm, Nomega, dcSpectrum, dcomega);
    if (iret != TRUE) {
      //Error Message will be added.
      return FALSE;
    }
    break;

  case SpectrumFD:
    iret = CalcSpectrumByFullDiag(X, Nomega,dcSpectrum,dcomega);
    break;

    // case CalcSpecByShiftedKlyrov will be added
  default:
    break;
  }

  if (iret != TRUE) {
    //Error Message will be added.
    return FALSE;
  }
  else {
    fprintf(stdoutMPI, "  End:  Calculating a spectrum.\n\n");
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_CalcSpectrumEnd, "a");
    iret = OutputSpectrum(X, Nomega, dcSpectrum, dcomega);
    return TRUE;
  }

}/*int CalcSpectrum*/

int GetExcitedState
(
 struct BindStruct *X,
 double complex *tmp_v0,/**< [out] Result v0 = H v1*/
 double complex *tmp_v1 /**< [in] v0 = H v1*/
 )
{
   if(X->Def.NSingleExcitationOperator > 0 && X->Def.NPairExcitationOperator > 0){
    fprintf(stderr, "Error: Both single and pair excitation operators exist.\n");
    return FALSE;
    }

    if(X->Def.NSingleExcitationOperator > 0){
  if(!GetSingleExcitedState(X,tmp_v0, tmp_v1)==TRUE){
    return FALSE;
  }
    }
else if(X->Def.NPairExcitationOperator >0){
  if(!GetPairExcitedState(X,tmp_v0, tmp_v1)==TRUE){
    return FALSE;
  }
    }

  return TRUE;
}

int SetOmega
(
 struct DefineList *X
){
  FILE *fp;
  char sdt[D_FileNameMax],ctmp[256];
  double domegaMax;
  double domegaMin;
  int istp=4;
  double E1, E2, E3, E4, Emax;
    long unsigned int iline_countMax=2;
    long unsigned int iline_count=2;


  if(X->iFlgSpecOmegaMax == TRUE && X->iFlgSpecOmegaMin == TRUE){
    return TRUE;
  }
  else{
    sprintf(sdt, cFileNameLanczosStep, X->CDataFileHead);
      childfopenMPI(sdt,"r", &fp);
    if(fp == NULL){
      fprintf(stdoutMPI, "Error: xx_Lanczos_Step.dat does not exist.\n");
      return FALSE;
    }
      fgetsMPI(ctmp, 256, fp); //1st line is skipped
      fgetsMPI(ctmp, 256, fp); //2nd line is skipped
      while(fgetsMPI(ctmp, 256, fp) != NULL){
          iline_count++;
      }
    iline_countMax=iline_count;
    iline_count=2;
    rewind(fp);
      fgetsMPI(ctmp, 256, fp); //1st line is skipped
      fgetsMPI(ctmp, 256, fp); //2nd line is skipped

      while(fgetsMPI(ctmp, 256, fp) != NULL) {
          sscanf(ctmp, "stp=%d %lf %lf %lf %lf %lf\n",
                 &istp,
                 &E1,
                 &E2,
                 &E3,
                 &E4,
                 &Emax);
          iline_count++;
          if(iline_count ==iline_countMax) break;
      }


      if(istp < 4){
      fprintf(stdoutMPI, "Error: Lanczos step must be greater than 4 for using spectrum calculation.\n");
      return FALSE;
    }
    //Read Lanczos_Step
    if(X->iFlgSpecOmegaMax == FALSE){
      X->dOmegaMax= Emax*(double)X->NsiteMPI;
    }
    if(X->iFlgSpecOmegaMin == FALSE){
      X->dOmegaMin= E1;
    }
  }

  return TRUE;
}
