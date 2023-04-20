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
#include "CalcSpectrum.h"
#include "CalcSpectrumByBiCG.h"
#include "CalcSpectrumByFullDiag.h"
#include "CalcTime.h"
#include "SingleEx.h"
#include "PairEx.h"
#include "wrapperMPI.h"
#include "FileIO.h"
#include "common/setmemory.h"
#include "readdef.h"
#include "sz.h"
#include "check.h"
#include "diagonalcalc.h"
/**
 * @file   CalcSpectrum.c
 * @version 1.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 *
 * @brief  File for givinvg functions of calculating spectrum
 *
 *
 */

 ///
 /// \brief Set target frequencies
 /// \param X [in, out] Struct to give and get the information of target frequencies.\n
 /// Output: dcOmegaMax, dcOmegaMin
 ///
 /// \retval FALSE Fail to set frequencies.
 /// \retval TRUE Success to set frequencies.
int SetOmega
(
  struct DefineList *X
) {
  FILE *fp;
  char sdt[D_FileNameMax], ctmp[256];
  int istp = 4;
  double E1, E2, E3, E4, Emax;
  long unsigned int iline_countMax = 2;
  long unsigned int iline_count = 2;


  if (X->iFlgSpecOmegaMax == TRUE && X->iFlgSpecOmegaMin == TRUE) {
    return TRUE;
  }
  else {
    if (X->iCalcType == Lanczos || X->iCalcType == FullDiag) {
      sprintf(sdt, cFileNameLanczosStep, X->CDataFileHead);
      childfopenMPI(sdt, "r", &fp);
      if (fp == NULL) {
        fprintf(stdoutMPI, "Error: xx_Lanczos_Step.dat does not exist.\n");
        return FALSE;
      }
      fgetsMPI(ctmp, 256, fp); //1st line is skipped
      fgetsMPI(ctmp, 256, fp); //2nd line is skipped
      while (fgetsMPI(ctmp, 256, fp) != NULL) {
        iline_count++;
      }
      iline_countMax = iline_count;
      iline_count = 2;
      rewind(fp);
      fgetsMPI(ctmp, 256, fp); //1st line is skipped
      fgetsMPI(ctmp, 256, fp); //2nd line is skipped

      while (fgetsMPI(ctmp, 256, fp) != NULL) {
        sscanf(ctmp, "stp=%d %lf %lf %lf %lf %lf\n",
          &istp,
          &E1,
          &E2,
          &E3,
          &E4,
          &Emax);
        iline_count++;
        if (iline_count == iline_countMax) break;
      }
      fclose(fp);
      if (istp < 4) {
        fprintf(stdoutMPI, "Error: Lanczos step must be greater than 4 for using spectrum calculation.\n");
        return FALSE;
      }
    }/*if (X->iCalcType == Lanczos || X->iCalcType == FullDiag)*/
    else
    {
      sprintf(sdt, cFileNameEnergy_Lanczos, X->CDataFileHead);
      childfopenMPI(sdt, "r", &fp);
      if (fp == NULL) {
        fprintf(stdoutMPI, "Error: xx_energy.dat does not exist.\n");
        return FALSE;
      }/*if (fp == NULL)*/
      fgetsMPI(ctmp, 256, fp); //1st line is skipped
      fgetsMPI(ctmp, 256, fp); //1st line is skipped
      sscanf(ctmp, "  Energy  %lf \n", &E1);
      Emax = LargeValue;
    }/**/
    //Read Lanczos_Step
    if (X->iFlgSpecOmegaMax == FALSE) {
      X->dcOmegaMax = Emax * (double)X->Nsite - E1;
    }
    if (X->iFlgSpecOmegaMin == FALSE) {
      X->dcOmegaMin = 0.0;
    }
  }/*Omegamax and omegamin is not specified in modpara*/

  return TRUE;
}
///
/// \brief Make the lists for the excited state; list_1, list_2_1 and list_2_2 (for canonical ensemble).
/// The original lists before the excitation are given by list_xxx_org
/// \param X [in, out] Struct to get and give information to make the lists for the excited state.\n
/// Output: iCalcModel (From HubbardNConserved to Hubbard), {Ne, Nup, Ndown, Nsite, Total2Sz} (update for MPI)
///
/// \param iFlgListModifed [out] If the list is modified due to the excitation, the value becomes TRUE(1), otherwise FALSE(0).
/// \retval -1 fail to make lists.
/// \retval 0  sucsess to make lists.
int MakeExcitedList(
  struct BindStruct *X,
  int *iFlgListModifed
) {
  long int j;
  *iFlgListModifed = FALSE;
  //To Get Original space
  if (check(X) == MPIFALSE) {
    FinalizeMPI();
    return -1;
  }

  X->Check.idim_maxOrg = X->Check.idim_max;
  X->Check.idim_maxMPIOrg = X->Check.idim_maxMPI;

  if (X->Def.NNSingleExcitationOperator > 0) {
    switch (X->Def.iCalcModel) {
    case HubbardGC:
      break;
    case HubbardNConserved:
    case KondoGC:
    case Hubbard:
    case Kondo:
      *iFlgListModifed = TRUE;
      break;
    case Spin:
    case SpinGC:
      return FALSE;
    }
  }
  else if (X->Def.NNPairExcitationOperator > 0) {
    switch (X->Def.iCalcModel) {
    case HubbardGC:
    case SpinGC:
    case HubbardNConserved:
      break;
    case KondoGC:
    case Hubbard:
    case Kondo:
    case Spin:
      if (X->Def.PairExcitationOperator[0][0][1] != X->Def.PairExcitationOperator[0][0][3]) {
        *iFlgListModifed = TRUE;
      }
      break;
    }
  }
  else {
    return FALSE;
  }

  if (*iFlgListModifed == TRUE) {
    if (GetlistSize(X) == TRUE) {
      list_1_org = lui_1d_allocate(X->Check.idim_max + 1);
#ifdef MPI
      unsigned long int MAXidim_max;
      MAXidim_max = MaxMPI_li(X->Check.idim_max);
      list_1buf_org = lui_1d_allocate(MAXidim_max + 1);
#endif // MPI
      list_2_1_org = lui_1d_allocate(X->Large.SizeOflist_2_1);
      list_2_2_org = lui_1d_allocate(X->Large.SizeOflist_2_2);
      if (list_1_org == NULL
        || list_2_1_org == NULL
        || list_2_2_org == NULL
        )
      {
        return -1;
      }
      for (j = 0; j < X->Large.SizeOflist_2_1; j++) {
        list_2_1_org[j] = 0;
      }
      for (j = 0; j < X->Large.SizeOflist_2_2; j++) {
        list_2_2_org[j] = 0;
      }

    }

    if (sz(X, list_1_org, list_2_1_org, list_2_2_org) != 0) {
      return FALSE;
    }

    if (X->Def.NNSingleExcitationOperator > 0) {
      switch (X->Def.iCalcModel) {
      case HubbardGC:
        break;
      case HubbardNConserved:
        if (X->Def.SingleExcitationOperator[0][0][2] == 1) { //cis
          X->Def.Ne = X->Def.NeMPI + 1;
        }
        else {
          X->Def.Ne = X->Def.NeMPI - 1;
        }
        break;
      case KondoGC:
      case Hubbard:
      case Kondo:
        if (X->Def.SingleExcitationOperator[0][0][2] == 1) { //cis
          X->Def.Ne = X->Def.NeMPI + 1;
          if (X->Def.SingleExcitationOperator[0][0][1] == 0) {//up
            X->Def.Nup = X->Def.NupOrg + 1;
            X->Def.Ndown = X->Def.NdownOrg;
          }
          else {//down
            X->Def.Nup = X->Def.NupOrg;
            X->Def.Ndown = X->Def.NdownOrg + 1;
          }
        }
        else {//ajt
          X->Def.Ne = X->Def.NeMPI - 1;
          if (X->Def.SingleExcitationOperator[0][0][1] == 0) {//up
            X->Def.Nup = X->Def.NupOrg - 1;
            X->Def.Ndown = X->Def.NdownOrg;

          }
          else {//down
            X->Def.Nup = X->Def.NupOrg;
            X->Def.Ndown = X->Def.NdownOrg - 1;
          }
        }
        break;
      case Spin:
      case SpinGC:
        return FALSE;
      }
    }
    else if (X->Def.NNPairExcitationOperator > 0) {
      X->Def.Ne = X->Def.NeMPI;
      switch (X->Def.iCalcModel) {
      case HubbardGC:
      case SpinGC:
      case HubbardNConserved:
        break;
      case KondoGC:
      case Hubbard:
      case Kondo:
        if (X->Def.PairExcitationOperator[0][0][1] != X->Def.PairExcitationOperator[0][0][3]) {
          if (X->Def.PairExcitationOperator[0][0][1] == 0) {//up
            X->Def.Nup = X->Def.NupOrg + 1;
            X->Def.Ndown = X->Def.NdownOrg - 1;
          }
          else {//down
            X->Def.Nup = X->Def.NupOrg - 1;
            X->Def.Ndown = X->Def.NdownOrg + 1;
          }
        }
        break;
      case Spin:
        if (X->Def.PairExcitationOperator[0][0][1] != X->Def.PairExcitationOperator[0][0][3]) {
          if (X->Def.iFlgGeneralSpin == FALSE) {
            if (X->Def.PairExcitationOperator[0][0][1] == 0) {//down
              X->Def.Nup = X->Def.NupOrg - 1;
              X->Def.Ndown = X->Def.NdownOrg + 1;
            }
            else {//up
              X->Def.Nup = X->Def.NupOrg + 1;
              X->Def.Ndown = X->Def.NdownOrg - 1;
            }
          }
          else {//for general spin
            X->Def.Total2Sz = X->Def.Total2SzMPI + 2 * (X->Def.PairExcitationOperator[0][0][1] - X->Def.PairExcitationOperator[0][0][3]);
          }
        }
        break;
      }
    }
    else {
      return FALSE;
    }
    //Update Infomation
    X->Def.Nsite = X->Def.NsiteMPI;

    if (check(X) == MPIFALSE) {
      FinalizeMPI();
      return FALSE;
    }
  }

  //set memory
  if (setmem_large(X) != 0) {
    fprintf(stdoutMPI, cErrLargeMem, iErrCodeMem);
    exitMPI(-1);
  }
  if (sz(X, list_1, list_2_1, list_2_2) != 0) {
    return FALSE;
  }
#ifdef MPI
  unsigned long int MAXidim_max, MAXidim_maxOrg;
  MAXidim_max = MaxMPI_li(X->Check.idim_max);
  MAXidim_maxOrg = MaxMPI_li(X->Check.idim_maxOrg);
  if (MAXidim_max < MAXidim_maxOrg) {
    free_cd_2d_allocate(v1buf);
    v1buf = cd_2d_allocate(MAXidim_maxOrg + 1, 1);
  }
#endif // MPI
  
  if (X->Def.iCalcModel == HubbardNConserved) {
    X->Def.iCalcModel = Hubbard;
  }

#ifdef _DEBUG
  if (*iFlgListModifed == TRUE) {
    for (j = 1; j <= X->Check.idim_maxOrg; j++) {
      fprintf(stdout, "Debug1: myrank=%d, list_1_org[ %ld] = %ld\n",
        myrank, j, list_1_org[j] + myrank * X->Def.OrgTpow[2 * X->Def.NsiteMPI - 1]);
    }

    for (j = 1; j <= X->Check.idim_max; j++) {
      fprintf(stdout, "Debug2: myrank=%d, list_1[ %ld] = %ld\n", myrank, j, list_1[j] + myrank * 64);
    }
  }
#endif
  return TRUE;
}
/// \brief Output spectrum.
/// \param X [in] Read information of the frequency origin.
/// \param Nomega [in] A total number of discrete frequencies.
/// \param dcSpectrum [in] Array of spectrum.
/// \param dcomega [in] Array of discrete frequencies.
/// \retval FALSE Fail to open the output file.
/// \retval TRUE Success to output the spectrum.
int OutputSpectrum(
  struct EDMainCalStruct *X,
  int nstate,
  int Nomega,
  int NdcSpectrum,
  double complex ***dcSpectrum,
  double complex **dcomega,
  double *energy0) 
{
  FILE *fp;
  char sdt[D_FileNameMax];
  int iomega, idcSpectrum, istate;

  //output spectrum
  for (istate = 0; istate < nstate; istate++) {
    sprintf(sdt, cFileNameCalcDynamicalGreen, X->Bind.Def.CDataFileHead, istate);
    if (childfopenMPI(sdt, "w", &fp) != 0) {
      return FALSE;
    }

    for (idcSpectrum = 0; idcSpectrum < NdcSpectrum; idcSpectrum++) {
      for (iomega = 0; iomega < Nomega; iomega++) {
        fprintf(fp, "%.10lf %.10lf %.10lf %.10lf \n",
          creal(dcomega[istate][iomega] - X->Bind.Def.dcOmegaOrg - energy0[istate]),
          cimag(dcomega[istate][iomega] - X->Bind.Def.dcOmegaOrg),
          creal(dcSpectrum[istate][iomega][idcSpectrum]), cimag(dcSpectrum[istate][iomega][idcSpectrum]));
      }/*for (i = 0; i < Nomega; i++)*/
      fprintf(fp, "\n");
    }

    fclose(fp);
  }
  return TRUE;
}/*int OutputSpectrum*/
/// \brief Parent function to calculate the excited state.
/// \param X [in] Struct to get number of excitation operators.
/// \param tmp_v0 [out] Result @f$ v_0 = H_{ex} v_1 @f$.
/// \param tmp_v1 [in] The original state before excitation  @f$ v_1 @f$.
/// \retval FALSE Fail to calculate the excited state.
/// \retval TRUE Success to calculate the excited state.
int GetExcitedState
(
  struct BindStruct *X,
  int nstate,
  double complex **tmp_v0,
  double complex **tmp_v1,
  int iEx
)
{
  if (X->Def.NNSingleExcitationOperator > 0) {
    if (GetSingleExcitedState(X, nstate, tmp_v0, tmp_v1, iEx) != TRUE) {
      return FALSE;
    }
  }
  else if (X->Def.NNPairExcitationOperator > 0) {
    if (GetPairExcitedState(X, nstate, tmp_v0, tmp_v1, iEx) != TRUE) {
      return FALSE;
    }
  }
  return TRUE;
}
/**
 * @brief A main function to calculate spectrum.
 *
 * @param X [in,out] CalcStruct list for getting and pushing calculation information \n
 * input: iFlgSpecOmegaOrg, dcOmegaMax, dcOmegaMin, iNOmega etc.\n
 * output: dcOmegaOrg, iFlagListModified.
 *
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
) {
  char sdt[D_FileNameMax];
  unsigned long int i;
  unsigned long int i_max = 0;
  int i_stp, NdcSpectrum;
  int iFlagListModified = FALSE;
  FILE *fp;
  double complex **v1Org, *vin; /**< Input vector to calculate spectrum function.*/

  //ToDo: Nomega should be given as a parameter
  int Nomega, nstate, istate;
  double complex OmegaMax, OmegaMin;
  double complex ***dcSpectrum;
  double complex **dcomega;
  double *energy0;
  size_t byte_size;

  if (X->Bind.Def.iCalcType == FullDiag) {
    nstate = X->Bind.Check.idim_max;
  }
  else if (X->Bind.Def.iCalcType == CG) {
    nstate = X->Bind.Def.k_exct;
  }
  else if (X->Bind.Def.iCalcType == TPQCalc || X->Bind.Def.iCalcType == cTPQ) {
    nstate = NumAve;
  }
  else {
    nstate = 1;
  }
  energy0 = d_1d_allocate(nstate);

  if (X->Bind.Def.iFlgCalcSpec == CALCSPEC_SCRATCH) {
    X->Bind.Def.Nsite = X->Bind.Def.NsiteMPI;
    X->Bind.Def.Total2Sz = X->Bind.Def.Total2SzMPI;
    X->Bind.Def.Ne = X->Bind.Def.NeMPI;
    X->Bind.Def.Nup = X->Bind.Def.NupMPI;
    X->Bind.Def.Ndown = X->Bind.Def.NdownMPI;
    if (GetlistSize(&(X->Bind)) == TRUE) {
      free_lui_1d_allocate(list_1);
      free_lui_1d_allocate(list_2_1);
      free_lui_1d_allocate(list_2_2);
    }
    free_d_1d_allocate(list_Diagonal);
    free_cd_2d_allocate(v0);
    v1Org = cd_2d_allocate(X->Bind.Check.idim_max + 1, nstate);
    for (i = 1; i <= X->Bind.Check.idim_max; i++) 
      for (istate =0;istate < nstate;istate ++)
        v1Org[i][istate] = v1[i][istate];
    free_cd_2d_allocate(v1);
#ifdef MPI
    free_lui_1d_allocate(list_1buf);
    free_cd_2d_allocate(v1buf);
#endif // MPI
    free_d_1d_allocate(X->Bind.Phys.num_down);
    free_d_1d_allocate(X->Bind.Phys.num_up);
    free_d_1d_allocate(X->Bind.Phys.num);
    free_d_1d_allocate(X->Bind.Phys.num2);
    free_d_1d_allocate(X->Bind.Phys.energy);
    free_d_1d_allocate(X->Bind.Phys.var);
    free_d_1d_allocate(X->Bind.Phys.doublon);
    free_d_1d_allocate(X->Bind.Phys.doublon2);
    free_d_1d_allocate(X->Bind.Phys.Sz);
    free_d_1d_allocate(X->Bind.Phys.Sz2);
    free_d_1d_allocate(X->Bind.Phys.s2);
  }/*if (X->Bind.Def.iFlgCalcSpec == CALCSPEC_SCRATCH)*/

  //set omega
  if (SetOmega(&(X->Bind.Def)) != TRUE) {
    fprintf(stderr, "Error: Fail to set Omega.\n");
    exitMPI(-1);
  }
  else {
    if (X->Bind.Def.iFlgSpecOmegaOrg == FALSE) {
      X->Bind.Def.dcOmegaOrg = I * (X->Bind.Def.dcOmegaMax - X->Bind.Def.dcOmegaMin) / (double)X->Bind.Def.iNOmega;
    }
  }
  //
  // Read energy origin of each state
  //
  if (X->Bind.Def.iCalcType == FullDiag) {
    strcpy(sdt, cFileNameEigenvalue_Lanczos);
    childfopenMPI(sdt, "r", &fp);
    for (istate = 0; istate < nstate; istate++) {
      fgetsMPI(sdt, D_FileNameMax, fp);
      sscanf(sdt, "%ld%lf", &i, &energy0[istate]);
    }
    fclose(fp);
  }
  else if (X->Bind.Def.iCalcType == CG) {
    sprintf(sdt, cFileNameEnergy_CG, X->Bind.Def.CDataFileHead);
    childfopenMPI(sdt, "r", &fp);
    for (istate = 0; istate < nstate; istate++) {
      fgetsMPI(sdt, D_FileNameMax, fp);
      fgetsMPI(sdt, D_FileNameMax, fp);
      sscanf(sdt, "  Energy  %lf", &energy0[istate]);
      fgetsMPI(sdt, D_FileNameMax, fp);
      fgetsMPI(sdt, D_FileNameMax, fp);
      fgetsMPI(sdt, D_FileNameMax, fp);
    }
    fclose(fp);
  }
  /*
   Set & malloc omega grid
  */
  Nomega = X->Bind.Def.iNOmega;
  dcomega = cd_2d_allocate(nstate, Nomega);
  OmegaMax = X->Bind.Def.dcOmegaMax + X->Bind.Def.dcOmegaOrg;
  OmegaMin = X->Bind.Def.dcOmegaMin + X->Bind.Def.dcOmegaOrg;
  for (istate = 0; istate < nstate; istate++) {
    for (i = 0; i < Nomega; i++) {
      dcomega[istate][i] = (OmegaMax - OmegaMin) / Nomega * i + OmegaMin + energy0[istate];
    }
  }

  fprintf(stdoutMPI, "\nFrequency range:\n");
  fprintf(stdoutMPI, "  Omega Max. : %15.5e %15.5e\n", creal(OmegaMax), cimag(OmegaMax));
  fprintf(stdoutMPI, "  Omega Min. : %15.5e %15.5e\n", creal(OmegaMin), cimag(OmegaMin));
  fprintf(stdoutMPI, "  Num. of Omega : %d\n", Nomega);

  if (X->Bind.Def.NNSingleExcitationOperator == 0){
    if (X->Bind.Def.NNPairExcitationOperator == 0) {
      fprintf(stderr, "Error: Any excitation operators are not defined.\n");
      exitMPI(-1);
    }
    else {
      NdcSpectrum = X->Bind.Def.NNPairExcitationOperator - 1;
    }
  }
  else if (X->Bind.Def.NNPairExcitationOperator == 0) {
    NdcSpectrum = X->Bind.Def.NNSingleExcitationOperator - 1;
  }
  else {
    fprintf(stderr, "Error: Both single and pair excitation operators exist.\n");
    exitMPI(-1);
  }
  dcSpectrum = cd_3d_allocate(nstate, Nomega, NdcSpectrum);

  //Make New Lists
  if (MakeExcitedList(&(X->Bind), &iFlagListModified) == FALSE) {
    return FALSE;
  }
  X->Bind.Def.iFlagListModified = iFlagListModified;
  
  //Make excited state
  StartTimer(6100);
  if (X->Bind.Def.iFlgCalcSpec == RECALC_NOT ||
    X->Bind.Def.iFlgCalcSpec == RECALC_OUTPUT_TMComponents_VEC ||
    (X->Bind.Def.iFlgCalcSpec == RECALC_INOUT_TMComponents_VEC && X->Bind.Def.iCalcType == CG)) {
    v1Org = cd_2d_allocate(X->Bind.Check.idim_maxOrg + 1, nstate);
    vin = cd_1d_allocate(X->Bind.Check.idim_maxOrg + 1);
    //input eigen vector
    StartTimer(6101);
    fprintf(stdoutMPI, "  Start: An Eigenvector is inputted in CalcSpectrum.\n");
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_InputEigenVectorStart, "a");
    for (istate = 0; istate < nstate; istate++) {
      sprintf(sdt, cFileNameInputEigen, nstate, istate, myrank);
      childfopenALL(sdt, "rb", &fp);
      if (fp == NULL) {
        fprintf(stderr, "Error: Inputvector file is not found.\n");
        exitMPI(-1);
      }
      byte_size = fread(&i_stp, sizeof(int), 1, fp);
      X->Bind.Large.itr = i_stp; //For TPQ
      byte_size = fread(&i_max, sizeof(long int), 1, fp);
      if (i_max != X->Bind.Check.idim_maxOrg) {
        fprintf(stderr, "Error: Invalid Inputvector file.\n");
        exitMPI(-1);
      }
      byte_size = fread(vin, sizeof(complex double), i_max, fp);
#pragma omp parallel for default(none) shared(v1Org,vin, i_max, istate), private(i)
      for (i = 1; i <= i_max; i++) {
        v1Org[i][istate] = vin[i];
      }
      fclose(fp);
    }/*for (ie = 0; ie < X->Def.k_exct; ie++)*/
    StopTimer(6101);
    if (byte_size == 0) printf("byte_size: %d \n", (int)byte_size);
 
    free_cd_1d_allocate(vin);
  }
  StopTimer(6100);

  diagonalcalc(&(X->Bind));

  int iret = TRUE;
  fprintf(stdoutMPI, "  Start: Calculating a spectrum.\n\n");
  TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_CalcSpectrumStart, "a");
  StartTimer(6200);
  switch (X->Bind.Def.iCalcType) {
  case CG:
    iret = CalcSpectrumByBiCG(X, v0, v1, Nomega, NdcSpectrum, dcSpectrum, dcomega, v1Org);
    if (iret != TRUE) {
      return FALSE;
    }
    break;
  case FullDiag:
    iret = CalcSpectrumByFullDiag(X, Nomega, NdcSpectrum, dcSpectrum, dcomega, v1Org);
    break;
  default:
    break;
  }
  StopTimer(6200);

  if (iret != TRUE) {
    fprintf(stderr, "  Error: The selected calculation type is not supported for calculating spectrum mode.\n");
    return FALSE;
  }

  fprintf(stdoutMPI, "  End:  Calculating a spectrum.\n\n");
  TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_CalcSpectrumEnd, "a");
  iret = OutputSpectrum(X, nstate, Nomega, NdcSpectrum, dcSpectrum, dcomega, energy0);
  free_cd_3d_allocate(dcSpectrum);
  free_cd_2d_allocate(dcomega);
  free_d_1d_allocate(energy0);
  return TRUE;

}/*int CalcSpectrum*/
