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
#include "CalcSpectrumByLanczos.h"
#include "CalcSpectrumByBiCG.h"
#include "CalcSpectrumByTPQ.h"
#include "CalcSpectrumByFullDiag.h"
#include "CalcTime.h"
#include "SingleEx.h"
#include "PairEx.h"
#include "ExcitationSectorShift.h"
#include "wrapperMPI.h"
#include "FileIO.h"
#include "./common/setmemory.h"
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

/// \brief Output spectrum.
/// \param X [in] Read information of the frequency origin.
/// \param Nomega [in] A total number of discrete frequencies.
/// \param dcSpectrum [in] Array of spectrum.
/// \param dcomega [in] Array of discrete frequencies.
/// \retval FALSE Fail to open the output file.
/// \retval TRUE Success to output the spectrum.
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
  if(childfopenMPI(sdt, "w", &fp)!=0){
    return FALSE;
  }

  for (i = 0; i < Nomega; i++) {
    fprintf(fp, "%.10lf %.10lf %.10lf %.10lf \n",
      creal(dcomega[i]-X->Bind.Def.dcOmegaOrg), cimag(dcomega[i]-X->Bind.Def.dcOmegaOrg),
      creal(dcSpectrum[i]), cimag(dcSpectrum[i]));
  }/*for (i = 0; i < Nomega; i++)*/

  fclose(fp);
  return TRUE;
}/*int OutputSpectrum*/

/// \brief Compute ONE dynamical spectrum from a single input eigenvector
///        |v1Org> at spectral shift OmegaOrg (= E_n).
///
/// Reentrant building block for the finite-temperature loop. It rebuilds the
/// frequency grid from OmegaOrg every call (so both the resolvent shift
/// (z-(H-E_n))^{-1} and the output frequency axis stay consistent with the
/// per-launch path), builds the ket (and optional bra) excited state from
/// v1Org, dispatches the solver, and returns the spectrum in dcSpectrum.
/// It does NOT write any output file; the caller is responsible for I/O and
/// for the excitation-list / sector setup (MakeExcitedList) and its teardown.
///
/// \param X          [in,out] Calculation struct. dcOmegaOrg is set to OmegaOrg.
/// \param v1Org      [in]  Input eigenvector in the original (pre-excitation) sector.
/// \param OmegaOrg   [in]  Spectral shift E_n for this eigenstate.
/// \param Nomega     [in]  Number of frequencies.
/// \param dcSpectrum [out] Spectrum array (Nomega).
/// \param dcomega    [out] Frequency grid (Nomega), rebuilt here from OmegaOrg.
/// \retval TRUE  success (including the zero-norm case, which returns a zero spectrum).
/// \retval FALSE solver failure.
static int CalcOneSpectrum(
  struct EDMainCalStruct *X,
  double complex *v1Org,
  double complex OmegaOrg,
  int Nomega,
  double complex *dcSpectrum,
  double complex *dcomega)
{
  unsigned long int i;
  double dnorm = 0.0;
  double complex *v0_Bra = NULL;
  double complex OmegaMax, OmegaMin;
  char sdt[D_FileNameMax];
  FILE *fp;
  int iret = TRUE;

  /* Rebuild the frequency grid from THIS eigenstate's shift so the resolvent
     shift and the written frequency axis match the per-launch behavior. */
  X->Bind.Def.dcOmegaOrg = OmegaOrg;
  OmegaMax = X->Bind.Def.dcOmegaMax + X->Bind.Def.dcOmegaOrg;
  OmegaMin = X->Bind.Def.dcOmegaMin + X->Bind.Def.dcOmegaOrg;
  for (i = 0; i < Nomega; i++) {
    dcomega[i] = (OmegaMax - OmegaMin) / Nomega * i + OmegaMin;
  }
  fprintf(stdoutMPI, "\nFrequency range:\n");
  fprintf(stdoutMPI, "  Omega Max. : %15.5e %15.5e\n", creal(OmegaMax), cimag(OmegaMax));
  fprintf(stdoutMPI, "  Omega Min. : %15.5e %15.5e\n", creal(OmegaMin), cimag(OmegaMin));
  fprintf(stdoutMPI, "  Num. of Omega : %d\n", Nomega);

  for (i = 0; i <= X->Bind.Check.idim_max; i++) v0[i] = 0;

  TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_CalcExcitedStateStart, "a");
  fprintf(stdoutMPI, "  Start: Calculating an excited vector.\n");
  StartTimer(6102);
  {
    ExcitationOperatorSet ketSet = {
      X->Bind.Def.NSingleExcitationOperator, X->Bind.Def.SingleExcitationOperator,
      X->Bind.Def.ParaSingleExcitationOperator, X->Bind.Def.NPairExcitationOperator,
      X->Bind.Def.PairExcitationOperator, X->Bind.Def.ParaPairExcitationOperator};
    if (GetExcitedState(&(X->Bind), &ketSet, v0, v1Org) != TRUE) {
      fprintf(stderr, "Error: failed to build the ket excited state A|phi>.\n");
      exitMPI(-1);
    }
  }
  StopTimer(6102);

  dnorm = NormMPI_dc(X->Bind.Check.idim_max, v0);
  if (fabs(dnorm) < pow(10.0, -15)) {
    fprintf(stderr, "Warning: Norm of an excited vector becomes 0.\n");
    for (i = 0; i < Nomega; i++) dcSpectrum[i] = 0;
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_CalcExcitedStateEnd, "a");
    return TRUE;
  }
#pragma omp parallel for default(none) private(i) shared(v1, v0) firstprivate(dnorm, X)
  for (i = 1; i <= X->Bind.Check.idim_max; i++) {
    v1[i] = v0[i] / dnorm;
  }

  /* Build the bra excited state B|phi> (un-normalized) for off-diagonal spectrum. */
  if (X->Bind.Def.NSingleExcitationOperatorBra > 0 || X->Bind.Def.NPairExcitationOperatorBra > 0) {
    ExcitationOperatorSet braSet = {
      X->Bind.Def.NSingleExcitationOperatorBra, X->Bind.Def.SingleExcitationOperatorBra,
      X->Bind.Def.ParaSingleExcitationOperatorBra, X->Bind.Def.NPairExcitationOperatorBra,
      X->Bind.Def.PairExcitationOperatorBra, X->Bind.Def.ParaPairExcitationOperatorBra};
    v0_Bra = cd_1d_allocate(X->Bind.Check.idim_max + 1);
    for (i = 0; i <= X->Bind.Check.idim_max; i++) v0_Bra[i] = 0;
    if (GetExcitedState(&(X->Bind), &braSet, v0_Bra, v1Org) != TRUE) {
      fprintf(stderr, "Error: failed to build the bra excited state B|phi>.\n");
      exitMPI(-1);
    }
    if (NormMPI_dc(X->Bind.Check.idim_max, v0_Bra) < pow(10.0, -15)) {
      fprintf(stderr, "Warning: Norm of the bra excited vector B|phi> is 0; the off-diagonal spectrum will be zero.\n");
    }
  }

  if (X->Bind.Def.iOutputExVec == 1) {
    sprintf(sdt, cFileNameOutputExcitedVec, X->Bind.Def.CDataFileHead, myrank);
    if (childfopenALL(sdt, "w", &fp) != 0) {
      if (v0_Bra != NULL) free_cd_1d_allocate(v0_Bra);
      return -1;
    }
    fprintf(fp, "%ld\n", X->Bind.Check.idim_max);
    for (i = 1; i <= X->Bind.Check.idim_max; i++) {
      fprintf(fp, "%.10lf, %.10lf\n", creal(v0[i]), cimag(v0[i]));
    }
    fclose(fp);
  }
  fprintf(stdoutMPI, "  End:   Calculating an excited vector.\n\n");
  TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_CalcExcitedStateEnd, "a");

  diagonalcalc(&(X->Bind));

  fprintf(stdoutMPI, "  Start: Calculating a spectrum.\n\n");
  TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_CalcSpectrumStart, "a");
  StartTimer(6200);
  switch (X->Bind.Def.iCalcType) {
    case Lanczos:
      iret = CalcSpectrumByLanczos(X, v1, dnorm, Nomega, dcSpectrum, dcomega);
      break;
    case CG:
      iret = CalcSpectrumByBiCG(X, v0, (v0_Bra != NULL) ? v0_Bra : v0, v1, vg, Nomega, dcSpectrum, dcomega);
      break;
    case TPQCalc:
      fprintf(stderr, "  Error: TPQ is not supported for calculating spectrum mode.\n");
      iret = FALSE;
      break;
    case FullDiag:
      iret = CalcSpectrumByFullDiag(X, Nomega, dcSpectrum, dcomega);
      break;
    default:
      break;
  }
  StopTimer(6200);
  if (v0_Bra != NULL) free_cd_1d_allocate(v0_Bra);
  /* c_CalcSpectrumEnd is logged once by the caller after the spectrum is written. */
  return iret;
}/*int CalcOneSpectrum*/

/// \brief Read one MPI-distributed eigenvector into v1Org (original, pre-excitation sector).
///
/// Reads the SpectrumVec file for the local MPI rank. When useIdxSuffix is TRUE the file name
/// is "<SpectrumVec>_<idx>_rank_<myrank>.dat" (used by the finite-T eigenstate loop, where
/// SpectrumVec is the common base, e.g. ".../zvo_eigenvec"); when FALSE it is the legacy
/// "<SpectrumVec>_rank_<myrank>.dat". GetFileNameByKW returns a pointer into a shared global
/// buffer, so we strcpy it to a local before appending (the loop calls this repeatedly).
///
/// \param X            [in,out] Calculation struct (Large.itr is set from the stored step).
/// \param idx          [in] Eigenstate index (used only when useIdxSuffix is TRUE).
/// \param useIdxSuffix [in] Whether to insert the "_<idx>" eigenstate suffix.
/// \param v1Org        [out] Eigenvector (length idim_maxOrg+1).
/// \retval TRUE  success.
/// \retval FALSE file missing or dimension mismatch.
static int ReadEigenVector(
  struct EDMainCalStruct *X,
  int idx,
  int useIdxSuffix,
  double complex *v1Org)
{
  char sdt[D_FileNameMax], base[D_FileNameMax];
  char *kw;
  FILE *fp;
  unsigned long int i_max = 0;
  int i_stp;
  size_t byte_size;

  int nfn;
  GetFileNameByKW(KWSpectrumVec, &kw);
  /* local copy (bounded): do NOT mutate the shared global file-name buffer.
     Fail explicitly on truncation rather than silently opening a different file. */
  if (snprintf(base, sizeof(base), "%s", kw) >= (int)sizeof(base)) {
    fprintf(stderr, "Error: SpectrumVec base name is too long for the file-name buffer.\n");
    return FALSE;
  }
  nfn = useIdxSuffix ? snprintf(sdt, sizeof(sdt), "%s_%d_rank_%d.dat", base, idx, myrank)
                     : snprintf(sdt, sizeof(sdt), "%s_rank_%d.dat", base, myrank);
  if (nfn >= (int)sizeof(sdt)) {
    fprintf(stderr, "Error: eigenvector file name is too long (base=%s).\n", base);
    return FALSE;
  }

  childfopenALL(sdt, "rb", &fp);
  if (fp == NULL) {
    fprintf(stderr, "Error: A file of Input vector (%s) does not exist.\n", sdt);
    return FALSE;
  }
  /* Validate every fread count. v1Org is reused across eigenstates in loop mode, so a
     truncated file must hard-fail here; otherwise the buffer would keep the previous
     state's tail and silently produce a mixed-state spectrum. */
  if (fread(&i_stp, sizeof(i_stp), 1, fp) != 1) {
    fprintf(stderr, "Error: failed to read the step header from Input vector (%s).\n", sdt);
    fclose(fp);
    return FALSE;
  }
  X->Bind.Large.itr = i_stp; /* For TPQ */
  if (fread(&i_max, sizeof(i_max), 1, fp) != 1) {
    fprintf(stderr, "Error: failed to read the dimension header from Input vector (%s).\n", sdt);
    fclose(fp);
    return FALSE;
  }
  if (i_max != X->Bind.Check.idim_maxOrg) {
    fprintf(stderr, "Error: myrank=%d, i_max=%ld\n", myrank, i_max);
    fprintf(stderr, "Error: A file of Input vector (%s) is incorrect.\n", sdt);
    fclose(fp);
    return FALSE;
  }
  byte_size = fread(v1Org, sizeof(complex double), i_max + 1, fp);
  fclose(fp);
  if (byte_size != (size_t)(i_max + 1)) {
    fprintf(stderr, "Error: truncated Input vector (%s): read %zu of %lu elements.\n",
            sdt, byte_size, (unsigned long)(i_max + 1));
    return FALSE;
  }
  return TRUE;
}/*int ReadEigenVector*/

/// \brief Read the idx-th eigen-energy E_idx from "<CDataFileHead>_energy.dat".
///
/// The energy file holds one "Energy <value>" line per computed eigenstate, in eigenstate
/// order; this returns the (idx)-th of them. Used as the per-state spectral shift OmegaOrg
/// in the finite-T loop (the same value DCore would otherwise pass as OmegaOrg per launch).
///
/// \param X   [in] Calculation struct (for CDataFileHead).
/// \param idx [in] Eigenstate index.
/// \param Ene [out] E_idx.
/// \retval TRUE  success.
/// \retval FALSE file missing or index not found.
static int ReadEigenEnergy(
  struct EDMainCalStruct *X,
  int idx,
  double *Ene)
{
  FILE *fp;
  char sdt[D_FileNameMax], ctmp[256], key[256];
  int count = 0;
  double val;

  if (snprintf(sdt, sizeof(sdt), cFileNameEnergy_Lanczos, X->Bind.Def.CDataFileHead) >= (int)sizeof(sdt)) {
    fprintf(stdoutMPI, "Error: energy file name is too long for the file-name buffer.\n");
    return FALSE;
  }
  childfopenMPI(sdt, "r", &fp);
  if (fp == NULL) {
    fprintf(stdoutMPI, "Error: %s does not exist.\n", sdt);
    return FALSE;
  }
  while (fgetsMPI(ctmp, 256, fp) != NULL) {
    if (sscanf(ctmp, "%255s %lf", key, &val) == 2 && strcmp(key, "Energy") == 0) {
      if (count == idx) {
        *Ene = val;
        fclose(fp);
        return TRUE;
      }
      count++;
    }
  }
  fclose(fp);
  fprintf(stdoutMPI, "Error: eigen-energy index %d not found in %s (found %d energies).\n",
          idx, sdt, count);
  return FALSE;
}/*int ReadEigenEnergy*/

/// \brief Write the spectrum for eigenstate idx to "<CDataFileHead>_DynamicalGreen_<idx>.dat".
///
/// Same format as OutputSpectrum, but with the per-eigenstate "_<idx>" suffix consumed by the
/// DCore finite-T reader (CalcSpectrumCore._read_spectrum).
static int OutputSpectrumIdx(
  struct EDMainCalStruct *X,
  int idx,
  int Nomega,
  double complex *dcSpectrum,
  double complex *dcomega)
{
  FILE *fp;
  char sdt[D_FileNameMax];
  int i;

  if (snprintf(sdt, sizeof(sdt), "%s_DynamicalGreen_%d.dat", X->Bind.Def.CDataFileHead, idx) >= (int)sizeof(sdt)) {
    fprintf(stderr, "Error: DynamicalGreen file name is too long for the file-name buffer.\n");
    return FALSE;
  }
  if (childfopenMPI(sdt, "w", &fp) != 0) {
    return FALSE;
  }
  for (i = 0; i < Nomega; i++) {
    fprintf(fp, "%.10lf %.10lf %.10lf %.10lf \n",
      creal(dcomega[i]-X->Bind.Def.dcOmegaOrg), cimag(dcomega[i]-X->Bind.Def.dcOmegaOrg),
      creal(dcSpectrum[i]), cimag(dcSpectrum[i]));
  }
  fclose(fp);
  return TRUE;
}/*int OutputSpectrumIdx*/

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
    unsigned long int i;
    int iFlagListModified = FALSE;
    double dnorm = 0.0;
    int iret = TRUE;
    int bAlreadyOutput = FALSE; /* TRUE once the eigenstate loop has written its own per-idx output */

    //ToDo: Nomega should be given as a parameter
    int Nomega;
    double complex OmegaMax, OmegaMin;
    double complex *dcSpectrum;
    double complex *dcomega;

    //set omega
    if (SetOmega(&(X->Bind.Def)) != TRUE) {
        fprintf(stderr, "Error: Fail to set Omega.\n");
        exitMPI(-1);
    } else {
        if (X->Bind.Def.iFlgSpecOmegaOrg == FALSE) {
            X->Bind.Def.dcOmegaOrg = I*(X->Bind.Def.dcOmegaMax - X->Bind.Def.dcOmegaMin) / (double) X->Bind.Def.iNOmega;
        }
    }
    /*
     Set & malloc omega grid
    */
    Nomega = X->Bind.Def.iNOmega;
    dcSpectrum = cd_1d_allocate(Nomega);
    dcomega = cd_1d_allocate(Nomega);
    OmegaMax = X->Bind.Def.dcOmegaMax + X->Bind.Def.dcOmegaOrg;
    OmegaMin = X->Bind.Def.dcOmegaMin + X->Bind.Def.dcOmegaOrg;
    for (i = 0; i < Nomega; i++) {
        dcomega[i] = (OmegaMax - OmegaMin) / Nomega * i + OmegaMin;
    }
    /* The frequency-range banner is printed inside CalcOneSpectrum, where the
       per-eigenstate OmegaOrg (and hence the actual grid) is known. */

  if (X->Bind.Def.NSingleExcitationOperator == 0 && X->Bind.Def.NPairExcitationOperator == 0) {
    fprintf(stderr, "Error: Any excitation operators are not defined.\n");
    exitMPI(-1);
  }
  /* The finite-T loop reuses one excited-vector buffer per eigenstate, so the
     single OutputExVec filename would be overwritten each iteration (only the last
     state would survive). Reject the combination until a per-state convention exists. */
  if (X->Bind.Def.iSpectrumLoopExct > 0 && X->Bind.Def.iOutputExVec == 1) {
    fprintf(stderr, "Error: OutputExVec=1 is not supported together with SpectrumLoopExct>0 "
                    "(the per-state excited vectors would overwrite each other).\n");
    exitMPI(-1);
  }
  /* The finite-T loop solves each eigenstate fresh (CalcSpec=Normal). The restart/save
     CalcSpec modes write single-state artifacts (TMComponents / recalcvec) to one fixed
     name, which the per-idx loop would overwrite, leaving only the last state's data.
     Require Normal until per-eigenstate restart filenames exist. */
  if (X->Bind.Def.iSpectrumLoopExct > 0 && X->Bind.Def.iFlgCalcSpec != RECALC_NOT) {
    fprintf(stderr, "Error: SpectrumLoopExct>0 requires CalcSpec=\"Normal\" "
                    "(restart/save modes would overwrite the per-state TMComponents/recalcvec files).\n");
    exitMPI(-1);
  }
  /* Off-diagonal (bra) excitation input validation. All checks are bra-gated:
     with no *Bra input this block is skipped and the diagonal path is unchanged. */
  if (X->Bind.Def.NSingleExcitationOperatorBra > 0 || X->Bind.Def.NPairExcitationOperatorBra > 0) {
    int ketSingle = (X->Bind.Def.NSingleExcitationOperator > 0);
    int ketPair   = (X->Bind.Def.NPairExcitationOperator > 0);
    int braSingle = (X->Bind.Def.NSingleExcitationOperatorBra > 0);
    int braPair   = (X->Bind.Def.NPairExcitationOperatorBra > 0);
    int isGeneralSpin = X->Bind.Def.iFlgGeneralSpin;
    int iCalcModel = X->Bind.Def.iCalcModel;
    SectorShift ketShift, braShift;
    if (X->Bind.Def.iCalcType != CG) {
      fprintf(stderr, "Error: off-diagonal (bra) spectrum requires method=\"CG\".\n");
      exitMPI(-1);
    }
    if (X->Bind.Def.iFlgCalcSpec != RECALC_NOT) {
      fprintf(stderr, "Error: off-diagonal (bra) spectrum currently requires CalcSpec=\"Normal\" (no restart/save); saved BiCG components carry no bra-operator metadata.\n");
      exitMPI(-1);
    }
    if (braSingle && braPair) {
      fprintf(stderr, "Error: SingleExcitationBra and PairExcitationBra cannot be used together.\n");
      exitMPI(-1);
    }
    if (ketSingle != braSingle || ketPair != braPair) {
      fprintf(stderr, "Error: ket and bra excitation operators must be the same type (both single or both pair).\n");
      exitMPI(-1);
    }
    if (ketSingle) {
      ketShift = GetExcitationOperatorSetShift(iCalcModel, isGeneralSpin, FALSE,
                   X->Bind.Def.SingleExcitationOperator, X->Bind.Def.NSingleExcitationOperator);
      braShift = GetExcitationOperatorSetShift(iCalcModel, isGeneralSpin, FALSE,
                   X->Bind.Def.SingleExcitationOperatorBra, X->Bind.Def.NSingleExcitationOperatorBra);
    } else {
      ketShift = GetExcitationOperatorSetShift(iCalcModel, isGeneralSpin, TRUE,
                   X->Bind.Def.PairExcitationOperator, X->Bind.Def.NPairExcitationOperator);
      braShift = GetExcitationOperatorSetShift(iCalcModel, isGeneralSpin, TRUE,
                   X->Bind.Def.PairExcitationOperatorBra, X->Bind.Def.NPairExcitationOperatorBra);
    }
    if (ketShift.valid == FALSE) {
      if (ketShift.reason == OFFDIAG_SHIFT_SET_INCONSISTENT)
        fprintf(stderr, "Error: the ket excitation operator set mixes operators with different Hilbert-sector shifts.\n");
      else
        fprintf(stderr, "Error: off-diagonal spectrum is not supported for this model / ket excitation operator.\n");
      exitMPI(-1);
    }
    if (braShift.valid == FALSE) {
      if (braShift.reason == OFFDIAG_SHIFT_SET_INCONSISTENT)
        fprintf(stderr, "Error: the bra excitation operator set mixes operators with different Hilbert-sector shifts.\n");
      else
        fprintf(stderr, "Error: off-diagonal spectrum is not supported for this model / bra excitation operator.\n");
      exitMPI(-1);
    }
    if (ketShift.dNe != braShift.dNe || ketShift.dNup != braShift.dNup ||
        ketShift.dNdown != braShift.dNdown || ketShift.dTotal2Sz != braShift.dTotal2Sz) {
      fprintf(stderr, "Error: ket and bra excitation operators map to different Hilbert sectors (sector mismatch).\n");
      exitMPI(-1);
    }
  }
  //Make New Lists
  if (MakeExcitedList(&(X->Bind), &iFlagListModified) == FALSE) {
    return FALSE;
  }
  X->Bind.Def.iFlagListModified = iFlagListModified;

    //Set Memory
    v1Org = cd_1d_allocate(X->Bind.Check.idim_maxOrg+1);
    for(i=0; i<X->Bind.Check.idim_maxOrg+1; i++){
      v1Org[i]=0;
    }
    
    //Make excited state
    StartTimer(6100);
    if (X->Bind.Def.iFlgCalcSpec == RECALC_NOT ||
        X->Bind.Def.iFlgCalcSpec == RECALC_OUTPUT_TMComponents_VEC ||
       (X->Bind.Def.iFlgCalcSpec == RECALC_INOUT_TMComponents_VEC && X->Bind.Def.iCalcType == CG)) {
        //input eigen vector
      StartTimer(6101);
      if (X->Bind.Def.iSpectrumLoopExct > 0) {
        /* Finite-temperature eigenstate loop. One HPhi launch covers all the
           thermally-relevant eigenstates: for each idx we read eigenvec_{idx} ONCE,
           build the excited state from it, solve with the per-state spectral shift
           OmegaOrg = E_idx, and write DynamicalGreen_{idx}.dat. This replaces DCore's
           per-(idx) relaunch, eliminating the per-state process-startup overhead. */
        int idx;
        int nloop = X->Bind.Def.iSpectrumLoopExct;
        /* Only the REAL part of the spectral shift is state-dependent (= E_idx).
           Preserve any configured imaginary broadening in dcOmegaOrg (OmegaIm /
           complex OmegaOrg) so the loop path matches the single-state path; with
           the DCore driver this imaginary part is 0, so the two are identical. */
        double dOmegaOrgIm = cimag(X->Bind.Def.dcOmegaOrg);
        double Elast = 0.0;
        fprintf(stdoutMPI, "  Start: finite-T spectrum loop over %d eigenstate(s).\n", nloop);
        TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_InputEigenVectorStart, "a");
        StopTimer(6100);
        /* Fail before writing any spectrum if SpectrumLoopExct exceeds the number of
           computed eigenstates (the energy file has one entry per eigenstate, matching
           the eigenvector files), so no partial DynamicalGreen_{idx}.dat is left behind. */
        if (ReadEigenEnergy(X, nloop - 1, &Elast) != TRUE) {
          fprintf(stderr, "Error: SpectrumLoopExct=%d exceeds the number of computed eigenstates.\n", nloop);
          iret = FALSE;
        }
        for (idx = 0; iret == TRUE && idx < nloop; idx++) {
          double Eidx = 0.0;
          if (ReadEigenVector(X, idx, TRUE, v1Org) != TRUE) { iret = FALSE; break; }
          if (ReadEigenEnergy(X, idx, &Eidx) != TRUE)        { iret = FALSE; break; }
          iret = CalcOneSpectrum(X, v1Org, Eidx + dOmegaOrgIm * I, Nomega, dcSpectrum, dcomega);
          if (iret != TRUE) break;
          if (OutputSpectrumIdx(X, idx, Nomega, dcSpectrum, dcomega) != TRUE) { iret = FALSE; break; }
        }
        StopTimer(6101);
        TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_InputEigenVectorEnd, "a");
        bAlreadyOutput = TRUE; /* the loop wrote one file per eigenstate; skip the single output below */
      } else {
        fprintf(stdoutMPI, "  Start: An Eigenvector is inputted in CalcSpectrum.\n");
        TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_InputEigenVectorStart, "a");
        if (ReadEigenVector(X, 0, FALSE, v1Org) != TRUE) {
          return -1;
        }
        StopTimer(6101);
        fprintf(stdoutMPI, "  End:   An Input vector is inputted in CalcSpectrum.\n\n");
        TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_InputEigenVectorEnd, "a");
        StopTimer(6100);
        /* Build the excited state from v1Org and solve, in one reentrant call. */
        iret = CalcOneSpectrum(X, v1Org, X->Bind.Def.dcOmegaOrg, Nomega, dcSpectrum, dcomega);
      }
  }
  else {
    /* Restart path (recalc from saved BiCG data): the excited state is NOT rebuilt
       and off-diagonal (bra) is unsupported here; solve directly. */
    StopTimer(6100);
    diagonalcalc(&(X->Bind));
    fprintf(stdoutMPI, "  Start: Calculating a spectrum.\n\n");
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_CalcSpectrumStart, "a");
    StartTimer(6200);
    switch (X->Bind.Def.iCalcType) {
      case Lanczos:
        iret = CalcSpectrumByLanczos(X, v1, dnorm, Nomega, dcSpectrum, dcomega);
        break;
      case CG:
        iret = CalcSpectrumByBiCG(X, v0, v0, v1, vg, Nomega, dcSpectrum, dcomega);
        break;
      case TPQCalc:
        fprintf(stderr, "  Error: TPQ is not supported for calculating spectrum mode.\n");
        iret = FALSE;
        break;
      case FullDiag:
        iret = CalcSpectrumByFullDiag(X, Nomega, dcSpectrum, dcomega);
        break;
      default:
        break;
    }
    StopTimer(6200);
  }

  /* Free the original-sector lists and the input eigenvector. */
  if (iFlagListModified == TRUE) {
    free(v1Org);
    free(list_1_org);
    free(list_2_1_org);
    free(list_2_2_org);
  }

  if (iret != TRUE) {
    /* The specific cause (unsupported solver type, failed read, etc.) was already
       reported where it occurred; keep this message neutral so it does not mask it. */
    fprintf(stderr, "  Error: spectrum calculation failed.\n");
    free_cd_1d_allocate(dcSpectrum);
    free_cd_1d_allocate(dcomega);
    return FALSE;
  }

  fprintf(stdoutMPI, "  End:  Calculating a spectrum.\n\n");
  TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_CalcSpectrumEnd, "a");
  /* The eigenstate loop already wrote one DynamicalGreen_{idx}.dat per state;
     only the single-eigenvector path needs the single OutputSpectrum here. */
  if (bAlreadyOutput == FALSE) {
    iret = OutputSpectrum(X, Nomega, dcSpectrum, dcomega);
  }
  free_cd_1d_allocate(dcSpectrum);
  free_cd_1d_allocate(dcomega);
  /* Mirror the loop path: report a failed write to the caller instead of success. */
  return (iret == TRUE) ? TRUE : FALSE;

}/*int CalcSpectrum*/

///
/// \brief Parent function to calculate the excited state.
/// \param X [in] Struct to get number of excitation operators.
/// \param tmp_v0 [out] Result @f$ v_0 = H_{ex} v_1 @f$.
/// \param tmp_v1 [in] The original state before excitation  @f$ v_1 @f$.
/// \retval FALSE Fail to calculate the excited state.
/// \retval TRUE Success to calculate the excited state.
int GetExcitedState
(
 struct BindStruct *X,
 const ExcitationOperatorSet *op,
 double complex *tmp_v0,
 double complex *tmp_v1
) {
  if (op->NSingle > 0 && op->NPair > 0) {
    fprintf(stderr, "Error: Both single and pair excitation operators exist.\n");
    return FALSE;
  }


  if (op->NSingle > 0) {
    if (GetSingleExcitedState(X, op->NSingle, op->Single, op->ParaSingle, tmp_v0, tmp_v1) != TRUE) {
      return FALSE;
    }
  } else if (op->NPair > 0) {
    if (GetPairExcitedState(X, op->NPair, op->Pair, op->ParaPair, tmp_v0, tmp_v1) != TRUE) {
      return FALSE;
    }
  } else {
    unsigned long int i_max, j;
    i_max = X->Check.idim_maxOrg;
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1)  \
  firstprivate(i_max) private(j)
    for (j = 1; j <= i_max; j++) tmp_v0[j] = tmp_v1[j];
  }

  return TRUE;
}

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
){
  FILE *fp;
  char sdt[D_FileNameMax],ctmp[256];
  int istp=4;
  double E1, E2, E3, E4, Emax;
    long unsigned int iline_countMax=2;
    long unsigned int iline_count=2;


  if(X->iFlgSpecOmegaMax == TRUE && X->iFlgSpecOmegaMin == TRUE){
    return TRUE;
  }
  else{
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
    if(X->iFlgSpecOmegaMax == FALSE){
      X->dcOmegaMax= Emax*(double)X->Nsite;
    }
    if(X->iFlgSpecOmegaMin == FALSE){
      X->dcOmegaMin= E1;
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

    if (X->Def.NSingleExcitationOperator > 0) {
        switch (X->Def.iCalcModel) {
            case HubbardGC:
                break;
            case HubbardNConserved:
            case Hubbard:
            case Kondo:
            case KondoGC:
            case tJ:
            case tJNConserved:
            case tJGC:
                *iFlgListModifed = TRUE;
                break;
            case Spin:
            case SpinGC:
                return FALSE;
        }
    } else if (X->Def.NPairExcitationOperator > 0) {
        switch (X->Def.iCalcModel) {
            case SpinGC:
            case HubbardNConserved:
            case HubbardGC:
            case KondoGC:
            case tJNConserved:
            case tJGC:
                break;
            case Hubbard:
            case Kondo:
            case tJ:
            case Spin:
                if (X->Def.PairExcitationOperator[0][1] != X->Def.PairExcitationOperator[0][3]) {
                    *iFlgListModifed = TRUE;
                }
                break;
        }
    } else {
        return FALSE;
    }

    if (*iFlgListModifed == TRUE) {
        if(GetlistSize(X)==TRUE) {
            list_1_org = lui_1d_allocate(X->Check.idim_max + 1);
#ifdef MPI
            list_1buf_org = lui_1d_allocate(X->Check.idim_maxMPI + 1);
            //lui_malloc1(list_1buf_org, X->Check.idim_maxMPI + 1);
#endif // MPI
            list_2_1_org = lui_1d_allocate(X->Large.SizeOflist_2_1);
            list_2_2_org = lui_1d_allocate(X->Large.SizeOflist_2_2);
            //lui_malloc1(list_2_1_org, X->Large.SizeOflist_2_1);
            //lui_malloc1(list_2_2_org, X->Large.SizeOflist_2_2);
            if(list_1_org==NULL
               || list_2_1_org==NULL
               || list_2_2_org==NULL
                    )
            {
                return -1;
            }
            for(j =0; j<X->Large.SizeOflist_2_1; j++){
                list_2_1_org[j]=0;
            }
            for(j =0; j<X->Large.SizeOflist_2_2; j++){
                list_2_2_org[j]=0;
            }

        }

        if (sz(X, list_1_org, list_2_1_org, list_2_2_org) != 0) {
            return FALSE;
        }

        if (X->Def.NSingleExcitationOperator > 0) {
            switch (X->Def.iCalcModel) {
                case HubbardGC:
                    break;
                case HubbardNConserved:
                case KondoNConserved:/*To be confirmed*/
                case tJNConserved:/*To be confirmed*/
                    if (X->Def.SingleExcitationOperator[0][2] == 1) { //cis
                        X->Def.Ne = X->Def.NeMPI + 1;
                    }
                    else{
                        X->Def.Ne = X->Def.NeMPI - 1;
                    }
                    break;
                case Hubbard:
                case Kondo:
                case KondoGC:
                case tJ:
                case tJGC:
                    if (X->Def.SingleExcitationOperator[0][2] == 1) { //cis
                        X->Def.Ne = X->Def.NeMPI + 1;
                        if (X->Def.SingleExcitationOperator[0][1] == 0) {//up
                            X->Def.Nup = X->Def.NupOrg + 1;
                            X->Def.Ndown=X->Def.NdownOrg;
                        } else {//down
                            X->Def.Nup=X->Def.NupOrg;
                            X->Def.Ndown = X->Def.NdownOrg + 1;
                        }
                    } else {//ajt
                        X->Def.Ne = X->Def.NeMPI - 1;
                        if (X->Def.SingleExcitationOperator[0][1] == 0) {//up
                            X->Def.Nup = X->Def.NupOrg - 1;
                            X->Def.Ndown=X->Def.NdownOrg;

                        } else {//down
                            X->Def.Nup=X->Def.NupOrg;
                            X->Def.Ndown = X->Def.NdownOrg - 1;
                        }
                    }
                    break;
                case Spin:
                case SpinGC:
                    return FALSE;
            }
        } else if (X->Def.NPairExcitationOperator > 0) {
            X->Def.Ne=X->Def.NeMPI;
            switch (X->Def.iCalcModel) {
                case SpinGC:
                case HubbardNConserved:
                case HubbardGC:
                case KondoNConserved:/*To be confirmed*/
                case tJNConserved:/*To be confirmed*/
                case tJGC:
                    break;
                case Hubbard:
                case Kondo:
                case KondoGC:
                case tJ:
                    if (X->Def.PairExcitationOperator[0][1] != X->Def.PairExcitationOperator[0][3]) {
                      if (X->Def.PairExcitationOperator[0][1] == 0) {//up
                        X->Def.Nup = X->Def.NupOrg + 1;
                        X->Def.Ndown = X->Def.NdownOrg - 1;
                      } else {//down
                        X->Def.Nup = X->Def.NupOrg - 1;
                        X->Def.Ndown = X->Def.NdownOrg + 1;
                      }
                    }
                    break;
              case Spin:
                if (X->Def.PairExcitationOperator[0][1] != X->Def.PairExcitationOperator[0][3]) {
                  if (X->Def.iFlgGeneralSpin == FALSE) {
                    if (X->Def.PairExcitationOperator[0][1] == 0) {//down
                      X->Def.Nup = X->Def.NupOrg - 1;
                      X->Def.Ndown = X->Def.NdownOrg + 1;
                    } else {//up
                      X->Def.Nup = X->Def.NupOrg + 1;
                      X->Def.Ndown = X->Def.NdownOrg - 1;
                    }
                  }
                  else{//for general spin
                      X->Def.Total2Sz = X->Def.Total2SzMPI+2*(X->Def.PairExcitationOperator[0][1]-X->Def.PairExcitationOperator[0][3]);
                  }
                }
                break;
            }
        } else {
            return FALSE;
        }
        //Update Infomation
        X->Def.Nsite=X->Def.NsiteMPI;

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

    if(X->Def.iCalcModel==HubbardNConserved){
        X->Def.iCalcModel=Hubbard;
    }

#ifdef _DEBUG
  if (*iFlgListModifed == TRUE) {
    for(j=1; j<=X->Check.idim_maxOrg; j++){
        fprintf(stdout, "Debug1: myrank=%d, list_1_org[ %ld] = %ld\n", myrank, j, list_1_org[j]+myrank*X->Def.OrgTpow[2*X->Def.NsiteMPI-1]);
    }

    for(j=1; j<=X->Check.idim_max; j++){
        fprintf(stdout, "Debug2: myrank=%d, list_1[ %ld] = %ld\n", myrank, j, list_1[j]+myrank* 64);
    }
    }
#endif

    return TRUE;
}
