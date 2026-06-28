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
    case CG: {
      /* Single bra: the namelist B|phi> if present, else the ket itself (diagonal G_AA). */
      double complex *braList[1] = { (v0_Bra != NULL) ? v0_Bra : v0 };
      iret = CalcSpectrumByBiCG(X, v0, braList, 1, v1, vg, Nomega, dcSpectrum, dcomega);
      break;
    }
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

/// \brief Write the spectrum for eigenstate idx (operator op, bra b) to a per-state output file.
///
/// Naming, by what is active:
///  - single operator, single bra (useOp==FALSE):   "<head>_DynamicalGreen_<idx>.dat" (Stage 1a).
///  - multiple operators (useOp==TRUE, useBra==FALSE): "<head>_DynamicalGreen_<idx>_<op>.dat".
///  - multiple bras (useBra==TRUE, op field forced on): "<head>_DynamicalGreen_<idx>_<op>_<bra>.dat".
/// The spectrum for bra b lives at dcSpectrum[i*stride + bra] (stride == nBra; Komega x(nBra,Nomega)
/// column-major layout), so this reads a strided slice. Same column format as OutputSpectrum;
/// consumed by the DCore finite-T reader (CalcSpectrumCore._read_spectrum).
static int OutputSpectrumIdx(
  struct EDMainCalStruct *X,
  int idx,
  int op,
  int useOp,
  int bra,
  int useBra,
  int stride,
  int Nomega,
  double complex *dcSpectrum,
  double complex *dcomega)
{
  FILE *fp;
  char sdt[D_FileNameMax];
  int i, nfn;

  if (useBra)
    nfn = snprintf(sdt, sizeof(sdt), "%s_DynamicalGreen_%d_%d_%d.dat", X->Bind.Def.CDataFileHead, idx, op, bra);
  else if (useOp)
    nfn = snprintf(sdt, sizeof(sdt), "%s_DynamicalGreen_%d_%d.dat", X->Bind.Def.CDataFileHead, idx, op);
  else
    nfn = snprintf(sdt, sizeof(sdt), "%s_DynamicalGreen_%d.dat", X->Bind.Def.CDataFileHead, idx);
  if (nfn >= (int)sizeof(sdt)) {
    fprintf(stderr, "Error: DynamicalGreen file name is too long for the file-name buffer.\n");
    return FALSE;
  }
  if (childfopenMPI(sdt, "w", &fp) != 0) {
    return FALSE;
  }
  for (i = 0; i < Nomega; i++) {
    fprintf(fp, "%.10lf %.10lf %.10lf %.10lf \n",
      creal(dcomega[i]-X->Bind.Def.dcOmegaOrg), cimag(dcomega[i]-X->Bind.Def.dcOmegaOrg),
      creal(dcSpectrum[i*stride + bra]), cimag(dcSpectrum[i*stride + bra]));
  }
  fclose(fp);
  return TRUE;
}/*int OutputSpectrumIdx*/

/// \brief Read one single-excitation operator set from an HPhi single_ex-format def file.
///
/// File layout (as written by the DCore driver): 5 header lines, the 2nd being "NSingle <N>",
/// followed by N lines "<site> <spin> <type> <re> <im>". Operators in one file are summed into
/// one excited state (e.g. c_i + i c_j). Allocates *pOps ([N][3]) and *pPara ([N]); the caller
/// frees them with free_i_2d_allocate / free_cd_1d_allocate. Used to load the operator sets
/// op=1.. for the multi-operator finite-T loop (set 0 comes from the namelist SingleExcitation).
///
/// \retval TRUE  success.
/// \retval FALSE file missing or malformed.
static int ReadSingleExcitationSet(
  const char *fname,
  unsigned int maxN,
  unsigned int *pN,
  int ***pOps,
  double complex **pPara)
{
  FILE *fp;
  char ctmp[256], kw[256];
  unsigned int N = 0, k;
  int site, spin, type;
  double re, im;
  int **ops;
  double complex *para;

  fp = fopenMPI(fname, "r");
  if (fp == NULL) {
    fprintf(stderr, "Error: single-excitation file (%s) does not exist.\n", fname);
    return FALSE;
  }
  /* Validate the full 5-line header. fgetsMPI returns NULL at EOF while leaving the
     previous buffer contents in place, so an unchecked read on a truncated file would
     silently reuse the prior line; require every header line to be present. */
  if (fgetsMPI(ctmp, 256, fp) == NULL) {                       /* line 1: separator */
    fprintf(stderr, "Error: %s is truncated (missing header).\n", fname);
    fclose(fp); return FALSE;
  }
  if (fgetsMPI(ctmp, 256, fp) == NULL ||                       /* line 2: "NSingle <N>" */
      sscanf(ctmp, "%255s %u", kw, &N) != 2 || strcmp(kw, "NSingle") != 0) {
    fprintf(stderr, "Error: malformed \"NSingle <N>\" header in %s.\n", fname);
    fclose(fp); return FALSE;
  }
  if (fgetsMPI(ctmp, 256, fp) == NULL ||                       /* lines 3-5: separators */
      fgetsMPI(ctmp, 256, fp) == NULL ||
      fgetsMPI(ctmp, 256, fp) == NULL) {
    fprintf(stderr, "Error: %s is truncated (incomplete header).\n", fname);
    fclose(fp); return FALSE;
  }
  if (N == 0) { fclose(fp); *pN = 0; *pOps = NULL; *pPara = NULL; return TRUE; }
  /* Bound NSingle (read from an external file) before allocating: a single-excitation set
     cannot have more than 2*Nsite distinct (site, spin) operators, so a larger value is a
     corrupt file. The codebase allocators abort on a failed malloc rather than returning a
     checkable status, so an absurd NSingle must be rejected here, not after allocation. */
  if (maxN > 0 && N > maxN) {
    fprintf(stderr, "Error: NSingle=%u in %s exceeds the maximum of %u operators.\n", N, fname, maxN);
    fclose(fp); return FALSE;
  }

  ops = i_2d_allocate(N, 3);
  para = cd_1d_allocate(N);
  for (k = 0; k < N; k++) {
    if (fgetsMPI(ctmp, 256, fp) == NULL ||
        sscanf(ctmp, "%d %d %d %lf %lf", &site, &spin, &type, &re, &im) != 5) {
      fprintf(stderr, "Error: malformed operator line %u in %s.\n", k, fname);
      fclose(fp);
      free_i_2d_allocate(ops);
      free_cd_1d_allocate(para);
      return FALSE;
    }
    ops[k][0] = site; ops[k][1] = spin; ops[k][2] = type;
    para[k] = re + I * im;
  }
  fclose(fp);
  *pN = N; *pOps = ops; *pPara = para;
  return TRUE;
}/*int ReadSingleExcitationSet*/

/// \brief Compute nBra dynamical spectra that share ONE ket BiCG solve (Stage 3 bra/ket reuse).
///
/// Generalises CalcOneSpectrum: the ket A|phi> is built from X->Def's single-excitation operator
/// (the current op set) and the resolvent (z-(H-E))^{-1}|A phi> is solved ONCE; it is then
/// projected onto nBra bras B_b|phi> via Komega's nl projections, giving
/// G_b(z) = <B_b phi|(z-(H-E))^{-1}|A phi> for every b in one solve. This cuts the BiCG count for
/// off-diagonal Green's functions from n_orb^2 to n_orb. The bras are un-normalized (the ket
/// carries the norm on both sides, as in the single-bra path) and live in the same excited sector
/// as the ket (one-body c_i/c_j on one spin), so the shared MakeExcitedList lists apply to both.
/// Writes nBra spectra into dcSpectrum with Komega's column-major layout dcSpectrum[i*nBra + b].
/// Does NOT write output; the caller does I/O and owns the excitation lists and v1Org/dcSpectrum.
///
/// \retval TRUE  success (including the zero-ket-norm case, which zeroes every spectrum).
/// \retval FALSE allocation or solver failure.
static int CalcSpectrumMultiBra(
  struct EDMainCalStruct *X,
  double complex *v1Org,
  double complex OmegaOrg,
  int Nomega,
  double complex *dcSpectrum,     /* [nBra*Nomega], column-major x(nBra,Nomega) */
  double complex *dcomega,
  int nBra,
  unsigned int *braSet_N,
  int ***braSet_ops,
  double complex **braSet_para)
{
  unsigned long int i;
  int b, iret = TRUE;
  double dnorm;
  double complex **braVecs, **braList;
  double complex OmegaMax, OmegaMin;

  /* Rebuild the frequency grid from THIS eigenstate's shift (same as CalcOneSpectrum). */
  X->Bind.Def.dcOmegaOrg = OmegaOrg;
  OmegaMax = X->Bind.Def.dcOmegaMax + X->Bind.Def.dcOmegaOrg;
  OmegaMin = X->Bind.Def.dcOmegaMin + X->Bind.Def.dcOmegaOrg;
  for (i = 0; i < (unsigned long int)Nomega; i++)
    dcomega[i] = (OmegaMax - OmegaMin) / Nomega * i + OmegaMin;
  fprintf(stdoutMPI, "\nFrequency range:\n");
  fprintf(stdoutMPI, "  Omega Max. : %15.5e %15.5e\n", creal(OmegaMax), cimag(OmegaMax));
  fprintf(stdoutMPI, "  Omega Min. : %15.5e %15.5e\n", creal(OmegaMin), cimag(OmegaMin));
  fprintf(stdoutMPI, "  Num. of Omega : %d\n", Nomega);

  /* Build the ket A|phi> from the current single-excitation op set (X->Def). */
  for (i = 0; i <= X->Bind.Check.idim_max; i++) v0[i] = 0;
  TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_CalcExcitedStateStart, "a");
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
  dnorm = NormMPI_dc(X->Bind.Check.idim_max, v0);
  if (fabs(dnorm) < pow(10.0, -15)) {
    fprintf(stderr, "Warning: Norm of an excited (ket) vector becomes 0.\n");
    for (i = 0; i < (unsigned long int)nBra * Nomega; i++) dcSpectrum[i] = 0;
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_CalcExcitedStateEnd, "a");
    return TRUE;
  }

  /* Build the nBra bras B_b|phi> (un-normalized, single-excitation only). */
  braVecs = (double complex **)malloc(nBra * sizeof(double complex *));
  braList = (double complex **)malloc(nBra * sizeof(double complex *));
  if (braVecs == NULL || braList == NULL) {
    fprintf(stderr, "Error: out of memory allocating bra vector tables (nBra=%d).\n", nBra);
    free(braVecs); free(braList);
    return FALSE;
  }
  for (b = 0; b < nBra; b++) { braVecs[b] = NULL; braList[b] = NULL; }
  for (b = 0; b < nBra; b++) {
    braVecs[b] = cd_1d_allocate(X->Bind.Check.idim_max + 1);
    for (i = 0; i <= X->Bind.Check.idim_max; i++) braVecs[b][i] = 0;
    {
      ExcitationOperatorSet braSet = {
        braSet_N[b], braSet_ops[b], braSet_para[b], 0, NULL, NULL};
      if (GetExcitedState(&(X->Bind), &braSet, braVecs[b], v1Org) != TRUE) {
        fprintf(stderr, "Error: failed to build bra excited state %d.\n", b);
        iret = FALSE; break;
      }
    }
    if (NormMPI_dc(X->Bind.Check.idim_max, braVecs[b]) < pow(10.0, -15))
      fprintf(stderr, "Warning: Norm of bra excited vector %d is 0; its spectrum will be zero.\n", b);
    braList[b] = braVecs[b];
  }
  TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_CalcExcitedStateEnd, "a");

  if (iret == TRUE) {
    diagonalcalc(&(X->Bind));
    StartTimer(6200);
    if (X->Bind.Def.iCalcType == CG) {
      iret = CalcSpectrumByBiCG(X, v0, braList, nBra, v1, vg, Nomega, dcSpectrum, dcomega);
    }
    else {
      fprintf(stderr, "Error: multi-bra spectrum (SpectrumNumBra>1) requires CalcType=CG (BiCG).\n");
      iret = FALSE;
    }
    StopTimer(6200);
  }

  for (b = 0; b < nBra; b++) if (braVecs[b] != NULL) free_cd_1d_allocate(braVecs[b]);
  free(braVecs); free(braList);
  return iret;
}/*int CalcSpectrumMultiBra*/

/// \brief Finite-temperature spectrum loop: eigenstates idx-outer, operator sets op-inner.
///
/// For each of the iSpectrumLoopExct eigenstates it reads the eigenvector ONCE and reuses it
/// across all iSpectrumNumOp operator sets (op-inner), so the eigenvector reads drop from
/// exct*nop to exct. Operator set 0 is the namelist SingleExcitation (already in X->Def);
/// sets 1.. are loaded from single_ex_<op>.def and must map to set 0's Hilbert sector (the
/// shared sector for which MakeExcitedList already built the excited lists). The per-state
/// spectral shift is E_idx (real) plus any configured imaginary broadening. Writes one
/// DynamicalGreen_<idx>[_<op>].dat per (eigenstate, operator). X->Def's single-excitation
/// pointers are always restored before return, and the operator sets loaded here (1..) are
/// freed; the caller still owns v1Org/dcSpectrum/dcomega and the excitation lists.
///
/// \retval TRUE  success.
/// \retval FALSE allocation, file-read, eigenstate-count, or sector-mismatch failure.
static int RunMultiOpFiniteTLoop(
  struct EDMainCalStruct *X,
  double complex *v1Org,
  int Nomega,
  double complex *dcSpectrum,
  double complex *dcomega)
{
  int idx, op, b, useOp, iret = TRUE;
  int nloop = X->Bind.Def.iSpectrumLoopExct;
  int nop = X->Bind.Def.iSpectrumNumOp;
  int nBra = X->Bind.Def.iSpectrumNumBra;
  /* Set 0 aliases the X->Def originals (owned by readdef); only sets 1.. are freed here. */
  unsigned int op0_N = X->Bind.Def.NSingleExcitationOperator;
  int **op0_ops = X->Bind.Def.SingleExcitationOperator;
  double complex *op0_para = X->Bind.Def.ParaSingleExcitationOperator;
  unsigned int *set_N;
  int ***set_ops;
  double complex **set_para;
  /* Bra operator sets (Stage 3 multi-bra). Set 0 is the namelist SingleExcitationBra; sets 1..
     are loaded from single_ex_bra_<b>.def. Only allocated/used when nBra > 1. */
  unsigned int *bra_N = NULL;
  int ***bra_ops = NULL;
  double complex **bra_para = NULL;
  /* Only the REAL part of the spectral shift is state-dependent (= E_idx). Preserve any
     configured imaginary broadening in dcOmegaOrg so the loop matches the single-state path;
     with the DCore driver this imaginary part is 0, so the two are identical. */
  double dOmegaOrgIm = cimag(X->Bind.Def.dcOmegaOrg);
  double Elast = 0.0;

  if (nop < 1) nop = 1;
  if (nBra < 1) nBra = 1;
  /* Multi-bra always carries the <op> field in the output name so file names stay unambiguous. */
  useOp = (nop > 1) || (nBra > 1);

  set_N    = (unsigned int *)    malloc(sizeof(unsigned int)     * nop);
  set_ops  = (int ***)           malloc(sizeof(int **)           * nop);
  set_para = (double complex **) malloc(sizeof(double complex *) * nop);
  if (set_N == NULL || set_ops == NULL || set_para == NULL) {
    fprintf(stderr, "Error: out of memory allocating operator-set tables (SpectrumNumOp=%d).\n", nop);
    free(set_N); free(set_ops); free(set_para);
    return FALSE;
  }
  set_N[0] = op0_N; set_ops[0] = op0_ops; set_para[0] = op0_para;
  for (op = 1; op < nop; op++) { set_N[op] = 0; set_ops[op] = NULL; set_para[op] = NULL; }

  /* Load operator sets 1..; for canonical models require every set to map to set 0's
     Hilbert sector. When sh0 is not valid the model has no particle-number sector for
     this excitation (grand canonical): every single-excitation operator stays in the
     same Hilbert space, so there is no cross-operator sector constraint to enforce and
     the sets are simply loaded. For canonical models sh0 is well defined and
     MakeExcitedList already built the excited lists for that single shared sector. */
  if (nop > 1) {
    SectorShift sh0 = GetExcitationOperatorSetShift(
      X->Bind.Def.iCalcModel, X->Bind.Def.iFlgGeneralSpin, FALSE, op0_ops, op0_N);
    int enforceSector = (sh0.valid == TRUE);
    for (op = 1; iret == TRUE && op < nop; op++) {
      char opfn[D_FileNameMax];
      /* A single-excitation set has at most 2*Nsite distinct (site, spin) operators. */
      unsigned int maxN = 2u * (unsigned int) X->Bind.Def.Nsite;
      snprintf(opfn, sizeof(opfn), "single_ex_%d.def", op);
      if (ReadSingleExcitationSet(opfn, maxN, &set_N[op], &set_ops[op], &set_para[op]) != TRUE) {
        iret = FALSE; break;
      }
      if (enforceSector) {
        SectorShift sh = GetExcitationOperatorSetShift(
          X->Bind.Def.iCalcModel, X->Bind.Def.iFlgGeneralSpin, FALSE, set_ops[op], set_N[op]);
        if (sh.valid == FALSE || sh.dNe != sh0.dNe || sh.dNup != sh0.dNup ||
            sh.dNdown != sh0.dNdown || sh.dTotal2Sz != sh0.dTotal2Sz) {
          fprintf(stderr, "Error: operator set %d maps to a different Hilbert sector than set 0; "
                          "all operators in one SpectrumNumOp run must share the sector.\n", op);
          iret = FALSE; break;
        }
      }
    }
  }

  /* Load the bra operator sets for Stage-3 multi-bra (one ket solve projected onto nBra bras).
     Set 0 is the namelist SingleExcitationBra; sets 1.. come from single_ex_bra_<b>.def. Every
     bra must map to the SAME excited sector as ket set 0 (one-body c_i/c_j on the same spin) —
     the single sector MakeExcitedList built the excited lists for. The check is skipped for
     grand-canonical models, where the shift is not well defined (as for the op sets). */
  if (iret == TRUE && nBra > 1) {
    bra_N    = (unsigned int *)    malloc(sizeof(unsigned int)     * nBra);
    bra_ops  = (int ***)           malloc(sizeof(int **)           * nBra);
    bra_para = (double complex **) malloc(sizeof(double complex *) * nBra);
    if (bra_N == NULL || bra_ops == NULL || bra_para == NULL) {
      fprintf(stderr, "Error: out of memory allocating bra-set tables (SpectrumNumBra=%d).\n", nBra);
      iret = FALSE;
    }
    else {
      SectorShift sh0 = GetExcitationOperatorSetShift(
        X->Bind.Def.iCalcModel, X->Bind.Def.iFlgGeneralSpin, FALSE, op0_ops, op0_N);
      int enforceSector = (sh0.valid == TRUE);
      unsigned int maxN = 2u * (unsigned int) X->Bind.Def.Nsite;
      bra_N[0]    = X->Bind.Def.NSingleExcitationOperatorBra;
      bra_ops[0]  = X->Bind.Def.SingleExcitationOperatorBra;
      bra_para[0] = X->Bind.Def.ParaSingleExcitationOperatorBra;
      for (b = 1; b < nBra; b++) { bra_N[b] = 0; bra_ops[b] = NULL; bra_para[b] = NULL; }
      if (bra_N[0] == 0) {
        fprintf(stderr, "Error: SpectrumNumBra>1 requires a namelist SingleExcitationBra (bra set 0).\n");
        iret = FALSE;
      }
      for (b = 0; iret == TRUE && b < nBra; b++) {
        if (b > 0) {
          char brafn[D_FileNameMax];
          snprintf(brafn, sizeof(brafn), "single_ex_bra_%d.def", b);
          if (ReadSingleExcitationSet(brafn, maxN, &bra_N[b], &bra_ops[b], &bra_para[b]) != TRUE) {
            iret = FALSE; break;
          }
        }
        if (enforceSector) {
          SectorShift sh = GetExcitationOperatorSetShift(
            X->Bind.Def.iCalcModel, X->Bind.Def.iFlgGeneralSpin, FALSE, bra_ops[b], bra_N[b]);
          if (sh.valid == FALSE || sh.dNe != sh0.dNe || sh.dNup != sh0.dNup ||
              sh.dNdown != sh0.dNdown || sh.dTotal2Sz != sh0.dTotal2Sz) {
            fprintf(stderr, "Error: bra set %d maps to a different excited sector than ket set 0; "
                            "all bras must share the ket's sector.\n", b);
            iret = FALSE; break;
          }
        }
      }
    }
  }

  /* Fail before writing any spectrum if SpectrumLoopExct exceeds the number of computed
     eigenstates (the energy file has one entry per eigenstate, matching the eigenvector
     files), so no partial DynamicalGreen output is left behind. */
  if (iret == TRUE && ReadEigenEnergy(X, nloop - 1, &Elast) != TRUE) {
    fprintf(stderr, "Error: SpectrumLoopExct=%d exceeds the number of computed eigenstates.\n", nloop);
    iret = FALSE;
  }

  for (idx = 0; iret == TRUE && idx < nloop; idx++) {
    double Eidx = 0.0;
    if (ReadEigenVector(X, idx, TRUE, v1Org) != TRUE) { iret = FALSE; break; } /* read ONCE per eigenstate */
    if (ReadEigenEnergy(X, idx, &Eidx) != TRUE)        { iret = FALSE; break; }
    for (op = 0; op < nop; op++) {
      /* Point X->Def at this operator set; the ket is built from X->Def. */
      X->Bind.Def.NSingleExcitationOperator    = set_N[op];
      X->Bind.Def.SingleExcitationOperator     = set_ops[op];
      X->Bind.Def.ParaSingleExcitationOperator = set_para[op];
      if (nBra > 1) {
        /* One ket solve, nBra bra projections (Stage 3). dcSpectrum holds x(nBra,Nomega). */
        iret = CalcSpectrumMultiBra(X, v1Org, Eidx + dOmegaOrgIm * I, Nomega,
                                    dcSpectrum, dcomega, nBra, bra_N, bra_ops, bra_para);
        if (iret != TRUE) break;
        for (b = 0; b < nBra; b++) {
          if (OutputSpectrumIdx(X, idx, op, useOp, b, TRUE, nBra, Nomega, dcSpectrum, dcomega) != TRUE) {
            iret = FALSE; break;
          }
        }
        if (iret != TRUE) break;
      }
      else {
        iret = CalcOneSpectrum(X, v1Org, Eidx + dOmegaOrgIm * I, Nomega, dcSpectrum, dcomega);
        if (iret != TRUE) break;
        if (OutputSpectrumIdx(X, idx, op, useOp, 0, FALSE, 1, Nomega, dcSpectrum, dcomega) != TRUE) { iret = FALSE; break; }
      }
    }
  }

  /* Restore set 0 so the caller's teardown frees the original arrays, and release sets 1.. */
  X->Bind.Def.NSingleExcitationOperator    = op0_N;
  X->Bind.Def.SingleExcitationOperator     = op0_ops;
  X->Bind.Def.ParaSingleExcitationOperator = op0_para;
  for (op = 1; op < nop; op++) {
    if (set_ops[op]  != NULL) free_i_2d_allocate(set_ops[op]);
    if (set_para[op] != NULL) free_cd_1d_allocate(set_para[op]);
  }
  free(set_N); free(set_ops); free(set_para);
  /* Release bra sets 1.. (set 0 aliases the X->Def namelist arrays, owned by readdef). */
  if (bra_N != NULL) {
    for (b = 1; b < nBra; b++) {
      if (bra_ops != NULL && bra_ops[b]  != NULL) free_i_2d_allocate(bra_ops[b]);
      if (bra_para != NULL && bra_para[b] != NULL) free_cd_1d_allocate(bra_para[b]);
    }
  }
  free(bra_N); free(bra_ops); free(bra_para);
  return iret;
}/*int RunMultiOpFiniteTLoop*/

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
    /* Multi-bra mode fills nBra spectra per BiCG solve (Komega x(nBra,Nomega)); size the
       buffer accordingly. Single-bra / single-state paths use the leading Nomega slice. */
    {
      int nBraAlloc = (X->Bind.Def.iSpectrumNumBra > 1) ? X->Bind.Def.iSpectrumNumBra : 1;
      dcSpectrum = cd_1d_allocate((unsigned long)nBraAlloc * Nomega);
    }
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
  /* Reject out-of-range Spectrum* keywords up front rather than silently normalizing them
     (later code clamps SpectrumNumOp/Bra<1 to 1 and treats SpectrumLoopExct<0 as legacy),
     so a misconfigured input fails loudly instead of running a different calculation. */
  if (X->Bind.Def.iSpectrumLoopExct < 0) {
    fprintf(stderr, "Error: SpectrumLoopExct must be >= 0 (got %d).\n", X->Bind.Def.iSpectrumLoopExct);
    exitMPI(-1);
  }
  if (X->Bind.Def.iSpectrumNumOp < 1) {
    fprintf(stderr, "Error: SpectrumNumOp must be >= 1 (got %d).\n", X->Bind.Def.iSpectrumNumOp);
    exitMPI(-1);
  }
  if (X->Bind.Def.iSpectrumNumBra < 1) {
    fprintf(stderr, "Error: SpectrumNumBra must be >= 1 (got %d).\n", X->Bind.Def.iSpectrumNumBra);
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
  /* Multiple operator sets per run (SpectrumNumOp>1) only take effect inside the finite-T
     eigenstate loop (op-inner reuse of the eigenvector). It swaps only the ket single
     -excitation set per operator, so it is incompatible with the bra (direct off-diagonal)
     path, which is not swapped. */
  if (X->Bind.Def.iSpectrumNumOp > 1) {
    if (X->Bind.Def.iSpectrumLoopExct == 0) {
      fprintf(stderr, "Error: SpectrumNumOp>1 requires SpectrumLoopExct>0 "
                      "(multiple operator sets are processed only inside the finite-T loop).\n");
      exitMPI(-1);
    }
    /* The multi-operator loop seeds set 0 from SingleExcitation and loads single_ex_<op>.def;
       it is single-excitation-specific (pair excitation is not handled). */
    if (X->Bind.Def.NSingleExcitationOperator == 0 || X->Bind.Def.NPairExcitationOperator > 0) {
      fprintf(stderr, "Error: SpectrumNumOp>1 is only supported with SingleExcitation operators "
                      "(not PairExcitation).\n");
      exitMPI(-1);
    }
    /* The single-bra off-diagonal path is not swapped per ket-op, so it is incompatible with
       SpectrumNumOp>1 — UNLESS multi-bra mode (SpectrumNumBra>1) is active, in which case the
       bras are handled per ket-op by the multi-bra projection and ket x bra compose. */
    if ((X->Bind.Def.NSingleExcitationOperatorBra > 0 || X->Bind.Def.NPairExcitationOperatorBra > 0)
        && X->Bind.Def.iSpectrumNumBra <= 1) {
      fprintf(stderr, "Error: SpectrumNumOp>1 is not supported with a single bra excitation operator "
                      "(use SpectrumNumBra>1 to project one ket solve onto multiple bras).\n");
      exitMPI(-1);
    }
  }
  /* Multi-bra (SpectrumNumBra>1): one ket BiCG solve projected onto nBra bras via Komega nl.
     Like SpectrumNumOp it only runs inside the finite-T loop, needs BiCG (CalcType CG), a
     single-excitation ket, and a namelist SingleExcitationBra as bra set 0. It composes with
     SpectrumNumOp (the full ket x bra grid is processed in one run). */
  if (X->Bind.Def.iSpectrumNumBra > 1) {
    if (X->Bind.Def.iSpectrumLoopExct == 0) {
      fprintf(stderr, "Error: SpectrumNumBra>1 requires SpectrumLoopExct>0 "
                      "(multiple bras are projected only inside the finite-T loop).\n");
      exitMPI(-1);
    }
    if (X->Bind.Def.iCalcType != CG) {
      fprintf(stderr, "Error: SpectrumNumBra>1 requires the BiCG solver (CalcType=CG).\n");
      exitMPI(-1);
    }
    if (X->Bind.Def.NSingleExcitationOperator == 0 || X->Bind.Def.NPairExcitationOperator > 0) {
      fprintf(stderr, "Error: SpectrumNumBra>1 is only supported with SingleExcitation kets "
                      "(not PairExcitation).\n");
      exitMPI(-1);
    }
    if (X->Bind.Def.NSingleExcitationOperatorBra == 0 || X->Bind.Def.NPairExcitationOperatorBra > 0) {
      fprintf(stderr, "Error: SpectrumNumBra>1 requires a namelist SingleExcitationBra (bra set 0) "
                      "and does not support PairExcitationBra.\n");
      exitMPI(-1);
    }
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
        /* Finite-temperature loop: one HPhi launch covers all the thermally-relevant
           eigenstates (idx-outer) and, for each, all operator sets (op-inner), reading
           every eigenvector once. Replaces DCore's per-(idx, operator) relaunch. */
        int nop_print = X->Bind.Def.iSpectrumNumOp < 1 ? 1 : X->Bind.Def.iSpectrumNumOp;
        fprintf(stdoutMPI, "  Start: finite-T spectrum loop over %d eigenstate(s) x %d operator(s).\n",
                X->Bind.Def.iSpectrumLoopExct, nop_print);
        TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_InputEigenVectorStart, "a");
        StopTimer(6100);
        iret = RunMultiOpFiniteTLoop(X, v1Org, Nomega, dcSpectrum, dcomega);
        StopTimer(6101);
        TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_InputEigenVectorEnd, "a");
        bAlreadyOutput = TRUE; /* the loop wrote one file per (eigenstate, operator); skip the single output below */
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
      case CG: {
        /* Restart path: diagonal only (bra == ket). */
        double complex *braList[1] = { v0 };
        iret = CalcSpectrumByBiCG(X, v0, braList, 1, v1, vg, Nomega, dcSpectrum, dcomega);
        break;
      }
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
