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
/**
 * @file CalcByTEM.c
 *
 * @brief Real-time evolution calculation using Suzuki-Trotter decomposition
 *
 * Implements time evolution of quantum states:
 *   |psi(t+dt)> = exp(-i H dt) |psi(t)>
 *
 * Method:
 * Uses Krylov subspace method (similar to Lanczos) for the matrix exponential.
 * The time evolution operator is approximated using a polynomial expansion
 * in the Krylov basis.
 *
 * Time-dependent Hamiltonian:
 * - Supports time-dependent transfer integrals (laser pulses, etc.)
 * - TETransfer and TEInterAll parameters read from input files
 * - Hamiltonian updated at each time step via MakeTEDTransfer/MakeTEDInterAll
 *
 * Workflow:
 * 1. Read initial state from file (must be provided)
 * 2. For each time step:
 *    a. Update time-dependent Hamiltonian if needed
 *    b. Apply exp(-i H dt) via Krylov approximation
 *    c. Every ExpecInterval steps, compute observables
 *
 * Output:
 * - Time-dependent observables (energy, correlations)
 * - Wavefunction at specified time points
 *
 * @author Kota Ido (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
#include "Common.h"
#include "readdef.h"
#include "FirstMultiply.h"
#include "Multiply.h"
#include "diagonalcalc.h"
#include "expec_energy_flct.h"
#include "expec_cisajs.h"
#include "expec_cisajscktaltdc.h"
#include "nbody_correlation.h"
#include "CalcByTEM.h"
#include "FileIO.h"
#include "wrapperMPI.h"
#include "HPhiTrans.h"

void MakeTEDTransfer(struct BindStruct *X, const int timeidx);
void MakeTEDInterAll(struct BindStruct *X, const int timeidx);

/**
 * @brief Main driver for real-time evolution calculation
 *
 * Evolves the initial state through NTETimeSteps time steps, computing
 * physical observables at intervals.
 *
 * Prerequisites:
 * - Initial wavefunction must be provided (iInputEigenVec = TRUE)
 * - NTETimeSteps must be larger than Lanczos_max
 *
 * Time evolution per step:
 * 1. If time-dependent terms exist, update Hamiltonian
 * 2. Apply exp(-i H dt) using Multiply() function
 * 3. Normalize if needed
 *
 * Observable output:
 * - SS_*.dat: Energy, spin correlations vs time
 * - Norm_*.dat: Norm evolution (should remain ~1)
 *
 * @param ExpecInterval Steps between observable calculations [in]
 * @param X Calculation parameters and state [in,out]
 *
 * @return 0 on success, -1 on error
 *
 * @author Kota Ido (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int CalcByTEM(
        const int ExpecInterval,
        struct EDMainCalStruct *X
) {
  char *defname;
  char sdt[D_FileNameMax];
  char sdt_phys[D_FileNameMax];
  char sdt_norm[D_FileNameMax];
  char sdt_flct[D_FileNameMax];
  int rand_i=0;
  int step_initial = 0;
  long int i_max = 0;
  FILE *fp;
  double Time = X->Bind.Def.Param.Tinit;
  double dt = ((X->Bind.Def.NLaser==0)? 0.0: X->Bind.Def.Param.TimeSlice);

  if(X->Bind.Def.NTETimeSteps < X->Bind.Def.Lanczos_max){
    fprintf(stdoutMPI, "Error: NTETimeSteps must be larger than Lanczos_max.\n");
    return -1;
  }
  step_spin = ExpecInterval;
  X->Bind.Def.St = 0;
  fprintf(stdoutMPI, "%s", cLogTEM_Start);
  if (X->Bind.Def.iInputEigenVec == FALSE) {
    fprintf(stderr, "Error: A file of Inputvector is not inputted.\n");
    return -1;
  } else {
    //input v1
    fprintf(stdoutMPI, "%s","An Initial Vector is inputted.\n");
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_InputEigenVectorStart, "a");
    GetFileNameByKW(KWSpectrumVec, &defname);
    strcat(defname, "_rank_%d.dat");
    sprintf(sdt, defname, myrank);
    childfopenALL(sdt, "rb", &fp);
    if (fp == NULL) {
      fprintf(stderr, "Error: A file of Inputvector does not exist.\n");
      fclose(fp);
      exitMPI(-1);
    }
    fread(&step_initial, sizeof(int), 1, fp);
    fread(&i_max, sizeof(long int), 1, fp);
    if (i_max != X->Bind.Check.idim_max) {
      fprintf(stderr, "Error: A file of Inputvector is incorrect.\n");
      fclose(fp);
      exitMPI(-1);
    }
    fread(v1, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
    fclose(fp);
    if (X->Bind.Def.iReStart == RESTART_NOT || X->Bind.Def.iReStart == RESTART_OUT) {
      step_initial = 0;
    }
  }

  if(X->Bind.Def.iOutputDataHead==1){
    sprintf(sdt_phys, "%s_%s", X->Bind.Def.CDataFileHead, cFileNameSS);
  }else{
    sprintf(sdt_phys, "%s", cFileNameSS);
  }
  if (childfopenMPI(sdt_phys, "w", &fp) != 0) {
    return -1;
  }
  fprintf(fp, "%s",cLogSS);
  fclose(fp);

  if(X->Bind.Def.iOutputDataHead==1){
    sprintf(sdt_norm, "%s_%s", X->Bind.Def.CDataFileHead, cFileNameNorm);
  }else{
    sprintf(sdt_norm, "%s", cFileNameNorm);
  }
  if (childfopenMPI(sdt_norm, "w", &fp) != 0) {
    return -1;
  }
  fprintf(fp, "%s",cLogNorm);
  fclose(fp);

  if(X->Bind.Def.iOutputDataHead==1){
    sprintf(sdt_flct, "%s_%s", X->Bind.Def.CDataFileHead, cFileNameFlct);
  }else{
    sprintf(sdt_flct, "%s", cFileNameFlct);
  }
  if (childfopenMPI(sdt_flct, "w", &fp) != 0) {
    return -1;
  }
  fprintf(fp, "%s",cLogFlct);
  fclose(fp);


  int iInterAllOffDiagonal_org = X->Bind.Def.NInterAll_OffDiagonal;
  int iTransfer_org = X->Bind.Def.EDNTransfer;
  int progress_interval = X->Bind.Def.Lanczos_max / 10;
  if (progress_interval < 1) progress_interval = 1;  /* avoid % 0 when Lanczos_max < 10 */
  for (step_i = step_initial; step_i < X->Bind.Def.Lanczos_max; step_i++) {
    X->Bind.Def.istep = step_i;

    //Reset total number of interactions (changed in MakeTED***function.)
    X->Bind.Def.EDNTransfer = iTransfer_org;
    X->Bind.Def.NInterAll_OffDiagonal = iInterAllOffDiagonal_org;

    if (step_i % progress_interval == 0) {
      fprintf(stdoutMPI, cLogTEStep, step_i, X->Bind.Def.Lanczos_max);
    }

    if(X->Bind.Def.NLaser !=0) {
      TransferWithPeierls(&(X->Bind), Time);
    }
    else {
      // common procedure
      Time = X->Bind.Def.TETime[step_i];
      if (step_i == 0) dt = 0.0;
      else {
        dt = X->Bind.Def.TETime[step_i] - X->Bind.Def.TETime[step_i - 1];
      }
      X->Bind.Def.Param.TimeSlice = dt;

      // Set interactions
      if(X->Bind.Def.NTETransferMax != 0 && X->Bind.Def.NTEInterAllMax!=0){
        fprintf(stdoutMPI,
                "Error: Time Evolution mode does not support TEOneBody and TETwoBody interactions at the same time. \n");
        return -1;
      }
      else if (X->Bind.Def.NTETransferMax > 0) { //One-Body type
        MakeTEDTransfer(&(X->Bind), step_i);
      }else if (X->Bind.Def.NTEInterAllMax > 0) { //Two-Body type
        MakeTEDInterAll(&(X->Bind), step_i);
      }
      //[e] Yoshimi
    }

    if(step_i == step_initial){
      TimeKeeperWithStep(&(X->Bind), cFileNameTEStep, cTEStep, "w", step_i);
    }
    else {
      TimeKeeperWithStep(&(X->Bind), cFileNameTEStep, cTEStep, "a", step_i);
    }
    MultiplyForTEM(&(X->Bind));
    //Add Diagonal Parts
    //Multiply Diagonal
    expec_energy_flct(&(X->Bind));

    if(X->Bind.Def.NLaser >0 ) Time+=dt;
    if (childfopenMPI(sdt_phys, "a", &fp) != 0) {
      return -1;
    }
    fprintf(fp, "%.16lf  %.16lf %.16lf %.16lf %.16lf %d\n", Time, X->Bind.Phys.energy, X->Bind.Phys.var,
            X->Bind.Phys.doublon, X->Bind.Phys.num, step_i);
    fclose(fp);

    if (childfopenMPI(sdt_norm, "a", &fp) != 0) {
      return -1;
    }
    fprintf(fp, "%.16lf %.16lf %d\n", Time, global_norm, step_i);
    fclose(fp);

    if (childfopenMPI(sdt_flct, "a", &fp) != 0) {
      return -1;
    }
    fprintf(fp, "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %d\n", Time,X->Bind.Phys.num,X->Bind.Phys.num2, X->Bind.Phys.doublon,X->Bind.Phys.doublon2, X->Bind.Phys.Sz,X->Bind.Phys.Sz2,step_i);
    fclose(fp);



    if (step_i % step_spin == 0) {
      if (expec_cisajs(&(X->Bind), v1) != 0) {
        return -1;
      }
      if (expec_cisajscktaltdc(&(X->Bind), v1) != 0) {
        return -1;
      }
      if (expec_nbodyg(&(X->Bind), v1) != 0) {
        return -1;
      }
    }
    if (X->Bind.Def.iOutputEigenVec == TRUE) {
      if (step_i % X->Bind.Def.Param.OutputInterval == 0) {
        sprintf(sdt, cFileNameOutputEigen, X->Bind.Def.CDataFileHead, step_i, myrank);
        if (childfopenALL(sdt, "wb", &fp) != 0) {
          fclose(fp);
          exitMPI(-1);
        }
        fwrite(&step_i, sizeof(step_i), 1, fp);
        fwrite(&X->Bind.Check.idim_max, sizeof(long int), 1, fp);
        fwrite(v1, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
        fclose(fp);
      }
    }
  }
  if (X->Bind.Def.iOutputEigenVec == TRUE) {
    sprintf(sdt, cFileNameOutputEigen, X->Bind.Def.CDataFileHead, rand_i, myrank);
    if (childfopenALL(sdt, "wb", &fp) != 0) {
      fclose(fp);
      exitMPI(-1);
    }
    fwrite(&step_i, sizeof(step_i), 1, fp);
    fwrite(&X->Bind.Check.idim_max, sizeof(long int), 1, fp);
    fwrite(v1, sizeof(complex double), X->Bind.Check.idim_max + 1, fp);
    fclose(fp);
  }

  fprintf(stdoutMPI, "%s",cLogTEM_End);
  return 0;
}

/// \brief Set transfer integrals at timeidx-th time
/// \param X struct for getting information of transfer integrals
/// \param timeidx index of time
void MakeTEDTransfer(struct BindStruct *X, const int timeidx) {
  int i,j;
  //Clear values
  for(i=0; i<X->Def.NTETransferMax ;i++) {
    for(j =0; j<4; j++) {
      X->Def.EDGeneralTransfer[i + X->Def.EDNTransfer][j] = 0;
    }
    X->Def.EDParaGeneralTransfer[i+X->Def.EDNTransfer]=0.0;
  }

  //Input values
  for(i=0; i<X->Def.NTETransfer[timeidx] ;i++){
    for(j =0; j<4; j++) {
      X->Def.EDGeneralTransfer[i + X->Def.EDNTransfer][j] = X->Def.TETransfer[timeidx][i][j];
    }
    X->Def.EDParaGeneralTransfer[i+X->Def.EDNTransfer]=X->Def.ParaTETransfer[timeidx][i];
  }
  X->Def.EDNTransfer += X->Def.NTETransfer[timeidx];
}

/// \brief Set interall interactions at timeidx-th time
/// \param X struct for getting information of interall interactions
/// \param timeidx index of time
void MakeTEDInterAll(struct BindStruct *X, const int timeidx) {
  int i, j;
  //Clear values
  for (i = 0; i < X->Def.NTEInterAllMax; i++) {
    for (j = 0; j < 8; j++) {
      X->Def.InterAll_OffDiagonal[i + X->Def.NInterAll_OffDiagonal][j] = 0;
    }
    X->Def.ParaInterAll_OffDiagonal[i + X->Def.NInterAll_OffDiagonal] = 0.0;
  }

  //Input values
  for (i = 0; i < X->Def.NTEInterAllOffDiagonal[timeidx]; i++) {
    for (j = 0; j < 8; j++) {
      X->Def.InterAll_OffDiagonal[i + X->Def.NInterAll_OffDiagonal][j] = X->Def.TEInterAllOffDiagonal[timeidx][i][j];
    }
    X->Def.ParaInterAll_OffDiagonal[i + X->Def.NInterAll_OffDiagonal] = X->Def.ParaTEInterAllOffDiagonal[timeidx][i];
  }
  X->Def.NInterAll_OffDiagonal += X->Def.NTEInterAllOffDiagonal[timeidx];
}
