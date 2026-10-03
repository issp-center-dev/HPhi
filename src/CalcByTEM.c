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
 * @brief Real-time evolution calculation using a Taylor polynomial
 *
 * Implements time evolution of quantum states:
 *   |psi(t+dt)> = exp(-i H dt) |psi(t)>
 *
 * Method:
 * Uses a truncated Taylor polynomial of the matrix exponential, followed
 * by normalization.
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
 *    b. Apply the Taylor polynomial of exp(-i H dt)
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
#include "anomalous_pair.h"
#include "green_output.h"
#include "CalcByTEM.h"
#include "FileIO.h"
#include "wrapperMPI.h"
#include "HPhiTrans.h"
#include <math.h>
#include <limits.h>
#include <inttypes.h>
#include "symmetry_checkpoint.h"
#include "symmetry_te.h"

void MakeTEDTransfer(struct BindStruct *X, const int timeidx);
void MakeTEDInterAll(struct BindStruct *X, const int timeidx);

/* A sector checkpoint is a state import, not a continuation of its solver clock. */
static int record_sector_te(const struct BindStruct *X, const struct SymmetryCheckpointInfo *source,
                             const struct SymmetryTEHamiltonian *view)
{
  char name[D_FileNameMax];
  FILE *fp = NULL;
  int error = 0;
  if (myrank == 0) {
    int n = snprintf(name, sizeof(name), "%s%ssymmetry_sector.dat",
                     X->Def.iOutputDataHead ? X->Def.CDataFileHead : "",
                     X->Def.iOutputDataHead ? "_" : "");
    error = n < 0 || (size_t)n + strlen(cParentOutputFolder) >= sizeof(name);
    if (!error) error = childfopenMPI(name, "a", &fp) != 0;
    if (!error) {
      fprintf(fp, "te_hamiltonian=%s\nte_steps=%u\nexpand_coef=%d\n"
              "source_method=%" PRIu64 "\nsource_state=%" PRIu64 "\nsource_step=%" PRIu64
              "\nsource_time=%.17g\nsource_hamiltonian_digest=%016" PRIx64 "\n",
              view ? "time_dependent" : "static", X->Def.Lanczos_max, X->Def.Param.ExpandCoef, source->source_method,
              source->state_index, source->step, source->time, source->hamiltonian_digest);
      for (unsigned int i = 0; i < X->Def.Lanczos_max; ++i) {
        fprintf(fp, "te_time_%u=%.17g\n", i, SymmetryTETime(&X->Def, i));
        if (view) fprintf(fp, "te_hamiltonian_%u=%016" PRIx64 "\n", i, view->digests[i]);
      }
      if (view) fprintf(fp, "te_integrator=right_endpoint_taylor\nte_plan_update=rebuild_with_fixed_basis\n");
      if (ferror(fp)) error = 1;
      if (fclose(fp) != 0) error = 1;
    }
  }
  if (SumMPI_i(error) != 0) {
    fprintf(stdoutMPI, "Error: failed to record symmetry TE schedule.\n");
    return -1;
  }
  return 0;
}

static int write_sector_te_vector(struct BindStruct *X, int step, double time,
                                  uint64_t state, int final)
{
  char name[D_FileNameMax];
  int n = final ? snprintf(name, sizeof(name), "%s_eigenvec_final_rank_%d.dat", X->Def.CDataFileHead, myrank)
                : snprintf(name, sizeof(name), "%s_eigenvec_%d_rank_%d.dat", X->Def.CDataFileHead, step, myrank);
  if (SumMPI_i(n < 0 || (size_t)n >= sizeof(name)) != 0) return -1;
  struct SymmetryCheckpointInfo info = {TimeEvolution, state, (uint64_t)step, time, 0};
  return WriteSymmetryCheckpoint(X, name, v1, &info);
}

static int write_raw_te_vector(struct BindStruct *X, int label, int next_step)
{
  char name[D_FileNameMax];
  FILE *fp = NULL;
  int length = snprintf(name, sizeof(name), "%s_eigenvec_%d_rank_%d.dat",
                        X->Def.CDataFileHead, label, myrank);
  int failed = length < 0 || (size_t)length + strlen(cParentOutputFolder) >= sizeof(name);
  if (SumMPI_i(failed) != 0) return -1;
  failed = childfopenALL(name, "wb", &fp) != 0;
  if (!failed) {
    failed = fwrite(&next_step, sizeof(next_step), 1, fp) != 1 ||
             fwrite(&X->Check.idim_max, sizeof(X->Check.idim_max), 1, fp) != 1 ||
             fwrite(v1, sizeof(*v1), X->Check.idim_max + 1, fp) != X->Check.idim_max + 1;
    if (fclose(fp) != 0) failed = 1;
  }
  if (SumMPI_i(failed) != 0) {
    fprintf(stdoutMPI, "Error: failed to write TE vector on one or more ranks.\n");
    return -1;
  }
  return 0;
}

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
static int calc_by_tem(
        const int ExpecInterval,
        struct EDMainCalStruct *X,
        struct SymmetryTEHamiltonian *view
) {
  struct SymmetryCheckpointInfo source = {0};
  char *defname = NULL;
  char sdt[D_FileNameMax];
  char sdt_phys[D_FileNameMax];
  char sdt_norm[D_FileNameMax];
  char sdt_flct[D_FileNameMax];
  int rand_i=0;
  int step_initial = 0;
  long int i_max = 0;
  FILE *fp = NULL;
  double Time = X->Bind.Def.Param.Tinit;
  double dt = ((X->Bind.Def.NLaser==0)? 0.0: X->Bind.Def.Param.TimeSlice);

  int invalid = X->Bind.Def.Param.ExpandCoef < 1 ||
                X->Bind.Def.Param.ExpandCoef == INT_MAX || ExpecInterval <= 0 ||
                (X->Bind.Def.iOutputEigenVec && X->Bind.Def.Param.OutputInterval <= 0);
  if (SumMPI_i(invalid) != 0) {
    fprintf(stdoutMPI, "Error: TE requires positive ExpandCoef, ExpecInterval and vector OutputInterval.\n");
    return -1;
  }
  if (!isfinite(Time) || !isfinite(dt) || dt < 0) {
    fprintf(stdoutMPI, "Error: TE requires finite initial time and nonnegative time step.\n");
    return -1;
  }
  if (!(view && X->Bind.Def.NLaser) && X->Bind.Def.NTETimeSteps < X->Bind.Def.Lanczos_max){
    fprintf(stdoutMPI, "Error: NTETimeSteps must be larger than Lanczos_max.\n");
    return -1;
  }
  if (X->Bind.Def.NLaser == 0) {
    invalid = X->Bind.Def.TETime == NULL;
    for (unsigned int i = 0; !invalid && i < X->Bind.Def.Lanczos_max; ++i) {
      double time = X->Bind.Def.TETime[i];
      invalid = !isfinite(time) ||
                (i > 0 && (time < X->Bind.Def.TETime[i-1] ||
                           !isfinite(time-X->Bind.Def.TETime[i-1])));
    }
    if (SumMPI_i(invalid) != 0) {
      fprintf(stdoutMPI, "Error: TE time grid must be finite and nondecreasing.\n");
      return -1;
    }
  }
  if (view && ValidateSymmetryTESchedule(&(X->Bind), view) != 0) return -1;
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
    invalid = GetFileNameByKW(KWSpectrumVec, &defname) != 0 || defname == NULL;
    int length = invalid ? -1 : snprintf(sdt, sizeof(sdt), "%s_rank_%d.dat", defname, myrank);
    if (SumMPI_i(length < 0 || (size_t)length + strlen(cParentOutputFolder) >= sizeof(sdt)) != 0) {
      fprintf(stdoutMPI, "Error: missing or overlong TE SpectrumVec filename.\n");
      return -1;
    }
    if (X->Bind.Def.iFlgSymmetryBasis) {
      if (ReadSymmetryCheckpoint(&(X->Bind), sdt, v1, &source) != 0) return -1;
      if (record_sector_te(&(X->Bind), &source, view) != 0) return -1;
    } else {
      invalid = childfopenALL(sdt, "rb", &fp) != 0;
      if (!invalid) {
        invalid = fread(&step_initial, sizeof(step_initial), 1, fp) != 1 ||
                  fread(&i_max, sizeof(i_max), 1, fp) != 1 ||
                  i_max < 0 || (unsigned long)i_max != X->Bind.Check.idim_max;
        if (!invalid) invalid = fread(v1, sizeof(*v1), X->Bind.Check.idim_max + 1, fp)
                                    != X->Bind.Check.idim_max + 1;
        if (!invalid && (fgetc(fp) != EOF || ferror(fp))) invalid = 1;
        if (fclose(fp) != 0) invalid = 1;
        fp = NULL;
      }
      if (SumMPI_i(invalid) != 0) {
        fprintf(stdoutMPI, "Error: missing, truncated or incompatible TE Inputvector file.\n");
        return -1;
      }
      double norm = 0;
      for (unsigned long i = 1; i <= X->Bind.Check.idim_max; ++i) {
        double re = creal(v1[i]), im = cimag(v1[i]);
        if (!isfinite(re) || !isfinite(im)) invalid = 1;
        norm += re*re + im*im;
      }
      norm = SumMPI_d(norm);
      if (SumMPI_i(invalid || !isfinite(norm) || norm <= 0) != 0) {
        fprintf(stdoutMPI, "Error: TE Inputvector has zero or non-finite global norm.\n");
        return -1;
      }
      if (X->Bind.Def.iReStart == RESTART_NOT || X->Bind.Def.iReStart == RESTART_OUT)
        step_initial = 0;
      int root_step = BcastMPI_i(0, step_initial);
      if (SumMPI_i(step_initial < 0 || (unsigned int)step_initial >= X->Bind.Def.Lanczos_max ||
                   step_initial != root_step) != 0) {
        fprintf(stdoutMPI, "Error: TE restart step must be consistent across ranks and within [0, Lanczos_max).\n");
        return -1;
      }
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

  if (X->Bind.Def.iReStart == RESTART_NOT || X->Bind.Def.iReStart == RESTART_OUT) {
    if (GreenOutputInitializeAggregateFiles(&(X->Bind)) != 0) {
      return -1;
    }
  }

  int iInterAllOffDiagonal_org = X->Bind.Def.NInterAll_OffDiagonal;
  int iTransfer_org = X->Bind.Def.EDNTransfer;
  if (step_initial > 0) {
    /* Raw TE files store the next step. Restore H|psi> at the preceding
     * time point: MultiplyForTEM consumes both psi (v1) and H psi (v0).
     * Without this, a restarted run silently loses its first Taylor term. */
    step_i = X->Bind.Def.istep = step_initial - 1;
    if (X->Bind.Def.NLaser > 0) {
      double previous_time = Time + (step_initial - 1) * dt;
      Time += step_initial * dt;
      if (!isfinite(previous_time) || !isfinite(Time)) return -1;
      TransferWithPeierls(&(X->Bind), previous_time);
    } else if (X->Bind.Def.NTETransferMax > 0) {
      MakeTEDTransfer(&(X->Bind), step_initial - 1);
    } else if (X->Bind.Def.NTEInterAllMax > 0) {
      MakeTEDInterAll(&(X->Bind), step_initial - 1);
    }
    for (unsigned long i = 1; i <= X->Bind.Check.idim_max; ++i) v0[i] = v1[i];
    if (expec_energy_flct(&(X->Bind)) != 0) return -1;
  }
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

    if (view) {
      Time = SymmetryTETime(&view->base, step_i);
      dt = step_i ? Time - SymmetryTETime(&view->base, step_i - 1) : 0;
      if (SelectSymmetryTEHamiltonian(&(X->Bind), view, step_i, Time) != 0 ||
          RebuildSymmetryTEPlan(&(X->Bind)) != 0) return -1;
      X->Bind.Def.Param.TimeSlice = dt;
      /* All Taylor powers use the SAME current H, including the first one. */
      for (unsigned long i = 1; i <= X->Bind.Check.idim_max; ++i) v0[i] = v1[i];
      if (expec_energy_flct(&(X->Bind)) != 0) return -1;
    }
    else if(X->Bind.Def.NLaser !=0) {
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
    if (MultiplyForTEM(&(X->Bind)) != 0) return -1;
    //Add Diagonal Parts
    //Multiply Diagonal
    if (expec_energy_flct(&(X->Bind)) != 0) return -1;

    if(!view && X->Bind.Def.NLaser >0 ) Time+=dt;
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
      if (expec_anomalousg(&(X->Bind), v1) != 0) {
        return -1;
      }
    }
    if (X->Bind.Def.iOutputEigenVec == TRUE) {
      if (step_i % X->Bind.Def.Param.OutputInterval == 0) {
        if (X->Bind.Def.iFlgSymmetryBasis) {
          if (write_sector_te_vector(&(X->Bind), step_i, Time, source.state_index, 0) != 0) return -1;
        } else if (write_raw_te_vector(&(X->Bind), step_i, step_i + 1) != 0) return -1;
      }
    }
  }
  if (X->Bind.Def.iOutputEigenVec == TRUE) {
    if (X->Bind.Def.iFlgSymmetryBasis) {
      if (write_sector_te_vector(&(X->Bind), step_i - 1, Time, source.state_index, 1) != 0) return -1;
    } else if (write_raw_te_vector(&(X->Bind), rand_i, step_i) != 0) return -1;
  }

  fprintf(stdoutMPI, "%s",cLogTEM_End);
  return 0;
}

int CalcByTEM(const int ExpecInterval, struct EDMainCalStruct *X)
{
  if (!SymmetryTEIsDynamic(&X->Bind.Def)) return calc_by_tem(ExpecInterval, X, NULL);
  struct SymmetryTEHamiltonian view = {0};
  view.base = X->Bind.Def;
  int result = calc_by_tem(ExpecInterval, X, &view);
  X->Bind.Def = view.base;
  FreeSymmetryTEHamiltonian(&view);
  return result;
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
