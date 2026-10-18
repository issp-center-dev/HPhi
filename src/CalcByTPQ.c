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
 * @file CalcByTPQ.c
 *
 * @brief Thermal Pure Quantum (TPQ) state calculation
 *
 * TPQ is a method to compute finite-temperature properties without
 * explicit thermal averaging. A single "thermal pure quantum state"
 * |psi_beta> represents the thermal ensemble at inverse temperature beta.
 *
 * Algorithm:
 * Starting from random state |psi_0>, apply imaginary time evolution:
 *   |psi_n> = (l - H/Ns)^n |psi_0>
 * where l is a large constant (LargeValue) and Ns is the number of sites.
 *
 * Physical quantities:
 *   \f$\langle O \rangle_\beta \approx \langle\psi_n|O|\psi_n\rangle / \langle\psi_n|\psi_n\rangle\f$
 *
 * The inverse temperature beta at step n is related to norm:
 *   beta ≈ 2n / (Ns * l)
 *
 * Advantages:
 * - No explicit diagonalization needed
 * - Memory: O(N) vs O(N^2) for full diagonalization
 * - Naturally parallelizable
 *
 * Statistical sampling:
 * - Multiple random initial states (NumAve samples)
 * - Average over samples for better statistics
 *
 * @version 0.1, 0.2
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
#include "FirstMultiply.h"
#include "Multiply.h"
#include "expec_energy_flct.h"
#include "expec_cisajs.h"
#include "expec_cisajscktaltdc.h"
#include "nbody_correlation.h"
#include "anomalous_pair.h"
#include "green_output.h"
#include "restart_checkpoint.h"
#include "restart_manifest.h"
#include "tpq_restart.h"
#include "CalcByTPQ.h"
#include "FileIO.h"
#include "wrapperMPI.h"
#include "CalcTime.h"
#ifdef HPHI_TEST_SYMMETRY_PROBE
#include "symmetry_sector_probe.h"
#include "symmetry_fulldiag_probe.h"
#endif

#ifdef HPHI_TEST_SYMMETRY_PROBE
static void tpq_report_agreed_failure(const char *stage)
{
    fprintf(stdout, "mTPQ agreed failure stage=%s rank=%d\n", stage, myrank);
    fflush(stdout);
}
#define HPHI_TPQ_REPORT_DRIVER_FAILURE() \
    tpq_report_agreed_failure("driver")
#else
#define HPHI_TPQ_REPORT_DRIVER_FAILURE() do { } while (0)
#endif

static int tpq_output_open(const char *name, const char *mode, FILE **fp)
{
    int failed = childfopenMPI(name, mode, fp) != 0 || *fp == NULL;
#ifdef HPHI_TEST_SYMMETRY_PROBE
    if (SymmetryFullDiagTestFail("tpq_output_open")) failed = 1;
#endif
    if (SumMPI_i(failed) == 0) return 0;
    if (*fp != NULL) {
        fclose(*fp);
        *fp = NULL;
    }
#ifdef HPHI_TEST_SYMMETRY_PROBE
    tpq_report_agreed_failure("output_open");
#endif
    fprintf(stdoutMPI,
            "Error: mTPQ output open failed on one or more MPI ranks.\n");
    return -1;
}

static int tpq_output_close(FILE **fp)
{
    int failed = *fp == NULL;
    if (*fp != NULL) {
        if (ferror(*fp) != 0) failed = 1;
        if (fclose(*fp) != 0) failed = 1;
        *fp = NULL;
    }
#ifdef HPHI_TEST_SYMMETRY_PROBE
    if (SymmetryFullDiagTestFail("tpq_output_close")) failed = 1;
#endif
    if (SumMPI_i(failed) == 0) return 0;
#ifdef HPHI_TEST_SYMMETRY_PROBE
    tpq_report_agreed_failure("output_close");
#endif
    fprintf(stdoutMPI,
            "Error: mTPQ output close failed on one or more MPI ranks.\n");
    return -1;
}

static int tpq_expectation_status(int status)
{
    /* Expectation helpers have left any ExpecLocal session before returning.
       Keep the injected rank-local status at that same driver boundary. */
#ifdef HPHI_TEST_SYMMETRY_PROBE
    if (SymmetryFullDiagTestFail("tpq_expec_status")) status = -1;
#endif
    return status;
}

/**
 * @brief Main driver for TPQ calculation
 *
 * Performs NumAve independent TPQ calculations with different random
 * initial states, computing physical quantities at each step.
 *
 * For each random sample:
 * 1. Generate random initial state |psi_0>
 * 2. Repeat for step = 0 to Lanczos_max:
 *    a. Apply (l - H/Ns) to |psi>
 *    b. Every ExpecInterval steps, compute observables:
 *       - Energy \f$\langle H\rangle\f$, variance, inverse temperature
 *       - Green's functions if requested
 *    c. Normalize |psi> to prevent overflow
 *
 * Output files:
 * - SS_rand*.dat, or SS_tpq.dat in aggregate mode: Energy,
 *   \f$\langle S^2\rangle\f$, etc. vs step
 * - Norm_rand*.dat, or Norm_tpq.dat in aggregate mode: Norm vs step
 *   (for beta calculation)
 * - Flct_rand*.dat, or Flct_tpq.dat in aggregate mode: Fluctuations
 *
 * @param NumAve Number of random samples to average [in]
 * @param ExpecInterval Steps between observable calculations [in]
 * @param X Calculation parameters and results [in,out]
 *
 * @return 0 on success, -1 on error
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int CalcByTPQ(
      const int NumAve,
      const int ExpecInterval,
      struct EDMainCalStruct *X
)
{
  char sdt_phys[D_FileNameMax];
  char sdt_norm[D_FileNameMax];
  char sdt_flct[D_FileNameMax];
  int rand_i, rand_max, iret;
  int step_iO=0;
  FILE *fp = NULL;
  double inv_temp, Ns;
  struct TimeKeepStruct tstruct;
  int green_output_initialized = 0;
  int tpq_data_output_initialized = 0;
  int tpq_data_output_aggregate = 0;
  struct TPQRestartPlan restart_plan;

  tstruct.tstart=time(NULL);
  
  rand_max = NumAve;
  step_spin = ExpecInterval;
  X->Bind.Def.St=0;
  fprintf(stdoutMPI, "%s", cLogTPQ_Start);
  tpq_data_output_aggregate = GreenOutputUsesTPQDataAggregate(&(X->Bind));
  if (TPQRestartPlanSets(&(X->Bind), rand_max, &restart_plan) != 0) return -1;
  if (TPQRestartRecoverOutputs(&(X->Bind), &restart_plan,
                               tpq_data_output_aggregate) != 0) {
    TPQRestartFreePlan(&restart_plan);
    return -1;
  }
#define TPQ_RETURN_ERROR() do { TPQRestartFreePlan(&restart_plan); return -1; } while (0)
#define TPQ_AGREE_OR_RETURN(status) do { \
        int tpq_local_failed = (status) != 0; \
        if (SumMPI_i(tpq_local_failed) != 0) { \
            /* The probe marker lets MPI tests prove every rank crossed the \
               agreement boundary before the top-level abort. */ \
            HPHI_TPQ_REPORT_DRIVER_FAILURE(); \
            fprintf(stdoutMPI, \
                    "Error: mTPQ driver operation failed on one or more MPI ranks.\n"); \
            TPQ_RETURN_ERROR(); \
        } \
    } while (0)
  for (rand_i = 0; rand_i<rand_max; rand_i++){
    int resumed_noop = 0;
    if(tpq_data_output_aggregate){
      TPQ_AGREE_OR_RETURN(GreenOutputTPQDataFileName(&(X->Bind), GreenOutputTPQDataSS, sdt_phys));
      TPQ_AGREE_OR_RETURN(GreenOutputTPQDataFileName(&(X->Bind), GreenOutputTPQDataNorm, sdt_norm));
      TPQ_AGREE_OR_RETURN(GreenOutputTPQDataFileName(&(X->Bind), GreenOutputTPQDataFlct, sdt_flct));
    }else if(X->Bind.Def.iOutputDataHead==1){
      int prefix_length;
      prefix_length = sprintf(sdt_phys, "%s_", X->Bind.Def.CDataFileHead);
      sprintf(sdt_phys + prefix_length, cFileNameSSRand, rand_i);
      prefix_length = sprintf(sdt_norm, "%s_", X->Bind.Def.CDataFileHead);
      sprintf(sdt_norm + prefix_length, cFileNameNormRand, rand_i);
      prefix_length = sprintf(sdt_flct, "%s_", X->Bind.Def.CDataFileHead);
      sprintf(sdt_flct + prefix_length, cFileNameFlctRand, rand_i);
    }else{
      sprintf(sdt_phys, cFileNameSSRand, rand_i);
      sprintf(sdt_norm, cFileNameNormRand, rand_i);
      sprintf(sdt_flct, cFileNameFlctRand, rand_i);
    }
    Ns = 1.0 * X->Bind.Def.NsiteMPI;
    fprintf(stdoutMPI, cLogTPQRand, rand_i+1, rand_max);
    iret=0;
    X->Bind.Def.irand=rand_i;

  //Make or Read initial vector
    if (restart_plan.sources[rand_i] != RESTART_SOURCE_NONE) {
      struct RestartVectorInfo info;
      StartTimer(3600);
      TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep,
                               cOutputVecStart, "a", rand_i, step_i);
      fprintf(stdoutMPI, "%s", cLogInputVecStart);
      if (restart_plan.sources[rand_i] == RESTART_SOURCE_V3) {
        const struct RestartManifestRow *row =
            RestartManifestFind(&restart_plan.manifest, (uint64_t)rand_i);
        if (RestartReadVector(&(X->Bind), row, RESTART_COMP_STATE,
                              RESTART_READ_CONTINUATION, v0, &info) != 0) {
          StopTimer(3600);
          TPQ_RETURN_ERROR();
        }
        step_i = (int)info.step;
        global_1st_norm = info.scalar_a;
      } else {
        int source_ranks = X->Bind.Def.RestartSourceRanks > 0
            ? X->Bind.Def.RestartSourceRanks : nproc;
        if (RestartReadLegacyRaw(&(X->Bind), cFileNameInputVector, rand_i,
                                 source_ranks, v0, &step_i) != 0 ||
            SumMPI_i(step_i != restart_plan.legacy_steps[rand_i]) != 0) {
          fprintf(stdoutMPI, "Error: legacy recovery payload read failed.\n");
          StopTimer(3600);
          TPQ_RETURN_ERROR();
        }
        /* Legacy tmpvec does not store the first norm. */
        global_1st_norm = 0.0;
      }
      TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep,
                               cOutputVecFinish, "a", rand_i, step_i);
      fprintf(stdoutMPI, "%s", cLogInputVecFinish);
      StopTimer(3600);
      X->Bind.Def.istep = step_i;
      StartTimer(3200);
      iret = expec_energy_flct(&(X->Bind));
      StopTimer(3200);
      TPQ_AGREE_OR_RETURN(iret);
      step_iO = step_i - 1;
      resumed_noop = step_i >= X->Bind.Def.Lanczos_max;
      if (resumed_noop)
        fprintf(stdoutMPI,
                "  mTPQ: no remaining steps (set %d, saved step %d, "
                "Lanczos_max %d).\n", rand_i, step_i,
                X->Bind.Def.Lanczos_max);
    }

    if (restart_plan.sources[rand_i] == RESTART_SOURCE_NONE) {
      if (X->Bind.Def.iReStart == RESTART_INOUT)
        fprintf(stdoutMPI, "Start to calculate in normal procedure.\n");
      StartTimer(3600);
      if (!restart_plan.any_resume && green_output_initialized == 0) {
        TPQ_AGREE_OR_RETURN(GreenOutputInitializeAggregateFiles(&(X->Bind)));
        green_output_initialized = 1;
      }
      if (tpq_data_output_aggregate) {
        if (!restart_plan.any_resume && tpq_data_output_initialized == 0) {
          TPQ_AGREE_OR_RETURN(GreenOutputInitializeTPQDataAggregateFiles(&(X->Bind)));
          tpq_data_output_initialized = 1;
        }
      } else {
        if (tpq_output_open(sdt_phys, "w", &fp) != 0) {
          TPQ_RETURN_ERROR();
        }
        fprintf(fp, "%s", cLogSSRand);
        if (tpq_output_close(&fp) != 0) TPQ_RETURN_ERROR();
// for norm
        if (tpq_output_open(sdt_norm, "w", &fp) != 0) {
          TPQ_RETURN_ERROR();
        }
        fprintf(fp, "%s", cLogNormRand);
        if (tpq_output_close(&fp) != 0) TPQ_RETURN_ERROR();
// for fluctuations
        if (tpq_output_open(sdt_flct, "w", &fp) != 0) {
          TPQ_RETURN_ERROR();
        }
        fprintf(fp, "%s", cLogFlctRand);
        if (tpq_output_close(&fp) != 0) TPQ_RETURN_ERROR();
      }

      StopTimer(3600);

      step_i = 0;

      StartTimer(3100);
      if(rand_i==0){
        TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cTPQStep, "w", rand_i, step_i);
      }
      else{
        TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cTPQStep, "a", rand_i, step_i);
      }
      /**@brief
      Initialize v1 and compute v0 = H*v1
      */
      iret = FirstMultiply(rand_i, &(X->Bind));
      if (SumMPI_i(iret != 0) != 0) {
        StopTimer(3100);
        TPQ_RETURN_ERROR();
      }
#ifdef HPHI_TEST_SYMMETRY_PROBE
      if (SymmetryProbeWriteVector(&(X->Bind), v1, "initial", rand_i, 0,
                                   global_1st_norm) != 0) {
        StopTimer(3100);
        TPQ_RETURN_ERROR();
      }
#endif
      inv_temp = 0.0;
      StopTimer(3100);
      if (tpq_output_open(sdt_phys, "a", &fp) != 0) {
        TPQ_RETURN_ERROR();
      }
      GreenOutputWriteTPQSSRow(fp, &(X->Bind), step_i, inv_temp);
      if (tpq_output_close(&fp) != 0) TPQ_RETURN_ERROR();
      // for norm
      if (tpq_output_open(sdt_norm, "a", &fp) != 0) {
        TPQ_RETURN_ERROR();
      }
      GreenOutputWriteTPQNormRow(fp, &(X->Bind), step_i, inv_temp, global_1st_norm, global_1st_norm);
      if (tpq_output_close(&fp) != 0) TPQ_RETURN_ERROR();
      /**@brief
      Compute expectation value at infinite temperature
      */
      X->Bind.Def.istep = 0;
      StartTimer(3300);
      iret=tpq_expectation_status(expec_cisajs(&(X->Bind), v1));
      StopTimer(3300);
      TPQ_AGREE_OR_RETURN(iret);

      StartTimer(3400);
      iret=expec_cisajscktaltdc(&(X->Bind), v1);
      StopTimer(3400);
      TPQ_AGREE_OR_RETURN(iret);
      iret=expec_nbodyg(&(X->Bind), v1);
      TPQ_AGREE_OR_RETURN(iret);
      iret=expec_anomalousg(&(X->Bind), v1);
      TPQ_AGREE_OR_RETURN(iret);

      /** @brief Compute v1=0, and compute v0 = H*v1 */
      StartTimer(3200);
      iret=expec_energy_flct(&(X->Bind)); //v0 = H*v1
      StopTimer(3200);
      TPQ_AGREE_OR_RETURN(iret);
      step_i += 1;
      inv_temp = (2.0 / Ns) / (LargeValue - X->Bind.Phys.energy / Ns);
      StartTimer(3600);
      if (tpq_output_open(sdt_phys, "a", &fp) != 0) {
        TPQ_RETURN_ERROR();
      }
      GreenOutputWriteTPQSSRow(fp, &(X->Bind), step_i, inv_temp);
      if (tpq_output_close(&fp) != 0) TPQ_RETURN_ERROR();
// for norm
      if (tpq_output_open(sdt_norm, "a", &fp) != 0) {
        TPQ_RETURN_ERROR();
      }
      GreenOutputWriteTPQNormRow(fp, &(X->Bind), step_i, inv_temp, global_norm, global_1st_norm);
      if (tpq_output_close(&fp) != 0) TPQ_RETURN_ERROR();
// for fluctuations
      if (tpq_output_open(sdt_flct, "a", &fp) != 0) {
        TPQ_RETURN_ERROR();
      }
      GreenOutputWriteTPQFlctRow(fp, &(X->Bind), step_i, inv_temp);
      if (tpq_output_close(&fp) != 0) TPQ_RETURN_ERROR();
//
      StopTimer(3600);
      step_i +=1;
      X->Bind.Def.istep = step_i;
      step_iO=0;
    }

    int progress_interval = (X->Bind.Def.Lanczos_max - step_iO) / 10;
    if (progress_interval < 1) progress_interval = 1;  /* avoid % 0 for small (Lanczos_max - step_iO) */
    for (step_i = X->Bind.Def.istep; step_i<X->Bind.Def.Lanczos_max; step_i++){
      X->Bind.Def.istep=step_i;
      if(step_i % progress_interval == 0){
        fprintf(stdoutMPI, cLogTPQStep, step_i, X->Bind.Def.Lanczos_max);
      }
      X->Bind.Def.istep=step_i;
      StartTimer(3600);
      TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cTPQStep, "a", rand_i, step_i);
      StopTimer(3600);
      StartTimer(3500);
      iret = Multiply(&(X->Bind));
      StopTimer(3500);
      TPQ_AGREE_OR_RETURN(iret);

      StartTimer(3200);
      iret=expec_energy_flct(&(X->Bind));
      StopTimer(3200);
      TPQ_AGREE_OR_RETURN(iret);

//
      inv_temp = (2.0*step_i / Ns) / (LargeValue - X->Bind.Phys.energy / Ns);

      StartTimer(3600);
      if(tpq_output_open(sdt_phys, "a", &fp)!=0){
        TPQ_RETURN_ERROR();
      }
      GreenOutputWriteTPQSSRow(fp, &(X->Bind), step_i, inv_temp);
// for
      if (tpq_output_close(&fp) != 0) TPQ_RETURN_ERROR();

      if(tpq_output_open(sdt_norm, "a", &fp)!=0){
        TPQ_RETURN_ERROR();
      }
      GreenOutputWriteTPQNormRow(fp, &(X->Bind), step_i, inv_temp, global_norm, global_1st_norm);
      if (tpq_output_close(&fp) != 0) TPQ_RETURN_ERROR();

// for fluctuations
      if (tpq_output_open(sdt_flct, "a", &fp) != 0) {
        TPQ_RETURN_ERROR();
      }
      GreenOutputWriteTPQFlctRow(fp, &(X->Bind), step_i, inv_temp);
      if (tpq_output_close(&fp) != 0) TPQ_RETURN_ERROR();
//
      StopTimer(3600);


      if (step_i%step_spin == 0){
        StartTimer(3300);
        iret=tpq_expectation_status(expec_cisajs(&(X->Bind),v1));
        StopTimer(3300);
        TPQ_AGREE_OR_RETURN(iret);

        StartTimer(3400);
        iret=expec_cisajscktaltdc(&(X->Bind), v1);
        StopTimer(3400);
        TPQ_AGREE_OR_RETURN(iret);
        iret=expec_nbodyg(&(X->Bind), v1);
        TPQ_AGREE_OR_RETURN(iret);
        iret=expec_anomalousg(&(X->Bind), v1);
        TPQ_AGREE_OR_RETURN(iret);
      }
      if ((X->Bind.Def.iReStart == RESTART_OUT ||
           X->Bind.Def.iReStart == RESTART_INOUT) &&
          X->Bind.Def.Param.CheckpointInterval > 0 &&
          step_i % X->Bind.Def.Param.CheckpointInterval == 0)
        if (TPQRestartWriteCheckpoint(&(X->Bind), rand_i, step_i + 1,
                                      global_1st_norm, 0.0, 0.0) != 0)
          TPQ_RETURN_ERROR();
    }

    if (!resumed_noop &&
        (X->Bind.Def.iReStart == RESTART_OUT ||
         X->Bind.Def.iReStart == RESTART_INOUT)) {
      TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecStart, "a", rand_i, step_i);
      fprintf(stdoutMPI, "%s", cLogOutputVecStart);
      if (TPQRestartWriteCheckpoint(&(X->Bind), rand_i, step_i,
                                    global_1st_norm, 0.0, 0.0) != 0)
        TPQ_RETURN_ERROR();
      TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecFinish, "a", rand_i, step_i);
      fprintf(stdoutMPI, "%s", cLogOutputVecFinish);
    }
  }
  fprintf(stdoutMPI, "%s", cLogTPQ_End);

  tstruct.tend=time(NULL);
  fprintf(stdoutMPI, cLogTPQEnd, (int)(tstruct.tend-tstruct.tstart));
  TPQRestartFreePlan(&restart_plan);
#undef TPQ_AGREE_OR_RETURN
#undef TPQ_RETURN_ERROR
  return TRUE;
}
