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
#include "CalcByTPQ.h"
#include "FileIO.h"
#include "wrapperMPI.h"
#include "CalcTime.h"

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
  char sdt[D_FileNameMax];
  char sdt_phys[D_FileNameMax];
  char sdt_norm[D_FileNameMax];
  char sdt_flct[D_FileNameMax];
  int rand_i, rand_max, iret;
  unsigned long int i_max;
  int step_iO=0;
  FILE *fp;
  double inv_temp, Ns;
  struct TimeKeepStruct tstruct;
  size_t byte_size;
  int green_output_initialized = 0;
  int tpq_data_output_initialized = 0;
  int tpq_data_output_aggregate = 0;

  tstruct.tstart=time(NULL);
  
  rand_max = NumAve;
  step_spin = ExpecInterval;
  X->Bind.Def.St=0;
  fprintf(stdoutMPI, "%s", cLogTPQ_Start);
  tpq_data_output_aggregate = GreenOutputUsesTPQDataAggregate(&(X->Bind));
  for (rand_i = 0; rand_i<rand_max; rand_i++){
    if(tpq_data_output_aggregate){
      if (GreenOutputTPQDataFileName(&(X->Bind), GreenOutputTPQDataSS, sdt_phys) != 0) {
        return -1;
      }
      if (GreenOutputTPQDataFileName(&(X->Bind), GreenOutputTPQDataNorm, sdt_norm) != 0) {
        return -1;
      }
      if (GreenOutputTPQDataFileName(&(X->Bind), GreenOutputTPQDataFlct, sdt_flct) != 0) {
        return -1;
      }
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
    if(X->Bind.Def.iReStart==RESTART_INOUT || X->Bind.Def.iReStart==RESTART_IN) {
      StartTimer(3600);
      TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecStart, "a", rand_i, step_i);
      fprintf(stdoutMPI, "%s", cLogInputVecStart);
      sprintf(sdt, cFileNameInputVector, rand_i, myrank);
      childfopenALL(sdt, "rb", &fp);
      if(fp==NULL){
        fprintf(stderr, "A file of Inputvector does not exist (rank %d).\n", myrank);
        iret=1;
      }
      /* All ranks must agree on the fallback: tmpvec files are per-rank and a
         partially missing set would otherwise split ranks across the restart
         and fresh-start branches, whose MPI collectives do not match. */
      iret = (int)MaxMPI_li((unsigned long int)iret);
      if(iret==1){
        if(fp != NULL) fclose(fp);
        fprintf(stdoutMPI, "Start to calculate in normal procedure.\n");
        StopTimer(3600);
      }else{
        byte_size = fread(&step_i, sizeof(step_i), 1, fp);
        byte_size = fread(&i_max, sizeof(long int), 1, fp);
        if(i_max != X->Bind.Check.idim_max){
          fprintf(stderr, "Error: A file of Inputvector is incorrect.\n");
          exitMPI(-1);
        }
        byte_size = fread(v0, sizeof(complex double), X->Bind.Check.idim_max+1, fp);
        TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecFinish, "a", rand_i, step_i);
        fprintf(stdoutMPI, "%s", cLogInputVecFinish);
        fclose(fp);
        StopTimer(3600);
        X->Bind.Def.istep=step_i;
        StartTimer(3200);
        iret=expec_energy_flct(&(X->Bind));
        StopTimer(3200);
        if(iret != 0) return -1;

        step_iO=step_i-1;
        if (byte_size == 0) printf("byte_size: %d \n", (int)byte_size);
      }
    }
    
    if(X->Bind.Def.iReStart==RESTART_NOT || X->Bind.Def.iReStart==RESTART_OUT || iret ==1) {
      StartTimer(3600);
      if (green_output_initialized == 0) {
        if (GreenOutputInitializeAggregateFiles(&(X->Bind)) != 0) {
          return -1;
        }
        green_output_initialized = 1;
      }
      if (tpq_data_output_aggregate) {
        if (tpq_data_output_initialized == 0) {
          if (GreenOutputInitializeTPQDataAggregateFiles(&(X->Bind)) != 0) {
            return -1;
          }
          tpq_data_output_initialized = 1;
        }
      } else {
        if (childfopenMPI(sdt_phys, "w", &fp) != 0) {
          return -1;
        }
        fprintf(fp, "%s", cLogSSRand);
        fclose(fp);
// for norm
        if (childfopenMPI(sdt_norm, "w", &fp) != 0) {
          return -1;
        }
        fprintf(fp, "%s", cLogNormRand);
        fclose(fp);
// for fluctuations
        if (childfopenMPI(sdt_flct, "w", &fp) != 0) {
          return -1;
        }
        fprintf(fp, "%s", cLogFlctRand);
        fclose(fp);
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
      FirstMultiply(rand_i, &(X->Bind));
      inv_temp = 0.0;
      StopTimer(3100);
      if (childfopenMPI(sdt_phys, "a", &fp) != 0) {
        return -1;
      }
      GreenOutputWriteTPQSSRow(fp, &(X->Bind), step_i, inv_temp);
      fclose(fp);
      // for norm
      if (childfopenMPI(sdt_norm, "a", &fp) != 0) {
        return -1;
      }
      GreenOutputWriteTPQNormRow(fp, &(X->Bind), step_i, inv_temp, global_1st_norm, global_1st_norm);
      fclose(fp);
      /**@brief
      Compute expectation value at infinite temperature
      */
      X->Bind.Def.istep = 0;
      StartTimer(3300);
      iret=expec_cisajs(&(X->Bind), v1);
      StopTimer(3300);
      if(iret !=0) return -1;

      StartTimer(3400);
      iret=expec_cisajscktaltdc(&(X->Bind), v1);
      StopTimer(3400);
      if(iret !=0) return -1;
      iret=expec_nbodyg(&(X->Bind), v1);
      if(iret !=0) return -1;
      iret=expec_anomalousg(&(X->Bind), v1);
      if(iret !=0) return -1;

      /** @brief Compute v1=0, and compute v0 = H*v1 */
      StartTimer(3200);
      iret=expec_energy_flct(&(X->Bind)); //v0 = H*v1
      StopTimer(3200);
      if(iret !=0) return -1;
      step_i += 1;
      inv_temp = (2.0 / Ns) / (LargeValue - X->Bind.Phys.energy / Ns);
      StartTimer(3600);
      if (childfopenMPI(sdt_phys, "a", &fp) != 0) {
        return -1;
      }
      GreenOutputWriteTPQSSRow(fp, &(X->Bind), step_i, inv_temp);
      fclose(fp);
// for norm
      if (childfopenMPI(sdt_norm, "a", &fp) != 0) {
        return -1;
      }
      GreenOutputWriteTPQNormRow(fp, &(X->Bind), step_i, inv_temp, global_norm, global_1st_norm);
      fclose(fp);
// for fluctuations
      if (childfopenMPI(sdt_flct, "a", &fp) != 0) {
        return -1;
      }
      GreenOutputWriteTPQFlctRow(fp, &(X->Bind), step_i, inv_temp);
      fclose(fp);
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
      Multiply(&(X->Bind));
      StopTimer(3500);

      StartTimer(3200);
      iret=expec_energy_flct(&(X->Bind));
      StopTimer(3200);
      if(iret !=0) return -1;

//
      inv_temp = (2.0*step_i / Ns) / (LargeValue - X->Bind.Phys.energy / Ns);

      StartTimer(3600);
      if(childfopenMPI(sdt_phys, "a", &fp)!=0){
        return FALSE;
      }
      GreenOutputWriteTPQSSRow(fp, &(X->Bind), step_i, inv_temp);
// for
      fclose(fp);

      if(childfopenMPI(sdt_norm, "a", &fp)!=0){
        return FALSE;
      }
      GreenOutputWriteTPQNormRow(fp, &(X->Bind), step_i, inv_temp, global_norm, global_1st_norm);
      fclose(fp);

// for fluctuations
      if (childfopenMPI(sdt_flct, "a", &fp) != 0) {
        return -1;
      }
      GreenOutputWriteTPQFlctRow(fp, &(X->Bind), step_i, inv_temp);
      fclose(fp);
//
      StopTimer(3600);


      if (step_i%step_spin == 0){
        StartTimer(3300);
        iret=expec_cisajs(&(X->Bind),v1);
        StopTimer(3300);
        if(iret !=0) return -1;

        StartTimer(3400);
        iret=expec_cisajscktaltdc(&(X->Bind), v1);
        StopTimer(3400);
        if(iret !=0) return -1;
        iret=expec_nbodyg(&(X->Bind), v1);
        if(iret !=0) return -1;
        iret=expec_anomalousg(&(X->Bind), v1);
        if(iret !=0) return -1;
      }
    }

    if(X->Bind.Def.iReStart== RESTART_OUT || X->Bind.Def.iReStart==RESTART_INOUT){
      TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecStart, "a", rand_i, step_i);
      fprintf(stdoutMPI, "%s", cLogOutputVecStart);
      sprintf(sdt, cFileNameOutputVector, rand_i, myrank);
      if(childfopenALL(sdt, "wb", &fp)!=0){
        exitMPI(-1);
      }
      fwrite(&step_i, sizeof(step_i), 1, fp);
      fwrite(&X->Bind.Check.idim_max, sizeof(X->Bind.Check.idim_max),1,fp);
      fwrite(v1, sizeof(complex double),X->Bind.Check.idim_max+1, fp);
      fclose(fp);
      TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecFinish, "a", rand_i, step_i);
      fprintf(stdoutMPI, "%s", cLogOutputVecFinish);
    }
  }
  fprintf(stdoutMPI, "%s", cLogTPQ_End);

  tstruct.tend=time(NULL);
  fprintf(stdoutMPI, cLogTPQEnd, (int)(tstruct.tend-tstruct.tstart));
  return TRUE;
}
