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
/*-------------------------------------------------------------*/
/*-------------------------------------------------------------
 * HPhi
 * timer program
 * "-lrt" option is needed for clock_gettime().
 *-------------------------------------------------------------
 * original code written by Satoshi Morita
 *-------------------------------------------------------------*/

#include "Common.h"
#include "FileIO.h"
#include "CalcTime.h"
#include "symmetry_basis.h"
#include "symmetry_matvec_plan.h"

#ifdef MPI
#include <mpi.h>
#endif

#ifdef MPI
static unsigned long long HashSymmetryBytes(unsigned long long hash,
                                            const void *data,
                                            size_t size)
{
  const unsigned char *bytes = (const unsigned char *)data;
  size_t i;
  for (i = 0U; i < size; i++) {
    hash ^= (unsigned long long)bytes[i];
    hash *= 1099511628211ULL;
  }
  return hash;
}

static unsigned long long SymmetryBasisDigest(
    const struct SymmetryBasisRuntime *sym)
{
  unsigned long long hash = 14695981039346656037ULL;
  unsigned long int index;
  hash = HashSymmetryBytes(hash, &sym->dim, sizeof(sym->dim));
  for (index = 1UL; index <= sym->dim; index++) {
    const struct SymmetryBasisVector *entry = &sym->basis[index];
    hash = HashSymmetryBytes(hash, &entry->rep_state,
                             sizeof(entry->rep_state));
    hash = HashSymmetryBytes(hash, &entry->orbit_size,
                             sizeof(entry->orbit_size));
    hash = HashSymmetryBytes(hash, &entry->stabilizer_size,
                             sizeof(entry->stabilizer_size));
    hash = HashSymmetryBytes(hash, &entry->norm, sizeof(entry->norm));
    hash = HashSymmetryBytes(hash, &entry->stabilizer_character_sum,
                             sizeof(entry->stabilizer_character_sum));
    hash = HashSymmetryBytes(hash, &entry->diagonal,
                             sizeof(entry->diagonal));
  }
  return hash;
}

static void OutputSymmetryRankStats(const struct BindStruct *X)
{
  static const int timer_ids[] = {
    1100, 1110, 1115, 1111, 1112, 1113, 1114,
    1101, 1120, 1121, 1122, 4113
  };
  static const char *work_keys[] = {
    "basis_raw_states",
    "basis_representative_candidates",
    "basis_compatible_survivors",
    "basis_transform_calls",
    "basis_orbit_metadata_calls",
    "basis_thread_count",
    "basis_thread_raw_states_max",
    "basis_thread_representative_candidates_max",
    "basis_thread_compatible_survivors_max",
    "basis_thread_transform_calls_max",
    "basis_gather_entries",
    "basis_gather_bytes",
    "plan_local_rows",
    "plan_local_nnz",
    "plan_row_nnz_max"
  };
  const size_t timer_count = sizeof(timer_ids) / sizeof(timer_ids[0]);
  const size_t work_count = sizeof(work_keys) / sizeof(work_keys[0]);
  const struct SymmetryMatvecPlan *plan;
  double timer_local[sizeof(timer_ids) / sizeof(timer_ids[0])];
  double timer_min[sizeof(timer_ids) / sizeof(timer_ids[0])];
  double timer_max[sizeof(timer_ids) / sizeof(timer_ids[0])];
  double timer_sum[sizeof(timer_ids) / sizeof(timer_ids[0])];
  unsigned long long work_local[sizeof(work_keys) / sizeof(work_keys[0])];
  unsigned long long work_min[sizeof(work_keys) / sizeof(work_keys[0])];
  unsigned long long work_max[sizeof(work_keys) / sizeof(work_keys[0])];
  unsigned long long work_sum[sizeof(work_keys) / sizeof(work_keys[0])];
  double row_mean_local;
  double row_mean_min;
  double row_mean_max;
  double row_mean_sum;
  unsigned long long basis_digest_local;
  unsigned long long basis_digest_min;
  unsigned long long basis_digest_max;
  char fileName[D_FileNameMax];
  FILE *fp;
  size_t i;

  if (X == NULL || X->Def.iFlgSymmetryBasis == FALSE || X->Sym == NULL ||
      X->Sym->enabled != TRUE) {
    return;
  }
  plan = X->Sym->matvec_plan;
  for (i = 0; i < timer_count; i++) timer_local[i] = Timer[timer_ids[i]];
  work_local[0] = X->Sym->basis_raw_states;
  work_local[1] = X->Sym->basis_representative_candidates;
  work_local[2] = X->Sym->basis_compatible_survivors;
  work_local[3] = X->Sym->basis_transform_calls;
  work_local[4] = X->Sym->basis_orbit_metadata_calls;
  work_local[5] = (unsigned long long)X->Sym->basis_thread_count;
  work_local[6] = X->Sym->basis_thread_raw_states_max;
  work_local[7] = X->Sym->basis_thread_representative_candidates_max;
  work_local[8] = X->Sym->basis_thread_compatible_survivors_max;
  work_local[9] = X->Sym->basis_thread_transform_calls_max;
  work_local[10] = X->Sym->basis_gather_entries;
  work_local[11] = X->Sym->basis_gather_bytes;
  work_local[12] = plan != NULL ? (unsigned long long)plan->local_dim : 0ULL;
  work_local[13] = plan != NULL ? (unsigned long long)plan->nnz : 0ULL;
  work_local[14] = plan != NULL ? (unsigned long long)plan->row_nnz_max : 0ULL;
  row_mean_local = plan != NULL && plan->local_dim > 0UL
                       ? (double)plan->nnz / (double)plan->local_dim
                       : 0.0;
  basis_digest_local = SymmetryBasisDigest(X->Sym);

  MPI_Allreduce(timer_local, timer_min, (int)timer_count, MPI_DOUBLE,
                MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(timer_local, timer_max, (int)timer_count, MPI_DOUBLE,
                MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(timer_local, timer_sum, (int)timer_count, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(work_local, work_min, (int)work_count, MPI_UNSIGNED_LONG_LONG,
                MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(work_local, work_max, (int)work_count, MPI_UNSIGNED_LONG_LONG,
                MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(work_local, work_sum, (int)work_count, MPI_UNSIGNED_LONG_LONG,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&row_mean_local, &row_mean_min, 1, MPI_DOUBLE,
                MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&row_mean_local, &row_mean_max, 1, MPI_DOUBLE,
                MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&row_mean_local, &row_mean_sum, 1, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&basis_digest_local, &basis_digest_min, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&basis_digest_local, &basis_digest_max, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);

  sprintf(fileName, "CalcTimerRankStats.dat");
  if (childfopenMPI(fileName, "w", &fp) != 0) return;
  fprintf(fp, "format=HPhiCalcTimerRankStats version=1 ranks=%d\n", nproc);
  for (i = 0; i < timer_count; i++) {
    fprintf(fp, "timer id=%d ranks=%d min=%.17g max=%.17g mean=%.17g\n",
            timer_ids[i], nproc, timer_min[i], timer_max[i],
            timer_sum[i] / (double)nproc);
  }
  for (i = 0; i < work_count; i++) {
    fprintf(fp, "work key=%s ranks=%d min=%llu max=%llu mean=%.17g\n",
            work_keys[i], nproc, work_min[i], work_max[i],
            (double)work_sum[i] / (double)nproc);
  }
  fprintf(fp,
          "work key=plan_row_nnz_mean ranks=%d min=%.17g max=%.17g mean=%.17g\n",
          nproc, row_mean_min, row_mean_max, row_mean_sum / (double)nproc);
  fprintf(fp,
          "basis_digest algorithm=fnv1a64-fields ranks=%d min=%016llx max=%016llx\n",
          nproc, basis_digest_min, basis_digest_max);
  fclose(fp);
}
#endif
/** 
 * 
 * @brief function for displaying elapse time
 * 
 * @version 2.0
 */

void StampTime(FILE *fp, char *str, int num){
#ifdef MPI
  char str1[256];
  sprintf(str1, "%-50s [%04d] %12.5lf\n", str, num, Timer[num]);
  fprintf(fp, "%s", str1);
#endif
}

/** 
 * 
 * @brief function for initializing Timer[]
 * 
 * @version 2.0
 */
void InitTimer() {
#ifdef MPI
  int i;
  int NTimer=10000;
  Timer       = (double*)malloc((NTimer)*sizeof(double));
  TimerStart  = (double*)malloc((NTimer)*sizeof(double));
  for(i=0;i<NTimer;i++) Timer[i]=0.0;
  for(i=0;i<NTimer;i++) TimerStart[i]=0.0;
#endif
  return;
}
/** 
 * 
 * @brief function for initializing elapse time [start]
 * 
 * @version 2.0
 */

void StartTimer(int n) {
#ifdef MPI
  TimerStart[n]=MPI_Wtime();
#endif
  return;
}
/** 
 * 
 * @brief function for calculating elapse time [elapse time=StartTimer-StopTimer]
 * 
 * @version 2.0
 */
void StopTimer(int n) {
#ifdef MPI
  Timer[n] += MPI_Wtime() - TimerStart[n];
#endif
  return;
}
/** 
 * 
 * @brief function for outputting elapse time for each function
 * 
 * @version 2.0
 */
void OutputTimer(struct BindStruct *X) {

#ifdef MPI
  char fileName[D_FileNameMax];
  FILE *fp;
  sprintf(fileName, "CalcTimer.dat"); //TBC
  childfopenMPI(fileName,"w", &fp);
  //fp = fopen(fileName, "w");
  //fp = childfopenMPI(fileName, "w");
  StampTime(fp, "All", 0);
  StampTime(fp, "  sz", 1000);
  StampTime(fp, "  symmetry basis build/activate", 1100);
  StampTime(fp, "    symmetry basis raw enumeration", 1110);
  StampTime(fp, "      symmetry basis MPI gather/reduction", 1115);
  StampTime(fp, "    symmetry basis sort/merge", 1111);
  StampTime(fp, "    symmetry representative hash build", 1112);
  StampTime(fp, "    symmetry diagonal materialization", 1113);
  StampTime(fp, "    symmetry dimension activation/validation", 1114);
  StampTime(fp, "  symmetry matvec plan build", 1101);
  StampTime(fp, "    symmetry plan count/prefix", 1120);
  StampTime(fp, "    symmetry plan storage allocation", 1121);
  StampTime(fp, "    symmetry plan fill", 1122);
  StampTime(fp, "  diagonalcalc", 2000);
  if(X->Def.iFlgCalcSpec == CALCSPEC_NOT){
    if(X->Def.iCalcType==TPQCalc || X->Def.iCalcType==cTPQ) {
      StampTime(fp, "  CalcByTPQ", 3000);
      StampTime(fp, "    FirstMultiply", 3100);
      StampTime(fp, "      rand   in FirstMultiply", 3101);
      StampTime(fp, "      mltply in FirstMultiply", 3102);
      StampTime(fp, "    expec_energy_flct        ", 3200);
      StampTime(fp, "      calc flctuation in expec_energy_flct ", 3201);
      StampTime(fp, "      mltply in expec_energy_flct ", 3202);
      StampTime(fp, "    expec_onebody            ", 3300);
      StampTime(fp, "    expec_twobody            ", 3400);
      StampTime(fp, "    Multiply                 ", 3500);
      StampTime(fp, "    FileIO                   ", 3600);
    }
    else if(X->Def.iCalcType==Lanczos){
      StampTime(fp, "  CalcByLanczos", 4000);
      StampTime(fp, "    LanczosEigenValue", 4100);
      StampTime(fp, "      mltply      in LanczosEigenValue", 4101);
      StampTime(fp, "      vec12       in LanczosEigenValue", 4102);
      StampTime(fp, "      DSEVvalue   in LanczosEigenValue", 4103);
      StampTime(fp, "      initial vector zero fill", 4110);
      StampTime(fp, "      initial vector random fill", 4111);
      StampTime(fp, "      initial vector local norm", 4112);
      StampTime(fp, "      initial vector MPI reduction", 4113);
      StampTime(fp, "      initial vector normalization", 4114);
      StampTime(fp, "    LanczosEigenVector", 4200);
      StampTime(fp, "      mltply      in LanczosEigenVector", 4201);
      StampTime(fp, "    expec_energy_flct", 4300);
      StampTime(fp, "      calc flctuation in expec_energy_flct ", 4301);
      StampTime(fp, "      mltply in expec_energy_flct ", 4302);
      StampTime(fp, "    CGEigenVector", 4400);
      StampTime(fp, "      mltply in CGEigenVector ", 4401);
      StampTime(fp, "    expec_onebody            ", 4500);
      StampTime(fp, "    expec_twobody            ", 4600);
      StampTime(fp, "    expec_TotalSz            ", 4700);
      StampTime(fp, "    FileIO                   ", 4800);
      StampTime(fp, "      Read Input Eigenvec ", 4801);    
    }
    else if(X->Def.iCalcType==FullDiag){
      StampTime(fp, "  CalcByFullDiag", 5000);
      StampTime(fp, "    MakeHam", 5100);
      StampTime(fp, "    LapackDiag", 5200);
      StampTime(fp, "    CalcPhys", 5300);
    StampTime(fp, "      calc flctuation in expec_energy_flct ", 5301);
    StampTime(fp, "      mltply in expec_energy_flct ", 5302);
        StampTime(fp, "    Output", 5400);
      StampTime(fp, "    OutputHam", 5500);
    }
  }
  else{ 
    StampTime(fp, "  CalcSpectrum by Lanczos method", 6000);
    StampTime(fp, "    Make excited state", 6100);
    StampTime(fp, "      Read origin state", 6101);
    StampTime(fp, "      Multiply excited operator", 6102);
    StampTime(fp, "    Calculate spectrum", 6200);
    if(X->Def.iCalcType==Lanczos){
      StampTime(fp, "      Read vector for recalculation", 6201);
      StampTime(fp, "      Read tridiagonal components for recalculation", 6202);
      StampTime(fp, "      Calculate tridiagonal components", 6203);
      StampTime(fp, "      Output tridiagonal components", 6204);
      StampTime(fp, "      Calculate spectrum by Lanczos method", 6205);
      StampTime(fp, "      Output vectors for recalculation", 6206);
    }
    else if(X->Def.iCalcType==FullDiag){
      StampTime(fp, "      MakeHam", 6301);
      StampTime(fp, "      lapackdiag", 6302);
      StampTime(fp, "      Calculate v1", 6303);
      StampTime(fp, "      Calculate spectrum", 6304);
    }
  }
  
  fprintf(fp,"================================================\n");
  
  StampTime(fp,"All mltply",1);
  StampTime(fp,"  symmetry input Allgatherv",1501);
  StampTime(fp,"  symmetry legacy beta scan",1502);
  StampTime(fp,"  symmetry local-row plan apply",1503);
  StampTime(fp,"  diagonal", 100);

  switch(X->Def.iCalcModel){
  case HubbardGC:
    StampTime(fp,"  HubbardGC", 200);
    StampTime(fp,"    trans    in HubbardGC", 210);
    StampTime(fp,"      double", 211);
    StampTime(fp,"      single", 212);
    StampTime(fp,"      inner", 213);
    StampTime(fp,"    interall in HubbardGC", 220);
    StampTime(fp,"      interPE", 221);
    StampTime(fp,"      inner", 222);
    StampTime(fp,"    pairhopp in HubbardGC", 230);
    StampTime(fp,"      interPE", 231);
    StampTime(fp,"      inner", 232);
    StampTime(fp,"    exchange in HubbardGC", 240);
    StampTime(fp,"      interPE", 241);
    StampTime(fp,"      inner", 242);
    break;
    
  case Hubbard:
  case tJ:
  case tJGC:
    StampTime(fp,"  Hubbard", 300);
    StampTime(fp,"    trans    in Hubbard", 310);
    StampTime(fp,"      double", 311);
    StampTime(fp,"      single", 312);
    StampTime(fp,"      inner", 313);
    StampTime(fp,"    interall in Hubbard", 320);
    StampTime(fp,"      interPE", 321);
    StampTime(fp,"      inner", 322);
    StampTime(fp,"    pairhopp in Hubbard", 330);
    StampTime(fp,"      interPE", 331);
    StampTime(fp,"      inner", 332);
    StampTime(fp,"    exchange in Hubbard", 340);
    StampTime(fp,"      interPE", 341);
    StampTime(fp,"      inner", 342);
    break;
    
  case Spin:
    fprintf(fp,"\n");
    StampTime(fp,"  Spin", 400);
    StampTime(fp,"    interall in Spin", 410);
    StampTime(fp,"      double", 411);
    StampTime(fp,"      single1", 412);
    StampTime(fp,"      single2", 413);
    StampTime(fp,"      inner", 414);
    StampTime(fp,"    exchange in Spin", 420);
    StampTime(fp,"      double", 421);
    StampTime(fp,"      single1", 422);
    StampTime(fp,"      single2", 423);
    StampTime(fp,"      inner", 424);
    break;
    
  case SpinGC:
    StampTime(fp,"  SpinGC", 500);
    StampTime(fp,"    trans    in SpinGC", 510);
    StampTime(fp,"      double", 511);
    StampTime(fp,"      inner", 512);
    StampTime(fp,"    interall in SpinGC", 520);
    StampTime(fp,"      double", 521);
    StampTime(fp,"      single", 522);
    StampTime(fp,"      inner", 523);
    StampTime(fp,"    exchange in SpinGC", 530);
    StampTime(fp,"      double", 531);
    StampTime(fp,"      single", 532);
    StampTime(fp,"      inner", 533);
    StampTime(fp,"    pairlift in SpinGC", 540);
    StampTime(fp,"      double", 541);
    StampTime(fp,"      single", 542);
    StampTime(fp,"      inner", 543);
    break;

  default:
    break;
  }
  fprintf(fp,"================================================\n");

  fclose(fp);
  OutputSymmetryRankStats(X);
  free(Timer);
  free(TimerStart);
#endif
  return;
}

/**
@page page_time Compute elapsed time for new functions

 Using StartTimer and StopTimer functions defined in time.c, we can measure the elapsed time for computation.

 1. Define an index and an output message by using StampTime function in time.c.
 For example, the index and the output message for the elapsed time of TPQ calculation is defined as follows.
 ```
       StampTime(fp, "  CalcByTPQ", 3000);
 ```

 2. Include CalcTime.h in the target source file.

 3. Set StartTimer and StopTimer functions in the region where you want to measure the time.
    It is noted that both functions must have the same index defined in time.c.
 For example, the elapsed time of TPQ calculation can be measured as follows.
 ```
       case TPQCalc:
        StartTimer(3000);
        if (CalcByTPQ(NumAve, X.Bind.Def.Param.ExpecInterval, &X) != TRUE) {
          FinalizeMPI();
          StopTimer(3000);
          return 0;
        }
        StopTimer(3000);
      break;

 ```

When above procedures were done, after calculation, you can see the elapsed time in CalcTimer.dat file as follows. The time unit is second.

```
 All                                                [0000]     37.94046
  sz                                               [1000]      0.00058
  diagonalcalc                                     [2000]      0.00046
  CalcByTPQ                                        [3000]     37.93129
```


*/
