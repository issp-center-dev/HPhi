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
#include "restart_checkpoint.h"
#include "restart_manifest.h"
#include "restart_output.h"
#ifdef HPHI_TEST_SYMMETRY_PROBE
#include "symmetry_sector_probe.h"
#endif

void MakeTEDTransfer(struct BindStruct *X, const int timeidx);
void MakeTEDInterAll(struct BindStruct *X, const int timeidx);

/* Driver helpers run outside ExpecLocal; agree even raw rank-local failures. */
#define TE_AGREE_OR_RETURN(call) do { int te_failed = (call) != 0; \
  if (SumMPI_i(te_failed) != 0) return -1; } while (0)

static int te_open(const char *name, const char *mode, FILE **fp)
{
  *fp = NULL;
  int failed = childfopenMPI(name, mode, fp) != 0;
  if (SumMPI_i(failed) != 0) {
    if (*fp) fclose(*fp);
    fprintf(stdoutMPI, "Error: failed to open TE output %s.\n", name);
    return -1;
  }
  return 0;
}

static double te_time_label(const struct BindStruct *X,
                            const struct SymmetryTEHamiltonian *view, int step, double dt)
{
  return view ? SymmetryTETime(&view->base, step - 1) :
         X->Def.NLaser ? X->Def.Param.Tinit + step * dt : X->Def.TETime[step - 1];
}

static int te_checkpoint_write(struct BindStruct *X, int next_step, double time)
{
  const uint64_t zero = 0;
  const struct RestartUnit unit = {1, &zero, 1};
  const struct RestartVectorInfo info = {TimeEvolution, 0, next_step, time, 1, 0, 0, 0, 0};
  double complex *const components[] = {v1};
  double complex *const *const vectors[] = {components};
  return RestartWriteUnit(X, &unit, vectors, &info, 0);
}

static int te_restore(struct BindStruct *X, const struct SymmetryTEHamiltonian *view,
                      double dt, int *step_initial, double *restored_time,
                      int *is_v3, int *legacy_step)
{
  struct RestartManifest manifest;
  enum RestartSource source;
  int exists;
  *is_v3 = 0; *legacy_step = 0;
  if (X->Def.iReStart != RESTART_INOUT && X->Def.iReStart != RESTART_IN) return 0;
  if (RestartManifestLoad(X, &manifest, &exists) != 0) return -1;
  /* Only the manifest selects continuation here. SpectrumVec is classified
   * by the import dispatcher below; its raw sentinel preserves legacy STEP. */
  if (RestartClassifySet(X, exists ? &manifest : NULL, 0, NULL,
                         &source, legacy_step) != 0) {
    RestartManifestFree(&manifest); return -1;
  }
  if (source != RESTART_SOURCE_V3) { RestartManifestFree(&manifest); return 0; }
  struct RestartVectorInfo info;
  int result = RestartReadVector(X, RestartManifestFind(&manifest, 0), RESTART_COMP_STATE,
                                RESTART_READ_CONTINUATION, v1, &info);
  RestartManifestFree(&manifest);
  if (result != 0) return -1;
  int step = (int)info.step; /* Reader validated 1 <= STEP <= INT_MAX-1. */
  int uses_table = !(view ? view->base.NLaser : X->Def.NLaser);
  if (SumMPI_i(uses_table && (unsigned int)step > X->Def.NTETimeSteps) != 0) {
    fprintf(stdoutMPI, "Error: time table shorter than the saved checkpoint (needs %d entries).\n", step);
    return -1;
  }
  double expected = te_time_label(X, view, step, dt);
  double tolerance = fmax(1e-15, 1e-12 * fabs(expected));
  if (SumMPI_i(!isfinite(expected) || fabs(info.time - expected) > tolerance) != 0) {
    fprintf(stdoutMPI, "Error: time grid changed: checkpoint time %.17g, grid %.17g at step %d.\n",
            info.time, expected, step);
    return -1;
  }
  *step_initial = step; *restored_time = info.time; *is_v3 = 1;
  return 0;
}

static int te_recover_green(struct BindStruct *X, const struct RestartSeriesRule *rule)
{
  const unsigned int counts[] = {X->Def.NCisAjt, X->Def.NCisAjtCkuAlvDC,
    X->Def.NTBody, X->Def.NFBody, X->Def.NSBody, X->Def.NNBodyG, X->Def.NAnomalousG};
  const struct RestartSeriesSpec spec = {"#", -1, 0, 1, rule, 0};
  for (int kind = GreenOutputOneBody; kind <= GreenOutputAnomalous; ++kind) {
    if (!counts[kind] || !GreenOutputKindUsesAggregate(X, (GreenOutputKind)kind)) continue;
    char name[D_FileNameMax];
    TE_AGREE_OR_RETURN(GreenOutputFileName(X, (GreenOutputKind)kind, name));
    if (RestartRecoverTimeSeries(X, name, &spec) != 0) return -1;
  }
  return 0;
}

/* Record the TE schedule and, for continuation, its restored source clock. */
static int record_sector_te(const struct BindStruct *X, const struct SymmetryCheckpointInfo *source,
                             const struct SymmetryTEHamiltonian *view, int continuation)
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
      fprintf(fp, "te_hamiltonian=%s\nte_steps=%u\nexpand_coef=%d\n",
              view ? "time_dependent" : "static", X->Def.Lanczos_max, X->Def.Param.ExpandCoef);
      /* Imports already record their complete provenance in the dispatcher. */
      if (continuation) fprintf(fp, "source_method=%" PRIu64 "\nsource_state=%" PRIu64
                                "\nsource_step=%" PRIu64 "\nsource_time=%.17g\ncontinuation=yes\n",
                                source->source_method, source->state_index, source->step, source->time);
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
                ((X->Bind.Def.iOutputEigenVec || X->Bind.Def.iReStart == RESTART_OUT ||
                  X->Bind.Def.iReStart == RESTART_INOUT) && X->Bind.Def.Param.OutputInterval <= 0);
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
  int is_v3 = 0, legacy_step = 0;
  double restored_time = 0;
  if (te_restore(&(X->Bind), view, dt, &step_initial, &restored_time, &is_v3, &legacy_step) != 0)
    return -1;
  if (!is_v3 && X->Bind.Def.iInputEigenVec == FALSE) {
    fprintf(stderr, "Error: A file of Inputvector is not inputted.\n");
    return -1;
  }
  if (!is_v3) {
    //input v1
    fprintf(stdoutMPI, "%s","An Initial Vector is inputted.\n");
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, c_InputEigenVectorStart, "a");
    invalid = GetFileNameByKW(KWSpectrumVec, &defname) != 0 || defname == NULL;
    int length = invalid ? -1 : snprintf(sdt, sizeof(sdt), "%s_rank_%d.dat", defname, myrank);
    if (SumMPI_i(length < 0 || (size_t)length + strlen(cParentOutputFolder) >= sizeof(sdt)) != 0) {
      fprintf(stdoutMPI, "Error: missing or overlong TE SpectrumVec filename.\n");
      return -1;
    }
    struct RestartVectorInfo provenance;
    enum RestartFileFormat format;
    int source_ranks;
    int imported = RestartImportVector(&(X->Bind), defname, v1, &provenance, &format, &source_ranks);
    if (imported == -1) return -1;
    if (imported == 0) {
      step_initial = 0;
      if (X->Bind.Def.iFlgSymmetryBasis) {
        /* Old sector OutputEigenVec retains the imported state label. v3
         * starts a new solver set at zero, independently of its provenance. */
        if (format != RESTART_FORMAT_V3) source.state_index = provenance.set;
#ifdef HPHI_TEST_SYMMETRY_PROBE
        double norm2 = 0;
        for (unsigned long i = 1; i <= X->Bind.Check.idim_max; ++i)
          norm2 += creal(conj(v1[i])*v1[i]);
        TE_AGREE_OR_RETURN(SymmetryProbeWriteVector(&(X->Bind), v1, "initial", 0, 0,
                                                    sqrt(SumMPI_d(norm2))));
#endif
        if (record_sector_te(&(X->Bind), NULL, view, 0) != 0) return -1;
      }
    } else { /* -2: same-rank legacy raw; keep ReStart=2/3 stored STEP. */
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
      int continuing = X->Bind.Def.iReStart == RESTART_INOUT || X->Bind.Def.iReStart == RESTART_IN;
      if (!continuing) step_initial = 0;
      int root_step = BcastMPI_i(0, step_initial);
      if (SumMPI_i(step_initial < (continuing ? 1 : 0) ||
                   (unsigned int)step_initial >= X->Bind.Def.Lanczos_max || step_initial != root_step) != 0) {
        fprintf(stdoutMPI, "Error: TE restart step must be consistent across ranks and within [0, Lanczos_max).\n");
        return -1;
      }
      if (continuing) legacy_step = step_initial;
    }
  }

  if (is_v3 && X->Bind.Def.iFlgSymmetryBasis) {
    source.source_method = TimeEvolution;
    source.step = step_initial; source.time = restored_time;
    if (record_sector_te(&(X->Bind), &source, view, 1) != 0) return -1;
  }
  int recover = is_v3 || (legacy_step > 0 &&
    (X->Bind.Def.iReStart == RESTART_INOUT || X->Bind.Def.iReStart == RESTART_IN));
  const struct RestartSeriesRule rule = {0, is_v3 ? step_initial : legacy_step};
  const struct RestartSeriesSpec series = {"#", -1, -1, 1, &rule, 0};

  if(X->Bind.Def.iOutputDataHead==1){
    sprintf(sdt_phys, "%s_%s", X->Bind.Def.CDataFileHead, cFileNameSS);
  }else{
    sprintf(sdt_phys, "%s", cFileNameSS);
  }
  if (recover) {
    if (RestartRecoverTimeSeries(&(X->Bind), sdt_phys, &series) != 0) return -1;
  } else {
    if (te_open(sdt_phys, "w", &fp) != 0) return -1;
    fprintf(fp, "%s", cLogSS);
    TE_AGREE_OR_RETURN(GreenOutputCloseChecked(fp));
  }

  if(X->Bind.Def.iOutputDataHead==1){
    sprintf(sdt_norm, "%s_%s", X->Bind.Def.CDataFileHead, cFileNameNorm);
  }else{
    sprintf(sdt_norm, "%s", cFileNameNorm);
  }
  if (recover) {
    if (RestartRecoverTimeSeries(&(X->Bind), sdt_norm, &series) != 0) return -1;
  } else {
    if (te_open(sdt_norm, "w", &fp) != 0) return -1;
    fprintf(fp, "%s", cLogNorm);
    TE_AGREE_OR_RETURN(GreenOutputCloseChecked(fp));
  }

  if(X->Bind.Def.iOutputDataHead==1){
    sprintf(sdt_flct, "%s_%s", X->Bind.Def.CDataFileHead, cFileNameFlct);
  }else{
    sprintf(sdt_flct, "%s", cFileNameFlct);
  }
  if (recover) {
    if (RestartRecoverTimeSeries(&(X->Bind), sdt_flct, &series) != 0) return -1;
  } else {
    if (te_open(sdt_flct, "w", &fp) != 0) return -1;
    fprintf(fp, "%s", cLogFlct);
    TE_AGREE_OR_RETURN(GreenOutputCloseChecked(fp));
  }

  if (recover) {
    if (te_recover_green(&(X->Bind), &rule) != 0) return -1;
  } else {
    TE_AGREE_OR_RETURN(GreenOutputInitializeAggregateFiles(&(X->Bind)));
  }
  int no_remaining = is_v3 && (unsigned int)step_initial >= X->Bind.Def.Lanczos_max;
  if (is_v3) Time = restored_time;
  if (no_remaining) fprintf(stdoutMPI, "TE: no remaining steps (saved step %d).\n", step_initial);

  int iInterAllOffDiagonal_org = X->Bind.Def.NInterAll_OffDiagonal;
  int iTransfer_org = X->Bind.Def.EDNTransfer;
  /* Raw no-op has only vector output, so no H psi reconstruction is needed.
   * In particular Laser has no saved-step table-length contract. Sector
   * output still needs its preceding Hamiltonian/digest selected. */
  if (step_initial > 0 && (!no_remaining || X->Bind.Def.iFlgSymmetryBasis)) {
    /* Raw TE files store the next step. Restore H|psi> at the preceding
     * time point: MultiplyForTEM consumes both psi (v1) and H psi (v0).
     * Without this, a restarted run silently loses its first Taylor term. */
    step_i = X->Bind.Def.istep = step_initial - 1;
    if (view) {
      if (SelectSymmetryTEHamiltonian(&(X->Bind), view, step_initial - 1, restored_time) != 0 ||
          RebuildSymmetryTEPlan(&(X->Bind)) != 0) return -1;
    } else if (X->Bind.Def.NLaser > 0) {
      double previous_time = X->Bind.Def.Param.Tinit + (step_initial - 1) * dt;
      Time = is_v3 ? restored_time : X->Bind.Def.Param.Tinit + step_initial * dt;
      if (!isfinite(previous_time) || !isfinite(Time)) return -1;
      TE_AGREE_OR_RETURN(TransferWithPeierls(&(X->Bind), previous_time));
    } else if (X->Bind.Def.NTETransferMax > 0) {
      MakeTEDTransfer(&(X->Bind), step_initial - 1);
    } else if (X->Bind.Def.NTEInterAllMax > 0) {
      MakeTEDInterAll(&(X->Bind), step_initial - 1);
    }
    for (unsigned long i = 1; i <= X->Bind.Check.idim_max; ++i) v0[i] = v1[i];
    TE_AGREE_OR_RETURN(expec_energy_flct(&(X->Bind)));
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
      TE_AGREE_OR_RETURN(expec_energy_flct(&(X->Bind)));
    }
    else if(X->Bind.Def.NLaser !=0) {
      TE_AGREE_OR_RETURN(TransferWithPeierls(&(X->Bind), Time));
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
    TE_AGREE_OR_RETURN(MultiplyForTEM(&(X->Bind)));
    //Add Diagonal Parts
    //Multiply Diagonal
    TE_AGREE_OR_RETURN(expec_energy_flct(&(X->Bind)));

    if(!view && X->Bind.Def.NLaser >0 ) Time+=dt;
    if (te_open(sdt_phys, "a", &fp) != 0) {
      return -1;
    }
    fprintf(fp, "%.16lf  %.16lf %.16lf %.16lf %.16lf %d\n", Time, X->Bind.Phys.energy, X->Bind.Phys.var,
            X->Bind.Phys.doublon, X->Bind.Phys.num, step_i);
    TE_AGREE_OR_RETURN(GreenOutputCloseChecked(fp));

    if (te_open(sdt_norm, "a", &fp) != 0) {
      return -1;
    }
    fprintf(fp, "%.16lf %.16lf %d\n", Time, global_norm, step_i);
    TE_AGREE_OR_RETURN(GreenOutputCloseChecked(fp));

    if (te_open(sdt_flct, "a", &fp) != 0) {
      return -1;
    }
    fprintf(fp, "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %d\n", Time,X->Bind.Phys.num,X->Bind.Phys.num2, X->Bind.Phys.doublon,X->Bind.Phys.doublon2, X->Bind.Phys.Sz,X->Bind.Phys.Sz2,step_i);
    TE_AGREE_OR_RETURN(GreenOutputCloseChecked(fp));



    if (step_i % step_spin == 0) {
      TE_AGREE_OR_RETURN(expec_cisajs(&(X->Bind), v1));
      TE_AGREE_OR_RETURN(expec_cisajscktaltdc(&(X->Bind), v1));
      TE_AGREE_OR_RETURN(expec_nbodyg(&(X->Bind), v1));
      TE_AGREE_OR_RETURN(expec_anomalousg(&(X->Bind), v1));
    }
    if (X->Bind.Def.iOutputEigenVec == TRUE) {
      if (step_i % X->Bind.Def.Param.OutputInterval == 0) {
        if (X->Bind.Def.iFlgSymmetryBasis) {
          if (write_sector_te_vector(&(X->Bind), step_i, Time, source.state_index, 0) != 0) return -1;
        } else if (write_raw_te_vector(&(X->Bind), step_i, step_i + 1) != 0) return -1;
      }
    }
    if ((X->Bind.Def.iReStart == RESTART_OUT || X->Bind.Def.iReStart == RESTART_INOUT) &&
        (step_i + 1) % X->Bind.Def.Param.OutputInterval == 0)
      if (te_checkpoint_write(&(X->Bind), step_i + 1, Time) != 0) return -1;
  }
  if (!no_remaining && (X->Bind.Def.iReStart == RESTART_OUT || X->Bind.Def.iReStart == RESTART_INOUT))
    if (te_checkpoint_write(&(X->Bind), step_i, Time) != 0) return -1;
  if (X->Bind.Def.iOutputEigenVec == TRUE) {
    if (X->Bind.Def.iFlgSymmetryBasis) {
      if (write_sector_te_vector(&(X->Bind), no_remaining ? step_initial : step_i - 1,
                                 no_remaining ? restored_time : Time, source.state_index, 1) != 0) return -1;
    } else if (write_raw_te_vector(&(X->Bind), no_remaining ? step_initial : rand_i, step_i) != 0) return -1;
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
