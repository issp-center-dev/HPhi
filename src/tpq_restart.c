#include <inttypes.h>
#include <stdlib.h>
#include <string.h>

#include "Common.h"
#include "green_output.h"
#include "restart_output.h"
#include "tpq_restart.h"
#include "wrapperMPI.h"
#ifdef HPHI_TEST_SYMMETRY_PROBE
#include "symmetry_fulldiag_probe.h"
#endif

static uint64_t keep_below(const struct TPQRestartPlan *plan, int set)
{
  if (plan->sources[set] == RESTART_SOURCE_V3) {
    const struct RestartManifestRow *row =
        RestartManifestFind(&plan->manifest, (uint64_t)set);
    return row == NULL ? 0 : row->step;
  }
  return (uint64_t)plan->legacy_steps[set];
}

int TPQRestartPlanSets(struct BindStruct *X, int numave,
                       struct TPQRestartPlan *plan)
{
  int failed;
  if (plan == NULL) return -1;
  memset(plan, 0, sizeof(*plan));
  plan->numave = numave;
  failed = X == NULL || numave <= 0;
  if (!failed) {
    plan->sources = calloc((size_t)numave, sizeof(*plan->sources));
    plan->legacy_steps = calloc((size_t)numave, sizeof(*plan->legacy_steps));
    failed = plan->sources == NULL || plan->legacy_steps == NULL;
  }
  if (SumMPI_i(failed) != 0) {
    fprintf(stdoutMPI, "Error: TPQ restart plan allocation failed.\n");
    TPQRestartFreePlan(plan);
    return -1;
  }
  if (X->Def.iReStart != RESTART_INOUT && X->Def.iReStart != RESTART_IN) {
    plan->boundaries_validated = 1;
    return 0;
  }
  if (RestartManifestLoad(X, &plan->manifest, &plan->manifest_exists) != 0)
    goto fail;
  if (plan->manifest_exists && plan->manifest.method != X->Def.iCalcType) {
    fprintf(stdoutMPI,
            "Error: restart manifest belongs to another method.\n");
    goto fail;
  }
  for (int set = 0; set < numave; ++set) {
    if (RestartClassifySet(X,
            plan->manifest_exists ? &plan->manifest : NULL,
            (uint64_t)set, cFileNameInputVector, &plan->sources[set],
            &plan->legacy_steps[set]) != 0)
      goto fail;
    if (plan->sources[set] == RESTART_SOURCE_NONE &&
        X->Def.iReStart == RESTART_IN) {
      fprintf(stdoutMPI,
              "Error: ReStart=3 requires a checkpoint for every set "
              "(set %d is missing).\n", set);
      goto fail;
    }
    plan->any_resume |= plan->sources[set] != RESTART_SOURCE_NONE;
  }
  /* Each set is an independent publication unit. Its headers have already
   * been checked against the runtime logical basis; compatible metadata
   * layouts may differ between sets after an interrupted ReStart=2. */
  plan->boundaries_validated = 1;
  return 0;
fail:
  TPQRestartFreePlan(plan);
  return -1;
}

static int recover_green_aggregates(struct BindStruct *X,
                                    const struct RestartSeriesRule *rules,
                                    int nrule)
{
  const unsigned int counts[] = {
    X->Def.NCisAjt, X->Def.NCisAjtCkuAlvDC, X->Def.NTBody,
    X->Def.NFBody, X->Def.NSBody, X->Def.NNBodyG, X->Def.NAnomalousG
  };
  struct RestartSeriesSpec spec = {"#", 0, 1, nrule, rules, 1};
  for (int kind = GreenOutputOneBody; kind <= GreenOutputAnomalous; ++kind) {
    if (counts[kind] == 0) continue;
    char name[D_FileNameMax];
    int failed = GreenOutputFileName(X, (GreenOutputKind)kind, name) != 0;
    if (SumMPI_i(failed) != 0) return -1;
    if (RestartRecoverTimeSeries(X, name, &spec) != 0) return -1;
  }
  return 0;
}

int TPQRestartRecoverOutputs(struct BindStruct *X,
                             const struct TPQRestartPlan *plan,
                             int aggregate)
{
  int failed = X == NULL || plan == NULL || !plan->boundaries_validated;
  if (SumMPI_i(failed) != 0) return -1;
  if (!plan->any_resume) return 0;
  struct RestartSeriesRule *rules =
      calloc((size_t)plan->numave, sizeof(*rules));
  failed = rules == NULL;
#ifdef HPHI_TEST_SYMMETRY_PROBE
  if (SymmetryFullDiagTestFail("restart_recover_alloc")) failed = 1;
#endif
  if (SumMPI_i(failed) != 0) {
    free(rules);
    fprintf(stdoutMPI, "Error: restart recovery allocation failed.\n");
    return -1;
  }
  int nrule = 0;
  for (int set = 0; set < plan->numave; ++set) {
    if (plan->sources[set] == RESTART_SOURCE_NONE) continue;
    rules[nrule].set = (uint64_t)set;
    rules[nrule].keep_below_step = keep_below(plan, set);
    ++nrule;
  }
  if (aggregate) {
    const GreenOutputTPQDataKind kinds[] = {
      GreenOutputTPQDataSS, GreenOutputTPQDataNorm, GreenOutputTPQDataFlct
    };
    struct RestartSeriesSpec spec = {"#", 0, 1, nrule, rules, 1};
    for (int kind = 0; kind < 3; ++kind) {
      char name[D_FileNameMax];
      failed = GreenOutputTPQDataFileName(X, kinds[kind], name) != 0;
      if (SumMPI_i(failed) != 0 ||
          RestartRecoverTimeSeries(X, name, &spec) != 0) {
        free(rules);
        return -1;
      }
    }
  } else {
    const char *patterns[] = {
      cFileNameSSRand, cFileNameNormRand, cFileNameFlctRand
    };
    for (int set = 0; set < plan->numave; ++set) {
      if (plan->sources[set] == RESTART_SOURCE_NONE) continue;
      struct RestartSeriesRule one = {(uint64_t)set, keep_below(plan, set)};
      struct RestartSeriesSpec spec = {"#", -1, -1, 1, &one, 0};
      for (int kind = 0; kind < 3; ++kind) {
        char name[D_FileNameMax];
        int n;
        if (X->Def.iOutputDataHead == 1) {
          n = snprintf(name, sizeof(name), "%s_", X->Def.CDataFileHead);
          if (n >= 0 && n < (int)sizeof(name))
            n += snprintf(name + n, sizeof(name) - (size_t)n,
                          patterns[kind], set);
        } else {
          n = snprintf(name, sizeof(name), patterns[kind], set);
        }
        failed = n < 0 || n >= (int)sizeof(name);
        if (SumMPI_i(failed) != 0 ||
            RestartRecoverTimeSeries(X, name, &spec) != 0) {
          free(rules);
          return -1;
        }
      }
    }
  }
  if (GreenOutputUsesAggregate(X) &&
      recover_green_aggregates(X, rules, nrule) != 0) {
    free(rules);
    return -1;
  }
  free(rules);
  return 0;
}

int TPQRestartWriteCheckpoint(struct BindStruct *X, int rand_i,
                              int next_step, double first_norm,
                              double inv_temp, double delta_tau)
{
  uint64_t set = (uint64_t)rand_i;
  double complex *components[1] = {v1};
  double complex **vectors[1] = {components};
  struct RestartVectorInfo info = {
    (uint64_t)X->Def.iCalcType, set, (uint64_t)next_step, 0.0,
    (uint64_t)NumAve, 0, first_norm, inv_temp, delta_tau
  };
  struct RestartUnit unit = {1, &set, 1};
  return RestartWriteUnit(X, &unit,
                          (double complex *const *const *)vectors, &info, 0);
}

void TPQRestartFreePlan(struct TPQRestartPlan *plan)
{
  if (plan == NULL) return;
  free(plan->sources);
  free(plan->legacy_steps);
  RestartManifestFree(&plan->manifest);
  memset(plan, 0, sizeof(*plan));
}
