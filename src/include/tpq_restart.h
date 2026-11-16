#ifndef HPHI_TPQ_RESTART_H
#define HPHI_TPQ_RESTART_H

#include "restart_checkpoint.h"
#include "restart_manifest.h"

struct BindStruct;

struct TPQRestartPlan {
  int numave;
  enum RestartSource *sources;
  int *legacy_steps;
  int any_resume;
  int boundaries_validated;
  struct RestartManifest manifest;
  int manifest_exists;
};

int TPQRestartPlanSets(struct BindStruct *X, int numave,
                       struct TPQRestartPlan *plan);
int TPQRestartRecoverOutputs(struct BindStruct *X,
                             const struct TPQRestartPlan *plan,
                             int aggregate);
int TPQRestartWriteCheckpoint(struct BindStruct *X, int rand_i,
                              int next_step, double first_norm,
                              double inv_temp, double delta_tau);
void TPQRestartFreePlan(struct TPQRestartPlan *plan);

#endif
