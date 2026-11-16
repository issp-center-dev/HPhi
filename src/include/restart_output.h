#ifndef HPHI_RESTART_OUTPUT_H
#define HPHI_RESTART_OUTPUT_H

#include <stddef.h>
#include <stdint.h>

struct BindStruct;

struct RestartSeriesRule {
  uint64_t set;
  uint64_t keep_below_step;
};

struct RestartSeriesSpec {
  const char *header_prefix;
  int set_column;
  int step_column;
  int nrule;
  const struct RestartSeriesRule *rules;
  int drop_unlisted_sets;
};

int RestartRecoverTimeSeriesFile(const char *path,
                                 const struct RestartSeriesSpec *spec,
                                 char *message, size_t size);

int RestartRecoverTimeSeries(const struct BindStruct *X, const char *relname,
                             const struct RestartSeriesSpec *spec);

#endif
