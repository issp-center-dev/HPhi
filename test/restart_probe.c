#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "restart_probe.h"
/* HPHI_TEST_RESTART_ABORT_AFTER="<set>:<step>" stops the run right after that checkpoint is published. */
int RestartTestAbortAfter(uint64_t set, uint64_t step)
{
  const char *spec = getenv("HPHI_TEST_RESTART_ABORT_AFTER");
  unsigned long long want_set, want_step;
  if (spec == NULL || sscanf(spec, "%llu:%llu", &want_set, &want_step) != 2) return 0;
  return set == want_set && step == want_step;
}
