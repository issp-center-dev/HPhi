#include <stdio.h>

#include "Common.h"
#include "CalcTime.h"
#include "restart_output.h"
#include "wrapperMPI.h"

int RestartRecoverTimeSeries(const struct BindStruct *X, const char *relname,
                             const struct RestartSeriesSpec *spec)
{
  char path[D_FileNameMax];
  char message[256] = "";
  int failed = 0;

  StartTimer(7120);
  if (myrank == 0) {
    int length = (X == NULL || relname == NULL || cParentOutputFolder == NULL)
                     ? -1
                     : snprintf(path, sizeof(path), "%s%s",
                                cParentOutputFolder, relname);
    failed = length < 0 || length >= (int)sizeof(path) ||
             RestartRecoverTimeSeriesFile(path, spec, message,
                                          sizeof(message)) != 0;
    if (failed)
      fprintf(stdoutMPI, "Error: output recovery failed: %s\n", message);
  }
  failed = SumMPI_i(failed) != 0;
  StopTimer(7120);
  return failed ? -1 : 0;
}
