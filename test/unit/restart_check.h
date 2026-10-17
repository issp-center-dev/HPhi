#ifndef HPHI_TEST_RESTART_CHECK_H
#define HPHI_TEST_RESTART_CHECK_H

#include <stdio.h>

/* Release builds define NDEBUG: never put the call under test inside assert(). */
#define CHECK(cond) do { \
  if (!(cond)) { \
    fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
    return 1; \
  } \
} while (0)

#endif
